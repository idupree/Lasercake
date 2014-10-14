/*

    Copyright Eli Dupree and Isaac Dupree, 2014

    This file is part of Lasercake.

    Lasercake is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Lasercake is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Lasercake.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <array>
#include <vector>

namespace impl {
  struct unsorted {};
  template<typename QueueEntry, typename ComesAfter>
  struct iterator_queue {
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, ComesAfter> data;
    QueueEntry& front() { return data.top(); }
    void pop() { return data.pop(); }
    void push(QueueEntry const& e) { return data.push(e); }
    bool empty()const { return data.empty(); }
  }
  template<typename QueueEntry>
  struct iterator_queue<QueueEntry, unsorted> {
    std::queue<QueueEntry> data;
    QueueEntry& front() { return data.front(); }
    void pop() { return data.pop(); }
    void push(QueueEntry const& e) { return data.push(e); }
    bool empty()const { return data.empty(); }
  }
}

// A map from KeyType to (finite union of regularized intervals, cf. Robust Regularized Set Operations on Polyhedra, Beat Bruderlin),
// with fast lookup by overlapping point, interval, etc.
template<typename KeyType, typename DomainType>
class bidirectional_interval_map {
public:
  struct interval_type {
    std::array<DomainType, 2> bounds;
  };
  struct value_type {
    interval_type interval;
    KeyType key;
  }
  struct tiebroken_domain_type {
    tiebroken_domain_type(DomainType value, uint64_t tiebreaker):value(value),tiebreaker(tiebreaker){}
    DomainType value;
    uint64_t tiebreaker;
  };
private:
  struct node {
    node(bool splits_based_on_which_bound):
      splits_based_on_which_bound(splits_based_on_which_bound),
      num_descendant_intervals(0),
      children({{nullptr,nullptr}}),
      values_here(){}
    bool splits_based_on_which_bound;
    size_t num_descendant_intervals;
    
    // only for non leaves
    tiebroken_domain_type child_separator;
    std::array<std::unique_ptr<node>, 2> children; // the first has the lower values, the second the higher values
    
    // only for leaves
    std::vector<value_type> values_here;
    void insert_value(value_type const& v) {
      for (value_type const& v2 : values_here) { assert(v2 != v); }
      values_here.push_back(v);
    }
    void erase_value(value_type const& v) {
      for (value_type& v2 : values_here) {
        if (v2 == v) {
          v2 = values_here.back();
          values_here.pop_back();
          return;
        }
      }
      assert(false);
    }
  };
  
public:
  
  template<typename ComesAfter = impl::unsorted>
  class iterator {
  public:
    iterator(impl::unsorted = impl::unsorted()):queue(),
      bound_minima_exist({{false,false}}),bound_maxima_exist({{false,false}}),min_length_exists(false){}
    iterator(std::enable_if_t<!std::is_same<ComesAfter, impl::unsorted>::value, ComesAfter const&> comes_after):comes_after(comes_after),queue(comes_after),
      bound_minima_exist({{false,false}}),bound_maxima_exist({{false,false}}),min_length_exists(false){}

    value_type& operator*() const {
      caller_error_if(queue.empty(), "can't dereference an empty iterator");
      return *queue.front().value;
    }
    value_type* operator->() const {
      caller_error_if(queue.empty(), "can't dereference an empty iterator");
      return queue.front().value;
    }
    value_type* operator++(int) {
      value_type* result = queue.front().value;
      ++*this;
      return result;
    }
    iterator& operator++() {
      queue.pop();
      advance_to_a_returnable();
      return *this;
    }

    bool operator==(iterator const& other) const {
      if (queue.empty()) { return other.queue.empty(); }
      if (other.queue.empty()) { return false; }
      return queue.front() == other.queue.front();
    }
    bool operator!=(iterator const& other) const { return !(*this == other); }
    
  private:
    struct queue_entry {
      node* n;
      value_type* value;
      std::array<interval_type, 2> possible_bounds;
    };
    impl::iterator_queue<queue_entry, ComesAfter> queue;
    ComesAfter comes_after;
    std::array<bool, 2> bound_minima_exist;
    std::array<bool, 2> bound_maxima_exist;
    bool min_length_exists;
    
    
    void queue_root(node* root) {
      queue_entry e;
      e.n = root;
      if (!e.n) return;
      e.possible_bounds[0][0] = neginf;
      e.possible_bounds[0][1] = inf;
      e.possible_bounds[1][0] = neginf;
      e.possible_bounds[1][1] = inf;
      queue.push(e);
    }
    void queue_node(queue_entry const& parent, bool use_second_child) {
      queue_entry e = parent;
      e.n = parent.n->children[use_second_child].get();
      if (!e.n) return;
      if (bound_minima_exist[parent.n->splits_based_on_which_bound] &&  use_second_child && (parent.n->child_separator.value < bound_minima[parent.n->splits_based_on_which_bound])) return;
      if (bound_maxima_exist[parent.n->splits_based_on_which_bound] && !use_second_child && (parent.n->child_separator.value > bound_maxima[parent.n->splits_based_on_which_bound])) return;
      e.possible_bounds[parent.n->splits_based_on_which_bound].bounds[!use_second_child] = parent.n->child_separator.value;
      if (min_length_exists && interval_type(e.possible_bounds[0].bounds[0], e.possible_bounds[1].bounds[1]).length_less_than(min_length)) return;
      queue.push(e);
    }
    void queue_value(node* n, value_type* value) {
      queue_entry e;
      e.n = n;
      e.value = value;
      std::array<DomainType, 2> const& bounds = value->interval.bounds;
      if (bound_minima_exist[0] && (bounds[0] < bound_minima[0])) return;
      if (bound_maxima_exist[0] && (bounds[0] > bound_maxima[0])) return;
      if (bound_minima_exist[1] && (bounds[1] < bound_minima[1])) return;
      if (bound_maxima_exist[1] && (bounds[1] > bound_maxima[1])) return;
      if (min_length_exists && value->interval.length_less_than(min_length)) return;
      queue.push(e);
    }
    void advance_to_a_returnable() {
      while (!queue.front().value) {
        queue_node(queue.front(), 0);
        queue_node(queue.front(), 1);
        for (value_type& v : queue.front().n->values_here) {
          queue_value(queue.front().n, &v);
        }
        queue.pop();
        if (queue.empty()) { return; }
      }
    }
    
    friend bidirectional_interval_map;
  };
  
private:
  class comes_after_by_bound_t {
  private:
    typedef iterator<comes_after_by_bound_t> iter_type;
  public:
    comes_after_by_bound_t(bidirectional_interval_map const* map, bool bound, bool descending):
      map(map),bound(bound),descending(descending){}
    
    // For a priority_queue, < is descending and > is ascending.
    // Also, a containing interval should always come before the contained interval, so we use
    // the *earliest* end of possible_bounds[bound].
    // So, for <, that's possible_bounds[bound][true] and for >, that's possible_bounds[bound][false].
    // a.possible_bounds[bound][]
    DomainType const& get_value(typename iter_type::queue_entry const& v) {
      return v.possible_bounds[bound][descending];
    }
    DomainType const& get_value(value_type const& v) {
      return v.interval.bounds[descending];
    }
    DomainType const& get_value(tiebroken_domain_type const& v) {
      return v.value;
    }
    DomainType const& get_tiebreaker(value_type const& v) {
      return map->get_tiebreaker(v.key);
    }
    DomainType const& get_tiebreaker(tiebroken_domain_type const& v) {
      return v.tiebreaker;
    }
    bool operator()(typename iter_type::queue_entry const& a, typename iter_type::queue_entry const& b) {
      if (get_value(a) != get_value(b)) {
        return (get_value(a) < get_value(b)) == descending;
      }
      if (a.value && b.value) {
        const uint64_t a_t = get_tiebreaker(a);
        const uint64_t b_t = get_tiebreaker(b);
        if (a_t != b_t) {
          return (a_t < b_t) == descending;
        }
      }
      else if (a.value || b.value) {
        return bool(b.value);
      }
      return false;
    }
    template<typename T1, typename T2>
    std::enable_if_t<!(std::is_same<T1, typename iter_type::queue_entry>::value && std::is_same<T2, typename iter_type::queue_entry>::value),
    bool> operator()(T1 const& a, T2 const& b) {
      if (get_value(a) != get_value(b)) {
        return (get_value(a) < get_value(b)) == descending;
      }
      const uint64_t a_t = get_tiebreaker(a);
      const uint64_t b_t = get_tiebreaker(b);
      if (a_t != b_t) {
        return (a_t < b_t) == descending;
      }
      return false;
    }
  private:
    bidirectional_interval_map const* map;
    bool bound;
    bool descending;
  }
    
  iterator& finish_iterator(iterator& i) {
    commit_changes();
    i.queue_root(root.get());
    i.advance_to_a_returnable();
  }
public:
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_containing(DomainType const& k, ComesAfter c = ComesAfter()) {
    iterator<ComesAfter> result(c);
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k;
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k;
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_overlapping(interval_type const& k, ComesAfter c = ComesAfter()) {
    iterator<ComesAfter> result(c);
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[1];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[0];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_subsuming(interval_type const& k, ComesAfter c = ComesAfter()) {
    iterator<ComesAfter> result(c);
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[0];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[1];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_subsumed_by(interval_type const& k, ComesAfter c = ComesAfter()) {
    iterator<ComesAfter> result(c);
    result.bound_minima_exist[0] = true;
    result.bound_minima[0] = k.bounds[0];
    result.bound_maxima_exist[1] = true;
    result.bound_maxima[1] = k.bounds[1];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_equal_to(interval_type const& k, ComesAfter c = ComesAfter()) {
    iterator<ComesAfter> result(c);
    result.bound_minima_exist[0] = true;
    result.bound_minima[0] = k.bounds[0];
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[0];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[1];
    result.bound_maxima_exist[1] = true;
    result.bound_maxima[1] = k.bounds[1];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> end()const {
    return iterator<ComesAfter>();
  }
  
  template<typename ComesAfter = impl::unsorted>
  boost::iterator_range<iterator<ComesAfter>> range_containing (   DomainType const& k, ComesAfter c = ComesAfter()) {
    return boost::iterator_range<iterator<ComesAfter>>(begin_containing (k), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  boost::iterator_range<iterator<ComesAfter>> range_overlapping(interval_type const& k, ComesAfter c = ComesAfter()) {
    return boost::iterator_range<iterator<ComesAfter>>(begin_overlapping(k), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  boost::iterator_range<iterator<ComesAfter>> range_subsuming  (interval_type const& k, ComesAfter c = ComesAfter()) {
    return boost::iterator_range<iterator<ComesAfter>>(begin_subsuming  (k), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  boost::iterator_range<iterator<ComesAfter>> range_subsumed_by(interval_type const& k, ComesAfter c = ComesAfter()) {
    return boost::iterator_range<iterator<ComesAfter>>(begin_subsumed_by(k), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  boost::iterator_range<iterator<ComesAfter>> range_equal_to   (interval_type const& k, ComesAfter c = ComesAfter()) {
    return boost::iterator_range<iterator<ComesAfter>>(begin_equal_to   (k), end<ComesAfter>());
  }
  
  comes_after_by_bound_t comes_after_by_bound(bool bound, bool descending)const {
    return comes_after_by_bound_t(this, bound, descending);
  }
private:
  
  
  // The tree maintains these invariants:
  // Nodes always alternate splits_based_on_which_bound.
  // No non-leaf directly contains any intervals.
  // 1) No leaf contains more than max_leaf_size intervals (it should be split).
  // 2) No non-leaf has <= max_leaf_size descendant values (it should be collapsed into a leaf).
  // 3) No non-leaf has more than twice as many descendants on one side as on the other.
  // (2) and (3) together bound the tree height at log_3(N).
  static const max_leaf_size = 6;
  
  static void collapse_node(node& n, std::unique_ptr<node>& descendant) {
    for (value_type& v : descendant.values_here) {
      n.insert_value(v);
    }
    if (descendant->children[0]) {
      collapse_node(n, descendant->children[0]);
      collapse_node(n, descendant->children[1]);
    }
    descendant = nullptr;
  }
  static void is_in_sorted_vector(value_type const& what, std::vector<value_type const*> const& vec, comes_after_by_bound_t const& comes_after) {
    if (vec.empty()) { return false; }
    size_t min = 0;
    size_t maxplus1 = vec.size();
    while (min+1 < maxplus1) {
      size_t mid = (min + maxplus1)>>1;
      if (comes_after(*vec[mid], what)) { maxplus1 = mid; }
      else                              { min      = mid; }
    }
    return *vec[min] == what;
  }
  // Insert M elements and erase P elements.
  // O(M+P+log(N)) amortized.
  static void insert_and_or_erase_impl(bool is_root, std::unique_ptr<node>& n, std::vector<value_type const*> const& values_to_insert, std::vector<value_type const*> const& values_to_erase) {
    if (!n) {
      assert(is_root);
      n.reset(new node(false));
    }
    n->num_descendant_intervals += values_to_insert.size();
    n->num_descendant_intervals -= values_to_erase .size();
    assert(is_root || (n->num_descendant_intervals >= min_leaf_size));
    if (n->num_descendant_intervals == 0) {
      assert(is_root);
      n = nullptr;
    }
    
    if (n->children[0]) { // i.e. "if non leaf"
      if (n->num_descendant_intervals <= max_leaf_size) {
        collapse_node(*n, n->children[0]);
        collapse_node(*n, n->children[1]);
      }
      else {
        std::array<std::vector<value_type const*>, 2> values_to_insert_by_child;
        std::array<std::vector<value_type const*>, 2> values_to_erase_by_child;
        for (auto& v : values_to_insert_by_child) { v.reserve(values_to_insert.size()); }
        for (auto& v : values_to_erase_by_child ) { v.reserve(values_to_erase .size()); }
        comes_after_by_bound_t is_higher = comes_after_by_bound(n->splits_based_on_which_bound, false);
        for (value_type const* v : values_to_insert) {
          values_to_insert_by_child[is_higher(*v, n->child_separator)].push_back(v);
        }
        for (value_type const* v : values_to_erase ) {
          values_to_erase_by_child [is_higher(*v, n->child_separator)].push_back(v);
        }
  #define NEW_NUM_DESCENDANTS(which_child) (n->children[which_child]->num_descendant_intervals \
    + values_to_insert_by_child[which_child].size() \
    - values_to_erase_by_child[which_child].size())
    
        for (uint32_t which_child = 0; which_child < 2; ++which_child) {
          if (NEW_NUM_DESCENDANTS(which_child) > 2*NEW_NUM_DESCENDANTS(!which_child)) {
            // reversed sense because std::sort wants a comes-before.
            // We want the *backs* of these vectors to be closest to the separator.
            std::sort(values_to_insert_by_child[0], comes_after_by_bound(n->splits_based_on_which_bound,  true));
            std::sort(values_to_insert_by_child[1], comes_after_by_bound(n->splits_based_on_which_bound, false));
            std::sort(values_to_erase_by_child[which_child], comes_after_by_bound(n->splits_based_on_which_bound, !which_child));
            // Regular sense, but we want the *front* of the thief iterator to be closest to the separator.
            comes_after_by_bound_t is_further_from_separator = comes_after_by_bound(n->splits_based_on_which_bound, !which_child);
            iterator<comes_after_by_bound_t> thief(is_further_from_separator);
            comes_after_by_bound_t is_closer_to_separator = comes_after_by_bound(n->splits_based_on_which_bound, which_child);
            std::vector<value_type const*> stolen;
            if (which_child) {
              thief.bound_minima_exist[0] = true;
              thief.bound_minima[0] = n->child_separator.value;
            }
            else {
              thief.bound_maxima_exist[1] = true;
              thief.bound_maxima[1] = n->child_separator.value;
            }
            thief.queue_root(n->children[which_child].get());
            thief.advance_to_a_returnable();
            while (NEW_NUM_DESCENDANTS(which_child) > NEW_NUM_DESCENDANTS(!which_child) + 1) {
              const value_type const* closest_to_insert = values_to_insert_by_child[which_child].back();
              const int64_t not_quite_as_close = (which_child ? 1 : -1);
              const DomainType closest_to_insert_bound = tiebroken_domain_type(
                closest_to_insert->interval.bounds[n->splits_based_on_which_bound],
                get_tiebreaker(closest_to_insert->key));
              const DomainType closest_to_steal_bound  = tiebroken_domain_type(
                            thief->interval.bounds[n->splits_based_on_which_bound],
                get_tiebreaker(            thief->key));
              if (is_further_from_separator(closest_to_insert_bound, closest_to_steal_bound)) {
                if ()
                if (!is_in_sorted_vector(*thief, values_to_erase_by_child[which_child], is_closer_to_separator)) {
                  stolen.push_back(&*thief);
                  n->child_separator = closest_to_steal_bound;
                  n->child_separator.tiebreaker += not_quite_as_close;
                }
                ++thief;
              }
              else {
                n->child_separator = closest_to_insert_bound;
                n->child_separator.tiebreaker += not_quite_as_close;
                values_to_insert_by_child[!which_child].push_back(closest_to_insert);
                values_to_insert_by_child[which_child].pop_back();
              }
            }
            for (value_type const* v : stolen) {
              values_to_erase_by_child [ which_child].push_back(v);
              values_to_insert_by_child[!which_child].push_back(v);
            }
          }
        }
        for (uint32_t which_child = 0; which_child < 2; ++which_child) {
          insert_and_or_erase_impl(false, n->children[which_child], values_to_insert_by_child[which_child], values_to_erase_by_child[which_child]);
        }
      }
    }
    else {
      for (value_type const* v : values_to_erase) {
        n->erase_value(*v);
      }
      for (value_type const* v : values_to_insert) {
        n->insert_value(*v);
      }
      assert (n->num_descendant_intervals == n->values_here.size());
      if (n->num_descendant_intervals > max_leaf_size) {
        n->children[0].reset(new node(!n->splits_based_on_which_bound));
        n->children[1].reset(new node(!n->splits_based_on_which_bound));
        
        // Ascending order. Reversed sense because std::sort wants a comes-before.
        std::sort(n->values_here, comes_after_by_bound(n->splits_based_on_which_bound, true));
        size_t median_idx = n->values_here.size() >> 1;
        n->child_separator = tiebroken_domain_type(
          n->values_here[median_idx].interval.bounds[n->splits_based_on_which_bound],
          get_tiebreaker(n->values_here[median_idx].key) - 1);
        std::array<std::vector<value_type*>, 2> values_to_insert_by_child;
        
        for (size_t i = 0; i < n->values_here.size(); ++i) {
          values_to_insert_by_child[i >= median_idx].push_back(&n->values_here[i]);
        }
        for (uint32_t which_child = 0; which_child < 2; ++which_child) {
          insert_and_or_erase_impl(false, n->children[which_child], values_to_insert_by_child[which_child], values_to_erase_by_child[which_child]);
        }
        n->values_here.clear();
      }
    }
  }
  
  void commit_changes() {
    std::vector<value_type const*> values_to_insert;
    std::vector<value_type const*> values_to_erase;
    
    for (value_type const& v : to_be_inserted) {
      values_to_insert.push_back(&v);
    }
    for (value_type const& v : to_be_erased) {
      values_to_erase.push_back(&v);
    }
    
    insert_and_or_erase_impl(true, root, values_to_insert, values_to_erase);
    
    to_be_inserted.clear();
    to_be_erased.clear();
  }
  void queue_erase(value_type const& v) {
    if (!to_be_inserted.erase(v)) { const auto p = to_be_erased.insert(v); assert(p.second); }
  }
  void queue_insert(value_type const& v) {
    if (!to_be_erased.erase(v)) { const auto p = to_be_inserted.insert(v); assert(p.second); }
  }
  void queue_imposition(value_type const& v, bool inserting) {
    caller_correct_if(v.interval.bounds[1] > v.interval.bounds[0], "intervals must have positive length");
    interval_type new_internal_interval = v.interval;
    
    key_metadata& m = metadata_by_key[v.key];
    if (!m.tiebreaker) {
      m.tiebreaker = next_tiebreaker;
      next_tiebreaker += 2;
    }
    
    bool already_filled;
    auto i = boost::prior(m.transitions.lower_bound(v.interval.bounds[0]));
    if (i == m.transitions.begin()) {
      already_filled = false;
    }
    else {
      const auto last_iter = boost::prior(i);
      already_filled = last_iter->second;
      if (last_iter->second) {
        queue_erase(last_iter->first, boost::next(last_iter)->first);
        if (inserting) {
          new_internal_interval.bounds[0] = last_iter->first;
        }
        else {
          queue_insert(v.key, last_iter->first, v.interval.bounds[0]);
        }
      }
    }
    if (already_filled != inserting) {
      m.transitions.insert(std::make_pair(v.interval.bounds[0], inserting));
      // If the insert fails, we need to delete instead, which we will do below.
    }
    
    while (true) {
      assert (i->second != already_filled);
      if ((i != m.transitions.end()) && (i->first <= v.interval.bounds[1])) {
        already_filled = i->second;
        if (i->second) {
          queue_erase(v.key, i->first, boost::next(i)->first);
        }
        if (!((i->first == v.interval.bounds[0]) && (i->second == inserting))) {
          m.transitions.erase(i++);
        }
        else {
          ++i;
        }
      }
      else {
        if (already_filled != inserting) {
          const auto p = m.transitions.insert(std::make_pair(v.interval.bounds[1], !inserting));
          assert (p.second);
          if (!inserting) {
            assert (i != m.transitions.end());
            assert (i->first != v.interval.bounds[1]);
            queue_insert(v.key, v.interval.bounds[1], i->first);
          }
        }
        else if (inserting) {
          new_internal_interval.bounds[1] = (i != m.transitions.end()) ? posinf : i->first;
          queue_insert(v.key, new_internal_interval);
        }
        break;
      }
    }
    if (m.transitions.empty()) {
      metadata_by_key.erase(v.key);
    }
  }
  uint64_t get_tiebreaker(KeyType k)const {
    auto m_iter = metadata_by_key.find(k);
    caller_correct_if(m_iter != metadata_by_key.end(), "getting a tiebreaker for a nonexistent key");
    return m_iter->second.tiebreaker;
  }
public:
  // TODO are these good names?
  void erase (value_type const& v) { queue_imposition(v, false); }
  void insert(value_type const& v) { queue_imposition(v, true ); }
  
  bidirectional_interval_map():root(nullptr),metadata_by_key(),next_tiebreaker(2){}
  
private:
  struct key_metadata {
    uint64_t tiebreaker;
    std::map<DomainType, bool> transitions;
  };
  
  std::unordered_set<value_type> to_be_erased;
  std::unordered_set<value_type> to_be_inserted;
  std::unique_ptr<node> root;
  std::unordered_map<KeyType, key_metadata> metadata_by_key;
  uint64_t next_tiebreaker;
  
};

namespace std {
  template<typename KeyType, typename DomainType>
  struct hash<typename bidirectional_interval_map<KeyType, DomainType>::value_type> {
    public:
    size_t operator()(typename bidirectional_interval_map<KeyType, DomainType>::value_type const& i)const {
      size_t seed = std::hash<KeyType>()(i.key);
      boost::hash_combine(seed, std::hash<DomainType>()(i.interval.bounds[0]));
      boost::hash_combine(seed, std::hash<DomainType>()(i.interval.bounds[1]));
    }
  };
}
