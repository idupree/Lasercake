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

#ifndef LASERCAKE_BIDIRECTIONAL_INTERVAL_MAP_HPP__
#define LASERCAKE_BIDIRECTIONAL_INTERVAL_MAP_HPP__

#include <array>
#include <vector>

namespace impl {
  struct unsorted {};
  template<typename QueueEntry, typename ComesAfter>
  struct iterator_queue {
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, ComesAfter> data;
    iterator_queue(ComesAfter const& comes_after):data(comes_after){}
    QueueEntry const& front()const { return data.top(); }
    void pop() { return data.pop(); }
    void push(QueueEntry const& e) { return data.push(e); }
    bool empty()const { return data.empty(); }
  };
  template<typename QueueEntry>
  struct iterator_queue<QueueEntry, unsorted> {
    std::queue<QueueEntry> data;
    QueueEntry const& front()const { return data.front(); }
    void pop() { return data.pop(); }
    void push(QueueEntry const& e) { return data.push(e); }
    bool empty()const { return data.empty(); }
  };
  
  template<typename DomainType>
  struct interval_type {
    interval_type(){}
    interval_type(DomainType a, DomainType b):bounds({{a,b}}){}
    std::array<DomainType, 2> bounds;
    bool operator==(interval_type const& o)const { return (bounds[0] == o.bounds[0]) && (bounds[1] == o.bounds[1]); }
    bool operator!=(interval_type const& o)const { return (bounds[0] != o.bounds[0]) || (bounds[1] != o.bounds[1]); }
  };
  template<typename KeyType, typename DomainType>
  struct keyed_interval {
    keyed_interval(KeyType key, DomainType a, DomainType b):key(key),interval(a,b){}
    KeyType key;
    interval_type<DomainType> interval;
    bool operator==(keyed_interval const& o)const { return (key == o.key) && (interval == o.interval); }
    bool operator!=(keyed_interval const& o)const { return (key != o.key) || (interval != o.interval); }
  };
}


// A map from KeyType to finite unions of regularized intervals
// (cf. Robust Regularized Set Operations on Polyhedra, Beat Bruderlin),
// (essentially, "the endpoints don't matter")
// with fast lookup by overlapping point, interval, etc.
template<typename KeyType, typename DomainType>
class bidirectional_interval_map {
public:
  typedef impl::interval_type<DomainType> interval_type;
  typedef impl::keyed_interval<KeyType, DomainType> value_type;
  struct tiebroken_domain_type {
    tiebroken_domain_type(){}
    tiebroken_domain_type(DomainType value, uint64_t tiebreaker):value(value),tiebreaker(tiebreaker){}
    DomainType value;
    uint64_t tiebreaker;
    bool operator==(tiebroken_domain_type const& o)const { return (value == o.value) && (tiebreaker == o.tiebreaker); }
    bool operator!=(tiebroken_domain_type const& o)const { return (value != o.value) || (tiebreaker != o.tiebreaker); }
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
      for (value_type const& v2 : values_here) { assert(v2 != v); if (v2.interval == v.interval) { std::cerr << "gosh"; } }
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
  
  struct iterator_queue_entry {
    iterator_queue_entry():n(nullptr),value(nullptr){}
    node const* n;
    value_type const* value;
    std::array<interval_type, 2> possible_bounds;
  };
public:
  template<typename ComesAfter = impl::unsorted>
  class iterator {
  private:
    struct contents {
      contents(impl::unsorted = impl::unsorted()):queue(),
        bound_minima_exist({{false,false}}),bound_maxima_exist({{false,false}}){}
      template<typename Hack = std::enable_if_t<!std::is_same<ComesAfter, impl::unsorted>::value>>
      contents(ComesAfter const& comes_after):comes_after(comes_after),queue(comes_after),
        bound_minima_exist({{false,false}}),bound_maxima_exist({{false,false}}){}
      
      ComesAfter comes_after;
      impl::iterator_queue<iterator_queue_entry, ComesAfter> queue;
      std::array<bool, 2> bound_minima_exist;
      std::array<bool, 2> bound_maxima_exist;
      std::array<DomainType, 2> bound_minima;
      std::array<DomainType, 2> bound_maxima;
    };
    
    std::shared_ptr<contents> c_;
  public:
    iterator():c_(nullptr){}
    iterator(ComesAfter comes_after):c_(new contents(comes_after)){}
    
    value_type const& operator*() const {
      caller_error_if(queue().empty(), "can't dereference an empty iterator");
      return *queue().front().value;
    }
    value_type const* operator->() const {
      caller_error_if(queue().empty(), "can't dereference an empty iterator");
      return queue().front().value;
    }
    value_type const* operator++(int) {
      value_type* result = queue().front().value;
      ++*this;
      return result;
    }
    iterator& operator++() {
      queue().pop();
      advance_to_a_returnable();
      return *this;
    }

    bool operator==(iterator const& other) const { return c_ == other.c_; }
    bool operator!=(iterator const& other) const { return !(*this == other); }
    
  private:
    impl::iterator_queue<iterator_queue_entry, ComesAfter> const& queue()const { return c_->queue; }
    impl::iterator_queue<iterator_queue_entry, ComesAfter>      & queue()      { return c_->queue; }
    
    void get_started(node const* root, interval_type const& observed_range) {
      iterator_queue_entry e;
      e.n = root;
      if (!e.n) {
        c_ = nullptr;
        return;
      }
      e.possible_bounds[0] = observed_range;
      e.possible_bounds[1] = observed_range;
      queue().push(e);
      advance_to_a_returnable();
    }
    void queue_node(iterator_queue_entry const& parent, bool use_second_child) {
      iterator_queue_entry e = parent;
      e.n = parent.n->children[use_second_child].get();
      if (!e.n) return;
      if (c_->bound_minima_exist[parent.n->splits_based_on_which_bound] && !use_second_child &&
        (parent.n->child_separator.value < c_->bound_minima[parent.n->splits_based_on_which_bound])) return;
      if (c_->bound_maxima_exist[parent.n->splits_based_on_which_bound] &&  use_second_child &&
        (parent.n->child_separator.value > c_->bound_maxima[parent.n->splits_based_on_which_bound])) return;
      e.possible_bounds[parent.n->splits_based_on_which_bound].bounds[!use_second_child] = parent.n->child_separator.value;
      queue().push(e);
    }
    void queue_value(node const* n, value_type const* value) {
      iterator_queue_entry e;
      e.n = n;
      e.value = value;
      std::array<DomainType, 2> const& bounds = value->interval.bounds;
      e.possible_bounds[0].bounds[0] = bounds[0];
      e.possible_bounds[0].bounds[1] = bounds[0];
      e.possible_bounds[1].bounds[0] = bounds[1];
      e.possible_bounds[1].bounds[1] = bounds[1];
      if (c_->bound_minima_exist[0] && (bounds[0] < c_->bound_minima[0])) return;
      if (c_->bound_maxima_exist[0] && (bounds[0] > c_->bound_maxima[0])) return;
      if (c_->bound_minima_exist[1] && (bounds[1] < c_->bound_minima[1])) return;
      if (c_->bound_maxima_exist[1] && (bounds[1] > c_->bound_maxima[1])) return;
      queue().push(e);
    }
    void advance_to_a_returnable() {
      if (queue().empty()) {
        c_ = nullptr;
        return;
      }
      while (!queue().front().value) {
        iterator_queue_entry e = queue().front();
        queue().pop();
        queue_node(e, 0);
        queue_node(e, 1);
        for (value_type const& v : e.n->values_here) {
          queue_value(e.n, &v);
        }
        if (queue().empty()) {
          c_ = nullptr;
          return;
        }
      }
    }
    
    friend bidirectional_interval_map;
  };
  template<typename ComesAfter>
  class iterator_range {
  public:
    iterator_range(iterator<ComesAfter> begin, iterator<ComesAfter> end):begin_(begin),end_(end){}
    iterator<ComesAfter> begin()const { return begin_; }
    iterator<ComesAfter> end()const { return end_; }
  private:
    iterator<ComesAfter> begin_;
    iterator<ComesAfter> end_;
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
    DomainType const& get_value( iterator_queue_entry const& v)const { return v.possible_bounds[bound].bounds[descending]; }
    DomainType const& get_value(           value_type const& v)const { return v.interval.bounds[bound]; }
    DomainType const& get_value(           value_type const* v)const { return get_value(*v); }
    DomainType const& get_value(tiebroken_domain_type const& v)const { return v.value; }
    uint64_t get_tiebreaker(           value_type const& v)const { return map->get_tiebreaker(v.key); }
    uint64_t get_tiebreaker(           value_type const* v)const { return get_tiebreaker(*v); }
    uint64_t get_tiebreaker(tiebroken_domain_type const& v)const { return v.tiebreaker; }
    bool operator()(iterator_queue_entry const& a, iterator_queue_entry const& b)const {
      if (get_value(a) != get_value(b)) {
        return (get_value(a) < get_value(b)) == descending;
      }
      if (a.value && b.value) {
        const uint64_t a_t = get_tiebreaker(a.value);
        const uint64_t b_t = get_tiebreaker(b.value);
        if (a_t != b_t) {
          return (a_t < b_t) == descending;
        }
      }
      else if (a.value || b.value) {
        // All nodes are processed before all values.
        // So, A has a value -> A comes after
        return bool(a.value);
      }
      return false;
    }
    template<typename T1, typename T2>
    std::enable_if_t<!(std::is_same<T1, iterator_queue_entry>::value && std::is_same<T2, iterator_queue_entry>::value),
    bool> operator()(T1 const& a, T2 const& b)const {
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
  };
    
  template<typename ComesAfter>
  iterator<ComesAfter>& finish_iterator(iterator<ComesAfter>& i)const {
    commit_changes();
    i.get_started(root.get(), observed_range);
    return i;
  }
public:
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_containing(DomainType const& k, ComesAfter c = ComesAfter())const {
    iterator<ComesAfter> result(c);
    result.c_->bound_maxima_exist[0] = true;
    result.c_->bound_maxima[0] = k;
    result.c_->bound_minima_exist[1] = true;
    result.c_->bound_minima[1] = k;
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_overlapping(interval_type const& k, ComesAfter c = ComesAfter())const {
    iterator<ComesAfter> result(c);
    result.c_->bound_maxima_exist[0] = true;
    result.c_->bound_maxima[0] = k.bounds[1];
    result.c_->bound_minima_exist[1] = true;
    result.c_->bound_minima[1] = k.bounds[0];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_subsuming(interval_type const& k, ComesAfter c = ComesAfter())const {
    iterator<ComesAfter> result(c);
    result.c_->bound_maxima_exist[0] = true;
    result.c_->bound_maxima[0] = k.bounds[0];
    result.c_->bound_minima_exist[1] = true;
    result.c_->bound_minima[1] = k.bounds[1];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_subsumed_by(interval_type const& k, ComesAfter c = ComesAfter())const {
    iterator<ComesAfter> result(c);
    result.c_->bound_minima_exist[0] = true;
    result.c_->bound_minima[0] = k.bounds[0];
    result.c_->bound_maxima_exist[1] = true;
    result.c_->bound_maxima[1] = k.bounds[1];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> begin_equal_to(interval_type const& k, ComesAfter c = ComesAfter())const {
    iterator<ComesAfter> result(c);
    result.c_->bound_minima_exist[0] = true;
    result.c_->bound_minima[0] = k.bounds[0];
    result.c_->bound_maxima_exist[0] = true;
    result.c_->bound_maxima[0] = k.bounds[0];
    result.c_->bound_minima_exist[1] = true;
    result.c_->bound_minima[1] = k.bounds[1];
    result.c_->bound_maxima_exist[1] = true;
    result.c_->bound_maxima[1] = k.bounds[1];
    return finish_iterator(result);
  }
  template<typename ComesAfter = impl::unsorted>
  iterator<ComesAfter> end()const {
    return iterator<ComesAfter>();
  }
  
  template<typename ComesAfter = impl::unsorted>
  iterator_range<ComesAfter> range_containing (   DomainType const& k, ComesAfter c = ComesAfter())const {
    return iterator_range<ComesAfter>(begin_containing (k, c), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  iterator_range<ComesAfter> range_overlapping(interval_type const& k, ComesAfter c = ComesAfter())const {
    return iterator_range<ComesAfter>(begin_overlapping(k, c), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  iterator_range<ComesAfter> range_subsuming  (interval_type const& k, ComesAfter c = ComesAfter())const {
    return iterator_range<ComesAfter>(begin_subsuming  (k, c), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  iterator_range<ComesAfter> range_subsumed_by(interval_type const& k, ComesAfter c = ComesAfter())const {
    return iterator_range<ComesAfter>(begin_subsumed_by(k, c), end<ComesAfter>());
  }
  template<typename ComesAfter = impl::unsorted>
  iterator_range<ComesAfter> range_equal_to   (interval_type const& k, ComesAfter c = ComesAfter())const {
    return iterator_range<ComesAfter>(begin_equal_to   (k, c), end<ComesAfter>());
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
  static const size_t max_leaf_size = 6;
  
  static bool is_in_unsorted_vector(value_type const& what, std::vector<value_type const*> const& vec) {
    for (auto const& v : vec) { if (*v == what) { return true; } }
    return false;
  }
  static bool is_in_unsorted_vector(value_type const& what, std::vector<value_type> const& vec) {
    for (auto const& v : vec) { if (v == what) { return true; } }
    return false;
  }
  static bool is_in_sorted_vector(value_type const& what, std::vector<value_type const*> const& vec, comes_after_by_bound_t const& comes_after) {
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
  bool value_in_tree(value_type const& v, node const* n)const {
    impl::unsorted f;
    iterator<impl::unsorted> result(f);
    result.c_->bound_minima_exist[0] = true;
    result.c_->bound_minima[0] = v.interval.bounds[0];
    result.c_->bound_maxima_exist[0] = true;
    result.c_->bound_maxima[0] = v.interval.bounds[0];
    result.c_->bound_minima_exist[1] = true;
    result.c_->bound_minima[1] = v.interval.bounds[1];
    result.c_->bound_maxima_exist[1] = true;
    result.c_->bound_maxima[1] = v.interval.bounds[1];
    result.get_started(n, observed_range);
      std::cerr << "A";
    for (; result != end(); ++result) {
      assert (result->interval == v.interval);
      std::cerr << "l";
      if (result->key == v.key) {
      std::cerr << "k";
        assert (*result == v);
        return true;
      }
    }
    return false;
  }
  void validate_tree()const {
    validate_tree(root);
  }
  void validate_tree(std::unique_ptr<node> const& n, std::vector<value_type const*> const* except = nullptr, std::vector<value_type const*> const* except2 = nullptr)const {
    if (n && n->children[0]) {
      for (uint32_t which_child = 0; which_child < 2; ++which_child) {
        validate_tree_recurse(n, n->children[which_child], which_child, except, except2);
        validate_tree(n->children[which_child], except, except2);
      }
    }
  }
  void validate_tree_recurse(std::unique_ptr<node> const& n, std::unique_ptr<node> const& n2, bool which_original_child, std::vector<value_type const*> const* except, std::vector<value_type const*> const* except2)const {
    if (n2->children[0]) {
      validate_tree_recurse(n, n2->children[0], which_original_child, except, except2);
      validate_tree_recurse(n, n2->children[1], which_original_child, except, except2);
    }
    else {
      for (auto const& v : n2->values_here) {
        comes_after_by_bound_t is_higher = comes_after_by_bound(n->splits_based_on_which_bound, false);
        if (!is_higher(v, n->child_separator) == which_original_child) {
          bool ok = false;
          if (except ) { for (auto const& v2 : *except ) { if (*v2 == v) { ok = true; }}}
          if (except2) { for (auto const& v2 : *except2) { if (*v2 == v) { ok = true; }}}
          assert (ok);
        }
      }
    }
  }
  static void collapse_node(std::unique_ptr<node>& n, std::unique_ptr<node>& descendant) {
    for (value_type const& v : descendant->values_here) {
      n->insert_value(v);
    }
    if (descendant->children[0]) {
      collapse_node(n, descendant->children[0]);
      collapse_node(n, descendant->children[1]);
    }
    descendant = nullptr;
  }
  
  // Insert M elements and erase P elements.
  // O((M+P)*log(N)) amortized.
  void insert_and_or_erase_impl(bool is_root, std::unique_ptr<node>& n,
      std::vector<value_type const*> const& values_to_insert, std::vector<value_type const*> const& values_to_erase)const {
    if (!n) {
      assert(is_root);
      n.reset(new node(false));
    }
    n->num_descendant_intervals += values_to_insert.size();
    n->num_descendant_intervals -= values_to_erase .size();
    if (n->num_descendant_intervals == 0) {
      assert(is_root);
      n = nullptr;
      return;
    }
    
    if (n->children[0]) { // i.e. "if non leaf"
      validate_tree(n, &values_to_erase);
      assert (n->children[1]);
      if (n->num_descendant_intervals <= max_leaf_size) {
      std::cerr << "E";
        for (value_type const* v : values_to_erase) {
          assert (value_in_tree(*v, n.get()));
        }
        collapse_node(n, n->children[0]);
        collapse_node(n, n->children[1]);
      std::cerr << "Q";
        for (value_type const* v : values_to_erase) {
          assert (value_in_tree(*v, n.get()));
        }
      std::cerr << "W";
      }
      else {
        std::array<std::vector<value_type const*>, 2> values_to_insert_by_child;
        std::array<std::vector<value_type const*>, 2> values_to_erase_by_child;
        for (auto& v : values_to_insert_by_child) { v.reserve(values_to_insert.size()); }
        for (auto& v : values_to_erase_by_child ) { v.reserve(values_to_erase .size()); }
        comes_after_by_bound_t is_higher = comes_after_by_bound(n->splits_based_on_which_bound, false);
        for (value_type const* v : values_to_insert) {
          assert (tiebroken_domain_type(
                  v->interval.bounds[n->splits_based_on_which_bound],
                  get_tiebreaker(v->key)) != n->child_separator);
          values_to_insert_by_child[is_higher(*v, n->child_separator)].push_back(v);
        }
        for (value_type const* v : values_to_erase) {
          assert (value_in_tree(*v, n.get()));
          assert (tiebroken_domain_type(
                  v->interval.bounds[n->splits_based_on_which_bound],
                  get_tiebreaker(v->key)) != n->child_separator);
          uint32_t which_child = is_higher(*v, n->child_separator);
          assert (!is_in_unsorted_vector(*v, values_to_erase_by_child[which_child]));
          values_to_erase_by_child[which_child].push_back(v);
          assert (value_in_tree(*v, n->children[which_child].get()));
        }
  #define NEW_NUM_DESCENDANTS(which_child) (n->children[which_child]->num_descendant_intervals \
    + values_to_insert_by_child[which_child].size() \
    - values_to_erase_by_child[which_child].size())
    
        std::vector<value_type> stolen;
        for (uint32_t which_child = 0; which_child < 2; ++which_child) {
          assert(n->children[which_child]->num_descendant_intervals + values_to_insert_by_child[which_child].size() >= values_to_erase_by_child[which_child].size());
          if (NEW_NUM_DESCENDANTS(which_child) > 2*NEW_NUM_DESCENDANTS(!which_child)) {
            // reversed sense because std::sort wants a comes-before.
            // We want the *backs* of these vectors to be closest to the separator.
            std::sort(values_to_insert_by_child[0].begin(), values_to_insert_by_child[0].end(),
                      comes_after_by_bound(n->splits_based_on_which_bound,  true));
            std::sort(values_to_insert_by_child[1].begin(), values_to_insert_by_child[1].end(),
                      comes_after_by_bound(n->splits_based_on_which_bound, false));
            std::sort(values_to_erase_by_child[which_child].begin(), values_to_erase_by_child[which_child].end(),
                      comes_after_by_bound(n->splits_based_on_which_bound, !which_child));
            // Regular sense, but we want the *front* of the thief iterator to be closest to the separator.
            comes_after_by_bound_t is_further_from_separator = comes_after_by_bound(n->splits_based_on_which_bound, !which_child);
            iterator<comes_after_by_bound_t> thief(is_further_from_separator);
            comes_after_by_bound_t is_closer_to_separator = comes_after_by_bound(n->splits_based_on_which_bound, which_child);
            thief.get_started(n->children[which_child].get(), observed_range);
            while (NEW_NUM_DESCENDANTS(which_child) > NEW_NUM_DESCENDANTS(!which_child) + 1 + 2*stolen.size()) {
              const int64_t not_quite_as_close = (which_child ? 1 : -1);
              tiebroken_domain_type closest_to_insert_bound;
              tiebroken_domain_type closest_to_steal_bound;
              bool use_thief  = (thief != end<comes_after_by_bound_t>());
              bool use_insert = (!values_to_insert_by_child[which_child].empty());
              value_type const* closest_to_insert;
              assert (use_thief || use_insert);
              if (use_thief) {
                closest_to_steal_bound = tiebroken_domain_type(
                              thief->interval.bounds[n->splits_based_on_which_bound],
                  get_tiebreaker(            thief->key));
              }
              if (use_insert) {
                closest_to_insert = values_to_insert_by_child[which_child].back();
                closest_to_insert_bound = tiebroken_domain_type(
                  closest_to_insert->interval.bounds[n->splits_based_on_which_bound],
                  get_tiebreaker(closest_to_insert->key));
              }
              if (use_thief && use_insert) {
                if (is_further_from_separator(closest_to_insert_bound, closest_to_steal_bound)) {
                  use_insert = false;
                }
                else {
                  use_thief = false;
                }
              }
              if (use_thief) {
                assert (closest_to_steal_bound != n->child_separator);
                assert (is_further_from_separator(closest_to_steal_bound, n->child_separator));
                if (!is_in_sorted_vector(*thief, values_to_erase_by_child[which_child], is_closer_to_separator)) {
                  assert (!is_in_unsorted_vector(*thief, values_to_erase_by_child[which_child]));
                  assert (!is_in_unsorted_vector(*thief, stolen));
                  stolen.push_back(*thief);
                  assert (value_in_tree(*thief, n->children[which_child].get()));
                  n->child_separator = closest_to_steal_bound;
                  n->child_separator.tiebreaker += not_quite_as_close;
                  assert (is_further_from_separator(n->child_separator, closest_to_steal_bound));
                  //validate_tree(n, &values_to_erase_by_child[which_child], &stolen);
                }
                ++thief;
              }
              else {
                n->child_separator = closest_to_insert_bound;
                n->child_separator.tiebreaker += not_quite_as_close;
                assert (is_further_from_separator(n->child_separator, closest_to_insert_bound));
                values_to_insert_by_child[!which_child].push_back(closest_to_insert);
                values_to_insert_by_child[which_child].pop_back();
                //validate_tree(n, &values_to_erase_by_child[which_child], &stolen);
              }
            }
            for (value_type const& v : stolen) {
              values_to_erase_by_child [ which_child].push_back(&v);
              values_to_insert_by_child[!which_child].push_back(&v);
            }
          }
        }
        validate_tree(n, &values_to_erase_by_child[0], &values_to_erase_by_child[1]);
        for (uint32_t which_child = 0; which_child < 2; ++which_child) {
          for (value_type const* v : values_to_insert_by_child[which_child]) {
            assert (is_higher(v, n->child_separator) == which_child);
          }
//           for (value_type const* v : values_to_erase_by_child[which_child]) {
//             assert (is_higher(v, n->child_separator) == which_child);
//           }
          insert_and_or_erase_impl(false, n->children[which_child], values_to_insert_by_child[which_child], values_to_erase_by_child[which_child]);
        }
        validate_tree(n, &values_to_erase_by_child[0], &values_to_erase_by_child[1]);
        validate_tree(n, &values_to_erase);
        validate_tree(n);
      }
    }
    
    if (!n->children[0]) {
      for (value_type const* v : values_to_erase) {
        n->erase_value(*v);
      }
      for (value_type const* v : values_to_insert) {
        n->insert_value(*v);
        assert (value_in_tree(*v, root.get()));
      }
      assert (n->num_descendant_intervals == n->values_here.size());
      if (n->num_descendant_intervals > max_leaf_size) {
        comes_after_by_bound_t is_higher = comes_after_by_bound(n->splits_based_on_which_bound, false);
        n->children[0].reset(new node(!n->splits_based_on_which_bound));
        n->children[1].reset(new node(!n->splits_based_on_which_bound));
        
        // Ascending order. Reversed sense because std::sort wants a comes-before.
        std::sort(n->values_here.begin(), n->values_here.end(), comes_after_by_bound(n->splits_based_on_which_bound, true));
        size_t median_idx = n->values_here.size() >> 1;
        n->child_separator = tiebroken_domain_type(
          n->values_here[median_idx].interval.bounds[n->splits_based_on_which_bound],
          get_tiebreaker(n->values_here[median_idx].key) - 1);
        std::array<std::vector<value_type const*>, 2> values_to_insert_by_child;
        std::vector<value_type const*> values_to_erase_by_child;
        
        for (size_t i = 0; i < n->values_here.size(); ++i) {
          values_to_insert_by_child[i >= median_idx].push_back(&n->values_here[i]);
          assert (is_higher(n->values_here[i], n->child_separator) == (i >= median_idx));
        }
        for (uint32_t which_child = 0; which_child < 2; ++which_child) {
          insert_and_or_erase_impl(false, n->children[which_child], values_to_insert_by_child[which_child], values_to_erase_by_child);
        }
        n->values_here.clear();
      }
      validate_tree(n);
    }
  }
  
  void commit_changes()const {
    std::vector<value_type const*> values_to_insert;
    std::vector<value_type const*> values_to_erase;
    
    for (value_type const& v : to_be_inserted) {
      values_to_insert.push_back(&v);
    }
    for (value_type const& v : to_be_erased) {
      values_to_erase.push_back(&v);
    }
    
    validate_tree();
    insert_and_or_erase_impl(true, root, values_to_insert, values_to_erase);
    validate_tree();
    
    to_be_inserted.clear();
    to_be_erased.clear();
    for (KeyType k : keys_maybe_gone) {
      auto m = metadata_by_key.find(k);
      if ((m != metadata_by_key.end()) && m->second.transitions.empty()) {
        metadata_by_key.erase(m);
      }
    }
    keys_maybe_gone.clear();
  }
  void queue_erase(value_type const& v) {
    if (!to_be_inserted.erase(v)) { const auto p = to_be_erased.insert(v); assert(p.second); }
  }
  void queue_insert(value_type const& v) {
    if (!to_be_erased.erase(v)) { const auto p = to_be_inserted.insert(v); assert(p.second); }
  }
  void queue_imposition(value_type const& v, bool inserting) {
    caller_correct_if(v.interval.bounds[1] > v.interval.bounds[0], "intervals must have positive length");
    if (metadata_by_key.empty() || (v.interval.bounds[0] < observed_range.bounds[0])) {
      observed_range.bounds[0] = v.interval.bounds[0];
    }
    if (metadata_by_key.empty() || (v.interval.bounds[1] > observed_range.bounds[1])) {
      observed_range.bounds[1] = v.interval.bounds[1];
    }
    value_type new_internal_value = v;
    
    key_metadata& m = metadata_by_key[v.key];
    if (!m.tiebreaker) {
      m.tiebreaker = next_tiebreaker;
      next_tiebreaker += 2;
    }
    
    bool b=false;
    bool e=false;
      std::cerr << inserting << "\n";
    for (auto const& i : m.transitions) {
      if (i.first >= v.interval.bounds[0] && !b) { std::cerr << "B"; b=true; }
      if (i.first == v.interval.bounds[0]) { std::cerr << "="; }
      if (i.first >= v.interval.bounds[1] && !e) { std::cerr << "E"; e=true; }
      if (i.first == v.interval.bounds[1]) { std::cerr << "="; }
      std::cerr << i.second;
    }
      std::cerr << "\n";
    
    bool already_filled;
    auto i = m.transitions.lower_bound(v.interval.bounds[0]);
    if (i == m.transitions.begin()) {
      already_filled = false;
    }
    else {
      const auto last_iter = boost::prior(i);
      already_filled = last_iter->second;
      if (last_iter->second) {
        queue_erase(value_type(v.key, last_iter->first, boost::next(last_iter)->first));
        if (inserting) {
          new_internal_value.interval.bounds[0] = last_iter->first;
        }
        else {
          queue_insert(value_type(v.key, last_iter->first, v.interval.bounds[0]));
        }
      }
    }
    if (already_filled != inserting) {
      m.transitions.insert(std::make_pair(v.interval.bounds[0], inserting));
      assert (key_active_before(v.key, v.interval.bounds[0]) == !inserting);
      // assert (key_active_after (v.key, v.interval.bounds[0]) ==  inserting);
      // If the insert fails, we need to delete instead, which we will do below.
    }
    else {
      assert (key_active_before(v.key, v.interval.bounds[0]) ==  inserting);
      //assert (key_active_after (v.key, v.interval.bounds[0]) ==  inserting);
      // Relying on the delete below
    }
    
    while (true) {
      assert (((i != m.transitions.end()) ? i->second : true) != already_filled);
      if ((i != m.transitions.end()) && (i->first <= v.interval.bounds[1])) {
        already_filled = i->second;
        if (i->second) {
          queue_erase(value_type(v.key, i->first, boost::next(i)->first));
        }
        if ((i->first == v.interval.bounds[0]) && (i->second == inserting)) {
          assert (false);
          assert (key_active_before(v.key, v.interval.bounds[0]) == !inserting);
          assert (key_active_after (v.key, v.interval.bounds[0]) ==  inserting);
          ++i;
        }
        else {
          m.transitions.erase(i++);
        }
      }
      else {
        if (already_filled != inserting) {
          const auto p = m.transitions.insert(std::make_pair(v.interval.bounds[1], !inserting));
          assert (p.second);
          if (!inserting) {
            assert (i != m.transitions.end());
            assert (i->first != v.interval.bounds[1]);
            queue_insert(value_type(v.key, v.interval.bounds[1], i->first));
          }
          assert (key_active_before(v.key, v.interval.bounds[1]) ==  inserting);
          assert (key_active_after (v.key, v.interval.bounds[1]) == !inserting);
        }
        else {
          assert (key_active_before(v.key, v.interval.bounds[1]) ==  inserting);
          assert (key_active_after (v.key, v.interval.bounds[1]) ==  inserting);
          if (inserting) {
            assert (i != m.transitions.end());
            new_internal_value.interval.bounds[1] = i->first;
          }
        }
        if (inserting) {
          queue_insert(new_internal_value);
        }
        break;
      }
    }
    if (m.transitions.empty()) {
      keys_maybe_gone.push_back(v.key);
    }
    assert (key_active_after (v.key, v.interval.bounds[0]) == inserting);
    assert (key_active_before(v.key, v.interval.bounds[1]) == inserting);
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
  template<typename T1, typename T2, typename...T>
  void erase (T1 data1, T2 data2, T... data){ erase (value_type(std::forward<T1>(data1), std::forward<T2>(data2), std::forward<T>(data)...)); }
  template<typename T1, typename T2, typename...T>
  void insert(T1 data1, T2 data2, T... data){ insert(value_type(std::forward<T1>(data1), std::forward<T2>(data2), std::forward<T>(data)...)); }
  
  bidirectional_interval_map():root(nullptr),metadata_by_key(),next_tiebreaker(2){}
  
  bool key_active_after(KeyType const& k, DomainType where) {
    auto m_iter = metadata_by_key.find(k);
    if (m_iter == metadata_by_key.end()) { return false; }
    auto t_iter = m_iter->second.transitions.upper_bound(where);
    if (t_iter == m_iter->second.transitions.end()) { return false; }
    return !t_iter->second;
  }
  bool key_active_before(KeyType const& k, DomainType where) {
    auto m_iter = metadata_by_key.find(k);
    if (m_iter == metadata_by_key.end()) { return false; }
    auto t_iter = m_iter->second.transitions.lower_bound(where);
    if (t_iter == m_iter->second.transitions.end()) { return false; }
    return !t_iter->second;
  }
private:
  struct key_metadata {
    uint64_t tiebreaker;
    std::map<DomainType, bool> transitions;
  };
  
  mutable std::unordered_set<value_type> to_be_erased;
  mutable std::unordered_set<value_type> to_be_inserted;
  mutable std::vector<KeyType> keys_maybe_gone;
  mutable std::unique_ptr<node> root;
  mutable std::unordered_map<KeyType, key_metadata> metadata_by_key;
  uint64_t next_tiebreaker;
  
  // hacks to init iterator possible_bounds:
  interval_type observed_range;
};

namespace std {
  template<typename KeyType, typename DomainType>
  struct hash<typename impl::keyed_interval<KeyType, DomainType>> {
    public:
    size_t operator()(typename impl::keyed_interval<KeyType, DomainType> const& i)const {
      size_t seed = std::hash<KeyType>()(i.key);
      seed ^= std::hash<DomainType>()(i.interval.bounds[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= std::hash<DomainType>()(i.interval.bounds[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
}

#endif
