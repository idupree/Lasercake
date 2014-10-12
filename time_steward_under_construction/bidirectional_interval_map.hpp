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

// A map from KeyType to (finite union of regularized intervals, cf. Robust Regularized Set Operations on Polyhedra, Beat Bruderlin),
// with fast lookup by overlapping point, interval, etc.
template<typename KeyType, typenameDomainType>
class bidirectional_interval_map {
public:
  typedef std::pair<interval, KeyType> value_type;
  struct tiebroken_domain_type {
    DomainType
    uint64_t tiebreaker;
  };
  struct interval {
    std::array<DomainType, 2> bounds;
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
    std::unordered_set<value_type> values_here;
  };
  
public:
  template<typename Ordering>
  class iterator {
  public:
    iterator():bound_minima_exist({{false,false}}),bound_maxima_exist({{false,false}}),min_length_exists(false){}

    value_type& operator*() const {
      caller_error_if(queue_.empty(), "can't dereference an empty iterator");
      return *queue.front().value;
    }
    value_type* operator->() const {
      caller_error_if(queue_.empty(), "can't dereference an empty iterator");
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
    std::array<bool, 2> bound_minima_exist;
    std::array<bool, 2> bound_maxima_exist;
    bool min_length_exists;
    
    struct queue_entry {
      node* n;
      value_type* value;
      std::array<interval, 2> possible_bounds;
    };
    std::priority_queue<queue_value_type_, std::vector<queue_value_type_>, reverse_first_ordering<CostOrdering, queue_value_type_> > queue_type_;
    
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
      if (bound_minima_exist[parent.n->splits_based_on_which_bound] &&  use_second_child && (parent.n->child_separator < bound_minima[parent.n->splits_based_on_which_bound])) return;
      if (bound_maxima_exist[parent.n->splits_based_on_which_bound] && !use_second_child && (parent.n->child_separator > bound_maxima[parent.n->splits_based_on_which_bound])) return;
      e.possible_bounds[parent.n->splits_based_on_which_bound].bounds[!use_second_child] = parent.n->child_separator;
      if (min_length_exists && interval(e.possible_bounds[0].bounds[0], e.possible_bounds[1].bounds[1]).length_less_than(min_length)) return;
      queue.push(e);
    }
    void queue_value(node* n, value_type* value) {
      queue_entry e;
      e.n = n;
      e.value = value;
      std::array<DomainType, 2> const& bounds = value->second.bounds;
      if (bound_minima_exist[0] && (bounds[0] < bound_minima[0])) return;
      if (bound_maxima_exist[0] && (bounds[0] > bound_maxima[0])) return;
      if (bound_minima_exist[1] && (bounds[1] < bound_minima[1])) return;
      if (bound_maxima_exist[1] && (bounds[1] > bound_maxima[1])) return;
      if (min_length_exists && value->second.length_less_than(min_length)) return;
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
    
    friend interval_tree;
  };
  
  iterator find_containing(KeyType const& k) {
    iterator result;
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k;
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k;
    result.queue_root(root.get());
    result.advance_to_a_returnable();
    return result;
  }
  iterator find_overlapping(interval const& k) {
    iterator result;
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[1];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[0];
    result.queue_root(root.get());
    result.advance_to_a_returnable();
    return result;
  }
  iterator find_subsuming(interval const& k) {
    iterator result;
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[0];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[1];
    result.queue_root(root.get());
    result.advance_to_a_returnable();
    return result;
  }
  iterator find_subsumed(interval const& k) {
    iterator result;
    result.bound_minima_exist[0] = true;
    result.bound_minima[0] = k.bounds[0];
    result.bound_maxima_exist[1] = true;
    result.bound_maxima[1] = k.bounds[1];
    result.queue_root(root.get());
    result.advance_to_a_returnable();
    return result;
  }
  iterator find_exact(interval const& k) {
    iterator result;
    result.bound_minima_exist[0] = true;
    result.bound_minima[0] = k.bounds[0];
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[0];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[1];
    result.bound_maxima_exist[1] = true;
    result.bound_maxima[1] = k.bounds[1];
    result.queue_root(root.get());
    result.advance_to_a_returnable();
    return result;
  }
  iterator end()const {
    return iterator();
  }
private:
  
  
  // The tree maintains these invariants:
  // Nodes always alternate splits_based_on_which_bound.
  // No non-leaf directly contains any intervals.
  // 1) No leaf contains more than max_leaf_size intervals.
  // 2) No leaf contains less than min_leaf_size intervals (except for when the root is just a leaf).
  // 3) No non-leaf has more than twice as many descendants on one side as on the other.
  // (2) and (3) together bound the tree height at log_3(N).
  static const min_leaf_size = 2;
  static const max_leaf_size = 6;
  static_assert (max_leaf_size+1 >= min_leaf_size * 2, "must be able to split a node that exceeds max size");
  
  static void collapse_node(node& n, std::unique_ptr<node>& descendant) {
    for (value_type& v : descendant.values_here) {
      n.values_here.insert(v);
    }
    if (descendant->children[0]) {
      collapse_node(n, descendant->children[0]);
      collapse_node(n, descendant->children[1]);
    }
    descendant = nullptr;
  }
  static void insert_and_or_erase_impl(bool is_root, std::unique_ptr<node>& n, std::vector<value_type*> const& values_to_insert, std::vector<value_type*> const& values_to_erase) {
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
        std::array<std::vector<value_type*>, 2> values_to_insert_by_child;
        std::array<std::vector<value_type*>, 2> values_to_erase_by_child;
        for (auto& v : values_to_insert_by_child) { v.reserve(values_to_insert.size()); }
        for (auto& v : values_to_erase_by_child ) { v.reserve(values_to_erase .size()); }
        for (value_type* v : values_to_insert) {
          values_to_insert_by_child[v.first.bounds[n->splits_based_on_which_bound] > n->child_separator].push_back(v);
        }
        for (value_type* v : values_to_erase ) {
          values_to_erase_by_child [v.first.bounds[n->splits_based_on_which_bound] > n->child_separator].push_back(v);
          if (                      v.first.bounds[n->splits_based_on_which_bound]== n->child_separator) {
            // for == values, we need to check both children.
            values_to_erase_by_child [true].push_back(v);
          }
        }
  #define NEW_NUM_DESCENDANTS(which_child) (n->children[which_child]->num_descendant_intervals \
    + values_to_insert_by_child[which_child].size() \
    - values_to_erase_by_child[which_child].size())
        std::sort(values_to_insert_by_child[0], low to high);
        std::sort(values_to_erase_by_child [0], low to high);
        std::sort(values_to_insert_by_child[1], high to low);
        std::sort(values_to_erase_by_child [1], high to low);
    
        for (which_child : true, false) {
          if (NEW_NUM_DESCENDANTS(which_child) > 2*NEW_NUM_DESCENDANTS(!which_child)) {
            iterator<> thief;
            std::vector<value_type*> stolen;
            if (which_child) {
              thief.bound_minima_exist[0] = true;
              thief.bound_minima[0] = n->child_separator;
              thief.comes_before = low to high
            }
            else {
              thief.bound_maxima_exist[1] = true;
              thief.bound_maxima[1] = n->child_separator;
              thief.comes_before = high to low
            }
            thief.queue_root(n->children[which_child].get());
            thief.advance_to_a_returnable();
            while (NEW_NUM_DESCENDANTS(which_child) > NEW_NUM_DESCENDANTS(!which_child) + 1) {
              if (thief.comes_before(                   thief->first.bounds[n->splits_based_on_which_bound]
                values_to_insert_by_child[which_child].back()->first.bounds[n->splits_based_on_which_bound])) {
                if (!values_to_erase.find(*thief)) {
                  stolen.push_back(*thief);
                  n->child_separator = thief->first.bounds[!which_child];
                }
                ++thief;
              }
              else {
                n->child_separator = values_to_insert_by_child[which_child].back()->first.bounds[n->splits_based_on_which_bound];
                values_to_insert_by_child[!which_child].push_back(values_to_insert_by_child[which_child].back());
                values_to_insert_by_child[which_child].pop_back();
              }
            }
            for (value_type* v : stolen) {
              values_to_erase_by_child [ which_child].push_back(v);
              values_to_insert_by_child[!which_child].push_back(v);
            }
          }
        }
        for (which_child : true, false) {
          insert_and_or_erase_impl(false, n->children[which_child], values_to_insert_by_child[which_child], values_to_erase_by_child[which_child]);
        }
      }
    }
    else {
      for (value_type* v : values_to_erase) {
        const auto p = n->values_here.erase(*v);
        caller_correct_if(p, "erasing an element that doesn't exist");
      }
      for (value_type* v : values_to_insert) {
        const auto p = n->values_here.insert(*v);
        assert(p.second);
      }
      if (n->num_descendant_intervals > max_leaf_size) {
        n->children[0].reset(new node(!n->splits_based_on_which_bound));
        n->children[1].reset(new node(!n->splits_based_on_which_bound));
        
        n->child_separator = median(n->values_here)
        std::array<std::vector<value_type*>, 2> values_to_insert_by_child;
        for (value_type* v : n->values_here) {
          values_to_insert_by_child[v->bounds[n->splits_based_on_which_bound] > n->child_separator].push_back(v);
        }
        for (which_child : true, false) {
          insert_and_or_erase_impl(false, n->children[which_child], values_to_insert_by_child[which_child], values_to_erase_by_child[which_child]);
        }
      }
    }
  }
  void queue_setting(value_type* v, bool inserting)
    interval internal_insertion = v->first;
    
    key_metadata& m = *metadata_by_key.insert(std::make_pair(v->second, key_metadata())).second;
    if (!m.tiebreaker) { m.tiebreaker = next_tiebreaker++; }
    bool already_filled = false;
    
    const auto p = m.transitions.insert(std::make_pair(v->first.bounds[0], inserting));
    auto i = boost::next(p.first);
    
    if (p.first == m.transitions.begin()) {
      if (!inserting) {
        internal_insertion.bounds[0] = neginf;
      }
    }
    else {
      const auto last_iter = boost::prior(p.first);
      internal_insertion.bounds[0] = last_iter->first;
      already_filled = last_iter->second;
    }
    if (already_filled == inserting) {
      assert(internal_insertion.bounds[0] != p.first);
      m.transitions.erase(p.first);
    }
    else {
      assert(p.first->second == inserting);
    }
    
    for (; i != m.transitions.end();) {
      assert (i->second != already_filled);
      if ((i != m.transitions.end()) && (i->first < v->first.bounds[1])) {
        m.transitions.erase(i++);
      }
      else {
        if (already_filled != inserting) {
          m.transitions.insert(std::make_pair(v->first.bounds[1], !inserting));
        }
        else {
          assert (i != m.transitions.end());
          internal_insertion.bounds[0] = i->first;
        }
        break;
      }
      already_filled = i->second;
    }
  }
  
public:
  // Insert M elements and erase P elements.
  // O(M+P+log(N)) amortized.
  void insert_and_or_erase(std::vector<value_type*> const& values_to_insert, std::vector<value_type*> const& values_to_erase) {
    std::vector<value_type*> internal_values_to_insert;
    std::vector<value_type*> internal_values_to_erase;
    insert_and_or_erase_impl(true, root, internal_values_to_insert, internal_values_to_erase);
  }
  void erase(value_type const& i) {
    std::vector<value_type*> values_to_insert;
    std::vector<value_type*> values_to_erase;
    values_to_erase.push_back(&i);
    insert_and_or_erase_impl(true, root, values_to_insert, values_to_erase);
  }
  void insert(value_type const& i) {
    std::vector<value_type*> values_to_insert;
    std::vector<value_type*> values_to_erase;
    values_to_insert.push_back(&i);
    insert_and_or_erase_impl(true, root, values_to_insert, values_to_erase);
  }
  
  bidirectional_interval_map():root(nullptr),metadata_by_key(),next_tiebreaker(1){}
  
private:
  struct key_metadata {
    uint64_t tiebreaker;
    std::map<DomainType, bool> transitions;
  };
  
  std::unique_ptr<node> root;
  std::unordered_map<KeyType, key_metadata> metadata_by_key;
  uint64_t next_tiebreaker;
  
};


