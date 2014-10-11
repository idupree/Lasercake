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

template<typename DomainType, typename MappedType>
class interval_tree {
public:
  typedef std::pair<interval, MappedType> value_type;
  struct interval {
    std::array<DomainType, 2> bounds;
  };
private:
  struct node {
    bool splits_based_on_which_bound;
    size_t num_descendant_intervals;
    
    // only for non leaves
    DomainType child_separator;
    std::array<std::unique_ptr<node>, 2> children; // the first has the lower values, the second the higher values
    
    // only for leaves
    std::unordered_map<interval, MappedType> intervals_here;
  };
  
  std::unique_ptr<node> root;
  
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
        for (value_type& v : queue.front().n->intervals_here) {
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
    return result;
  }
  iterator find_overlapping(interval const& k) {
    iterator result;
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[1];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[0];
    result.queue_root(root.get());
    return result;
  }
  iterator find_subsuming(interval const& k) {
    iterator result;
    result.bound_maxima_exist[0] = true;
    result.bound_maxima[0] = k.bounds[0];
    result.bound_minima_exist[1] = true;
    result.bound_minima[1] = k.bounds[1];
    result.queue_root(root.get());
    return result;
  }
  iterator find_subsumed(interval const& k) {
    iterator result;
    result.bound_minima_exist[0] = true;
    result.bound_minima[0] = k.bounds[0];
    result.bound_maxima_exist[1] = true;
    result.bound_maxima[1] = k.bounds[1];
    result.queue_root(root.get());
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
    return result;
  }
  iterator end()const {
    return iterator();
  }
  
  void erase(iterator i) {
    i.queue.front().n
  }
  
  
  
};


