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

#ifndef LASERCAKE_MANUAL_ORDERER_HPP__
#define LASERCAKE_MANUAL_ORDERER_HPP__

#include <list>

#ifndef WORST_CASE_MANUAL_ORDERER
template<typename ValueType> class manual_orderer;
namespace manual_orderer_impl {
typedef uint64_t idx_type;
typedef uint32_t num_bits_type;
const idx_type no_idx = std::numeric_limits<idx_type>::max();

template<typename ValueType>
struct entry {
  template <class... Args>
  entry(idx_type idx, Args&&... args):idx(idx),ref_count(0),contents(std::forward<Args>(args)...){}
  idx_type idx;
  size_t ref_count;
  ValueType contents;
};

// an entry_ref is a non-threadsafe reference-counted smart pointer
// that is LessThanComparable by a idx it stores next to the contents,
// and EqualityComparable by the pointer.  You must not compare entry_refs
// created by different manual_orderer instances.
// 
// manual_orderer is the only code allowed to access or change the idx;
// it ensures that idx is equal iff pointers are equal, and that it only
// changes idx in ways that are not observable using the less-than operations.
// As an exception, entries can have idx be no_idx, in which case
// less-than-comparison ops throw an exception and idx might later
// be updated to a real idx.
template<typename ValueType>
struct entry_ref {
public:
  entry_ref():data(nullptr){}
  entry_ref(entry_ref const& o):data(o.data) { inc_ref(); }
  bool operator< (entry_ref const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx <  o.data->idx;
  }
  bool operator> (entry_ref const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx >  o.data->idx;
  }
  bool operator<=(entry_ref const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx <= o.data->idx;
  }
  bool operator>=(entry_ref const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx >= o.data->idx;
  }
  bool operator==(entry_ref const& o)const { return data == o.data; }
  bool operator!=(entry_ref const& o)const { return data != o.data; }
  explicit operator bool()const { return bool(data); }
  ValueType& operator*()const { return data->contents; }
  ValueType* operator->()const { return &data->contents; }
  entry_ref& operator=(entry_ref const& o) {
    o.inc_ref(); // do this before dec_ref in case lhs and rhs are the same object
    dec_ref();
    data = o.data;
    return *this;
  }
  ~entry_ref() { dec_ref(); }
private:
  void inc_ref()const {
    if (data) { ++data->ref_count; }
  }
  void dec_ref() {
    if (data && (--data->ref_count == 0)) {
      delete data;
      data = nullptr;
    }
  }
  explicit entry_ref(entry<ValueType>* data):data(data){ inc_ref(); }
  entry<ValueType>* data;
  friend class manual_orderer<ValueType>;
  friend struct std::hash<entry_ref>;
};
} // end namespace manual_orderer_impl

//#define AUDIT_ORDERED_STUFF


// manual_orderer is a data structure that lets you place
// arbitrary objects in an O(1)-comparable order, but you have to refer
// to them by entry_ref instead of the object itself.
//
// Implementation limitation: if you put too many more than 2**32 items in,
// this data structure will break (currently, by being buggy, not even
// by throwing an exception).
//
// Issue: the contained data aren't destroyed unless you both drop all
// its entry_refs and destroy the manual_orderer.
template<typename ValueType>
class manual_orderer {
  typedef manual_orderer_impl::entry<ValueType> entry;
  typedef manual_orderer_impl::idx_type idx_type;
  typedef manual_orderer_impl::num_bits_type num_bits_type;
public:
  typedef manual_orderer_impl::entry_ref<ValueType> entry_ref;

  // create a ValueType in the heap and an entry_ref refcounted pointer
  // to it.  It's not in the ordering until you put it there using one of
  // this manual_orderer's put_*() methods.
  template <class... Args>
  entry_ref construct(Args&&... args) {
    // TODO: better allocator
    return entry_ref(new entry(manual_orderer_impl::no_idx, std::forward<Args>(args)...));
  }

  // TODO: allow moving (maybe), and deletion

  // put_only puts the first object in the ordering; you can't use it
  // after it's been called once
  void put_only(entry_ref moving) {
    place_at(0, moving);
  }

  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just prior to relative_to.
  void put_before(entry_ref moving, entry_ref relative_to) {
#ifdef AUDIT_ORDERED_STUFF
    std::map<idx_type, entry_ref> old_data = data;
#endif
    const idx_type idx = relative_to.data->idx;
    std::pair<idx_type,idx_type> p = make_room_for_split(idx & ~3, 0, idx & 3);
    move_entry(idx, p.second);
    place_at(p.first, moving);
#ifdef AUDIT_ORDERED_STUFF
    assert(relative_to > moving);
    entry_ref prev;
    for (auto& f : old_data) {
      if (prev) assert(prev < f.second);
      if (f.second < relative_to) assert(f.second < moving);
      if (f.second > relative_to) assert(f.second > moving);
      prev = f.second;
      //std::cerr << f.second.data->idx << ", ";
    }
#endif
  }
  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just after relative_to.
  void put_after(entry_ref moving, entry_ref relative_to) {
#ifdef AUDIT_ORDERED_STUFF
    std::map<idx_type, entry_ref> old_data = data;
#endif
    const idx_type idx = relative_to.data->idx;
    std::pair<idx_type,idx_type> p = make_room_for_split(idx & ~3, 0, idx & 3);
    move_entry(idx, p.first);
    place_at(p.second, moving);
#ifdef AUDIT_ORDERED_STUFF
    assert(relative_to < moving);
    entry_ref prev;
    for (auto& f : old_data) {
      if (prev) assert(prev < f.second);
      if (f.second < relative_to) assert(f.second < moving);
      if (f.second > relative_to) assert(f.second > moving);
      prev = f.second;
      //std::cerr << f.second.data->idx << ", ";
    }
#endif
  }
private:
  // TODO: a custom hashtable for this so we don't need the extra layer of pointers
#ifdef AUDIT_ORDERED_STUFF
  std::map<idx_type, entry_ref> data;
#else
  std::unordered_map<idx_type, entry_ref> data;
#endif
  bool idx_exists(idx_type idx)const { return data.find(idx) != data.end(); }
  
  // The data forms - mathematically, but not in memory - a lenient B-tree with 2,3,or 4 children for each node.
  // When it overflows to 5, it splits into a 3-node and a 2-node.
  // An index's Nth bit-pair indicates which child it is at the Nth level of the tree (from the leaves).
  // Since a node's first child always exists, parent nodes (which have no data)
  // are represented by the ID of the first descendant leaf, plus their size/level.
  std::pair<idx_type,idx_type> make_room_for_split(idx_type prefix, num_bits_type child_shift, idx_type which_child) {
    const idx_type child_size = 1 << child_shift;
    if (idx_exists(prefix+3*child_size)) {
      // split this node to make room for splitting children
      const num_bits_type parent_shift = child_shift + 2;
      const idx_type parent_size = 1 << parent_shift;
      const idx_type next_which_child_mask = parent_size * 0x3;
      std::pair<idx_type,idx_type> p = make_room_for_split(prefix & ~next_which_child_mask, parent_shift, (prefix & next_which_child_mask) >> parent_shift);
      assert(p.first >= prefix);
      assert(p.second >= p.first + parent_size);
      assert(!idx_exists(p.second));
      
      // Don't move which_child, because the caller will move it.
      // Don't put anything into result.first or result.second.
      if (which_child != 3) { move_subtree_if_necessary(
        prefix + 3*child_size, child_shift, /*(which_child < 3) ?*/p.second + 1*child_size/*:  p.second            */); }
      if (which_child != 2) { move_subtree_if_necessary(
        prefix + 2*child_size, child_shift,   (which_child < 2) ?  p.second                 :  p.first + 2*child_size); }
      if (which_child != 1) { move_subtree_if_necessary(
        prefix + 1*child_size, child_shift,   (which_child < 1) ?  p.first  + 2*child_size  :  p.first + 1*child_size); }
      if (which_child != 0) { move_subtree_if_necessary(
        prefix               , child_shift, /*(which_child < 0) ?  p.first  + 1*child_size  :*/p.first               ); }
      
      if (which_child == 0) {
        return std::pair<idx_type,idx_type>(p.first               , p.first  + 1*child_size);
      }
      else if (which_child == 1) {
        return std::pair<idx_type,idx_type>(p.first + 1*child_size, p.first  + 2*child_size);
      }
      else if (which_child == 2) {
        return std::pair<idx_type,idx_type>(p.first + 2*child_size, p.second              );
      }
      else if (which_child == 3) {
        return std::pair<idx_type,idx_type>(p.second,               p.second + 1*child_size);
      }
      else {
        assert(false);
      }
    }

    // Make result.second empty. Don't move which_child, because the caller will handle it.
    assert (which_child < 3);
    if ((which_child < 2) && idx_exists(prefix + 2*child_size)) { move_subtree(prefix + 2*child_size, child_shift, prefix + 3*child_size); }
    if  (which_child < 1)                                       { move_subtree(prefix + 1*child_size, child_shift, prefix + 2*child_size); }
    return std::pair<idx_type,idx_type>(prefix + which_child*child_size, prefix + (which_child+1)*child_size);
  }
  void move_subtree_if_necessary(idx_type prefix, num_bits_type node_shift, idx_type new_prefix) {
    if (prefix != new_prefix) { move_subtree(prefix, node_shift, new_prefix); }
  }
  void move_subtree(idx_type prefix, num_bits_type node_shift, idx_type new_prefix) {
    if (node_shift == 0) {
      move_entry(prefix, new_prefix);
    }
    else {
      const num_bits_type child_shift = node_shift-2;
      const idx_type child_size = 1 << child_shift;
      for (int i = 3; i >= 0; --i) {
        const idx_type child_prefix = prefix + child_size*i;
        if (idx_exists(child_prefix)) {
          move_subtree(child_prefix, child_shift, new_prefix + child_size*i);
        }
      }
    }
  }
  void move_entry(idx_type idx, idx_type new_idx) {
    if (new_idx != idx) {
      auto i = data.find(idx);
      if (i != data.end()) {
        i->second.data->idx = new_idx;
        auto p = data.insert(std::make_pair(new_idx, i->second));
        assert(p.second);
        data.erase(i);
      }
    }
  }
  void place_at(idx_type idx, entry_ref moving) {
    assert(moving.data->idx == manual_orderer_impl::no_idx);
    moving.data->idx = idx;
    auto p = data.insert(std::make_pair(idx, moving));
    assert(p.second);
  }
};

namespace std {
  template<typename ValueType>
  struct hash<typename manual_orderer_impl::entry_ref<ValueType>> {
    public:
    size_t operator()(manual_orderer_impl::entry_ref<ValueType> const& i)const {
      // Nondeterministic - but the ordering of STL hash tables is already nondeterministic.
      return std::hash<decltype(i.data)>()(i.data);
    }
  };
}
#else

namespace manual_orderer_impl {
// O(log n) amortized, maybe O(log^2 n) worst case?

typedef uint64_t idx_type;
typedef uint32_t num_bits_type;

const num_bits_type level_size_shift = 4;
const uint64_t level_size = 1 << level_size_shift;
const num_bits_type supply_size_shift = 2;
const uint64_t supply_size = 1 << supply_size_shift;
const uint32_t desired_supplies_per_level = level_size - supply_size;
const uint32_t max_supplies_per_level = level_size;
// In the time it takes us to use up a supply, move far enough to get all the way from its refill source
// even if it's the first supply among the supplies_per_level supplies drawing from the same source.
// max_gap_speed can be lenient
const uint32_t max_gap_speed = divide(max_supplies_per_level * level_size, supply_size, rounding_strategy<round_up, negative_is_forbidden>());

typedef uint64_t ordered_category;
ordered_category new_unique_ordered_category() {
  static ordered_category a = 0;
  return ++a;
}

struct node;
struct entry_or_node_base {
  entry_or_node_base():parent(nullptr){}
  node* parent;
};

struct entry_base : public entry_or_node_base {
  entry_base():entry_or_node_base(),idx(0),category(new_unique_ordered_category()),ref_count(0){}
  idx_type idx;
  ordered_category category;
  size_t ref_count;
  entry_base* prev();
  entry_base* next();
  void put(entry_base* o, bool after);
  void erase_from_parent();
};

struct node : public entry_or_node_base {
  node(entry_base* sole_child):
    supply_desired_size(supply_size),
    supply_current_size(uint64_t(-1)),
    supply_reserved(supply_desired_size),
    next_refill_progress_emptying_supply(-1),
    children(),
    supply_lower_bound(sole_child),
    next_refill_lower_bound(sole_child)
  {
    children.push_back(sole_child);
  }
private:
  node(uint64_t size, entry_base* end):
    supply_desired_size(size),
    supply_current_size(uint64_t(-1)),
    supply_reserved(supply_desired_size),
    next_refill_progress_emptying_supply(-1),
    children(),
    supply_lower_bound(end),
    next_refill_lower_bound(end)
  {}
  uint64_t supply_desired_size;
  uint64_t supply_current_size;
  uint64_t supply_reserved;
  int64_t next_refill_progress_emptying_supply;
public:
  std::list<entry_or_node_base*> children; // TODO better data structure
private:
  entry_base* supply_lower_bound;
  entry_base* next_refill_lower_bound;
  
public:
  void erase(entry_base* a) {
    validate();
    assert (is_bottom_level());
    bool found = false;
    for (auto i = children.begin(); i != children.end(); ) {
      entry_base* e = (entry_base*)*i;
      if (found) { --e->idx; }
      if (e == a) {
        children.erase(i++);
        found = true;
      }
      else {
        ++i;
      }
    }
    assert (found);
    if (supply_lower_bound == a) {
      supply_lower_bound = (entry_base*)children.back();
    }
    if (next_refill_lower_bound == a) {
      next_refill_lower_bound = (entry_base*)children.back();
    }
    claim_from_supply(1, false);
    reserve_one_from_supply(false);
    validate();
  }
  void insert(entry_base* a, entry_base* existing_child, bool after) {
    validate();
    assert (is_bottom_level());
    a->idx = existing_child->idx;
    reserve_one_from_supply();
    claim_from_supply(1);
    bool found = false;
    for (auto i = children.begin(); i != children.end(); ++i) {
      entry_base* e = (entry_base*)*i;
      if (found) { ++e->idx; }
      if (e == existing_child) {
        a->parent = this;
        children.insert(after ? boost::next(i) : i, a);
        found = true;
        if (!after) { ++e->idx; }
      }
    }
    assert (found);
    if (supply_lower_bound->idx < ((entry_base*)children.back())->idx) {
      assert (children.back() == a);
      if (next_refill_lower_bound == supply_lower_bound) {
        next_refill_lower_bound = a;
      }
      supply_lower_bound = a;
    }
    validate();
  }
  entry_or_node_base* prev_sibling(entry_or_node_base* existing_child) {
    if (children.front() == existing_child) {
      if (parent) {
        entry_or_node_base* prev = parent->prev_sibling(this);
        if (prev) {
          return ((node*)prev)->children.back();
        }
      }
      return nullptr;
    }
    else {
      for (auto i = children.begin(); i != children.end(); ++i) {
        if (*i == existing_child) {
          return *boost::prior(i);
        }
      }
      assert (false);
    }
  }
  entry_or_node_base* next_sibling(entry_or_node_base* existing_child) {
    if (children.back() == existing_child) {
      if (parent) {
        entry_or_node_base* next = parent->next_sibling(this);
        if (next) {
          return ((node*)next)->children.front();
        }
      }
      return nullptr;
    }
    else {
      for (auto i = children.begin(); i != children.end(); ++i) {
        if (*i == existing_child) {
          return (node*)*boost::next(i);
        }
      }
      assert (false);
    }
  }
private:
  bool is_bottom_level()const { return supply_desired_size == supply_size; }
  void require_parent() {
    if (!parent) {
      parent = new node(supply_desired_size << level_size_shift, supply_lower_bound);
      parent->children.push_back(this);
      supply_current_size = supply_desired_size;
      supply_reserved = supply_desired_size;
    }
  }
  node* require_next_sibling(node* existing_child) {
    if (children.back() == existing_child) {
      assert (children.size() <= max_supplies_per_level);
      if (children.size() >= desired_supplies_per_level) {
        require_parent();
        return (node*)parent->require_next_sibling(this)->children.front();
      }
      else {
        node* new_sibling = new node(existing_child->supply_desired_size, existing_child->supply_lower_bound);
        new_sibling->parent = this;
        children.push_back(new_sibling);
        return new_sibling;
      }
    }
    else {
      for (auto i = children.begin(); i != children.end(); ++i) {
        if (*i == existing_child) {
          return (node*)*boost::next(i);
        }
      }
      assert (false);
    }
  }
  void claim_from_supply(uint64_t amount, bool reserving = true) {
    if (reserving) {
      assert (supply_reserved >= amount);
      supply_reserved -= amount;
      supply_current_size -= amount;
    }
    else {
      supply_reserved += amount;
      supply_current_size += amount;
    }
  }
  void join_refill_to_supply(bool reserving = true) {
    if (reserving) {
      assert (next_refill_lower_bound == supply_lower_bound);
      assert (next_refill_progress_emptying_supply == -1);
      parent->claim_from_supply(supply_desired_size);
      supply_current_size += supply_desired_size;
      next_refill_lower_bound = parent->supply_lower_bound;
      next_refill_progress_emptying_supply = supply_desired_size;
    }
    else {
      assert (next_refill_lower_bound == parent->supply_lower_bound);
      assert (next_refill_progress_emptying_supply == supply_desired_size);
      assert (supply_current_size >= supply_desired_size);
      parent->claim_from_supply(supply_desired_size, false);
      supply_current_size -= supply_desired_size;
      next_refill_lower_bound = supply_lower_bound;
      next_refill_progress_emptying_supply = -1;
    }
  }
  int64_t progress_emptying_supply() {
    if (supply_desired_size < (1ULL<<level_size_shift)) {
      return supply_lower_bound->idx & (supply_desired_size - 1);
    }
    return supply_lower_bound->idx & (supply_desired_size - 1) & ~((supply_desired_size>>level_size_shift) - 1);
  }
  void roll_refill_once(bool reserving = true) {
    if (reserving) {
      assert (next_refill_lower_bound->idx >= supply_lower_bound->idx);
      if (next_refill_lower_bound == supply_lower_bound) {
        if (next_refill_progress_emptying_supply != supply_desired_size) { assert (next_refill_progress_emptying_supply == progress_emptying_supply()); }
        assert (supply_lower_bound);
        entry_base* prev =  supply_lower_bound->prev();
        assert (prev);
        supply_lower_bound = prev;
        const int64_t new_progress_emptying_supply = progress_emptying_supply();
        if (new_progress_emptying_supply > next_refill_progress_emptying_supply) {
          next_refill_progress_emptying_supply = -1;
        }
        else {
          next_refill_progress_emptying_supply = new_progress_emptying_supply;
        }
        
        const idx_type last_child_boundary_idx = is_bottom_level() ?
          ((entry_base*)children.back())->idx : ((node*)children.back())->supply_lower_bound->idx;
        if (last_child_boundary_idx > supply_lower_bound->idx) {
          node* new_child_parent = parent->require_next_sibling(this);
          children.back()->parent = new_child_parent;
          if (new_child_parent->children.empty()) {
            new_child_parent->supply_lower_bound = (entry_base*)children.back();
            new_child_parent->next_refill_lower_bound = (entry_base*)children.back();
          }
          new_child_parent->children.push_front(children.back());
          children.pop_back();
          new_child_parent->validate();
        }
      }
      next_refill_lower_bound->idx += supply_desired_size;
      entry_base* next = next_refill_lower_bound->next();
      if (next) { assert (next_refill_lower_bound->idx < next->idx); }
      assert (next_refill_lower_bound);
      entry_base* prev =  next_refill_lower_bound->prev();
      assert (prev);
      next_refill_lower_bound = prev;
      assert (next_refill_lower_bound->idx >= supply_lower_bound->idx);
    }
    else {
      if (next_refill_progress_emptying_supply != supply_desired_size) {
        if (next_refill_progress_emptying_supply != -1) { assert (next_refill_progress_emptying_supply == progress_emptying_supply()); }
        supply_lower_bound = supply_lower_bound->next();
        const int64_t new_progress_emptying_supply = progress_emptying_supply();
        if (new_progress_emptying_supply < next_refill_progress_emptying_supply) {
          next_refill_progress_emptying_supply = supply_desired_size;
        }
        else {
          next_refill_progress_emptying_supply = new_progress_emptying_supply;
        }
        
        node* possible_new_child_parent = (node*)parent->next_sibling(this);
        if (possible_new_child_parent) {
          possible_new_child_parent->validate();
          const idx_type first_child_boundary_idx = is_bottom_level() ?
            ((entry_base*)possible_new_child_parent->children.front())->idx : ((node*)possible_new_child_parent->children.front())->supply_lower_bound->idx;
          if (first_child_boundary_idx <= supply_lower_bound->idx) {
            possible_new_child_parent->children.front()->parent = this;
            children.push_back(possible_new_child_parent->children.front());
            possible_new_child_parent->children.pop_front();
            if (possible_new_child_parent->children.empty()) {
              assert (possible_new_child_parent == parent->children.back());
              delete possible_new_child_parent;
              parent->children.pop_back();
              assert (this == parent->children.back());
              assert (parent->children.size() < desired_supplies_per_level);
            }
            else {
              possible_new_child_parent->validate();
            }
          }
        }
      }
      next_refill_lower_bound = next_refill_lower_bound->next();
      next_refill_lower_bound->idx -= supply_desired_size;
      assert (next_refill_lower_bound->idx > next_refill_lower_bound->prev()->idx);
    }
  }
  void reserve_one_from_supply(bool reserving = true) {
    validate();
    
    if (reserving) {
      if (children.size() >= desired_supplies_per_level) {
        require_parent();
      }
      ++supply_reserved;
        
      if (parent) {
        parent->reserve_one_from_supply();
        if (supply_current_size < supply_reserved) {
          join_refill_to_supply();
        }
        for (uint32_t i = 0; (i < max_gap_speed) && (next_refill_progress_emptying_supply >= 0); ++i) {
          roll_refill_once();
        }
      }
    }
    else {
      --supply_reserved;
      
      if (parent) {
        for (uint32_t i = 0; (i < max_gap_speed) && (next_refill_lower_bound != parent->supply_lower_bound); ++i) {
          roll_refill_once(false);
        }
        if (supply_current_size-supply_reserved > supply_desired_size) {
          join_refill_to_supply(false);
        }
        parent->reserve_one_from_supply(false);
      
        if ((children.size() < desired_supplies_per_level) && (parent->parent == nullptr) && (parent->children.front() == this)) {
          assert (parent->children.size() == 1);
          delete parent;
          parent = nullptr;
        }
      }
    }
    
    validate();
  }
  void validate() {
    assert (supply_current_size >= supply_reserved);
    assert (children.size() <= max_supplies_per_level);
    if (parent) { assert (parent->supply_lower_bound->idx >= supply_lower_bound->idx); }
    assert (supply_lower_bound->idx <= next_refill_lower_bound->idx);
    for (auto child : children) { assert(child->parent == this); }
    assert (!children.empty());
    if (is_bottom_level()) {
      idx_type i = 0;
      for (auto child : children) {
        assert(child == children.front() || ((entry_base*)child)->idx > i);
        i = ((entry_base*)child)->idx;
      }
      assert (supply_lower_bound->idx >= i);
    }
  }
};

entry_base* entry_base::prev() {
  return (entry_base*)parent->prev_sibling(this);
}
entry_base* entry_base::next() {
  return (entry_base*)parent->next_sibling(this);
}
void entry_base::put(entry_base* o, bool after) {
  erase_from_parent();
  category = o->category;
  if (!o->parent) {
    o->parent = new node(o);
  }
  o->parent->insert(this, o, after);
}
void entry_base::erase_from_parent() {
  if (parent) {
    parent->erase(this);
    if ((parent->children.size() == 1) && (parent->parent == nullptr)) {
      parent->children.front()->parent = nullptr;
      delete parent;
    }
    parent = nullptr;
  }
}
  
// a manually_orderable is a non-threadsafe reference-counted smart pointer
// that is LessThanComparable by a idx it stores next to the contents,
// and EqualityComparable by the pointer.
// You must not compare entry_refs for which ... TODO explain

template<typename ValueType>
struct entry : public entry_base {
  template <class... Args>
  entry(Args&&... args):entry_base(),contents(std::forward<Args>(args)...){}
  ValueType contents;
};

struct construct {};
template<typename ValueType>
struct manually_orderable {
public:
  template <class... Args>
  manually_orderable(construct, Args&&... args):data(new entry<ValueType>(std::forward<Args>(args)...)){ inc_ref(); }
  manually_orderable():data(nullptr){}
  manually_orderable(manually_orderable const& o):data(o.data) { inc_ref(); }
  template<typename OtherValueType>
  friend struct manually_orderable;
  
  template<typename OtherValueType>
  void put_before(manually_orderable<OtherValueType> const& o) {
    data->put(o.data, false);
  }
  template<typename OtherValueType>
  void put_after(manually_orderable<OtherValueType> const& o) {
    data->put(o.data, true);
  }
  
  template<typename OtherValueType>
  bool operator< (manually_orderable<OtherValueType> const& o)const {
    caller_correct_if(data->category == o.data->category, "comparing items not ordered relative to each other");
    return data->idx <  o.data->idx;
  }
  template<typename OtherValueType>
  bool operator> (manually_orderable<OtherValueType> const& o)const {
    caller_correct_if(data->category == o.data->category, "comparing items not ordered relative to each other");
    return data->idx >  o.data->idx;
  }
  template<typename OtherValueType>
  bool operator<=(manually_orderable<OtherValueType> const& o)const {
    caller_correct_if(data->category == o.data->category, "comparing items not ordered relative to each other");
    return data->idx <= o.data->idx;
  }
  template<typename OtherValueType>
  bool operator>=(manually_orderable<OtherValueType> const& o)const {
    caller_correct_if(data->category == o.data->category, "comparing items not ordered relative to each other");
    return data->idx >= o.data->idx;
  }
  bool operator==(manually_orderable const& o)const { return data == o.data; }
  bool operator!=(manually_orderable const& o)const { return data != o.data; }
  explicit operator bool()const { return bool(data); }
  ValueType& operator*()const { return data->contents; }
  ValueType* operator->()const { return &data->contents; }
  manually_orderable& operator=(manually_orderable const& o) {
    o.inc_ref(); // do this before dec_ref in case lhs and rhs are the same object
    dec_ref();
    data = o.data;
    return *this;
  }
  ~manually_orderable() { dec_ref(); }
private:
  void inc_ref()const {
    if (data) { ++data->ref_count; }
  }
  void dec_ref() {
    if (data && (--data->ref_count == 0)) {
      data->erase_from_parent();
      delete data;
      data = nullptr;
    }
  }
  explicit manually_orderable(entry<ValueType>* data):data(data){ inc_ref(); }
  entry<ValueType>* data;
  friend struct std::hash<manually_orderable>;
};

} // end namespace manual_orderer_impl

namespace std {
  template<typename ValueType>
  struct hash<manual_orderer_impl::manually_orderable<ValueType>> {
    public:
    size_t operator()(manual_orderer_impl::manually_orderable<ValueType> const& i)const {
      // Nondeterministic - but the ordering of STL hash tables is already nondeterministic.
      return std::hash<decltype(i.data)>()(i.data);
    }
  };
}

// manual_orderer is a data structure that lets you place
// arbitrary objects in an O(1)-comparable order, but you have to refer
// to them by manually_orderable instead of the object itself.
template<typename ValueType>
class manual_orderer {
public:
  typedef manual_orderer_impl::manually_orderable<ValueType> manually_orderable;
  typedef manually_orderable entry_ref;

  // create a ValueType in the heap and an manually_orderable refcounted pointer
  // to it.  It's not in the ordering until you put it there using one of
  // this manual_orderer's put_*() methods.
  template <class... Args>
  manually_orderable construct(Args&&... args) {
    // TODO: better allocator
    return manually_orderable(manual_orderer_impl::construct(), std::forward<Args>(args)...);
  }

  // put_only puts the first object in the ordering; you can't use it
  // after it's been called once
  void put_only(manually_orderable) {}

  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just prior to relative_to.
  void put_before(manually_orderable moving, manually_orderable relative_to) {
    moving.put_before(relative_to);
  }
  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just after relative_to.
  void put_after(manually_orderable moving, manually_orderable relative_to) {
    moving.put_after(relative_to);
  }
};
#endif

#endif
