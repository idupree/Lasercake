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

// O(log n) amortized, maybe O(log^2 n) worst case?


const num_bits_type level_size_shift = 4;
const uint64_t level_size = 1 << level_size_shift;
const num_bits_type supply_size_shift = 2;
const uint64_t supply_size = 1 << supply_size_shift;
const uint32_t desired_supplies_per_level = (level_size-supply_size) / supply_size;
const uint32_t max_supplies_per_level = level_size / supply_size;
// In the time it takes us to use up a supply, move far enough to get all the way from its refill source
// even if it's the first supply among the supplies_per_level supplies drawing from the same source.
// max_gap_speed can be lenient
const uint32_t max_gap_speed = divide(max_supplies_per_level * level_size, supply_size, rounding_strategy<round_up, negative_is_forbidden>);

typedef uint64_t ordered_category;
ordered_category new_unique_ordered_category() {
  static ordered_category a = 0;
  return ++ordered_category;
}

struct entry_base {
  entry_base():idx(0),category(new_unique_ordered_category()),ref_count(0),parent(nullptr){}
  idx_type idx;
  ordered_category category;
  size_t ref_count;
  node* parent;
  entry_base* prev() {
    return parent->prev_sibling(this);
  }
  entry_base* next() {
    return parent->next_sibling(this);
  }
  void put(entry_base* o, bool after) {
    if (parent) {
      parent->erase(this);
    }
    category = o->category;
    o->parent->insert(a, o, after);
  }
};

struct node {
  uint64_t supply_desired_size;
  uint64_t supply_current_size;
  uint64_t supply_reserved;
  int64_t next_refill_progress_emptying_supply;
  node* parent;
  typedef (node* or entry_base*) child;
  std::deque<child> children;
  entry_base* supply_lower_bound;
  entry_base* next_refill_lower_bound;
  
  void erase(entry_base* a) {
    assert (is_bottom_level());
    bool found = false;
    for (auto i = children.begin(); i != children.end(); ) {
      if (found) { --i->idx; }
      if (*i == a) {
        children.erase(i++);
        found = true;
      }
      else {
        ++i;
      }
    }
    reserve_one_from_supply(false);
  }
  void insert(entry_base* a, entry_base* existing_child, bool after) {
    assert (is_bottom_level());
    a->idx = existing_child->idx;
    reserve_one_from_supply();
    bool found = false;
    for (auto i = children.begin(); i != children.end(); ++i) {
      if (found) { ++i->idx; }
      if (*i == existing_child) {
        children.insert(after ? boost::next(i) : i, a);
        found = true;
        if (!after) { ++i->idx; }
      }
    }
  }
  
  child prev_sibling(child existing_child) {
    if (children.front() == existing_child) {
      if (parent) {
        return parent->prev_sibling(this)->children.back();
      }
      else {
        return nullptr;
      }
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
  child next_sibling(child existing_child) {
    if (children.back() == existing_child) {
      if (parent) {
        return parent->prev_sibling(this)->children.front();
      }
      else {
        return nullptr;
      }
    }
    else {
      for (auto i = children.begin(); i != children.end(); ++i) {
        if (*i == existing_child) {
          return *boost::next(i);
        }
      }
      assert (false);
    }
  }
  void require_parent() {
    if (!parent) {
      parent = new node();
      parent->supply_lower_bound = supply_lower_bound;
      parent->supply_desired_size = supply_desired_size << level_size_shift;
      parent->supply_current_size = parent->supply_desired_size;
      parent->supply_reserved = 0;
      parent->children.push_back(this);
    }
  }
  node* require_next_sibling(node* existing_child) {
    if (children.back() == existing_child) {
      assert (children.size() <= max_supplies_per_level);
      if (children.size() >= desired_supplies_per_level) {
        require_parent();
        return parent->next_sibling(this)->children.front();
      }
      else {
        node* new_sibling = new node();
        new_sibling->supply_lower_bound = existing_child->supply_lower_bound;
        new_sibling->supply_desired_size = existing_child->supply_desired_size;
        new_sibling->supply_current_size = new_sibling->supply_desired_size;
        new_sibling->supply_reserved = 0;
        children.push_back();
      }
    }
    else {
      for (auto i = children.begin(); i != children.end(); ++i) {
        if (*i == existing_child) {
          return *boost::next(i);
        }
      }
      assert (false);
    }
  }
  void claim_from_supply(uint64_t amount, bool reserving = true) {
    if (reserving) {
      assert (reserved >= amount);
      reserved -= amount;
      current_size -= amount;
    }
    else {
      reserved += amount;
      current_size += amount;
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
    return supply_lower_bound->idx & (supply_desired_size - 1) & ~((supply_desired_size>>level_size_shift) - 1);
  }
  void roll_refill_once(bool reserving = true) {
    if (reserving) {
      if (next_refill_lower_bound == supply_lower_bound) {
        if (next_refill_progress_emptying_supply != supply_desired_size) { assert (next_refill_progress_emptying_supply == progress_emptying_supply()); }
        supply_lower_bound = supply_lower_bound->prev();
        const int64_t new_progress_emptying_supply = progress_emptying_supply();
        if (new_progress_emptying_supply > next_refill_progress_emptying_supply) {
          next_refill_progress_emptying_supply = -1;
        }
        else {
          next_refill_progress_emptying_supply = new_progress_emptying_supply;
        }
        
        if (is_bottom_level()) {
          if (children.back()->idx > supply_lower_bound->idx) {
            node* new_child_parent = parent->require_next_sibling(this);
            children.back()->parent = new_child_parent;
            new_child_parent->children.push_front(children.back());
            children.pop_back();
          }
        }
        else {
          if (children.back()->supply_lower_bound->idx > supply_lower_bound->idx) {
            node* new_child_parent = parent->require_next_sibling(this);
            children.back()->parent = new_child_parent;
            new_child_parent->children.push_front(children.back());
            children.pop_back();
          }
        }
      }
      next_refill_lower_bound->idx += supply_desired_size;
      assert (next_refill_lower_bound->idx < next_refill_lower_bound->next()->idx);
      next_refill_lower_bound = next_refill_lower_bound->prev();
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
        
        node* possible_new_child_parent = parent->next_sibling(this);
        if (is_bottom_level()) {
          if (possible_new_child_parent->children.front()->idx <= supply_lower_bound->idx) {
            possible_new_child_parent->children.front()->parent = this;
            children.push_back(possible_new_child_parent->children.front());
            possible_new_child_parent->children.pop_front();
          }
        }
        else {
          if (children.front()->supply_lower_bound->idx <= supply_lower_bound->idx) {
            possible_new_child_parent->children.front()->parent = this;
            children.push_back(possible_new_child_parent->children.front());
            possible_new_child_parent->children.pop_front();
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
        if (current_size < supply_reserved) {
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
        if (current_size-supply_reserved > supply_desired_size) {
          join_refill_to_supply(false);
        }
        parent->reserve_one_from_supply(false);
      
        if (children.size() < desired_supplies_per_level) {
          assert (parent->parent == nullptr);
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
    if (parent) { assert (parent->supply_lower_bound >= supply_lower_bound); }
  }
};

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

template<typename ValueType>
struct manually_orderable {
public:
  template <class... Args>
  manually_orderable(Args&&... args):data(new entry<ValueType>(std::forward<Args>(args)...)){}
  manually_orderable():data(nullptr){}
  manually_orderable(manually_orderable const& o):data(o.data) { inc_ref(); }
  template<typename OtherValueType>
  friend struct manually_orderable;
  
  template<typename OtherValueType>
  void put_before(manually_orderable<OtherValueType> const& o) {
    data->put_before(o.data);
  }
  template<typename OtherValueType>
  void put_after(manually_orderable<OtherValueType> const& o) {
    data->put_after(o.data);
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
      delete data;
      data = nullptr;
    }
  }
  explicit manually_orderable(entry<ValueType>* data):data(data){ inc_ref(); }
  entry<ValueType>* data;
  friend class manual_orderer<ValueType>;
  friend struct std::hash<manually_orderable>;
};

namespace std {
  template<typename ValueType>
  struct hash<typename manually_orderable<ValueType>> {
    public:
    size_t operator()(manually_orderable<ValueType> const& i)const {
      // Nondeterministic - but the ordering of STL hash tables is already nondeterministic.
      return std::hash<decltype(i.data)>()(i.data);
    }
  };
}

#endif
