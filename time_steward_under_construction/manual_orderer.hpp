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
#define NEW_AMORTIZED_MANUAL_ORDERER

#ifndef WORST_CASE_MANUAL_ORDERER
#ifndef NEW_AMORTIZED_MANUAL_ORDERER
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

const uint64_t no_idx = std::numeric_limits<uint64_t>::max();
class manual_orderer_base;
struct entry_base {
  entry_base():idx(no_idx),ref_count(0),prev(nullptr),next(nullptr),owner(nullptr){}
  uint64_t idx;
  size_t ref_count;
  entry_base* prev;
  entry_base* next;
  manual_orderer_base* owner;
};

template<typename ValueType>
struct entry : public entry_base {
  template <class... Args>
  entry(Args&&... args):entry_base(),contents(std::forward<Args>(args)...){}
  ValueType contents;
};

template<typename ValueType>
struct NoDecRefCallback { void operator()(ValueType&, size_t)const{} };

template<typename ValueType, class DecRefCallback = NoDecRefCallback<ValueType>> class manual_orderer;

template<typename ValueType, class DecRefCallback = NoDecRefCallback<ValueType>>
struct manually_orderable {
public:
  template <class... Args>
  static manually_orderable construct(Args&&... args) {
    return manually_orderable(new entry<ValueType>(std::forward<Args>(args)...));
  }
  manually_orderable():data(nullptr){}
  manually_orderable(manually_orderable const& o):data(o.data) { inc_ref(); }
  
  bool operator< (manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx <  o.data->idx;
  }
  bool operator> (manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx >  o.data->idx;
  }
  bool operator<=(manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx <= o.data->idx;
  }
  bool operator>=(manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
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
    if (data) {
      DecRefCallback()(data->contents, data->ref_count);
      if (data) {
        --data->ref_count;
        if (data->ref_count == 0) {
          if (data->owner) {
            data->owner->erase(data);
          }
          entry<ValueType>* was_data = data;
          data = nullptr;
          delete was_data;
        }
      }
    }
  }
  explicit manually_orderable(entry<ValueType>* data):data(data){ inc_ref(); }
  entry<ValueType>* data;
  friend struct std::hash<manually_orderable>;
  friend class manual_orderer<ValueType, DecRefCallback>;
};

const uint32_t max_children_per_block_shift = 4;
const uint64_t max_children_per_block = 1ULL << max_children_per_block_shift;
const uint64_t min_children_per_block = 5;

constexpr inline uint64_t block_size(uint32_t level) {
  return 1ULL << (level*max_children_per_block_shift);
}
constexpr inline uint64_t position_in_block_mask(uint32_t level) {
  return block_size(level)-1;
}
constexpr inline uint64_t block_start_mask(uint32_t level) {
  return ~position_in_block_mask(level);
}
constexpr inline uint64_t block_start(uint64_t idx, uint32_t level) {
  return idx & block_start_mask(level);
}
constexpr inline uint64_t next_block_start(uint64_t idx, uint32_t level) {
  return (idx & block_start_mask(level)) + block_size(level);
}
constexpr inline uint64_t position_in_block(uint64_t idx, uint32_t level) {
  return idx & position_in_block_mask(level);
}
constexpr inline uint64_t which_child_is_block(uint64_t idx, uint32_t level) {
  return (block_start(idx, level) & position_in_block_mask(level+1)) / block_size(level);
}

class manual_orderer_base {
protected:
  std::unordered_map<uint64_t, entry_base*> data;
  entry_base* last_entry;
  
  entry_base* get(uint64_t idx) {
    auto i = data.find(idx);
    if (i == data.end()) { return nullptr; }
    return i->second;
  }
  void set(uint64_t idx, entry_base* e) {
    if (e) {
      data[idx] = e;
    }
    else {
      data.erase(idx);
    }
  }
  
  entry_base* first_in_block(uint64_t idx, uint32_t level) {
    return get(block_start(idx, level));
  }
  entry_base* last_in_block(uint64_t idx, uint32_t level) {
    uint64_t after_this = next_block_start(idx, level);
    auto a = get(after_this);
    while (!a) {
      if (last_entry->idx < after_this) { return last_entry; }
      ++level;
      after_this = next_block_start(idx, level);
      a = get(after_this);
    }
    assert (a->prev->idx < after_this);
    return a->prev;
  }
  bool block_is_full(uint64_t idx, uint32_t level) {
    uint64_t after_this = next_block_start(idx, level);
    bool result = last_in_block(idx, level)->idx >= after_this - block_size(level-1);
    assert (result == (num_children_in_block(idx, level) == max_children_per_block));
    return result;
  }
  uint64_t num_children_in_block(uint64_t idx, uint32_t level, bool lenient = false) {
    entry_base* l = last_in_block(idx, level);
    if (l->idx < block_start(idx, level)) { return 0; }
    uint64_t num_children = which_child_is_block(l->idx, level-1) + 1;
    assert (num_children <= max_children_per_block);
    if (l->next) {
      assert (num_children >= min_children_per_block - lenient);
    }
    return num_children;
  }
  
  entry_base* move_entries_up(entry_base* first, uint64_t last, uint64_t dist) {
    auto i = first;
    for (; i && i->idx >= last; i = i->prev) {
      set(i->idx, nullptr);
      i->idx += dist;
      assert (!i->next || i->next->idx > i->idx);
      set(i->idx, i);
    }
    return i;
  }
  entry_base* move_entries_down(entry_base* first, uint64_t last, uint64_t dist) {
    auto i = first;
    for (; i && i->idx <= last; i = i->next) {
      set(i->idx, nullptr);
      i->idx -= dist;
      assert (!i->prev || i->prev->idx < i->idx);
      set(i->idx, i);
    }
    return i;
  }
  
  ~manual_orderer_base() {
    for (auto p : data) {
      entry_base* e = p.second;
      e->idx = no_idx;
      e->prev = nullptr;
      e->next = nullptr;
      e->owner = nullptr;
    }
  }
  
  void validate() {
    return;
    for (auto p : data) {
      entry_base* e = p.second;
      if (e->prev) {
        assert (e->prev->next == e);
        assert (e->prev->idx < e->idx);
      }
      if (e->next) {
        assert (e->next->prev == e);
        assert (e->next->idx > e->idx);
      }
      else {
        assert (e == last_entry);
      }
      assert (e->idx != no_idx);
      assert (e->owner == this);
    }
    
    for (uint32_t level = 1; level < 15; ++level) {
      for (uint64_t start = 0; start <= last_entry->idx; start += block_size(level)) {
        bool f = true;
        for (uint64_t cstart = start; cstart < start + block_size(level); cstart += block_size(level-1)) {
          if (!get(cstart)) { f = false; }
          else { assert(f); }
        }
      }
    }
  }

public:
  void insert(entry_base* inserted_entry, entry_base* existing_entry, bool after) {
    validate();
    uint32_t nonfull_level = 1;
    while (block_is_full(existing_entry->idx, nonfull_level)) {
      ++nonfull_level;
      caller_correct_if(nonfull_level*max_children_per_block_shift < 64, "manual_orderer overflowed");
    }
    
    for (uint32_t level = nonfull_level; level != 0; --level) {
      uint64_t shove_dist = block_size(level-1);
      uint64_t split_dist = shove_dist >> 1;
      uint64_t last_shoved = next_block_start(existing_entry->idx, level-1);
      uint64_t last_split = last_shoved - split_dist;
      if (level == 1 && !after) {
        --last_shoved;
      }
      auto i = last_in_block(existing_entry->idx, level);
      i = move_entries_up(i, last_shoved, shove_dist);
      move_entries_up(i, last_split, split_dist);
    }
    
    if (after) {
      set(existing_entry->idx+1, inserted_entry);
      inserted_entry->idx = existing_entry->idx+1;
      inserted_entry->prev = existing_entry;
      inserted_entry->next = existing_entry->next;
      if (existing_entry->next) { existing_entry->next->prev = inserted_entry; }
      existing_entry->next = inserted_entry;
      if (last_entry == existing_entry) { last_entry = inserted_entry; }
      
      assert (inserted_entry->idx > existing_entry->idx);
      assert (!inserted_entry->next || inserted_entry->idx < inserted_entry->next->idx);
    }
    else {
      set(existing_entry->idx-1, inserted_entry);
      inserted_entry->idx = existing_entry->idx-1;
      inserted_entry->prev = existing_entry->prev;
      inserted_entry->next = existing_entry;
      if (existing_entry->prev) { existing_entry->prev->next = inserted_entry; }
      existing_entry->prev = inserted_entry;
      
      assert (inserted_entry->idx < existing_entry->idx);
      assert (!inserted_entry->prev || inserted_entry->idx > inserted_entry->prev->idx);
    }
    inserted_entry->owner = this;
    for (uint32_t level = nonfull_level; level != 0; --level) {
      assert (first_in_block(existing_entry->idx, level));
    }
    validate();
  }
  
  void erase(entry_base* erased_entry) {
    validate();
    uint64_t idx = erased_entry->idx;
    if (last_entry == erased_entry) { last_entry = erased_entry->prev; }
    if (erased_entry->prev) { erased_entry->prev->next = erased_entry->next; }
    if (erased_entry->next) { erased_entry->next->prev = erased_entry->prev; }
    erased_entry->idx = no_idx;
    erased_entry->prev = nullptr;
    erased_entry->next = nullptr;
    erased_entry->owner = nullptr;
    set(idx, nullptr);
    
    for (uint32_t level = 1; ; ++level) {
      move_entries_down(first_in_block(idx+block_size(level-1), level-1), next_block_start(idx, level)-1, block_size(level-1));
      
      uint64_t num_children = num_children_in_block(idx, level, true);
      if (num_children >= min_children_per_block) {
        break;
      }
      else {
        uint64_t which_child_are_we = which_child_is_block(idx, level);
        uint64_t prev_num_children = 0;
        uint64_t next_num_children = 0;
        uint64_t somewhere_in_next = idx + block_size(level);
        uint64_t somewhere_in_prev = idx - block_size(level);
        if (which_child_are_we > 0) {
          prev_num_children = num_children_in_block(somewhere_in_prev, level);
        }
        if (which_child_are_we < max_children_per_block-1) {
          next_num_children = num_children_in_block(somewhere_in_next, level);
        }
        assert (prev_num_children || next_num_children);
        bool prefer_prev = prev_num_children && ((!next_num_children) || (prev_num_children < next_num_children));
        if (prefer_prev) {
          if (prev_num_children + num_children <= max_children_per_block) {
            uint64_t retract_dist = block_size(level) - prev_num_children*block_size(level-1);
            uint64_t last_retracted = next_block_start(idx, level)-1;
            move_entries_down(first_in_block(idx, level), last_retracted, retract_dist);
            assert (first_in_block(somewhere_in_prev, level));
          }
          else {
            uint64_t transferred_children = (prev_num_children - num_children) >> 1;
            assert (prev_num_children - transferred_children >= max_children_per_block/2);
            assert (num_children + transferred_children >= max_children_per_block/2);
            auto i = last_in_block(idx, level);
            i = move_entries_up(i, block_start(idx, level), transferred_children*block_size(level-1));
            uint64_t steal_end = block_start(somewhere_in_prev, level) + (prev_num_children-transferred_children)*block_size(level-1);
            move_entries_up(i, steal_end, block_start(idx, level) - steal_end);
            assert (first_in_block(idx, level));
            assert (first_in_block(somewhere_in_prev, level));
            break;
          }
        }
        else {
          if (next_num_children + num_children <= max_children_per_block) {
            uint64_t retract_dist = block_size(level) - num_children*block_size(level-1);
            uint64_t last_retracted = next_block_start(somewhere_in_next, level)-1;
            move_entries_down(first_in_block(somewhere_in_next, level), last_retracted, retract_dist);
            assert (first_in_block(idx, level));
            idx = somewhere_in_next;
          }
          else {
            uint64_t transferred_children = (next_num_children - num_children) >> 1;
            assert (next_num_children - transferred_children >= max_children_per_block/2);
            assert (num_children + transferred_children >= max_children_per_block/2);
            uint64_t steal_end = block_start(somewhere_in_next, level) + transferred_children*block_size(level-1) - 1;
            auto i = first_in_block(somewhere_in_next, level);
            i = move_entries_down(i, steal_end, block_size(level) - num_children*block_size(level-1));
            move_entries_down(i, next_block_start(somewhere_in_next, level)-1, transferred_children*block_size(level-1));
            assert (first_in_block(idx, level));
            assert (first_in_block(somewhere_in_next, level));
            break;
          }
        }
      }
    }
    validate();
  }
};


// manual_orderer is a data structure that lets you place
// arbitrary objects in an O(1)-comparable order, but you have to refer
// to them by manually_orderable instead of the object itself.
template<typename ValueType, class DecRefCallback>
class manual_orderer : public manual_orderer_base {
public:
  typedef manually_orderable<ValueType, DecRefCallback> entry_ref;

  // put_only puts the first object in the ordering; you can't use it
  // after it's been called once
  void put_only(entry_ref m) {
    m.data->prev = nullptr;
    m.data->next = nullptr;
    m.data->owner = this;
    m.data->idx = 0;
    last_entry = m.data;
    set(0, m.data);
  }

  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just prior to relative_to.
  void put_before(entry_ref moving, entry_ref relative_to) {
    insert(moving.data, relative_to.data, false);
  }
  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just after relative_to.
  void put_after(entry_ref moving, entry_ref relative_to) {
    insert(moving.data, relative_to.data, true);
  }
};
  
  
  
  
  /*
  
uint64_t supply_size(uint32_t level) {
  return 1ULL << (level*2);
}
uint64_t block_size(uint32_t level) {
  uint64_t result = 4;
  uint32_t level2 = level*2;
  for (uint32_t i = 0; i < level2; i += 2) {
    result <<= 2;
    result += (1ULL << i);
  }
  return result;
}
uint64_t block_idx(uint64_t idx, uint32_t level, bool use_next) {
  uint64_t result = 0;
  uint32_t end = level*2 - 2;
  uint64_t block_size_ = max_block_size;
  for (uint32_t i = max_level*2; ; i -= 2) {
    assert (block_size_ = block_size(i));
    uint64_t position_in_block = idx % block_size_;
    result += (idx-position_in_block);
    idx = position_in_block;
    if (i == end) { break; }
    block_size_ -= (1ULL << i);
    block_size_ >>= 2;
  }
  if (use_next) { result += block_size_; }
  return result;
}

void insert(entry_base* inserted_entry, entry_base* existing_entry, bool after) {
  const uint64_t idx = existing_entry->idx;
  std::array<uint64_t, max_possible_levels> block_sizes;
  std::array<uint64_t, max_possible_levels> our_block_starts;
  std::array<uint64_t, max_possible_levels> next_block_starts;
  uint64_t block_size_ = 4;
  for (uint32_t level = 0; ; ++level) {
    assert (block_size_ = block_size(level));
    block_sizes[level] = block_size_;
    if (level == max_level) { break; }
    block_size_ <<= 2;
    block_size_ += supply_size(level);
  }
  uint64_t block_start = 0;
  for (uint32_t level = max_level; ; --level) {
    uint64_t position_in_block = idx % block_sizes[level];
    block_start += (idx-position_in_block);
    idx = position_in_block;
    our_block_starts[level] = block_start;
    assert (block_start == block_idx(idx, level, true));
    next_block_starts[level] = block_start + block_sizes[level];
    if (level == 0) { break; }
    block_sizes[level] -= supply_size(level);
    block_sizes[level] >>= 2;
  }
  
  uint32_t nonfull_level = 0;
  while (true) { 
    auto a = get(next_block_starts[level]);
    auto last_entry_maybe_in_supply;
    if (a) {
      last_entry_maybe_in_supply = a->prev;
    }
    else {
      assert (last->idx < next_block_starts[level]);
      last_entry_maybe_in_supply = last;
    }
    if (last_entry_maybe_in_supply->idx >= next_block_starts[level]-supply_size(level)) {
      ++nonfull_level;
    }
    else { break; }
  }
  
  for (uint32_t level = nonfull_level; level != uint32_t(-1); --level) {
    uint64_t move_dist;
    uint64_t last_moved;
    if (level == 0) {
      move_dist = 1;
      last_moved = idx + after;
    }
    else {
      move_dist = block_sizes[level-1];
      last_moved = next_block_starts[level-1] - supply_size(level-1);
    }
    for (auto i = get(next_block_starts[level])->prev; i && i->idx >= last_moved; i = i->prev) {
      set(i->idx, nullptr);
      i->idx += move_dist;
      set(i->idx, i);
    }
  }
  
  if (after) {
    set(idx+1, inserted_entry);
    inserted_entry->prev = existing_entry;
    inserted_entry->next = existing_entry->next;
    existing_entry->next->prev = inserted_entry;
    existing_entry->next = inserted_entry;
  }
  else {
    set(idx, inserted_entry);
    inserted_entry->prev = existing_entry->prev;
    inserted_entry->next = existing_entry;
    existing_entry->prev->next = inserted_entry;
    existing_entry->prev = inserted_entry;
  }
}

*/
} // end namespace manual_orderer_impl

namespace std {
  template<typename ValueType, class DecRefCallback>
  struct hash<typename manual_orderer_impl::manually_orderable<ValueType, DecRefCallback>> {
    public:
    size_t operator()(manual_orderer_impl::manually_orderable<ValueType, DecRefCallback> const& i)const {
      // Nondeterministic - but the ordering of STL hash tables is already nondeterministic.
      return std::hash<decltype(i.data)>()(i.data);
    }
  };
}

using manual_orderer_impl::manual_orderer;
#endif
#else


/*

node implements:


push_front [O(level) worst case] insert at the beginning; works only (half cap) times
pop_back [O(level) worst case]
insert_after: O(level^2) worst case insert anywhere, as long as entries < cap

if we can get push_front and pop_back to be O(1) then the whole thing becomes O(log n).
The reason to seek this particular improvement:
This algorithm already achieves O(log n) *index changes* per insert.
The pop_back-into-push_front operation only changes one index.


struct upper_level_node {
  bool can_push_front() {
    return !(children[0] && !children[0]->can_push_front());
  }
  void push_front(entry e) {
    if (!children[1]) {
      create children[1];
    }
    if (children[1]->can_push_front()) {
      children[1]->push_front(e);
    }
    else {
      if (!children[0]) {
        create children[0];
      }
      assert (children[0]->can_push_front());
      children[0]->push_front(e);
    }
  }
  entry pop_back() {
    entry result = children[last_child]->pop_back();
    if (children[last_child]->empty()) {
      destroy children[last_child];
      --last_child;
    }
    return result;
  }
  void insert_after(inserted, relative_to) {
    child_idx which = TODO;
    if ((which < 2) && (children[0]->entries + children[1]->entries >= entry_cap(level-1))) {
      first of the last two children.push_front(last of the first two children.pop_back());
    }
    children[which].insert_after(inserted, relative_to);
  }
};

struct bottom_level_node {
  bool can_push_front() {
    return !entries[0];
  }
  void push_front(entry e) {
    --first_entry;
    entries[first_entry] = e;
    e.idx = entry_0_idx + first_entry;
  }
  entry pop_back() {
    entry result = entries[last_entry];
    entries[last_entry] = null;
    --last_entry;
    return result;
  }
  void insert_after(inserted, relative_to) {
    child_idx which = relative_to->idx - entry_0_idx;
    for (child_idx i = last_entry; i > which; --i) {
      ++entries[i]->idx;
      entries[i+1] = entries[i];
    }
    ++last_entry;
    inserted->idx = relative_to->idx + 1;
    entries[which+1] = inserted;
  }
};

 */


namespace manual_orderer_impl {

typedef uint64_t idx_type;
const idx_type no_idx = std::numeric_limits<idx_type>::max();
struct entry_base {
  entry_base():idx(no_idx),ref_count(0){}
  idx_type idx;
  size_t ref_count;
};

template<typename ValueType>
struct entry : public entry_base {
  template <class... Args>
  entry(Args&&... args):entry_base(),contents(std::forward<Args>(args)...){}
  ValueType contents;
};


// TODO abstract more of this (e.g. references to "64")
const uint32_t bottom_level_node_bits = 4;
const idx_type bottom_level_node_size = (idx_type(1) << bottom_level_node_bits);
const uint32_t bottom_level = 1;
constexpr idx_type entry_cap(uint32_t level) {
  return bottom_level_node_size << (level - bottom_level);
}
constexpr idx_type node_size(uint32_t level) {
  return bottom_level_node_size << ((level - bottom_level) * 2);
}
const uint32_t max_level = bottom_level + ((64 - bottom_level_node_bits) >> 1);

struct node_base {
  uint32_t level;
  idx_type entries;
  idx_type beginning;
  // default-construct root node
  node_base():level(max_level),entries(0),beginning(0) {}
  node_base(node_base* parent, child_idx which_child):level(parent->level - 1),entries(0),beginning(parent->beginning + node_size(level) * which_child) {}
  bool empty() {
    return entries == 0;
  }
};

struct upper_level_node : public node_base {
  std::array<node_base*, 4> children;
  child_idx first_child;
  child_idx last_child;
  
  bottom_level_node(node_base* parent, child_idx which_child):node_base(parent, which_child),first_child(2),last_child(4){}
  bool can_push_front() {
    if (!children[0]) { return true; }
    if (level - 1 == bottom_level) { return ((bottom_level_node*)children[0])->can_push_front(); }
    else {                           return (( upper_level_node*)children[0])->can_push_front(); }
  }
  void create_child(child_idx which) {
    assert (!children[which]);
    if (level - 1 == bottom_level) { children[which] = new bottom_level_node(this, which); }
    else {                           children[which] = new  upper_level_node(this, which); }
  }
  void destroy_child(child_idx which) {
    assert (children[which]);
    assert (children[which]->empty());
    if (level - 1 == bottom_level) {}
    else {
      for (child_idx i = 0; i < 4; ++i) {
        assert (!(( upper_level_node*)children[which])->children(i));
      }
    }
    delete children[which];
    children[which] = nullptr;
  }
  idx_type first_half_entries() {
    return (children[0] && children[0]->entries) + (children[1] && children[1]->entries);
  }
  void push_front(entry_base* e) {
    ++entries;
    if (first_child > 1) {
      create_child(1);
    }
    if (children[1]->can_push_front()) {
      children[1]->push_front(e);
    }
    else {
      if (!children[0]) {
        create_child(0);
      }
      assert (children[0]->can_push_front());
      children[0]->push_front(e);
    }
  }
  entry_base* pop_back() {
    --entries;
    entry_base* result = children[last_child]->pop_back();
    if (children[last_child]->empty()) {
      destroy_child(last_child);
      --last_child;
    }
    return result;
  }
  void insert(entry_base* inserted_entry, entry_base* existing_entry, bool after) {
    ++entries;
    const child_idx which = which_child(level, existing_entry->idx - beginning);
    if ((which < 2) && (first_half_entries() >= entry_cap(level-1))) {
      const child_idx min_right_child = children[2] ? 2 : 3;
      const child_idx max_left_child = children[1] ? 1 : 0;
      min_right_child.push_front(max_left_child.pop_back());
    }
    if (level - 1 == bottom_level) { return ((bottom_level_node*)children[which])->insert_after(inserted_entry, existing_entry, after); }
    else {                           return (( upper_level_node*)children[which])->insert_after(inserted_entry, existing_entry, after); }
  }
};

struct bottom_level_node : public node_base {
  std::array<entry_base*, bottom_level_node_size> entries;
  child_idx first_entry;
  child_idx last_entry;
  
  bottom_level_node(node_base* parent, child_idx which_child):node_base(parent, which_child),first_entry((bottom_level_node_size>>1)+1),last_entry(bottom_level_node_size>>1){}
  
  bool can_push_front() {
    return !entries[0];
  }
  void push_front(entry_base* e) {
    ++entries;
    --first_entry;
    entries[first_entry] = e;
    e->idx = beginning + first_entry;
  }
  entry_base* pop_back() {
    --entries;
    entry_base* result = entries[last_entry];
    entries[last_entry] = nullptr;
    --last_entry;
    return result;
  }
  void insert(entry_base* inserted_entry, entry_base* existing_entry, bool after) {
    ++entries;
    child_idx which = existing_entry->idx - beginning;
    for (child_idx i = last_entry; i > which; --i) {
      ++entries[i]->idx;
      entries[i+1] = entries[i];
    }
    ++last_entry;
    assert (last_entry < bottom_level_node_size);
    if (after) {
      inserted_entry->idx = existing_entry->idx + 1;
      entries[which+1] = inserted_entry;
    }
    else {
      inserted_entry->idx = existing_entry->idx;
      existing_entry->idx = existing_entry->idx + 1;
      entries[which] = inserted_entry;
      entries[which+1] = existing_entry;
    }
  }
};

template<typename ValueType>
struct NoDecRefCallback { void operator()(ValueType&, size_t)const{} };

template<typename ValueType, class DecRefCallback = NoDecRefCallback<ValueType>> class manual_orderer;

template<typename ValueType, class DecRefCallback = NoDecRefCallback<ValueType>>
struct manually_orderable {
public:
  template <class... Args>
  static manually_orderable construct(Args&&... args) {
    return manually_orderable(new entry<ValueType>(std::forward<Args>(args)...));
  }
  manually_orderable():data(nullptr){}
  manually_orderable(manually_orderable const& o):data(o.data) { inc_ref(); }
  
  bool operator< (manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx <  o.data->idx;
  }
  bool operator> (manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx >  o.data->idx;
  }
  bool operator<=(manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
    return data->idx <= o.data->idx;
  }
  bool operator>=(manually_orderable const& o)const {
    caller_correct_if(data->idx != no_idx && o.data->idx != no_idx, "comparing not-yet-ordered item(s)");
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
    if (data) {
      DecRefCallback()(data->contents, data->ref_count);
      if (data) {
        --data->ref_count;
      }
    }
  }
  explicit manually_orderable(entry<ValueType>* data):data(data){ inc_ref(); }
  entry<ValueType>* data;
  friend struct std::hash<manually_orderable>;
  friend class manual_orderer<ValueType, DecRefCallback>;
};


// manual_orderer is a data structure that lets you place
// arbitrary objects in an O(1)-comparable order, but you have to refer
// to them by manually_orderable instead of the object itself.
template<typename ValueType, class DecRefCallback>
class manual_orderer : public manual_orderer_base {
private:
  upper_level_node* root;
public:
  typedef manually_orderable<ValueType, DecRefCallback> entry_ref;

  // put_only puts the first object in the ordering; you can't use it
  // after it's been called once
  void put_only(entry_ref m) {
    root.push_front(m.data);
  }

  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just prior to relative_to.
  void put_before(entry_ref moving, entry_ref relative_to) {
    root.insert(moving.data, relative_to.data, false);
  }
  // relative_to must already have been put in the ordering.
  // Puts "moving" in the ordering just after relative_to.
  void put_after(entry_ref moving, entry_ref relative_to) {
    root.insert(moving.data, relative_to.data, true);
  }
};
} // end namespace manual_orderer_impl

namespace std {
  template<typename ValueType, class DecRefCallback>
  struct hash<typename manual_orderer_impl::manually_orderable<ValueType, DecRefCallback>> {
    public:
    size_t operator()(manual_orderer_impl::manually_orderable<ValueType, DecRefCallback> const& i)const {
      // Nondeterministic - but the ordering of STL hash tables is already nondeterministic.
      return std::hash<decltype(i.data)>()(i.data);
    }
  };
}

using manual_orderer_impl::manual_orderer;





































namespace manual_orderer_impl {
  
typedef uint64_t idx_type;
const uint32_t node_max_entries_per_level_shift = 4;
constexpr inline idx_type node_max_entries(uint32_t level) {
  return 1ULL << (level*node_max_entries_per_level_shift);
}
constexpr inline idx_type node_idxspace_size(uint32_t level) {
  return 1ULL << (level*(node_max_entries_per_level_shift+1));
}
constexpr inline idx_type node_entries_to_start_splitting_threshold(uint32_t level) {
  return (node_max_entries(1)-1) * node_max_entries(level-1);
}
  
struct wafwaff {
  bottom_level_node* node;
  size_t* which_child;
  entry_base* get() {
    return node->children[which_child];
  }
  void set(entry_base* e) {
    node->children[which_child] = e;
  }
  void inc() {
    ++which_child;
    if (which_child >= node->num_children) {
      node = node->next;
      which_child = 0;
    }
  }
  operator bool() const {
    return bool(node);
  }
}

struct manual_orderer { 
  wafwaff gc_walker_junk_begin;
  wafwaff gc_walker_junk_end;
  uint32_t max_level;
  
  void insert(entry_base* inserted_entry, entry_base* existing_entry, bool after) {
    std::vector<node_base*> path;
    std::vector<size_t> child_indices;
    path.push_back(&root);
    for (size_t i = 0; i < max_level-1; ++i) {
      auto a = path.back()->which_child_contains(existing_entry);
      path.push_back(a.first);
      child_indices.push_back(a.second);
    }
    bottom_level_node* b = (bottom_level_node*)path.back();
    for (size_t i = 0; i < path.size(); ++i) {
      node_base* n = path[i];
      ++n->num_entries;
    }
    for (size_t i = 0; i < path.size()-1; ++i) {
      upper_level_node* u = (upper_level_node*)path[i];
      u->advance_child_resupplies();
    }
    for (size_t i = 1; i < path.size()-1; ++i) {
      upper_level_node* u = (upper_level_node*)path[i];
      u->advance_resupply();
      u->advance_resupply();
    }
    b->advance_resupply();
    b->advance_resupply();
    for (size_t i = 0; i < path.size()-1; ++i) {
      upper_level_node* u = (upper_level_node*)path[i];
      uint32_t level = path.size()-i;
      if (u->num_entries > node_max_entries(level)) {
        u->split((upper_level_node*)path[i-1], this);
      }
    }
    if (b->num_entries > node_max_entries(1)) {
      b->split((upper_level_node*)path[path.size()-2], this);
      
      if (b->entries[0]->idx > existing_entry->idx) {
        b = b->prev;
      }
    }
    b->insert(inserted_entry, existing_entry, after);
  }
  void erase(entry_base* erased_entry) {
    assert (erased_entry->is_junk());
    gc_walk();
    gc_walk();
    gc_drop();
    gc_drop();
  }
  void gc_walk() {
    if (!gc_walker_junk_end) {
      gc_walker_junk_begin.node = first;
      gc_walker_junk_end.node = first;
      gc_walker_junk_begin.which_child = 0;
      gc_walker_junk_end.which_child = 0;
    }
    entry_base* next = gc_walker_junk_end.get();
    if (!next->is_junk()) {
      entry_base* start = gc_walker_junk_begin.get();
      idx_type temp = next->idx;
      next->idx = start->idx;
      start->idx = next;
      gc_walker_junk_begin.set(next);
      gc_walker_junk_end.set(start);
      gc_walker_junk_begin.inc();
    }
    gc_walker_junk_end.inc();
  }
};

idx_type entry_capacity_provided_to_each_new_sibling_by_refill(uint32_t level) {
  return divide(max_dist_to_parent_in_entries(level), refill_entries_moved_per_insert, rounding_strategy<round_up, negative_is_forbidden>());
}
idx_type ideal_refill_size(uint32_t level) {
  return entry_capacity_provided_to_each_new_sibling_by_refill(level)*reserved_per_insert(level);
}
idx_type worst_case_inserts_to_receive_refill_from_source(uint32_t level) {
  return divide(max_dist_to_parent_in_entries(level), refill_entries_moved_per_insert, rounding_strategy<round_up, negative_is_forbidden>());
}

struct node_base {
  bottom_level_node* resupply_next_moved;
  idx_type num_entries;
  
  void advance_resupply_impl(idx_type dist) {
    for (size_t i = 0; i < resupply_next_moved->num_children; ++i) {
      resupply_next_moved->children[i]->idx += dist;
    }
    resupply_next_moved = resupply_next_moved->prev;
  }
  idx_type remaining_capacity() {
    return num_entries * node_max_entries(level);
  }
  idx_type max_dist_to_refill_in_entries() {
    return refill_entries_moved_per_insert * remaining_capacity();
  }
    return resupply_next_moved ? ideal_refill_size(level) : 0;
  }
  idx_type refill_size() {
    return resupply_next_moved ? ideal_refill_size(level) : 0;
  }
  idx_type min_inserts_till_next_refill_receipt() { // "next refill" not including the current approaching one if any
    return remaining_capacity() + (resupply_next_moved ? entry_capacity_provided_to_each_new_sibling_by_refill(level) : 0);
  }
  idx_type min_inserts_till_next_refill_request() {
    return min_inserts_till_next_refill_receipt() - worst_case_inserts_to_receive_refill_from_source(level);
  }
  void validate() {
    // Each node must...
    // ...have enough unreserved supply to last until it fills up
    assert (unreserved_supply() == remaining_capacity()*reserved_per_insert(level));
    // ...have a refill coming if it couldn't get here from the source by the time the node is full
    if (entries_to_parent_supply > max_dist_to_refill_in_entries()) {
      assert (resupply_next_moved);
    }
    // ...have enough reserved from parent supply to provide for the next refill-pair in time
    assert (reserved_from_parent_supply() + reserved_per_insert(level+1)*min_inserts_till_next_refill_request() == ideal_refill_size(level));
    
#if 0
    // Each node must...
    // if a refill is coming,
    if (refill_size()>0) {
      // ...have enough unreserved supply to last until the refill arrives
      assert (unreserved_supply()*refill_entries_moved_per_insert >= reserved_per_insert(level)*refill_dist_to_parent_in_entries());
      // ...have half its refill be big enough to repeat the feat
      assert (refill_size()*refill_entries_moved_per_insert >= 2*reserved_per_insert(level)*max_dist_to_parent_in_entries(level));
      // (thus ideal_refill_size(level) = 2*reserved_per_insert(level)*max_dist_to_parent_in_entries(level) / refill_entries_moved_per_insert)
      // ...run a tight ship
      assert ((unreserved_supply() - reserved_per_insert(level))*refill_entries_moved_per_insert < reserved_per_insert(level)*refill_dist_to_parent_in_entries());
      assert (refill_size() == ideal_refill_size(level));
    }
    else {
      // ...have enough unreserved supply to last until *a* refill arrives
      assert (unreserved_supply()*refill_entries_moved_per_insert >= reserved_per_insert(level)*max_dist_to_parent_in_entries(level));
      // ...run a tight ship
      assert ((unreserved_supply() - reserved_per_insert(level))*refill_entries_moved_per_insert >= reserved_per_insert(level)*max_dist_to_parent_in_entries(level));
    }
    uint64_t inserts_till_next_receipt = (unreserved_supply()+refill_size()) / reserved_per_insert(level);
    uint64_t inserts_till_next_request = inserts_till_next_receipt - (max_dist_to_parent_in_entries(level)/refill_entries_moved_per_insert);
    assert (inserts_till_next_request <= inserts_till_next_receipt); // should be a repeat of the above
    // ...have enough reserved from parent supply to make the next refill in time
    assert (reserved_from_parent_supply() + inserts_till_next_request*reserved_per_insert(level+1) >= ideal_refill_size(level));
    // ...run a tight ship
    assert (reserved_from_parent_supply() + (inserts_till_next_request-1)*reserved_per_insert(level+1) < ideal_refill_size(level));
    
    // in summary, force:
#endif
  }
};
struct bottom_level_node : public node_base {
  bottom_level_node* prev;
  bottom_level_node* next;
  std::array<entry_base*, node_max_entries(1)> entries;
  
  void advance_resupply() {
    if (resupply_next_moved != this) {
      advance_resupply_impl(node_max_entries(1));
    }
  }
};
struct upper_level_node : public node_base {
  uint32_t num_children;
  // TODO: maybe children should be a std::array<bottom_level_node, node_max_entries(1)*3>*
  //                              (or std::array< upper_level_node, node_max_entries(1)*3>*)
  // so that the children are in the same memory region.
  std::array<node_base*, node_max_entries(1)*3> children;
  bottom_level_node* last_descendant;
  bottom_level_node* split_stopper;
  
  void advance_resupply() {
    idx_type threshold = node_entries_to_start_splitting_threshold(level)
    if (num_entries < threshold) {
      if (resupply_next_moved != last_descendant) {
        advance_resupply_impl(node_max_entries(level));
      }
    }
    else {
      if (num_entries == threshold) {
        split_stopper = current_best_split_stopper();
      }
      if (resupply_next_moved != split_stopper) {
        advance_resupply_impl(node_max_entries(level) / 2);
      }
    }
  }
};



















struct supply {
  entry_base* first_served_entry;
  entry_base* last_served_entry;
  supply* middle_served_supply;
  uint64_t steps_to_parent;
  uint64_t num_served_entries;
  uint64_t current_size;
  uint64_t amount_reserved;
  uint64_t amount_reserved_from_parent_supply_to_replace_resupplies;
  
  uint64_t remaining_capacity()const {
    return max_served_entries(level) - num_served_entries;
  }
  uint64_t entries_above_half()const {
    const uint64_t half = max_served_entries(level)>>1;
    if (num_served_entries < half) { return 0; }
    return num_served_entries - half;
  }
  
  int64_t resupply_steps_after_this(resupply const& r)const {
    return steps_to_parent - r.steps_to_parent;
  }
  int64_t resupply_max_steps_after_this(bool is_middle)const {
    return remaining_capacity()*max_resupply_steps_per_insert - (is_middle ? middle_served_supply->steps_to_parent : 0);
  }
  int64_t resupply_max_steps_to_parent()const {
    return entries_above_half()*max_resupply_steps_per_insert;
  }
  uint64_t resupply_reserve_needed()const {
    const uint64_t split_claim_size = supply_target_size(level)*4;
    const uint64_t leeway = remaining_capacity()*max_reserved_per_insert;
    if (split_claim_size > leeway) {
      return split_claim_size - leeway;
    }
    return 0;
  }
  void validate_resupply(bool is_middle)const {
    resupply const& r = is_middle ? middle_resupply : end_resupply;
    assert (resupply_steps_after_this(r) <= resupply_max_steps_after_this(is_middle));
    assert (r.steps_to_parent <= resupply_max_steps_to_parent());
    assert (r.amount_reserved_from_parent_supply_for_replacement >= resupply_reserve_needed());
  }
  
  void validate()const {
    assert (num_served_entries <= max_served_entries(level));
    assert (served_supplies.size() >= 2);
    assert (served_supplies.size() <= max_served_supplies);
    assert (served_supplies.front().first_served_entry == first_served_entry);
    assert (served_supplies.back().last_served_entry == last_served_entry);
    uint64_t confirm_amount_reserved = 0;
    uint64_t confirm_num_served_entries = 0;
    for (auto i : served_supplies) {
      assert (i.level == level - 1);
      confirm_amount_reserved += i.amount_reserved_from_parent_supply_to_replace_resupplies;
      confirm_num_served_entries += i.num_served_entries;
    }
    assert (confirm_num_served_entries == num_served_entries);
    assert (confirm_amount_reserved == amount_reserved);
    assert (current_size >= amount_reserved);
    if (last_served_entry->next()) {
      const uint64_t observed_gap_size = last_served_entry->next()->idx - last_served_entry->idx;
      assert (observed_gap_size > current_size);
    }
    validate_resupply(true);
    validate_resupply(false);
  }
  
  void step_resupply_upward(resupply& r) {
    r.lower_end->idx -= supply_target_size(level);
    assert (r.lower_end->idx > r.lower_end->prev()->idx);
    r.lower_end = r.lower_end.next();
    --r.steps_to_parent;
  }
  void step_resupply_downward(resupply& r) {
    r.lower_end->idx += supply_target_size(level);
    assert (r.lower_end->idx < r.lower_end->next()->idx);
    r.lower_end = r.lower_end.prev();
    ++r.steps_to_parent;
  }
  void step_resupply_as_needed(resupply& r, int64_t max_steps, bool is_middle) {
    int64_t forced_steps_upward = r.steps_to_parent - resupply_max_steps_to_parent();
    if (forced_steps_upward > 0) {
      assert (forced_steps_upward <= max_steps);
      while (forced_steps_upward > 0) {
        step_resupply_upward(r);
        --forced_steps_upward;
      }
    }
    else if (forced_steps_upward < 0) {
      // middle_resupply moves downward impatiently, to generate the leeway for
      // its destination to jump down without stretching the max speed.
      // end_resupply only moves as needed.
      const int64_t forced_steps_downward = resupply_steps_after_this(r) - resupply_max_steps_after_this(is_middle);
      assert (forced_steps_downward <= max_steps);
      assert (forced_steps_downward + forced_steps_upward <= 0);
      int64_t steps_downward = forced_steps_downward;
      if (is_middle) {
        const int64_t desired_steps_downward = r.steps_to_parent - middle_served_supply->steps_to_parent;
        steps_downward = std::min(std::min(desired_steps_downward, -forced_steps_upward), max_steps);
        assert (steps_downward >= forced_steps_downward);
      }
      while (steps_downward > 0) {
        step_resupply_once(r);
        --steps_downward;
      }
    }
  }
  void reserve_for_resupplies_as_needed() {
    resupply& r = is_middle ? middle_resupply : end_resupply;
    const int64_t reserve_needed = resupply_reserve_needed();
    const int64_t change = reserve_needed - int64_t(amount_reserved_from_parent_supply_to_replace_resupplies);
    assert (change <=  max_reserved_per_insert);
    assert (change >= -max_reserved_per_insert);
    amount_reserved_from_parent_supply_to_replace_resupplies = change;
    parent->amount_reserved += change;
  }
  void update_resupply_based_on_insertion(entry_base* inserted_entry, entry_base* existing_entry, bool after, bool is_middle) {
    resupply& r = is_middle ? middle_resupply : end_resupply;
    if (existing_entry->idx > r.lower_end->idx) {
      ++r.steps_to_parent;
    }
    if (existing_entry->idx == r.lower_end->idx && after) {
      r.lower_end = inserted_entry;
    }
    step_resupply_as_necessary(r, 1, is_middle);
  }
  void update_resupply_based_on_erasure(entry_base* erased_entry, bool is_middle) {
    resupply& r = is_middle ? middle_resupply : end_resupply;
    if (erased_entry->idx > r.lower_end->idx) {
      --r.steps_to_parent;
    }
    if (erased_entry->idx == r.lower_end->idx && after) {
      r.lower_end = erased_entry.prev();
    }
    step_resupply_as_necessary(r, 1, is_middle);
  }
  void split() {
    assert (num_served_entries == max_served_entries());
    assert (amount_reserved_from_parent_supply_to_replace_resupplies == supply_target_size(level)*4);
    
    supply* s1 = new supply();
      
    s1->last_served_entry = s0->last_served_entry;
    s0->last_served_entry = s0->middle_resupply.lower_end;
    s1->first_served_entry = s0->last_served_entry->next();
    s1->num_served_entries = s0->middle_served_supply->steps_to_parent;
    s0->num_served_entries -= s1->num_served_entries;
      
    parent->amount_reserved -= supply_target_size(s0->level)*4;
    parent->current_size -= supply_target_size(s0->level)*4;
    s0->middle_resupply.lower_end = parent->last_entry;
    s0->middle_resupply.steps_to_parent = 0;
    s0->   end_resupply.lower_end = parent->last_entry;
    s0->   end_resupply.steps_to_parent = 0;
    s1->middle_resupply.lower_end = parent->last_entry;
    s1->middle_resupply.steps_to_parent = 0;
    s1->   end_resupply.lower_end = parent->last_entry;
    s1->   end_resupply.steps_to_parent = 0;
    
    s0->current_size = supply_target_size(s0->level);
    s1->current_size = supply_target_size(s0->level);
    s0->amount_reserved = 0;
    s1->amount_reserved = 0;
    for (auto i : served_supplies) {
      if (i.last_served_entry->idx < s1->first_served_entry->idx) {
        s0->amount_reserved += i.amount_reserved_from_parent_supply_to_replace_resupplies;
      }
      else {
        s1->amount_reserved += i.amount_reserved_from_parent_supply_to_replace_resupplies;
      }
    }
    s0->amount_reserved_from_parent_supply_to_replace_resupplies = 0;
    s1->amount_reserved_from_parent_supply_to_replace_resupplies = 0;
    s0->update_middle_served_supply();
    s1->update_middle_served_supply();
  }
  void join() {
    assert (num_served_entries + old_sibling->num_served_entries == max_served_entries(level));
    assert (s0->amount_reserved_from_parent_supply_to_replace_resupplies == 0);
    assert (s1->amount_reserved_from_parent_supply_to_replace_resupplies == 0);
    assert (s0->middle_resupply.lower_end == parent->last_entry);
    assert (s0->middle_resupply.steps_to_parent == 0);
    assert (s0->   end_resupply.lower_end == parent->last_entry);
    assert (s0->   end_resupply.steps_to_parent == 0);
    assert (s1->middle_resupply.lower_end == parent->last_entry);
    assert (s1->middle_resupply.steps_to_parent == 0);
    assert (s1->   end_resupply.lower_end == parent->last_entry);
    assert (s1->   end_resupply.steps_to_parent == 0);
    assert (s0->current_size == supply_target_size(s0->level));
    assert (s1->current_size == supply_target_size(s0->level));
    
    parent->amount_reserved += supply_target_size(s0->level)*4;
    parent->current_size += supply_target_size(s0->level)*4;
    s0->amount_reserved_from_parent_supply_to_replace_resupplies = supply_target_size(s0->level)*4;
    
    s0->amount_reserved += s1->amount_reserved;
    s0->num_served_entries = max_served_entries(s0->level);
    s0->middle_resupply.lower_end = s0->last_served_entry;
    s0->middle_resupply.steps_to_parent = s0->steps_to_parent;
    s0->   end_resupply.lower_end = s1->last_served_entry;
    s0->   end_resupply.steps_to_parent = s1->steps_to_parent;
    s0->last_served_entry = s1->last_served_entry;
    s0->current_size = 0;
    s0->update_middle_served_supply();
  }
  


struct supply {
  entry_base* first_served_entry;
  entry_base* last_served_entry;
  supply* middle_served_supply;
  uint64_t num_served_entries;
  uint64_t supply_current_size;
  uint64_t supply_free_size;
  supply* resupply_source;
  uint64_t num_entries_between_this_and_resupply_source;
  
  struct resupply {
    entry_base* lower_end;
    int64_t steps_before_parent_supply;
    int64_t steps_after_this_node_supply;
  };
  resupply end_resupply;
  resupply middle_resupply;
  
  
  uint64_t max_served_entries()const {
    return ;
  }
  uint64_t remaining_capacity()const {
    return max_served_entries() - num_served_entries;
  }
  int64_t resupply_max_steps_to_dst()const {
    return remaining_capacity()*max_resupply_steps_per_insert;
  }
  void validate()const {
    assert (entries < max_entries());
    if (next_sibling) {
      assert (next_sibling->prev_sibling == this);
      assert (next_sibling->first_entry == last_entry->next());
    }
    if (prev_sibling) {
      assert (prev_sibling->next_sibling == this);
      assert (prev_sibling->last_entry == first_entry->prev());
    }
    if (resupply_source_node) {
      assert (resupply_source_node->last_entry->idx > last_entry->idx);
      
    }
    if (last_entry->next()) {
      const uint64_t observed_gap_size = last_entry->next()->idx - last_entry->idx;
      assert (observed_gap_size >= supply_current_size);
    }
    assert (supply_current_size >= supply_free_size);
    
    if (is_bottom_level()) {
      assert (supply_free_size >= remaining_capacity());
    }
    else {
      assert (supply_free_size >= remaining_capacity()*2);
    }
    assert (
      middle_resupply.steps_after_this_node_supply + middle_resupply.steps_before_parent_supply == 
         end_resupply.steps_after_this_node_supply +    end_resupply.steps_before_parent_supply
    );
    assert (middle_resupply.steps_before_parent_supply >= 0);
    assert (   end_resupply.steps_before_parent_supply >= 0);
    assert (middle_resupply.steps_after_this_node_supply <= resupply_max_steps_after_this_node_supply( true));
    assert (   end_resupply.steps_after_this_node_supply <= resupply_max_steps_after_this_node_supply(false));
    int64_t total_close_siblings = 1;
    for (node* sib = next_sibling; sib && sib->resupply_source_node == resupply_source_node; sib = sib->next_sibling) {
      ++total_close_siblings;
    }
    for (node* sib = prev_sibling; sib && sib->resupply_source_node == resupply_source_node; sib = sib->prev_sibling) {
      ++total_close_siblings;
    }
    assert (total_close_siblings <= max_close_siblings);
  }
  void step_resupply_once(resupply& r) {
    r.lower_end->idx += resupply_size();
    assert (r.lower_end->idx < r.lower_end->next()->idx);
    r.lower_end = r.lower_end.prev();
    --r.steps_after_this_node_supply;
    ++r.steps_before_parent_supply;
  }
  void unstep_resupply_once(resupply& r) {
    r.lower_end->idx -= resupply_size();
    assert (r.lower_end->idx > r.lower_end->prev()->idx);
    r.lower_end = r.lower_end.next();
    ++r.steps_after_this_node_supply;
    --r.steps_before_parent_supply;
  }
  void step_resupply_as_necessary(resupply& r, int64_t max_steps, bool is_middle) {
    if (is_middle) {
      int64_t unsteps = -int64_t(entries >> 1) - r.steps_after_this_node_supply;
      if (unsteps > 0) {
        if (unsteps > max_steps) { unsteps = max_steps; }
        while (unsteps > 0) {
          unstep_resupply_once(r);
          --unsteps;
        }
        return;
      }
    }
    int64_t steps = r.steps_after_this_node_supply - resupply_max_steps_after_this_node_supply(is_middle);
    assert (steps <= max_steps);
    while (steps > 0) {
      step_resupply_once(r);
      --steps;
    }
  }
  void update_resupply_based_on_insertion(resupply& r, entry_base* inserted_entry, entry_base* existing_entry, bool after, bool is_middle) {
    if (last_entry->idx < existing_entry->idx && existing_entry->idx <= r.lower_end->idx) {
      ++r.steps_after_this_node_supply;
    }
    if (last_entry->idx >= existing_entry->idx && existing_entry->idx > r.lower_end->idx) {
      --r.steps_after_this_node_supply;
    }
    if (existing_entry->idx > r.lower_end->idx) {
      ++r.steps_before_parent_supply;
    }
    if (existing_entry->idx == r.lower_end->idx && after) {
      r.lower_end = inserted_entry;
    }
    step_resupply_as_necessary(r, 1, is_middle);
  }
  void update_resupply_based_on_erasure(resupply& r, entry_base* erased_entry, bool is_middle) {
    if (last_entry->idx < erased_entry->idx && erased_entry->idx <= r.lower_end->idx) {
      --r.steps_after_this_node_supply;
    }
    if (last_entry->idx >= erased_entry->idx && erased_entry->idx > r.lower_end->idx) {
      ++r.steps_after_this_node_supply;
    }
    if (erased_entry->idx > r.lower_end->idx) {
      --r.steps_before_parent_supply;
    }
    if (erased_entry->idx == r.lower_end->idx && after) {
      r.lower_end = erased_entry.prev();
    }
    step_resupply_as_necessary(r, 1, is_middle);
  }
  bool was_inserted(entry_base* inserted_entry, entry_base* existing_entry, bool after) {
    validate();
    assert (existing_entry->idx >= first_entry->idx);
    assert (existing_entry->idx <= last_entry->idx);
    
    if (inserted_entry->idx > last_entry->idx) {
      last_entry = inserted_entry;
    }
    ++entries;
    supply_free_size -= is_bottom_level() ? 1 : 2;
    
    assert (entries <= max_entries());
    update_resupply_based_on_insertion(middle_resupply, inserted_entry, existing_entry, after,  true);
    update_resupply_based_on_insertion(   end_resupply, inserted_entry, existing_entry, after, false);
    for (node* sib = next_sibling; sib && sib->resupply_source_node == resupply_source_node; sib = sib->next_sibling) {
      update_resupply_based_on_insertion(sib->middle_resupply, inserted_entry, existing_entry, after,  true);
      update_resupply_based_on_insertion(sib->   end_resupply, inserted_entry, existing_entry, after, false);
    }
    for (node* sib = prev_sibling; sib && sib->resupply_source_node == resupply_source_node; sib = sib->prev_sibling) {
      update_resupply_based_on_insertion(sib->middle_resupply, inserted_entry, existing_entry, after,  true);
      update_resupply_based_on_insertion(sib->   end_resupply, inserted_entry, existing_entry, after, false);
    }
    step_resupply_as_necessary(middle_resupply, max_resupply_steps_per_insert,  true);
    step_resupply_as_necessary(   end_resupply, max_resupply_steps_per_insert, false);
    
    node* bigger_node = (resupply_source_node->first_entry->idx <= existing_entry->idx) ? resupply_source_node : resupply_source_node->prev_sibling;
    bool parent_split = bigger_node->was_inserted(inserted_entry, existing_entry, after);
    if (parent_split) {
      for (node* sib = next_sibling; sib && sib->resupply_source_node == bigger_node; sib = sib->next_sibling) {
        if (sib->last_entry->idx < bigger_node->prev_sibling->last_entry->idx) {
          sib->resupply_source_node = bigger_node->prev_sibling;
        }
      }
      for (node* sib = prev_sibling; sib && (sib->resupply_source_node == resupply_source_node || sib->resupply_source_node == bigger_node); sib = sib->prev_sibling) {
        if (sib->resupply_source_node == bigger_node && sib->last_entry->idx < bigger_node->prev_sibling->last_entry->idx) {
          sib->resupply_source_node = bigger_node->prev_sibling;
        }
      }
      if (resupply_source_node == bigger_node && last_entry->idx < bigger_node->prev_sibling->last_entry->idx) {
        resupply_source_node = bigger_node->prev_sibling;
      }
    }
    if (entries == max_entries) {
      assert (end_resupply.steps_after_this_node_supply == 0);
      node* new_sibling = new node();
      
      new_sibling->first_entry = first_entry;
      new_sibling->last_entry = middle_resupply.lower_end;
      first_entry = middle_resupply.lower_end->next();
      if (new_sibling->last_entry->idx < resupply_source_node->prev_sibling->last_entry->idx) {
        new_sibling->resupply_source_node = resupply_source_node->prev_sibling;
      }
      else {
        new_sibling->resupply_source_node = resupply_source_node;
      }
      new_sibling->next_sibling = next_sibling;
      next_sibling->prev_sibling = new_sibling;
      next_sibling = new_sibling;
      new_sibling->prev_sibling = this;
      
      new_sibling->middle_resupply.lower_end = new_sibling->resupply_source_node->last_entry;
      new_sibling->   end_resupply.lower_end = new_sibling->resupply_source_node->last_entry;
      new_sibling->resupply_source_node->supply_current_size -= resupply_size();
      new_sibling->resupply_source_node->supply_current_size -= resupply_size();
      new_sibling->middle_resupply.steps_before_parent_supply =
      new_sibling->   end_resupply.steps_before_parent_supply = 0;
      new_sibling->middle_resupply.steps_after_this_node_supply =
      new_sibling->   end_resupply.steps_after_this_node_supply = ;
      
      middle_resupply.lower_end = resupply_source_node->last_entry;
         end_resupply.lower_end = resupply_source_node->last_entry;
      resupply_source_node->supply_current_size -= resupply_size();
      resupply_source_node->supply_current_size -= resupply_size();
      
      new_sibling->validate();
      validate();
      return true;
    }
    validate();
    return false;
  }
  bool will_be_erased(entry_base* inserted_entry, entry_base* erased_entry) {
    validate();
    assert (erased_entry->idx >= first_entry->idx);
    assert (erased_entry->idx <= last_entry->idx);
    
    if (erased_entry == last_entry) {
      last_entry = last_entry->prev();
    }
    --entries;
    supply_free_size -= is_bottom_level() ? 1 : 2;
    
    update_resupply_based_on_erasure(middle_resupply, erased_entry,  true);
    update_resupply_based_on_erasure(   end_resupply, erased_entry, false);
    for (node* sib = next_sibling; sib && sib->resupply_source_node == resupply_source_node; sib = sib->next_sibling) {
      update_resupply_based_on_erasure(sib->middle_resupply, erased_entry,  true);
      update_resupply_based_on_erasure(sib->   end_resupply, erased_entry, false);
    }
    for (node* sib = prev_sibling; sib && sib->resupply_source_node == resupply_source_node; sib = sib->prev_sibling) {
      update_resupply_based_on_erasure(sib->middle_resupply, erased_entry,  true);
      update_resupply_based_on_erasure(sib->   end_resupply, erased_entry, false);
    }
    
    validate();
  }
};




  
  

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
