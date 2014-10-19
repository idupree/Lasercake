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


const num_bits_type level_size_shift = 4;
const uint64_t level_size = 1 << level_size_shift;
const num_bits_type supply_size_shift = 1;
const uint64_t supply_size = 1 << supply_size_shift;
const uint32_t supplies_per_level = (level_size-supply_size) / supply_size;
// In the time it takes us to use up a supply, move far enough to get all the way from its refill source
// even if it's the first supply among the supplies_per_level supplies drawing from the same source.
// max_gap_speed can be lenient, since it only affects the worst-case time and not the amortized time.
const uint32_t max_gap_speed = divide(supplies_per_level * level_size, supply_size, rounding_strategy<round_up, negative_is_forbidden>);
constexpr uint64_t entries_beneath_supply(num_bits_type level) {
  if (level == 0) {
    return 1;
  }
  return entries_beneath_supply(level-1) * supplies_per_level;
}
constexpr num_bits_type calc_max_level() {
  for (num_bits_type level = 0; ; ++level) {
    if (entries_beneath_supply(level) >= uint64_t(-1) / supplies_per_level) {
      return level;
    }
  }
}
const num_bits_type max_level = calc_max_level();
const num_bits_type entries_beneath_max_level = entries_beneath_supply(max_level);

struct supply {
  uint64_t max_size;
  uint64_t current_size;
  uint64_t reserved;
  int64_t next_refill_progress_emptying_self;
  supply* refill_source;
  entry* lower_bound;
  entry* next_refill_lower_bound;
  
  void claim(uint64_t amount) {
    assert (reserved >= amount);
    reserved -= amount;
    current_size -= amount;
  }
  void reserve_one() {
    if (!refill_source) {
      // We're the near-infinite supply at the right. If we run out, it's over.
      ++reserved;
      caller_correct_if(current_size >= reserved, "manual_orderer overflowed");
      
      uint64_t entries = current_size-reserved;
      uint64_t entries_this_level = entries_beneath_max_level;
      for (int i = max_level; i >= 0; --i) {
        while (entries >= entries_this_level) {
          entries -= entries_this_level;
        }
        if (entries == 0) {
          refill_source = new supply();
          refill_source->lower_bound = lower_bound;
          refill_source->max_size = max_size;
          refill_source->current_size = current_size - (1ULL<<i);
          refill_source->reserved = reserved;
          max_size = (1ULL<<i);
          current_size = max_size;
          reserved = 0;
        }
        entries_this_level = entries_this_level / supplies_per_level;
      }
    }
    else {
      refill_source->reserve_one();
      ++reserved;
      if (current_size < reserved) {
        assert (next_refill_lower_bound == lower_bound);
        assert (next_refill_progress_emptying_self == -1);
        refill_source->claim(max_size);
        current_size += max_size;
        next_refill_lower_bound = refill_source->lower_bound;
        next_refill_progress_emptying_self = max_size;
      }
      for (uint32_t i = 0; (i < max_gap_speed) && (next_refill_progress_emptying_self >= 0); ++i) {
        if (next_refill_lower_bound == lower_bound) {
          lower_bound = lower_bound->prev;
          int64_t new_progress_emptying_self = lower_bound->idx & (max_size - 1);
          assert (new_progress_emptying_self != next_refill_progress_emptying_self);
          if (new_progress_emptying_self > next_refill_progress_emptying_self) {
            next_refill_progress_emptying_self = -1;
          }
          else {
            next_refill_progress_emptying_self = new_progress_emptying_self;
          }
        }
        next_refill_lower_bound->idx += max_size;
        next_refill_lower_bound = next_refill_lower_bound->prev;
      }
      assert (current_size >= reserved);
    }
  }
};


num_bits_type bits_per_level = 4;
num_children_t gap_size_in_children = 2;
num_children_t hypothetical_children_per_node = (1<<bits_per_level)
num_children_t max_children_per_node = hypothetical_children_per_node - gap_size_in_children;
num_children_t node_last_child = max_children_per_node-1;
// In the time it takes us to fill the front gap, move the back gap all the way to the front.
// So, in gap_size_in_children*child_size, move by (hypothetical_children_per_node-gap_size_in_children)*child_size.
 max_gap_speed = ((hypothetical_children_per_node+gap_size_in_children-1)/gap_size_in_children) - 1;
struct node {
  leaf_ref insert_spot;
  leaf_ref gap_move_spot;
  leaf_ref shed_spot;
  uint64_t num_descendants;
  uint64_t first_hypothetical_idx;
  
  bool 
  node* parent;
  num_children_t which_child_of_parent;
  std::array<std::unique_ptr<node> or entry_ref, max_children_per_node> children;
  num_children_t last_empty_child
  // Doesn't go into earlier nodes. Returns a null ref for the first leaf in this node.
  leaf_ref prev_leaf() {
  }
  void insert_first (input) {
    node* node_to_insert_into = this;
    while (above bottom level) {
      num_children_t child_to_insert_into = 0;
      for (num_children_t i = 0; i < max_children_per_node; ++i) {
        if (node_to_insert_into->children[i+1] && node_to_insert_into->children[i+1]->saturated()) {
          break;
        }
        child_to_insert_into = i;
      }
      if (!node_to_insert_into->children[child_to_insert_into]) {
        node_to_insert_into->children[child_to_insert_into].reset(new node(node_to_insert_into));
      }
      node_to_insert_into = node_to_insert_into->children[child_to_insert_into].get();
    }
    for (num_children_t i = 0; i < max_children_per_node; ++i) {
      if (node_to_insert_into->children[i+1]) {
        node_to_insert_into->children[i] = input;
        while (node_to_insert_into) {
          ++node_to_insert_into.num_descendants;
          node_to_insert_into = node_to_insert_into->parent;
        }
        break;
      }
    }
    if (saturated()) {
      for (int i = 0; (i < max_gap_speed) && (gap_move_spot  first_hypothetical_idx); ++i) {
        gap_move_spot->idx += this->gap_length();
        while (true) {
          if (insert_spot.child > 0) {
            --gap_move_spot.child;
          }
          else {
            gap_move_spot.n->first_hypothetical_idx += this->gap_length();
            if (insert_spot.n == this) {
              nowhere_left_to_insert = true;
              break;
            }
            insert_spot.child = insert_spot.n->which_child_of_parent();
            insert_spot.n = insert_spot.n->parent;
          }
          else { break; }
        }
        leaf_ref last = prev_leaf(gap_move_spot);
        if (last) {
          gap_move_spot = last;
        }
        else {
          gap_move_spot = insert_spot;
        }
      }
      if (nowhere_left_to_insert) {
        assert (gap_move_spot == insert_spot)
        for (num_children_t i = node_last_child; i >= max_children_per_node-gap_size_in_children; --i) {
          assert (children[i]->empty());
        }
        for (num_children_t i = node_last_child; i >= gap_size_in_children; --i) {
          children[i].swap(children[i-1]);
        }
        gap_move_spot = prev_leaf(shed_spot);
        insert_spot = end of new empty child, children[gap_size_in_children-1];
      }
      insert_spot.n->children[insert_spot.child] = input;
      shed_last();
    }
    
    
    
    bool nowhere_left_to_insert = false;
    assert(insert_spot.n->children[insert_spot.child]);
    while (true) {
      --insert_spot.child;
      if (insert_spot.child < gap_size_in_children) {
        insert_spot.n->saturated = true;
        if (insert_spot.n == this) {
          nowhere_left_to_insert = true;
          break;
        }
        insert_spot.child = insert_spot.n->which_child_of_parent();
        insert_spot.n = insert_spot.n->parent;
      }
      else { break; }
    }
    if (!nowhere_left_to_insert) {
      assert (!insert_spot.n->children[insert_spot.child]);
      while (above bottom level) {
        insert_spot.n->children[insert_spot.child].reset(new node(insert_spot.n));
        insert_spot.n = insert_spot.n->children[insert_spot.child].get();
        insert_spot.child = node_last_child;
      }
    }
    if (!saturated) { assert(!nowhere_left_to_insert); }
    if (saturated) {
      for (int i = 0; (i < max_gap_speed) && (gap_move_spot > insert_spot); ++i) {
        gap_move_spot->idx += this->gap_length();
        leaf_ref last = prev_leaf(gap_move_spot);
        if (last) {
          gap_move_spot = last;
        }
        else {
          gap_move_spot = insert_spot;
        }
      }
      if (nowhere_left_to_insert) {
        assert (gap_move_spot == insert_spot)
        for (num_children_t i = node_last_child; i >= max_children_per_node-gap_size_in_children; --i) {
          assert (children[i]->empty());
        }
        for (num_children_t i = node_last_child; i >= gap_size_in_children; --i) {
          children[i].swap(children[i-1]);
        }
        gap_move_spot = prev_leaf(shed_spot);
        insert_spot = end of new empty child, children[gap_size_in_children-1];
      }
      insert_spot.n->children[insert_spot.child] = input;
      shed_last();
    }
  }
  void shed(input) {
    if (input > parent->shed_spot) {
      parent->shed(input);
    }
    else {
      parent->children[which_child_of_parent()+1].insert(input);
    }
  }
  void update_shed_spot() {
    if (!shed_spot.n->children[shed_spot.child]) {
      // The only way to get rid of it is if a 
    }
  }
  void shed_last() {
    entry_ref input = shed_spot.n->children[shed_spot.child];
    shed_spot.n->children[shed_spot.child] = nullptr;
    while (shed_spot.n->empty()) {
      assert (shed_spot.n != this);
      shed_spot.child = shed_spot.n->which_child_of_parent();
      shed_spot.n = shed_spot.n->parent;
      shed_spot.n->children[shed_spot.child] = nullptr;
    }
    assert (shed_spot.child > 0);
    --shed_spot.child;
    shed(input);
  }
};


#endif
