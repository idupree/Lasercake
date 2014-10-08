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
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <inttypes.h>
#include <assert.h>
#include <limits>
#include <iostream>
#include <boost/next_prior.hpp>
#include <memory>

/*
 
 I expect that this explanatory comment can be much improved.

 A time_steward sits outside the flow of time, looking down upon it.
 It wields impossible power over a universe of Entities.
 Its power comes in the form of FiatEvents, which are an
 event at a point in time imposed upon the world, such as, in a
 platform game, the player character jumping.

 There are also ConsequentialEvents, such as when that
 jumping player character smacks their head into a ceiling and
 thus stops moving upwards.

 In between those two Events, there is a state of
 continuous change.  It has no effect other than helping
 the (timesteward??) determine when the next ConsequentialEvent
 is, and letting the client determine a state of the world
 at a time in between events.



 Nothing else changes.

 


in a rocket
 game, the player character turning on their engines.

 
 
 A single time_steward manages one universe of Entities and their changes over time,
   1) such that the changes are discrete changes and occur at specific instants,
   2) such that, for a given time and Entity, the time and nature of its upcoming changes can be determined solely
        from the states of a small finite number of other entities,
   3) for the purpose of making a deterministic simulation with efficient access and mutation of the world state at arbitrary times.
 
 (3) allows various things:
   it can be used to smooth over networked games' lag (continue simulating locally assuming no network input,
     then go back and add the network input at the correct time as it arrives, and recompute only what is actually affected by that input).
   it can be used to smooth over momentary complicated situations (continuously simulate several seconds ahead,
     and add the user's input as it comes, similar to the above), and allows for parallelism (do multiple changes and overwrite any that happen to conflict).
   it allows us to give the user an interface into the world where ze can look back in time to any moment in the history,
     and perhaps try changing something, creating a new branch.
     In this manner, we can create a sandbox game with zero risk from momentary mistakes, since the user can scroll back in time whenever something goes wrong.
   The time_steward also has the advantage that neither Entities nor durations of time consume computer resources by themselves:
     only the Events need processing. Hence, any Entity that is seeing little action automatically yields processor time
     to other parts of the simulation - and if all Entities are fairly inactive, the time_steward can quickly pass a very large amount of time
     (as opposed to, for instance, having simulation steps that each require a minimal update).
     
 The rest of this comment essentially describes the technical meanings of (1) and (2) - 
   the restrictions we make on the universe in order to make these features possible.
 In this comment, we will discuss two types of things: Data, which exists in the program in a literal sense;
   and Facts, in the sense of mathematical facts, which are true but are not necessarily represented in the same way they are stated.
 The existence of "Entities" and "Events" is a mere Fact; when discussing Data, we will more specifically refer to EntityData, EntityIDs, and such.
     
 From (1), we have the Fact of an "EntityThroughoutTime", which is a partitioning of the time line into intervals, plus
   a mapping from each interval to a (constant throughout that interval) EntityData.
   This EntityData may describe a trajectory, so that the Entity represents a continuously moving object.
   
 (2) is never a hard limit; the time_steward will still function even if we must examine *every* entity to determine *any* entity's next change.
   Localizing the change function is simply a way to derive performance advantages.
 
 As a loosely described example, imagine a space full of moving circles that bounce off of each other.
   The Entities are the circles, and their EntityData is a radius and the coefficients of an affine function from time to circle center location.
   The Events are the moments of collision,
   and the time and nature of a Event can be derived from the EntityData of the nearby circles.
     (How we determine which circles to call 'nearby' is a complex matter better discussed elsewhere.)
   Suppose we have a group of circles that are all 'nearby' each other.
   For each pair, if those two are on a collision course for each other, it is a Fact that they will collide at a specific time if they remain on those courses.
   For various reasons, we choose to require the client code to represent that Fact as another Entity, which will exist in that form
     from the moment the two circles have assumed their current trajectories until the moment one of them diverts.
     To wit:

 Some Entities are persistent; others are temporary EventAnticipationEntities whose EntityData is simply a time and a description of Event.
   If a EventAnticipationEntity continues to exist until its time, then the Event happens.

 Fact:
   A time_steward (with given parameters) is defined entirely by
     (A) a single EntityData that is the initial state of a single object (called the "global object"), and
     (B) a collection of Events imposed by fiat, their existences and data contingent on nothing within the universe.
   We call A and B together a "Specification".
   It is a Fact that a Specification defines the entire history/future of a universe. Hence, given any Specification, there are
     a fixed (probably infinite) collection of Events that will happen.
     Data-wise, we only ever know of a finite collection of Events that might happen, *some* of which we know to actually happen.
   There are exactly two kinds of Events: Those imposed by fiat, which we will call FiatEvents,
     and those specified by EventAnticipationEntities, which we will call ConsequentialEvents.
   Other than the initial state of the global object, every EntityData is created by exactly one Event.
 
 A time_steward's external interface:
   Construct/initialize time_steward with a given global object EntityData.
   Add or remove one of the FiatEvents.
   Query the EntityData of an entity with a specific EntityID infinitesimally before or after a specific time.
     (To avoid confusion about the state at boundary moments, we never discuss "state AT a time", only infinitesimally before or after that time.)
   A few functions to manage resource usage.
 
 In code structure, we have:
   "External client code", such as UI code, creates a time_steward, and calls the external interface functions described above,
     passing in EventData, which are function objects.
   To answer the external client code's queries, the time_steward may call the EventData.
     We will refer to the EventData of ConsequentialEvents as "internal client code",
     and to the EventData of FiatEvents as "liminal client code",
       because it is subject to some of the restrictions of both internal and external client code.
     (These two kinds of Events are the only kinds of Events that exist.)
 
 (The global object exists to have a fixed EntityID as a point of reference for the client code.
   External and liminal client code aren't allowed to remember any other EntityIDs,
   because if a change is made in the past, those Entities might never exist.
   Instead, that code must ask for EntityData corresponding to the global object EntityID,
   and that EntityData may contain other EntityIDs.)
 
 A EventData is a function object that receives a single argument: A time_steward Accessor object tied to a specific time.
   It then uses the Accessor
     to create Entities, for which the Accessor will provide new EntityIDs,
     and/or
     to modify some Entities, which, internally, means marking that their current EntityData ends at this time
       and the modified EntityData begins at this time. Since EntityData is always copied, any large Entity should either
       use persistent data structures or be broken up into multiple Entities.
   EventData can only base its behavior on information from the Accessor and its own fields (information it was constructed with).
   The Accessor allows it to query the EntityData of any EntityID infinitesimally-before its time,
     and it tracks which EntityIDs are so queried, so that, even if external client code changes the past, the time_steward might not
     need to call this EventData again - if the change in the past doesn't affect any of the queried EntityData,
     then the time_steward knows that this EventData will give the same result as before, so the time_steward
     can keep using that result.
   A FiatEvent, as mentioned above, cannot know any EntityIDs except the global object ID, except through the Accessor.
   A ConsequentialEvent may know any EntityIDs that were known by the earlier Event that created this one's EventAnticipationEntity.
     (Those EntityIDs were meaningful at the time, and if a past change undoes the creation of one of them,
     then the earlier Event would never have known it either. If internal or liminal code changes an Entity in a way
     that would make this Event consider it invalid, that code is responsible for cleaning up the references to that Entity.
     Events cannot destroy an Entity per se, but they can replace it with a null EntityData and remove all references to it.)
 
 How does the time_steward work internally?
   If the time_steward knows (Fact-wise) all the EntityData at a given time A, then it knows all the EventAnticipationEntities.
     Thus it knows both the corresponding ConsequentialEvents and all the FiatEvents after A.
   Among those, one comes first, at a time B.
   No other Event exists that could change the Entities from (their state at A) before B.
   So the time_steward effectively "knows" the EntityData at all times in (A, B).
   Then it can apply the Event at B to those EntityData, which yields the state immediately after B - from B to whenever the next change is.
   This is part of how the process works Data-wise, too. There is a moving time to which the time_steward is updated,
     which is defined to be the time of the first Event that the time_steward knows about but hasn't processed yet - hereafter ExposedEvents.
   Hence, when someone inserts a FiatEvent far in the past, that FiatEvent becomes the first ExposedEvent,
     and hence, the time_steward knows that it "doesn't know for sure" any EntityData after that point.
     
 The time_steward may also have some "provisional knowledge" - EntityData etc. that it's computed beyond the updated-to-time.
   This can happen either directly (the time_steward processes a Event that's not the first),
     or indirectly (a FiatEvent is added prior to the first ExposedEvent).
   Provisional knowledge may be invalidated by changes that come later code-wise but earlier in-world.
   To see how we do this, consider a digraph where the nodes are EntityDataIntervals and Events.
     A Event that creates a new EntityData (and thus begins an EntityDataInterval) flows to that EntityDataInterval.
     A Event that *references* an Entity flows *from* the EntityDataInterval into which it's referencing.
   Whenever an EntityData is changed at a certain time, we undo all changes downstream of it in the digraph.
   This isn't quite accurate - if an entity is changed in the middle of an interval, we only need to invalidate references that
     came *after* that middle time. But the basic idea is there.
   
 The above assumes that Event times are unique (and totally ordered).
   We choose to enforce that fact because it makes work simpler for the client code.
   In order to do this, we consider two kinds of times: BaseTimes, which client code gives us;
     and ExtendedTimes, which are a BaseTime plus a distinguisher we create to order them even if their BaseTime is the same.
   This isn't a trivial task, because the ordering we create must be the same regardless of the order we processed the Events
     that brought us to this situation.
   To accomplish this, we use hashing. Each time and each entity gets a 128-bit unique ID using siphash.
     The global object EntityID is 0.
     Any other EntityID is hash(TimeID of when it was created, index among Entities created at that time).
     The TimeID of a ConsequentialEvent is hash(its BaseTime, the EntityID of its EventAnticipationEntity).
     The TimeID of a          FiatEvent is hash(its BaseTime, a distinguisher provided by the external client code).
 
 Some possibilities for improvement:
 - Have a way to discard some of the old data to save memory (at the price of needing extra processing time if/when we go back.)
 - Make time_steward itself a persistent data structure.
 */


// hackly copied from /usr/include/c++/4.9.0/bits/unique_ptr.h
// so we don't need c++1y libs just c++11
namespace futurestd {
  using namespace std;
  template<typename _Tp>
    struct _MakeUniq
    { typedef unique_ptr<_Tp> __single_object; };

  template<typename _Tp>
    struct _MakeUniq<_Tp[]>
    { typedef unique_ptr<_Tp[]> __array; };

  template<typename _Tp, size_t _Bound>
    struct _MakeUniq<_Tp[_Bound]>
    { struct __invalid_type { }; };

  /// std::make_unique for single objects
  template<typename _Tp, typename... _Args>
    inline typename _MakeUniq<_Tp>::__single_object
    make_unique(_Args&&... __args)
    { return unique_ptr<_Tp>(new _Tp(std::forward<_Args>(__args)...)); }

  /// std::make_unique for arrays of unknown bound
  template<typename _Tp>
    inline typename _MakeUniq<_Tp>::__array
    make_unique(size_t __num)
    { return unique_ptr<_Tp>(new typename remove_extent<_Tp>::type[__num]()); }

  /// Disable std::make_unique for arrays of known bound
  template<typename _Tp, typename... _Args>
    inline typename _MakeUniq<_Tp>::__invalid_type
    make_unique(_Args&&...) = delete;
}
using futurestd::make_unique;



uint64_t siphash24(const void *src,
                   unsigned long src_sz,
                   const char key[16]);

// TODO: move all the data structures stuff before namespace time_steward_system to a better file
template<typename ...T> using persistent_map = std::map<T...>;
template<typename ...T> using persistent_set = std::set<T...>;

namespace ordered_stuff_system {
typedef uint64_t idx_type;
typedef uint32_t num_bits_type;
const idx_type no_idx = std::numeric_limits<idx_type>::max();

template<typename ValueType> class ordered_stuff;
namespace impl {
template<typename ValueType>
struct entry {
  template <class... Args>
  entry(idx_type idx, Args&&... args):idx(idx),ref_count(0),contents(std::forward<Args>(args)...){}
  idx_type idx;
  size_t ref_count;
  ValueType contents;
};
template<typename ValueType>
struct entry_ref {
public:
  entry_ref():data(nullptr){}
  bool operator< (entry_ref const& o)const { return data->idx <  o.data->idx; }
  bool operator> (entry_ref const& o)const { return data->idx >  o.data->idx; }
  bool operator<=(entry_ref const& o)const { return data->idx <= o.data->idx; }
  bool operator>=(entry_ref const& o)const { return data->idx >= o.data->idx; }
  bool operator==(entry_ref const& o)const { return data == o.data; }
  bool operator!=(entry_ref const& o)const { return data != o.data; }
  explicit operator bool()const { return bool(data); }
  ValueType& operator*()const { return data->contents; }
  ValueType* operator->()const { return &data->contents; }
  entry_ref& operator=(entry_ref const& o){ dec_ref(); data = o.data; if (data) { ++data->ref_count; } }
private:
  void dec_ref() {
    if (--data->ref_count == 0) {
      delete data;
    }
  }
  entry_ref(entry<ValueType>* data):data(data){ ++data->ref_count; }
  entry<ValueType>* data;
  friend class ordered_stuff<ValueType>;
  friend struct std::hash<entry_ref>;
};
}

//#define AUDIT_ORDERED_STUFF

template<typename ValueType>
class ordered_stuff {
  typedef impl::entry<ValueType> entry;
public:
  typedef impl::entry_ref<ValueType> entry_ref;
  template <class... Args>
  entry_ref construct(Args&&... args) {
    // TODO: don't leak memory (reference count?)
    // TODO: better allocator
    return entry_ref(new entry(no_idx, std::forward<Args>(args)...));
  }

  // TODO: allow moving and deletion
  void put_only(entry_ref moving) {
    place_at(0, moving);
  }
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
    assert(moving.data->idx == no_idx);
    moving.data->idx = idx;
    auto p = data.insert(std::make_pair(idx, moving));
    assert(p.second);
  }
};
} // namespace ordered_stuff_system
using ordered_stuff_system::ordered_stuff;
namespace std {
  template<typename ValueType>
  struct hash<typename ordered_stuff_system::impl::entry_ref<ValueType>> {
    public:
    size_t operator()(ordered_stuff_system::impl::entry_ref<ValueType> const& i)const {
      // Nondeterministic - but the ordering of STL hash tables is already nondeterministic.
      return std::hash<decltype(i.data)>()(i.data);
    }
  };
}

namespace time_steward_system {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
// 128 bits of random stuff (presumed unique to this library)
const char siphash_key[16] = { 0xb7, 0xac, 0x3d, 0xf8, 0xc3, 0xa2, 0x8c, 0xd9, 0x3a, 0x10, 0x91, 0x68, 0x09, 0x02, 0x74, 0x0e };
#pragma GCC diagnostic pop

// siphash_id takes anything that's hashable and turns it into a
// fixed-size, statistically almost certainly unique, equality-comparable
// and orderable data structure (i.e. 128 bits of data).  The ==
// is true iff the things are the same; the < is an arbitrary total ordering.
namespace impl {
  template<typename... Args> struct combining_helper;
}
namespace persistent_siphash_id_trie_system {
typedef uint32_t num_bits_type;
template<typename MappedType = void, num_bits_type bits_per_level = 4>
class persistent_siphash_id_trie;
}
class siphash_id {
  template<typename... Args> friend struct impl::combining_helper;
  template<typename MappedType, persistent_siphash_id_trie_system::num_bits_type bits_per_level>
  friend class persistent_siphash_id_trie_system::persistent_siphash_id_trie;
  typedef std::array<uint64_t,2> data_type;
public:
  // Implementation note:
  // A reply from the SipHash developers about our idea to get a
  // large enough result from SipHash to be statistically unique
  // by hashing each of N and N+1 and keeping both results:
  //
  // Hi Eli,
  //
  // the development of a 128-bit SipHash is under discussion, see e.g.
  // this Twitter thread:
  //https://twitter.com/Kryptoblog/status/397285808714960897
  //
  // The trick you use seems to work, but it'd perhaps be more elegant to
  // just hash the same info under two distinct keys, or (say) under key K
  // and key (K xor 0xffff..ffff). That's still twice as slow as a single
  // pass though.
  //
  // What you should not do obviously is "hash the hash" to get another 64
  // bits, as any collision would propagate.
  //
  // When/if a SipHash-128 is developed, we'll advertize it on Twitter.
  //
  // Best,
  //
  // JP
  
public:
  // Accepts any combination of siphash_ids and types castable to uint64_t.
  // Combines all the entropy from all the arguments.
  template<typename... Args>
  static siphash_id combining(Args... args);
  
  // null siphash
  constexpr siphash_id() : data_({{ 0, 0 }}) {}
  constexpr static siphash_id null() {
    return siphash_id();
  }
  // least and greatest nonnull siphash
  constexpr static siphash_id least() {
    return siphash_id(data_type{{ 1, 0 }});
  }
  constexpr static siphash_id greatest() {
    return siphash_id(data_type{{ uint64_t(0)-1, uint64_t(0)-1 }});
  }

  bool operator<(siphash_id const& o)const {
    if (data_[1] < o.data_[1]) return true;
    if (data_[1] > o.data_[1]) return false;
    return data_[0] < o.data_[0];
  }
  bool operator>(siphash_id const& o)const { return o < *this; }
  bool operator<=(siphash_id const& o)const { return !(o < *this); }
  bool operator>=(siphash_id const& o)const { return !(o > *this); }
  bool operator==(siphash_id const& o)const {
    return (data_[0] == o.data_[0]) && (data_[1] == o.data_[1]);
  }
  bool operator!=(siphash_id const& o)const {
    return (data_[0] != o.data_[0]) || (data_[1] != o.data_[1]);
  }
  explicit operator bool()const {
    return (*this) != null();
  }
  
  friend inline std::ostream& operator<<(std::ostream& os, siphash_id const& s) {
    os << s.data_[0] << "," << s.data_[1];
    return os;
  }
private:
  constexpr siphash_id(data_type data):data_(data){}
  data_type data_;
  friend struct std::hash<siphash_id>;
};
namespace impl {
  template<typename... Args> struct combining_helper;
  template<> struct combining_helper<> {
    static const size_t next_idx = 0;
    static inline void enter_arg(uint64_t[]){}
  };
  template<typename... Tail> struct combining_helper<siphash_id, Tail...> {
    static const size_t next_idx = combining_helper<Tail...>::next_idx + 2;
    static inline void enter_arg(uint64_t in[], siphash_id head, Tail... tail) {
      combining_helper<Tail...>::enter_arg(in, tail...);
      in[combining_helper<Tail...>::next_idx] = head.data_[0];
      in[combining_helper<Tail...>::next_idx+1] = head.data_[1];
    }
  };
  template<typename Head, typename... Tail> struct combining_helper<Head, Tail...> {
    static const size_t next_idx = combining_helper<Tail...>::next_idx + 1;
    static inline void enter_arg(uint64_t in[], Head head, Tail... tail) {
      combining_helper<Tail...>::enter_arg(in, tail...);
      in[combining_helper<Tail...>::next_idx] = uint64_t(head);
    }
  };
  template<typename Head, typename... Tail> struct combining_helper<Head&, Tail...> : public combining_helper<Head, Tail...> {};
}
template<typename... Args>
siphash_id siphash_id::combining(Args... args) {
  siphash_id result;
  uint64_t in[impl::combining_helper<Args...>::next_idx];
  impl::combining_helper<Args...>::enter_arg(in, args...);
  result.data_[0] = siphash24((void*)in, sizeof(in), siphash_key);
  ++in[0];
  result.data_[1] = siphash24((void*)in, sizeof(in), siphash_key);
  return result;
}

} //end namespace time_steward_system
namespace std {
  template<>
  struct hash<time_steward_system::siphash_id> {
    public:
    size_t operator()(time_steward_system::siphash_id const& i)const {
      // Since the siphash is already a hash, just return some of it.
      return i.data_[0];
    }
  };
}
namespace time_steward_system {
  

namespace persistent_siphash_id_trie_system {
template<typename KeyType, typename MappedType, bool IsMap> struct value_typer;
template<typename KeyType, typename MappedType> struct value_typer<KeyType, MappedType, false> {
  typedef KeyType type;
  static KeyType get_key(type const& t) { return t; }
};
template<typename KeyType, typename MappedType> struct value_typer<KeyType, MappedType, true> {
  typedef std::pair<KeyType, MappedType> type;
  static KeyType get_key(type const& t) { return t.first; }
};
template<typename MappedType, num_bits_type bits_per_level>
class persistent_siphash_id_trie {
public:
  static_assert((128 % bits_per_level) == 0, "I don't think everything works right if 128 isn't a multiple of bits_per_level");
  typedef siphash_id key_type;
  typedef MappedType mapped_type;
  typedef value_typer<key_type, mapped_type, !std::is_void<mapped_type>::value> value_typer_t;
  typedef typename value_typer_t::type value_type;
private:  
  typedef typename boost::uint_t<(1<<bits_per_level)>::fast one_bit_per_child_t;
  typedef size_t ref_count_t;
  static const ref_count_t is_leaf_bit = 1ULL<<63;
  static const ref_count_t ref_count_mask = is_leaf_bit-1;
  struct child_header {
    child_header(ref_count_t r):ref_count_and_is_leaf(r){}
    ref_count_t ref_count_and_is_leaf;
  };
  struct leaf_t : public child_header {
    template <class... Args>
    leaf_t(Args&&... args):child_header(is_leaf_bit),value(std::forward<Args>(args)...){}
    value_type value;
    key_type id() { return value_typer_t::get_key(value); }
  };
  struct node_header : public child_header {
    node_header():child_header(0),children_exist(0){}
    one_bit_per_child_t children_exist;
  };
  class child_ptr {
  public:
    child_ptr(){}
    child_ptr(decltype(nullptr)):data(nullptr){}
    explicit child_ptr(char* data):data(data){}
    inline bool is_leaf()const { return (((child_header*)(data))->ref_count_and_is_leaf) & is_leaf_bit; }
    leaf_t& leaf() { return (*(leaf_t*)(data)); }
    leaf_t const& leaf()const { return (*(leaf_t const*)(data)); }
    node_header& node() { return (*(node_header*)(data)); }
    node_header const& node()const { return (*(node_header const*)(data)); }
    child_ptr& child_by_idx(num_bits_type i)const { return (*(child_ptr*)(data + sizeof(node_header) + i*sizeof(child_ptr))); }
    child_ptr& child_by_bit(num_bits_type b)const { return child_by_idx(popcount(node().children_exist & (b-1))); }
    num_bits_type num_children()const { return popcount(node().children_exist); }
    explicit operator bool()const { return bool(data); }
    bool operator==(child_ptr o)const { return data==o.data; }
    bool operator!=(child_ptr o)const { return data!=o.data; }
    void inc_ref()const {
      if (data) {
        ++(((child_header*)(data))->ref_count_and_is_leaf);
      }
    }
    void dec_ref()const {
      if (data && ((--(((child_header*)(data))->ref_count_and_is_leaf)) & ref_count_mask) == 0) {
        if (is_leaf()) {
          leaf().~leaf_t();
        }
        else {
          const num_bits_type num_children_ = num_children();
          for (num_bits_type i = 0; i < num_children_; ++i) {
            child_by_idx(i).dec_ref();
            child_by_idx(i).~child_ptr();
          }
          node().~node_header();
        }
        delete[] data;
      }
    }
  private:
    char* data;
  };
  template <class... Args>
  static child_ptr allocate_leaf(Args&&... args) {
    char* ptr = new char[sizeof(leaf_t)];
    new(ptr) leaf_t(std::forward<Args>(args)...);
    return child_ptr(ptr);
  }
  static child_ptr allocate_node(num_bits_type num_children) {
    char* ptr = new char[sizeof(node_header) + num_children*sizeof(child_ptr)];
    new(ptr) node_header();
    for (num_bits_type i = 0; i < num_children; ++i) {
      new(ptr + sizeof(node_header) + i*sizeof(child_ptr)) child_ptr();
    }
    return child_ptr(ptr);
  }
  child_ptr root;
  
  static inline num_bits_type which_child(siphash_id id, num_bits_type which_bits) {
    return (id.data_[which_bits>=64]>>which_bits) & ((1<<bits_per_level)-1);
  }
  
  persistent_siphash_id_trie(child_ptr root):root(root){ root.inc_ref(); }
public:
  persistent_siphash_id_trie():root(nullptr){}
  // Not set equality. Notably, though, trie.erase (something not in the trie) == trie.
  bool operator==(persistent_siphash_id_trie o)const { return root==o.root; }
  bool operator!=(persistent_siphash_id_trie o)const { return root!=o.root; }
  bool empty()const { return !bool(root); }
  
  persistent_siphash_id_trie(persistent_siphash_id_trie const& o):root(o.root){ root.inc_ref(); }
  persistent_siphash_id_trie& operator=(persistent_siphash_id_trie const& o) {
    o.root.inc_ref();
    root.dec_ref();
    root = o.root;
    return *this;
  }
  ~persistent_siphash_id_trie() { root.dec_ref(); }
    
  class iterator : public boost::iterator_facade<iterator, const value_type, boost::bidirectional_traversal_tag> {
    public:
      //iterator() : path_len(0) {}
      explicit operator bool()const { return id && back(); }
    private:
      explicit iterator(child_ptr root):path_len(1) { path[0] = root; }
      friend class boost::iterator_core_access;
      friend class persistent_siphash_id_trie;
      
      // State: We can be
      // 1) an iterator to an element (path is the path ending in the leaf, id redundant (== this->first or *this))
      // 2) an iterator to an empty child (path is the path ending in nullptr, possibly with a leaf before it, id is needed to specify which child)
      // 3) the end iterator (path is the singleton root, id null)
      // --end is the last element, ++end is the first element
      // The empty set has one iterator which is its own next and prior.
      void increment() {
        if (!path[0]) return;
        if (id) {
          num_bits_type child = which_child(id, back_bits());
          
          while (!back() || back().is_leaf() || !(back().node().children_exist & ~((1<<(child+1))-1))) {
            if (back() && !back().is_leaf()) { assert(popcount(back().node().children_exist & ((1U<<child)-1)) + 1 == back().num_children()); }
            bool check = bool(back());
            if (path_len == 1) {
              id = siphash_id::null();
              return;
            }
            --path_len;
            child = which_child(id, back_bits());
            if (check) { assert(back().node().children_exist & (1<<child)); }
          }
          
          num_bits_type idx = popcount(back().node().children_exist & ((1U<<child)-1)) + 1;
          assert (idx < back().num_children());
          assert (idx >= 0);
          path[path_len] = back().child_by_idx(idx);
          ++path_len;
        }
        
        while (!back().is_leaf()) {
          path[path_len] = back().child_by_idx(0);
          ++path_len;
        }
        id = back().leaf().id();
      }
      void decrement() {
        if (!path[0]) return;
        if (id) {
          num_bits_type child = which_child(id, back_bits());
          
          while (!back() || back().is_leaf() || !(back().node().children_exist & ((1<<child)-1))) {
            if (back() && !back().is_leaf()) { assert(popcount(back().node().children_exist & ((1U<<child)-1)) == 0); }
            bool check = bool(back());
            if (path_len == 1) {
              id = siphash_id::null();
              return;
            }
            --path_len;
            child = which_child(id, back_bits());
            if (check) { assert(back().node().children_exist & (1<<child)); }
          }
          
          num_bits_type idx = popcount(back().node().children_exist & ((1U<<child)-1)) + 1;
          assert (idx < back().num_children());
          assert (idx >= 0);
          path[path_len] = back().child_by_idx();
          ++path_len;
        }
        while (!back().is_leaf()) {
          path[path_len] = back().child_by_idx(back().num_children()-1);
          ++path_len;
        }
        id = back().leaf().id();
      }
      bool equal(iterator const& other)const { return id == other.id; }
      value_type const& dereference()const { return back().leaf().value; }
      
      child_ptr& back() { return path[path_len-1]; }
      child_ptr const& back()const { return path[path_len-1]; }
      num_bits_type bits(num_bits_type path_idx)const { return 128 - (path_idx+1)*bits_per_level; }
      num_bits_type back_bits()const { return 128 - path_len*bits_per_level; }
      
      siphash_id id;
      num_bits_type path_len;
      std::array<child_ptr, 128/bits_per_level> path;
  };
  iterator begin()const { return ++iterator(root); }
  iterator end()const { return iterator(root); }
  
  template<typename Hack = void>
  /*std::enable_if_t<std::is_same<key_type, value_type>::value, persistent_siphash_id_trie>*/ persistent_siphash_id_trie insert(key_type k)const {
    return insert(find(k));
  }
  template<typename Hack = void>
  static /*std::enable_if_t<std::is_same<key_type, value_type>::value, persistent_siphash_id_trie>*/ persistent_siphash_id_trie insert(iterator const& i) {
    persistent_siphash_id_trie result = set_leaf(i, allocate_leaf(i.id));
    assert(result.find(i.id).back().leaf().id() == i.id);
    return result;
  }
  template<class... Args>
  persistent_siphash_id_trie emplace(key_type k, Args&&... args)const {
    return emplace(find(k), std::forward<Args>(args)...);
  }
  template<class... Args>
  static persistent_siphash_id_trie emplace(iterator const& i, Args&&... args) {
    persistent_siphash_id_trie result = set_leaf(i, allocate_leaf(i.id, std::forward<Args>(args)...));
    assert(result.find(i.id).back().leaf().id() == i.id);
    return result;
  }
  persistent_siphash_id_trie erase(key_type k)const {
    return erase(find(k));
  }
  static persistent_siphash_id_trie erase(iterator const& i) {
//     std::vector<child_ptr> old_leaves;
//     for (iterator j = persistent_siphash_id_trie(i.path[0]).begin(); j != persistent_siphash_id_trie(i.path[0]).end(); ++j) {
//       if (j.id != i.id) old_leaves.push_back(j.back());
//     }
    persistent_siphash_id_trie result = unset_leaf(i);
    assert(!result.find(i.id));
//     size_t idx = 0;
//     for (iterator j = result.begin(); j != result.end(); ++j) {
//       assert(old_leaves[idx++] == j.back());
//     }
    return result;
  }
  iterator find(key_type k)const {
    caller_correct_if(bool(k), "You can't find() for the null ID");
    iterator i(root);
    i.id = k;
    while (i.back()) {
      if (i.back().is_leaf()) {
        if (i.back().leaf().id() == k) {
          return i;
        }
        else {
          i.path[i.path_len++] = nullptr;
        }
      }
      else {
        const num_bits_type child_bit = 1<<which_child(i.id, i.back_bits());
        if (i.back().node().children_exist & child_bit) {
          i.path[i.path_len] = i.back().child_by_bit(child_bit);
          ++i.path_len;
        }
        else {
          i.path[i.path_len++] = nullptr;
        }
      }
    }
    return i;
  }
  iterator lower_bound(key_type k)const {
    iterator i = find(k);
    if (i) return i;
    return ++i;
  }
  iterator upper_bound(key_type k)const {
    return ++find(k);
  }
private:
  static persistent_siphash_id_trie unset_leaf(iterator const& i) {
    if (i.back()) {
      assert(i.back().is_leaf());
      return walk_change_upwards(i, i.path_len - 1, nullptr);
    }
    else {
      return i.path[0];
    }
  }
  static persistent_siphash_id_trie set_leaf(iterator const& i, child_ptr new_leaf) {
    child_ptr walker;
    num_bits_type replacing_idx = i.path_len - 1;
    if (i.back()) {
      assert(i.back().is_leaf());
      walker = new_leaf;
    }
    else if (i.path_len < 2) {
      return new_leaf;
    }
    else {
      child_ptr parent = i.path[i.path_len - 2];
      if (parent.is_leaf()) {
        --replacing_idx;
        num_bits_type which_bits = i.bits(replacing_idx);
        const siphash_id k1 = parent.leaf().id();
        num_bits_type child0 = which_child(i.id, which_bits);
        num_bits_type child1 = which_child(k1, which_bits);
        child_ptr split_at;
        if (child0 == child1) {
          split_at = walker = allocate_node(1);
          while (true) {
            num_bits_type old_child = child0;
            child_ptr old_split_at = split_at;
            which_bits = which_bits - bits_per_level;
            child0 = which_child(i.id, which_bits);
            child1 = which_child(k1, which_bits);
            old_split_at.node().children_exist = 1<<old_child;
            if (child0 == child1) {
              old_split_at.child_by_idx(0) = split_at = allocate_node(1);
              split_at.inc_ref();
              assert(old_split_at.num_children() == 1);
              assert(old_split_at.child_by_idx(0));
            }
            else {
              old_split_at.child_by_idx(0) = split_at = allocate_node(2);
              split_at.inc_ref();
              assert(old_split_at.num_children() == 1);
              assert(old_split_at.child_by_idx(0));
              break;
            }
          }
        }
        else {
          split_at = walker = allocate_node(2);
        }
        split_at.node().children_exist = (1<<child0) | (1<<child1);
        split_at.child_by_idx(child1 > child0) = parent;
        split_at.child_by_idx(child0 > child1) = new_leaf;
        assert(split_at.num_children() == 2);
        assert(split_at.child_by_idx(0));
        assert(split_at.child_by_idx(1));
        parent.inc_ref();
        new_leaf.inc_ref();
      }
      else {
        walker = new_leaf;
      }
    }
    return walk_change_upwards(i, replacing_idx, walker);
  }
    
  static persistent_siphash_id_trie walk_change_upwards(iterator const& i, num_bits_type replacing_idx, child_ptr walker) {
    child_ptr old_child = i.path[replacing_idx];
    --replacing_idx;
    for (; replacing_idx != num_bits_type(-1); --replacing_idx) {
      child_ptr old_parent = i.path[replacing_idx];
      num_bits_type which_bits = i.bits(replacing_idx);
      assert(!old_parent.is_leaf());
      
      const num_bits_type child_idx = which_child(i.id, which_bits);
      assert (walker != old_child);
      const num_bits_type new_num_children = old_parent.num_children() + (walker && !old_child) - (old_child && !walker);
      if ((new_num_children > 1) || (walker && !walker.is_leaf())) {
        child_ptr walker_parent = allocate_node(new_num_children);
        num_bits_type old_idx = 0;
        num_bits_type new_idx = 0;
        for (num_bits_type i = 0; i < (1<<bits_per_level); ++i) {
          if (i == child_idx) {
            if (old_child) { ++old_idx; }
            if (walker) {
              walker_parent.node().children_exist |= 1 << i;
              walker_parent.child_by_idx(new_idx++) = walker;
              walker.inc_ref();
            }
          }
          else if (old_parent.node().children_exist & (1<<i)) {
            walker_parent.node().children_exist |= 1 << i;
            child_ptr old_child_ptr = old_parent.child_by_idx(old_idx++);
            walker_parent.child_by_idx(new_idx++) = old_child_ptr;
            old_child_ptr.inc_ref();
          }
        }
        assert(new_idx == new_num_children);
        assert(walker_parent.num_children() == new_num_children);
        walker = walker_parent;
      }
      else if ((new_num_children == 1) && (!walker)) {
        assert(old_parent.num_children() == 2);
        num_bits_type other_child_bit = old_parent.node().children_exist & ~(1<<child_idx);
        child_ptr other_child = old_parent.child_by_bit(other_child_bit); // TODO optimize
        if (other_child.is_leaf()) {
          walker = other_child;
        }
        else {
          walker = allocate_node(1);
          walker.node().children_exist = other_child_bit;
          walker.child_by_idx(0) = other_child;
          other_child.inc_ref();
        }
      }
      
      old_child = old_parent;
    }
    
    return walker;
  }
};
} // namespace persistent_siphash_id_trie_system

typedef persistent_siphash_id_trie_system::persistent_siphash_id_trie<> persistent_id_set;
template<typename T>
using persistent_id_map = persistent_siphash_id_trie_system::persistent_siphash_id_trie<T>;

template<typename IntType>
struct combine_hash_with_integer {
  siphash_id operator()(siphash_id const& hash, IntType i)const noexcept {
    return siphash_id::combining(hash, static_cast<uint64_t>(i));
  }
};

template<typename TimeType = int64_t, TimeType Never = std::numeric_limits<TimeType>::min(), TimeType MinTime = Never+1, TimeType MaxTime = std::numeric_limits<TimeType>::max(), class CombineHashWithTimeTypeFuncType = combine_hash_with_integer<TimeType>>
struct time_type_info {
  static_assert(Never < MinTime, "Never must be less than MinTime");
  typedef TimeType time_type;
  static const TimeType never = Never;
  static const TimeType min_time = MinTime;
  static const TimeType max_time = MaxTime;
  typedef CombineHashWithTimeTypeFuncType combine_hash_with_time_type_func_type;
};

typedef siphash_id entity_id;
typedef siphash_id trigger_id;
const entity_id global_object_id = siphash_id::least();

typedef size_t field_base_id;
template<typename T>
using optional = boost::optional<T>;
using boost::none;

namespace fields_list_impl {
  class fields_list_nature_base {};
  class empty_fields_list_nature : public fields_list_nature_base {};
  class nonempty_fields_list_nature : public fields_list_nature_base {};
  class generic_input_nature : public fields_list_nature_base {};
  class fields_list_entry_nature : public fields_list_nature_base {};
  
  struct two { char c[2]; };
  template<class T> std::enable_if_t<std::is_base_of<fields_list_nature_base, typename T::fields_list_nature>::value, char> test(int);
  template<class T> two test(...);
  template<class T, bool is_specified> struct get_fields_list_nature_impl { typedef typename T::fields_list_nature type; };
  template<class T> struct get_fields_list_nature_impl<T, false> { typedef generic_input_nature type; };
  template<class T> struct get_fields_list_nature { typedef typename get_fields_list_nature_impl<T, sizeof(test<T>(0))==1>::type type; };
  
  template<typename T>
  struct default_field_traits {
    static inline T null() { return T(); }
    static inline bool is_null(T const& t) { return bool(t); }
  };
  template<typename T, typename Traits>
  struct null_const_ref {
    static inline T const& get() {
      static const T value = Traits::null();
      return value;
    }
  };
  
  template<typename IdentifyingType, typename InnerDataType = IdentifyingType, typename InnerDataTraits = default_field_traits<InnerDataType>, bool PerID = false> class fields_list_entry;
  template<typename IdentifyingType, typename InnerDataType, typename InnerDataTraits>
  class fields_list_entry<IdentifyingType, InnerDataType, InnerDataTraits, false> {
  public:
    typedef fields_list_entry_nature fields_list_nature;
    typedef IdentifyingType identifying_type;
    typedef InnerDataType inner_data_type;
    typedef InnerDataTraits traits;
    static const bool per_id = false;
  };
  template<typename IdentifyingType, typename InnerDataType, typename InnerDataTraits>
  class fields_list_entry<IdentifyingType, InnerDataType, InnerDataTraits, true> {
  public:
    typedef fields_list_entry_nature fields_list_nature;
    typedef IdentifyingType identifying_type;
    typedef InnerDataType inner_data_type;
    typedef InnerDataTraits traits;
    static const bool per_id = true;
  };
  template<typename FieldsListEntry, template<typename> class Data, bool Persistent, bool PerID = FieldsListEntry::per_id> struct fields_list_entry_data; 
  template<typename FieldsListEntry, template<typename> class Data, bool Persistent>
  struct fields_list_entry_data<FieldsListEntry, Data, Persistent, false> {
    typedef Data<typename FieldsListEntry::inner_data_type> type;
    type      & get()      { return t; }
    type const& get()const { return t; }
    type const* find()const { return &t; }
    template<class... Args>
    type& set(type& t, Args&&... args){ return t = type(std::forward<Args>(args)...); }
  private:
    type t;
  };
  template<typename FieldsListEntry, template<typename> class Data>
  struct fields_list_entry_data<FieldsListEntry, Data, true, true> {
    typedef Data<typename FieldsListEntry::inner_data_type> inner_type;
    typedef persistent_id_map<inner_type> type;
    inner_type const& get(siphash_id id)const {
      auto i = t.find(id);
      if (i) { return i->second; }
      return null_const_ref<typename FieldsListEntry::inner_data_type, typename FieldsListEntry::traits>::get();
    }
    template<class... Args>
    inner_type const& set(siphash_id id, Args&&... args) {
      auto i = t.find(id);
      t = t.emplace(id, std::forward<Args>(args)...);
      return i->second;
    }
  private:
    type t;
  };
  template<typename FieldsListEntry, template<typename> class Data>
  struct fields_list_entry_data<FieldsListEntry, Data, false, true> {
    typedef Data<typename FieldsListEntry::inner_data_type> inner_type;
    typedef std::map<siphash_id, inner_type> type;
    inner_type& get(siphash_id id) {
      return t[id];
    }
    inner_type const* find(siphash_id id)const {
      auto i = t.find(id);
      if (i != t.end()) { return &i->second; }
      return nullptr;
    }
  private:
    type t;
  };
  
  template<typename Input, typename InputNature> struct make_fields_list_entry;
  template<typename Input> struct make_fields_list_entry<Input, generic_input_nature> { typedef fields_list_entry<Input, optional<Input>> type; };
  template<typename Input> struct make_fields_list_entry<Input, fields_list_entry_nature> { typedef Input type; };
  
  template<typename FieldsList, typename FieldID> struct get_field_entry;
  template<typename FieldsList, typename FieldID, bool Same> struct get_field_entry_impl;
  template<typename FieldsList, typename FieldID> struct get_field_entry_impl<FieldsList, FieldID, true> { typedef typename FieldsList::head type; };
  template<typename FieldsList, typename FieldID> struct get_field_entry_impl<FieldsList, FieldID, false> { typedef typename get_field_entry<typename FieldsList::tail, FieldID>::type type; };
  template<typename FieldsList, typename FieldID> struct get_field_entry {
    typedef typename get_field_entry_impl<FieldsList, FieldID, std::is_same<FieldID, typename FieldsList::head::identifying_type>::value>::type type;
  };
  template<typename FieldsList, typename FieldID>
  using field_entry = typename get_field_entry<FieldsList, FieldID>::type;
  template<typename FieldsList, typename FieldID>
  using field_data = typename field_entry<FieldsList, FieldID>::inner_data_type;
  
  template<typename ...Input> struct make_fields_list;
  template<typename HeadNature, typename Head, typename ...Tail> struct make_fields_list_contents;
  
  class empty_fields_list {
  public:
    typedef empty_fields_list_nature fields_list_nature;
    static const size_t size = 0;
    //template<typename T> static constexpr field_base_id idx_of() { static_assert(false, "No field of that type exists"); }
    // TODO find a way to make invalid-field errors more readable
  };
  template<typename Head, typename ...Tail>
  class nonempty_fields_list {
  public:
    typedef nonempty_fields_list_nature fields_list_nature;
    typedef typename make_fields_list_entry<Head, typename get_fields_list_nature<Head>::type>::type head;
    typedef typename make_fields_list<Tail...>::type tail;
    static const field_base_id idx = tail::size;
    static const size_t size = tail::size + 1;
    template<typename T> static constexpr std::enable_if_t< std::is_same<T, typename head::identifying_type>::value, field_base_id> idx_of() { return idx; }
    template<typename T> static constexpr std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value, field_base_id> idx_of() { return tail::template idx_of<T>(); }
    template<typename T> static inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value, field_data<nonempty_fields_list, T> const&>
    get_null_const_ref() { return null_const_ref<typename head::inner_data_type, typename head::traits>::get(); }
    template<typename T> static inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value, field_data<nonempty_fields_list, T> const&>
    get_null_const_ref() { return tail::template get_null_const_ref<T>(); }
  };
  
  template<typename HeadNature, typename Head, typename ...Tail>
  struct make_fields_list_contents                                             { typedef nonempty_fields_list<Head, Tail...> type; };
  template<typename Head, typename ...Tail>
  struct make_fields_list_contents<   empty_fields_list_nature, Head, Tail...> { typedef    empty_fields_list                type; };
  template<typename Head, typename ...Tail>
  struct make_fields_list_contents<nonempty_fields_list_nature, Head, Tail...> {
    typedef typename make_fields_list<typename Head::head, typename Head::tail, Tail...>::type type;
  };
  
  template<>
  struct make_fields_list<> { typedef empty_fields_list type; };
  template<typename Head, typename ...Tail>
  struct make_fields_list<Head, Tail...> { typedef typename make_fields_list_contents<typename get_fields_list_nature<Head>::type, Head, Tail...>::type type; };

  template<typename T> using identity = T;
  template<class FieldsList, bool Persistent, template<typename> class repeated = identity>
  class foreach_field {
    typedef typename FieldsList::head head;
    fields_list_entry_data<head, repeated, Persistent> head_;
    foreach_field<typename FieldsList::tail, Persistent, repeated> tail_;
  public:
    template<typename T, typename... Ext> inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value,
    repeated<field_data<FieldsList, T>>      &> get(Ext... ext)      { return head_.get(ext...); }
    template<typename T, typename... Ext> inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value,
    repeated<field_data<FieldsList, T>>      &> get(Ext... ext)      { return tail_.template get<T>(ext...); }
    template<typename T, typename... Ext> inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value,
    repeated<field_data<FieldsList, T>> const*> find(Ext... ext)const { return head_.find(ext...); }
    template<typename T, typename... Ext> inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value,
    repeated<field_data<FieldsList, T>> const*> find(Ext... ext)const { return tail_.template find<T>(ext...); }
    template<typename T, typename... Args> inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value && !field_entry<FieldsList, T>::per_id,
    repeated<field_data<FieldsList, T>> const&> set(                  Args&&... args) { return head_.set(head_,       std::forward<Args>(args)...); }
    template<typename T, typename... Args> inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value && !field_entry<FieldsList, T>::per_id,
    repeated<field_data<FieldsList, T>> const&> set(                  Args&&... args) { return tail_.template set<T, Args...>(       std::forward<Args>(args)...); }
    template<typename T, typename... Args> inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value &&  field_entry<FieldsList, T>::per_id,
    repeated<field_data<FieldsList, T>> const&> set(siphash_id which, Args&&... args) { return head_.set(head_, which, std::forward<Args>(args)...); }
    template<typename T, typename... Args> inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value &&  field_entry<FieldsList, T>::per_id,
    repeated<field_data<FieldsList, T>> const&> set(siphash_id which, Args&&... args) { return tail_.template set<T, Args...>(which, std::forward<Args>(args)...); }
  };
  template<bool Persistent, template<typename> class repeated>
  class foreach_field<empty_fields_list, Persistent, repeated> {};
  
  template<class FieldsList, template<typename> class repeated = identity>
  class foreach_field_base {
    typedef typename FieldsList::head head;
    repeated<typename head::inner_data_type> head_;
    foreach_field_base<typename FieldsList::tail, repeated> tail_;
  public:
    template<typename T> inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value,
    field_data<FieldsList, T>      &> get()      { return head_; }
    template<typename T> inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value,
    field_data<FieldsList, T>      &> get()      { return tail_.template get<T>(); }
    template<typename T> inline std::enable_if_t< std::is_same<T, typename head::identifying_type>::value,
    field_data<FieldsList, T> const&> get()const { return head_; }
    template<typename T> inline std::enable_if_t<!std::is_same<T, typename head::identifying_type>::value,
    field_data<FieldsList, T> const&> get()const { return tail_.template get<T>(); }
  };
  template<template<typename> class repeated>
  class foreach_field_base<empty_fields_list, repeated> {};
  
  template<class FieldsList, typename repeated>
  class foreach_field_base_array : public std::array<repeated, FieldsList::size> {
  public:
    template<typename T> inline repeated& get() { return (*this)[FieldsList::template idx_of<T>()]; }
    template<typename T> inline repeated const& get()const { return (*this)[FieldsList::template idx_of<T>()]; }
  };
  
  template<class FieldsList, class FuncClass>
  using func_type = decltype(&FuncClass::template func<typename FieldsList::head::identifying_type>);
  template<class FieldsList, class FuncClass>
  class function_array : public foreach_field_base_array<FieldsList, func_type<FieldsList, FuncClass>> {
  public:
    function_array() { add_fptr<FieldsList>(); }
  private:
    template<class SubList> inline std::enable_if_t<!std::is_same<SubList, empty_fields_list>::value> add_fptr() {
      typedef typename SubList::head::identifying_type T;
      foreach_field_base_array<FieldsList, func_type<FieldsList, FuncClass>>::template get<T>() = &FuncClass::template func<T>;
      add_fptr<typename SubList::tail>();
    }
    template<class SubList> inline std::enable_if_t< std::is_same<SubList, empty_fields_list>::value> add_fptr() {}
  };
}
template<typename ...Input> using fields_list = typename fields_list_impl::make_fields_list<Input...>::type;
template<typename IdentifyingType, typename InnerDataType = IdentifyingType, typename InnerDataTraits = fields_list_impl::default_field_traits<InnerDataType>>
using field        = fields_list_impl::fields_list_entry<IdentifyingType, optional<InnerDataType>, InnerDataTraits, false>;
template<typename IdentifyingType, typename InnerDataType = IdentifyingType, typename InnerDataTraits = fields_list_impl::default_field_traits<InnerDataType>>
using field_per_id = fields_list_impl::fields_list_entry<IdentifyingType, optional<InnerDataType>, InnerDataTraits, true>;
template<typename FieldsList, typename FieldID>
using field_data = fields_list_impl::field_data<FieldsList, FieldID>;

struct field_id {
  field_id(field_base_id base):base(base),which(){}
  field_id(field_base_id base, siphash_id which):base(base),which(which){}
  field_base_id base;
  siphash_id which;
  bool operator==(field_id const& o)const { return (base == o.base) && (which == o.which); }
};

struct entity_field_id {
  entity_field_id(entity_id e, field_id f):e(e),f(f){}
  entity_id e;
  field_id f;
  bool operator==(entity_field_id const& o)const { return (e == o.e) && (f == o.f); }
};

} //end namespace time_steward_system
namespace std {
  template<>
  struct hash<time_steward_system::entity_field_id> {
    public:
    size_t operator()(time_steward_system::entity_field_id const& i)const {
      return std::hash<time_steward_system::siphash_id>()(i.e) ^ std::hash<time_steward_system::siphash_id>()(i.f.which) ^ i.f.base;
    }
  };
}
namespace time_steward_system {

template<typename TimeSteward>
class time_steward_accessor {
private:
  typedef typename TimeSteward::entity_fields entity_fields;
  typedef typename TimeSteward::extended_time extended_time;
  typedef typename TimeSteward::time_type time_type;
  typedef typename TimeSteward::event event;
  typedef typename TimeSteward::trigger trigger;
  
  struct field_metadata {
    bool accessed_preexisting_state;
    bool ever_modified;
    bool accessed(time_steward_accessor const& acc, entity_field_id const& id) {
      // If we already overwrote a field,
      // then this is not accessing the field's state *before* this time.
      if (!ever_modified && !accessed_preexisting_state) {
        accessed_preexisting_state = true;
        acc.entity_fields_preexisting_state_accessed.push_back(id);
        return true;
      }
      return false;
    }
    void modified(time_steward_accessor const& acc, entity_field_id const& id) {
      if (!ever_modified) {
        ever_modified = true;
        acc.entity_fields_modified.push_back(id);
      }
    }
  };
  
  template<typename FieldData>
  struct field_info {
    FieldData value;
    field_metadata metadata;
  };
    
  struct entity_info {
    entity_id id;
    fields_list_impl::foreach_field<entity_fields, false, field_info> fields;
  };
  
public:
  struct entity_ref {
  public:
    entity_ref():data(nullptr){}
    explicit operator bool()const { return bool(data); }
    entity_id id()const { return data->id; }
  private:
    entity_info* data;
    entity_ref(entity_info* data):data(data){}
    friend class time_steward_accessor;
  };
  
private:
  template<typename FieldID, typename... Ext>
  inline field_data<entity_fields, FieldID>& get_impl(entity_ref e, bool modified, Ext... ext)const {
    const field_id idx(entity_fields::template idx_of<FieldID>(), ext...);
    auto& f = e.data->fields.template get<FieldID>(ext...);
    auto& result = f.value;
    const entity_field_id id(e.id(), idx);
    if (f.metadata.accessed(*this, id)) {
      result = ts_->template get_provisional_entity_field_before<FieldID>(e.id(), time_, ext...);
    }
    if (modified) { f.metadata.modified(*this, id); }
    return result;
  }
  template<typename FieldID, typename... Ext>
  inline field_data<entity_fields, FieldID>& set_impl(entity_ref e, field_data<entity_fields, FieldID> new_contents, Ext... ext) {
    const field_id idx(entity_fields::template idx_of<FieldID>(), ext...);
    auto& f = e.data->fields.template get<FieldID>(ext...);
    auto& result = f.value;
    const entity_field_id id(e.id(), idx);
    f.metadata.modified(*this, id);
    result = new_contents;
    return result;
  }
  
public:
  inline entity_ref get(entity_id id)const {
    entity_info& e = entities[id];
    e.id = id;
    return entity_ref(&e);
  }
  
  template<typename FieldID, typename... Ext>
  inline field_data<entity_fields, FieldID> const& get    (entity_ref e, Ext... ext)const { return get_impl<FieldID>(e, false, ext...); }
  template<typename FieldID, typename... Ext>
  inline field_data<entity_fields, FieldID>      & get_mut(entity_ref e, Ext... ext)      { return get_impl<FieldID>(e,  true, ext...); }
  template<typename FieldID>
  inline field_data<entity_fields, FieldID>& set(entity_ref e,                   field_data<entity_fields, FieldID> new_contents) { return set_impl<FieldID>(e, new_contents       ); }
  template<typename FieldID>
  inline field_data<entity_fields, FieldID>& set(entity_ref e, siphash_id which, field_data<entity_fields, FieldID> new_contents) { return set_impl<FieldID>(e, new_contents, which); }
  
  void anticipate_event(time_type when, std::shared_ptr<const event> e) {
    assert(is_trigger_); // TODO an exception
    assert(when >= time_->base_time); // TODO an exception
    const extended_time ext_when = (when == time_->base_time) ?
      /*TimeSteward::*/ts_->make_event_time(time_, create_id()) :
      /*TimeSteward::*/ts_->make_event_time(when , create_id());
    assert(ext_when > time_);
    new_upcoming_events.push_back(std::make_pair(ext_when, e));
  };
  void set_trigger(trigger_id id, std::shared_ptr<const trigger> t) {
    trigger_changes[id] = t;
  };
  trigger_id set_trigger(std::shared_ptr<const trigger> t) {
    trigger_id id = create_id();
    set_trigger(id, t);
    return id;
  };
  time_type now()const {
    return time_->base_time;
  }
  entity_ref create_entity() {
    // All entity IDs must be unique.
    // This function is one of the two ways that new entity IDs are created.
    // (The other is the fixed entity ID of the global object.)
    // As long as the other requirement is preserved - that
    // no two events are allowed to happen at extended_times with the same ID -
    // this will result in all unique entity IDs.
    return get(create_id());
  }
    
  time_steward_accessor(time_steward_accessor const&) = delete;
  time_steward_accessor(time_steward_accessor&&) = delete;
  time_steward_accessor& operator=(time_steward_accessor const&) = delete;
  time_steward_accessor& operator=(time_steward_accessor&&) = delete;
private:
  TimeSteward const* ts_;
  extended_time time_;
  bool is_trigger_;
  size_t ids_created_;
  mutable std::unordered_map<entity_id, entity_info> entities;
  
  std::unordered_map<trigger_id, std::shared_ptr<const trigger>> trigger_changes;
  std::vector<std::pair<extended_time, std::shared_ptr<const event>>> new_upcoming_events;
  mutable std::vector<entity_field_id> entity_fields_modified; // hack - we COULD make this non-mutable, but get_impl is easier this way
  mutable std::vector<entity_field_id> entity_fields_preexisting_state_accessed;

  siphash_id create_id() { return siphash_id::combining(time_->id, ids_created_++); }
    
  time_steward_accessor(TimeSteward const* ts, extended_time time, bool is_trigger):ts_(ts),time_(time),is_trigger_(is_trigger),ids_created_(0){}
  void process_event(event const* e) {
    (*e)(this);
  }
  friend TimeSteward;
};


template<class FieldsList/*, class Event = */, class TimeTypeInfo = time_type_info<>>
class time_steward {
public:
  typedef FieldsList entity_fields;
  typedef time_steward_accessor<time_steward> accessor;
  friend class time_steward_accessor<time_steward>;
  
  typedef typename TimeTypeInfo::time_type time_type;
  static const time_type min_time = TimeTypeInfo::min_time;
  static const time_type max_time = TimeTypeInfo::max_time;
  static const time_type never = TimeTypeInfo::never;
  
  // TODO can the action/event/trigger stuff be template parameters rather than always being virtual classes to be subclassed?
  class event {
  public:
    virtual void operator()(accessor* accessor)const = 0;
    virtual ~event(){}
  };
  typedef event trigger;
private:
  typedef typename TimeTypeInfo::combine_hash_with_time_type_func_type combine_hash_with_time_type_func_type;
  
  static const uint32_t not_a_trigger = std::numeric_limits<uint32_t>::max();
  struct extended_time_metadata {
    typedef typename ordered_stuff<extended_time_metadata>::entry_ref extended_time;
    
    extended_time_metadata(time_type base_time):base_time(base_time){}
    extended_time_metadata(time_type base_time, siphash_id id):
      base_time(base_time),id(id),trigger_iteration(not_a_trigger){}
    extended_time_metadata(time_type base_time, siphash_id id, uint32_t trigger_iteration, extended_time closest_non_trigger_ancestor):
      base_time(base_time),id(id),trigger_iteration(trigger_iteration),closest_non_trigger_ancestor(closest_non_trigger_ancestor){}
    
    struct sort_sibling_extended_times_by_id {
      bool operator()(extended_time a, extended_time b)const {
        if (a->trigger_iteration != b->trigger_iteration) { return a->trigger_iteration < b->trigger_iteration; }
        return a->id < b->id;
      }
    };
    
    // TODO: can we optimize for size in common cases?
    time_type base_time;
    siphash_id id; // unused for base_time_roots
    uint32_t trigger_iteration; // constant not_a_trigger for all non-triggers
    extended_time closest_non_trigger_ancestor; // unused for all non-triggers
    std::set<extended_time, sort_sibling_extended_times_by_id> children; // usually empty
    
    bool operator==(extended_time_metadata const& o)const {
      return id == o.id;
    }
  };
  
  typedef typename extended_time_metadata::extended_time extended_time;
  struct sort_extended_times_by_base_time {
    bool operator()(extended_time a, extended_time b)const {
      return a->base_time < b->base_time;
    }
  };
  // statics must also be declared outside the class like other globals
  //   but we can't figure out how to do that WRT templating (id U+cErimeRjHMjQ)
  /*static*/ mutable ordered_stuff<extended_time_metadata> all_extended_times;
  /*static*/ mutable std::set<extended_time, sort_extended_times_by_base_time> base_time_roots;
  
  /*static*/ extended_time get_base_time_root(time_type t)const {
    if (base_time_roots.empty()) {
      // create a sentinel with no children at the end
      // time_type() is a hack - without it, compiler tries to pass max_time as reference and gets undefined reference
      const extended_time sentinel = all_extended_times.construct(time_type(max_time));
      all_extended_times.put_only(sentinel);
      base_time_roots.insert(sentinel);
    }
    const extended_time result = all_extended_times.construct(t);
    const auto i = base_time_roots.lower_bound(result);
    assert(i != base_time_roots.end());
    if ((*i)->base_time == t) { return *i; }
    all_extended_times.put_before(result, *i);
    base_time_roots.emplace_hint(i, result);
    return result;
  }
  /*static*/ extended_time get_base_time_end(time_type t)const {
    extended_time root = get_base_time_root(t);
    if (root->children.empty()) {
      // create a sentinel with no children at the end - TODO reduce duplicate code ID GszZcyzY/wwUQg
      const extended_time sentinel = all_extended_times.construct(t, siphash_id::greatest());
      all_extended_times.put_after(sentinel, root);
      root->children.insert(sentinel);
      return sentinel;
    }
    return *boost::prior(root->children.end());
  }
  template <class... Args>
  /*static*/ extended_time make_extended_time_impl(extended_time parent, Args&&... args)const {
    if (parent->children.empty()) {
      // create a sentinel with no children at the end - TODO reduce duplicate code ID GszZcyzY/wwUQg
      const extended_time sentinel = all_extended_times.construct(parent->base_time, siphash_id::greatest());
      all_extended_times.put_after(sentinel, parent);
      parent->children.insert(sentinel);
    }
    const extended_time result = all_extended_times.construct(parent->base_time, std::forward<Args>(args)...);
    const auto i = parent->children.lower_bound(result);
    assert(i != parent->children.end());
    if (**i == *result) { return *i; }
    all_extended_times.put_before(result, *i);
    parent->children.emplace_hint(i, result);
    return result;
  }
  /*static*/ extended_time make_event_time(time_type t, siphash_id id)const {
    return make_extended_time_impl(get_base_time_root(t), id);
  }
  /*static*/ extended_time make_event_time(extended_time t, siphash_id id)const {
    return make_extended_time_impl(t->closest_non_trigger_ancestor ? t->closest_non_trigger_ancestor : t, id);
  }
  /*static*/ extended_time trigger_call_time(trigger_id tid, extended_time when_triggered)const {
    const extended_time t = when_triggered->closest_non_trigger_ancestor ? when_triggered->closest_non_trigger_ancestor : when_triggered;
    uint32_t trigger_iteration = (when_triggered->trigger_iteration == not_a_trigger) ? 0 : when_triggered->trigger_iteration;
    
    siphash_id id = siphash_id::combining(t->id, tid, trigger_iteration);
    if ((when_triggered->trigger_iteration == trigger_iteration) && (when_triggered->id >= id)) {
      ++trigger_iteration;
      id = siphash_id::combining(t->id, tid, trigger_iteration);
    }
    
    const extended_time result = make_extended_time_impl(t, id, trigger_iteration, t);
    assert (result > when_triggered);
    return result;
  }
  
  struct trigger_call_info {
    trigger_call_info():field_changes_andor_creations_triggering_this(0){}
    uint32_t field_changes_andor_creations_triggering_this;
    std::set<extended_time> anticipated_events;
  };
  struct event_pile_info {
    event_pile_info(std::shared_ptr<const event> instigating_event, trigger_id tid = trigger_id::null()) :
      instigating_event (instigating_event),
      tid(tid),
      creation_cut_off (0),
      has_been_executed (false)
      {}
    std::shared_ptr<const event> instigating_event;
    trigger_id tid;
    std::unordered_set<entity_field_id> entity_fields_pile_accessed;
    std::unordered_set<entity_field_id> entity_fields_pile_modified;
    std::unordered_set<trigger_id> triggers_changed;
    bool creation_cut_off;
    bool has_been_executed;
    bool should_be_executed()const { return instigating_event && !creation_cut_off; }
  };
  
  struct field_metadata_throughout_time {
    std::map<extended_time, persistent_id_set> triggers_pointing_at_this_changes;
    std::set<extended_time> event_piles_which_accessed_this;
  };
  template<typename FieldData>
  using field_throughout_time = std::map<extended_time, FieldData>;
  struct entity_throughout_time_info {
    fields_list_impl::foreach_field<entity_fields, false, field_throughout_time> fields;
  };
  typedef std::map<extended_time, trigger_call_info> trigger_throughout_time_info;
  
  void update_trigger(bool undoing, trigger_id tid, extended_time field_change_andor_creation_time, bool force_trigger_change = false, std::shared_ptr<const trigger> new_trigger = nullptr) {
    trigger_throughout_time_info& trigger_info = triggers[tid];
    const auto next_call_iter = trigger_info.upper_bound(field_change_andor_creation_time); // upper_bound: if we're IN a trigger call, we still want to trigger it again
    if (next_call_iter != trigger_info.begin()) {
      const auto cut_off_call_iter = boost::prior(next_call_iter);
      const auto end_iter2 = (next_call_iter == trigger_info.end()) ?
        cut_off_call_iter->second.anticipated_events.end() : 
        cut_off_call_iter->second.anticipated_events.lower_bound(next_call_iter->first); // upper/lower shouldn't matter: anticipated events are never triggers
      for (auto i = cut_off_call_iter->second.anticipated_events.upper_bound(field_change_andor_creation_time); i != end_iter2; ++i) {
        const auto pile_iter = event_piles.find(*i);
        if (undoing) {
          pile_iter->second.creation_cut_off = false;
          if (pile_iter->second.should_be_executed()) { event_piles_not_correctly_executed.insert(pile_iter->first); }
        }
        else {
          pile_iter->second.creation_cut_off = true;
          if (pile_iter->second.has_been_executed) { event_piles_not_correctly_executed.insert(pile_iter->first); }
          else                                     { event_piles_not_correctly_executed.erase (pile_iter->first); }
        }
      }
      if (!undoing && !force_trigger_change) { new_trigger = event_piles.find(cut_off_call_iter->first)->second.instigating_event; }
    }
    
    const extended_time new_trigger_call_time = trigger_call_time(tid, field_change_andor_creation_time);
    if (undoing) {
      const auto i = trigger_info.find(new_trigger_call_time); assert(i != trigger_info.end());
      --i->second.field_changes_andor_creations_triggering_this;
      if (i->second.field_changes_andor_creations_triggering_this == 0) {
        trigger_info.erase(i);
        erase_instigating_event(new_trigger_call_time);
      }
    }
    else {
      // TODO minor: fix duplicate code (id TG0+Xm06xcLdjQ)
      const auto p = trigger_info.insert(std::make_pair(new_trigger_call_time, trigger_call_info()));
      if (p.second && new_trigger) {
        insert_instigating_event(new_trigger_call_time, event_pile_info(new_trigger, tid));
      }
      ++p.first->second.field_changes_andor_creations_triggering_this;
    }
  }
  void update_triggers(bool undoing, entity_field_id id, extended_time field_change_time) {
    auto& metadata = get_field_metadata_throughout_time(id);
    const auto next_triggers_iter = metadata.triggers_pointing_at_this_changes.lower_bound(field_change_time);
    if (next_triggers_iter == metadata.triggers_pointing_at_this_changes.begin()) return;
    const auto triggers_iter = boost::prior(next_triggers_iter);
    for (auto tid : triggers_iter->second) {
      update_trigger(undoing, tid, field_change_time);
    }
  }
  
  struct copy_field_change_from_accessor {
    template<typename FieldID>
    static std::enable_if_t<!fields_list_impl::field_entry<entity_fields, FieldID>::per_id>
    impl(time_steward* ts, extended_time time, decltype(accessor::entity_info::fields)& src, entity_id id, siphash_id) {
      const auto p = ts->get_field_throughout_time<FieldID>(id       ).insert(std::make_pair(time, src.template get<FieldID>(     ).value));
      assert(p.second);
    }
    template<typename FieldID>
    static std::enable_if_t< fields_list_impl::field_entry<entity_fields, FieldID>::per_id>
    impl(time_steward* ts, extended_time time, decltype(accessor::entity_info::fields)& src, entity_id id, siphash_id which) {
      const auto p = ts->get_field_throughout_time<FieldID>(id, which).insert(std::make_pair(time, src.template get<FieldID>(which).value));
      assert(p.second);
    }
    template<typename FieldID>
    static void func(time_steward* ts, extended_time time, decltype(accessor::entity_info::fields)& src, entity_id id, siphash_id which = siphash_id::null()) {
      impl<FieldID>(ts, time, src, id, which);
    }
  };
  struct undo_field_change {
    template<typename FieldID>
    static std::enable_if_t<!fields_list_impl::field_entry<entity_fields, FieldID>::per_id>
    impl(time_steward* ts, extended_time time, entity_id id, siphash_id) {
      const auto p = ts->get_field_throughout_time<FieldID>(id).erase(time);
      assert(p);
    }
    template<typename FieldID>
    static std::enable_if_t< fields_list_impl::field_entry<entity_fields, FieldID>::per_id>
    impl(time_steward* ts, extended_time time, entity_id id, siphash_id which) {
      const auto p = ts->get_field_throughout_time<FieldID>(id, which).erase(time);
      assert(p);
    }
    template<typename FieldID>
    static void func(time_steward* ts, extended_time time, entity_id id, siphash_id which = siphash_id::null()) {
      impl<FieldID>(ts, time, id, which);
    }
  };
  
  // statics must also be declared outside the class like other globals
  //   but we can't figure out how to do that WRT templating (id U+cErimeRjHMjQ)
  /*static*/ const fields_list_impl::function_array<entity_fields, copy_field_change_from_accessor> copy_field_change_from_accessor_funcs;
  /*static*/ const fields_list_impl::function_array<entity_fields, undo_field_change> undo_field_change_funcs;
  
  // The events map doesn't need to be ordered (even though it has a meaningful order)
  // because it's just for looking up events whose times we know.
  typedef std::unordered_map<extended_time, event_pile_info> event_piles_map;
  typedef std::unordered_map<entity_id, entity_throughout_time_info> entities_map;
  typedef std::unordered_map<entity_field_id, field_metadata_throughout_time> field_metadata_map;
  entities_map entities;
  field_metadata_map field_metadata;
  template<typename FieldID, typename... Ext>
  field_throughout_time<field_data<entity_fields, FieldID>>& get_field_throughout_time(entity_id const& id, Ext... ext) {
    return entities[id].fields.template get<FieldID>(ext...);
  }
  field_metadata_throughout_time& get_field_metadata_throughout_time(entity_field_id const& id) {
    return field_metadata[id];
  }
public:
  time_steward() {}
  
  void insert_fiat_event(time_type time, uint64_t distinguisher, std::shared_ptr<const event> e) {
    // This function is one of exactly two places where
    // (persistent) extended_times are created.
    // Uniqueness justification:
    // TODO
    const extended_time t = make_event_time(time, combine_hash_with_time_type_func_type()(siphash_id::combining(distinguisher), time));
    // TODO throw an exception if the user inserts two events at the same time with the same distinguisher
    insert_instigating_event(t, event_pile_info(e));
  }
  void erase_fiat_event(time_type time, uint64_t distinguisher) {
    const extended_time t = make_event_time(time, combine_hash_with_time_type_func_type()(siphash_id::combining(distinguisher), time));
    erase_instigating_event(t);
  }
  
  // For external client code to examine the state - 
  //   e.g., for display.
  // Output of these should never find their way into internal client code,
  //   and it's a little strange (but ok) to have them affect FiatEvents.
  /*template<typename T>
  entity_ref<T const> get_entity_data_after(entity_id<T> const& id, time_type const& time) {
    // Note: For collision safety, these extended_times must not persist.
    return dynamic_pointer_cast<T const>(entity_ref<entity const>(
      get_actual_entity_data_before(id, extended_time(time, extended_time::last))));
  }*/
  std::unique_ptr<accessor> accessor_after(time_type const& time) {
    const extended_time end_time = get_base_time_end(time);
    update_through_time(end_time);
    return std::unique_ptr<accessor>(new accessor(this, end_time, false));
  }
  
  // Some functions for external client code to regulate how much processing the time_steward
  //   does in advance of when it's needed.
  void update_until_time(time_type const& time) {
    update_until_time(get_base_time_root(time));
  }
  void update_through_time(time_type const& time) {
    update_through_time(get_base_time_end(time));
  }
private:
  event_piles_map event_piles;
  std::set<extended_time> event_piles_not_correctly_executed;
  std::unordered_map<trigger_id, trigger_throughout_time_info> triggers;
  
  void update_through_time(extended_time time) {
    while (!is_updated_through(time)) execute_event_pile(*event_piles_not_correctly_executed.begin());
  }
  void update_until_time  (extended_time time) {
    while (!is_updated_until  (time)) execute_event_pile(*event_piles_not_correctly_executed.begin());
  }
  bool is_updated_until  (extended_time time)const {
    return event_piles_not_correctly_executed.empty() || *event_piles_not_correctly_executed.begin() >= time;
  }
  bool is_updated_through(extended_time time)const {
    return event_piles_not_correctly_executed.empty() || *event_piles_not_correctly_executed.begin() >  time;
  }
  /*entity const* get_actual_entity_data_before(entity_id<entity> const& id, extended_time time) {
    update_until_time(time);
    return get_provisional_entity_data_before(id, time);
  }*/
  template<typename FieldID, typename... Ext>
  field_data<entity_fields, FieldID> const& get_provisional_entity_field_before(entity_id id, extended_time time, Ext... ext)const {
    const auto entity_stream_iter = entities.find(id);
    if (entity_stream_iter == entities.end()) { return entity_fields::template get_null_const_ref<FieldID>(); }
    
    entity_throughout_time_info const& e = entity_stream_iter->second;
    field_throughout_time<field_data<entity_fields, FieldID>> const* f = e.fields.template find<FieldID>(ext...);
    if (!f) { return entity_fields::template get_null_const_ref<FieldID>(); }
    
    const auto next_change_iter = f->upper_bound(time);
    if (next_change_iter == f->begin()) { return entity_fields::template get_null_const_ref<FieldID>(); }
    return boost::prior(next_change_iter)->second;
  }
  
  void insert_instigating_event(extended_time time, event_pile_info const& e) {
    const auto p1 = event_piles.insert(std::make_pair(time, e));
    assert(p1.second);
    const auto p2 = event_piles_not_correctly_executed.insert(time);
    assert(*p2.first == time);
    assert(p2.second);
  }
  void erase_instigating_event(extended_time time) {
    const auto event_pile_iter = event_piles.find(time);
    if (event_pile_iter != event_piles.end()) {
      event_pile_info& pile_info = event_pile_iter->second;
      
      pile_info.instigating_event = nullptr;
      if (pile_info.has_been_executed) {
        event_piles_not_correctly_executed.insert(time);
      }
      else {
        event_piles_not_correctly_executed.erase(time);
        event_piles.erase(event_pile_iter);
      }
    }
  }
  // One of the essential, nonintuitive strengths of this system is that
  // execute_event_pile() is a safe operation.
  // If you execute a change that won't actually happen, it will just end
  // up getting invalidated later.
  void execute_event_pile(extended_time time) {
    const auto event_pile_iter = event_piles.find(time);
    assert(event_pile_iter != event_piles.end());
    event_pile_info& pile_info = event_pile_iter->second;
    
    // If the event already has been carried out in a different way,
    // we need to undo that before proceeding.
    const bool event_pile_deleted_for_being_out_of_date = unexecute_event_pile(time);
    if (event_pile_deleted_for_being_out_of_date) return;
    
    if (pile_info.should_be_executed()) {
      // Let's be very explicit about how long the accessor, which is a bit of a hack, is allowed to exist.
      accessor a(this, time, bool(pile_info.tid));
      a.process_event(pile_info.instigating_event.get());
      
      for (entity_field_id const& id : a.entity_fields_preexisting_state_accessed) {
        pile_info.entity_fields_pile_accessed.insert(id);
        auto& metadata = get_field_metadata_throughout_time(id);
        const auto p = metadata.event_piles_which_accessed_this.insert(time);
        assert(p.second);
        if (pile_info.tid) {
          persistent_id_set old_triggers;
          const auto next_triggers_iter = metadata.triggers_pointing_at_this_changes.lower_bound(time); // lower/upper bound shouldn't matter because there shouldn't be an entry now.
          if (next_triggers_iter != metadata.triggers_pointing_at_this_changes.begin()) {
            old_triggers = boost::prior(next_triggers_iter)->second;
          }
          // Hack sZInnn3FkZy0yg: insert even if there's no change, so that the below code can catch it.
          metadata.triggers_pointing_at_this_changes.insert(std::make_pair(time, old_triggers.insert(pile_info.tid)));
        }
      }
      for (entity_field_id const& id : a.entity_fields_modified) {
        pile_info.entity_fields_pile_modified.insert(id);
        (*copy_field_change_from_accessor_funcs[id.f.base])(this, time, a.entities.find(id.e)->second.fields, id.e, id.f.which);
        update_triggers(false, id, time);
      }
      for (std::pair<trigger_id, std::shared_ptr<const trigger>> const& t : a.trigger_changes) {
        pile_info.triggers_changed.insert(t.first);
        update_trigger(false, t.first, time, true, t.second);
      }
      if (pile_info.tid) {
        auto& trigger_info = triggers.find(pile_info.tid)->second;
        const auto call_iter = trigger_info.find(time);
        for (std::pair<extended_time, std::shared_ptr<const event>> const& ev : a.new_upcoming_events) {
          insert_instigating_event(ev.first, event_pile_info(ev.second));
          call_iter->second.anticipated_events.insert(ev.first);
        }
        if (call_iter != trigger_info.begin()) {
          const auto last_call_iter = boost::prior(call_iter);
          const extended_time last_call_time = last_call_iter->first;
          event_pile_info const& last_call_pile_info = event_piles.find(last_call_time)->second;
          for (entity_field_id accessed_by_last_call : last_call_pile_info.entity_fields_pile_accessed) {
            auto& metadata = get_field_metadata_throughout_time(accessed_by_last_call);
            const auto triggers_iter = metadata.triggers_pointing_at_this_changes.find(last_call_time);
            assert (triggers_iter != metadata.triggers_pointing_at_this_changes.end());
            const persistent_id_set old_triggers = triggers_iter->second;
            const auto new_triggers = old_triggers.erase(pile_info.tid);
            assert (new_triggers != old_triggers);
            // Hack sZInnn3FkZy0yg: we inserted above, so this insert fails.
            metadata.triggers_pointing_at_this_changes.insert(std::make_pair(time, new_triggers));
          }
        }
      }
      
      // a is destroyed here
    }
    
    pile_info.has_been_executed = true;
    event_piles_not_correctly_executed.erase(time);
  }
  
  // In essence, this just forgets the consequences of how a
  // particular event happened, without forgetting the fact
  // that it's going to happen.
  // So it too is a "safe" operation (no amount of calling it
  // with bogus values will change us from valid to invalid),
  // because if you undo anything, it's slated to be redone.
  bool unexecute_event_pile(extended_time time) {
    const auto event_pile_iter = event_piles.find(time);
    assert (event_pile_iter != event_piles.end());
    event_pile_info& pile_info = event_pile_iter->second;
    if (!pile_info.has_been_executed) return false;
    
    for (entity_field_id const& id : pile_info.entity_fields_pile_accessed) {
      auto& metadata = get_field_metadata_throughout_time(id);
      const auto j = metadata.event_piles_which_accessed_this.erase(time);
      assert(j);
      if (pile_info.tid) {
        metadata.triggers_pointing_at_this_changes.erase(time);
      }
    }
    for (entity_field_id const& id : pile_info.entity_fields_pile_modified) {
      (*undo_field_change_funcs[id.f.base])(this, time, id.e, id.f.which);
      update_triggers(true, id, time);
    }
    for (trigger_id const& tid : pile_info.triggers_changed) {
      update_trigger(true, tid, time);
    }
    if (pile_info.tid) {
      auto& trigger_info = triggers.find(pile_info.tid)->second;
      auto call_iter = trigger_info.find(time);
      for (extended_time t : call_iter->second.anticipated_events) {
        erase_instigating_event(t);
      }
      call_iter->second.anticipated_events.clear();
      if (call_iter != trigger_info.begin()) {
        for (entity_field_id former_trigger_access : event_piles.find(boost::prior(call_iter)->first)->second.entity_fields_pile_accessed) {
          auto& metadata = get_field_metadata_throughout_time(former_trigger_access);
          metadata.triggers_pointing_at_this_changes.erase(time);
        }
      }
    }
    pile_info.entity_fields_pile_accessed.clear();
    pile_info.entity_fields_pile_modified.clear();
    pile_info.triggers_changed.clear();
    
    pile_info.has_been_executed = false;
    if (pile_info.should_be_executed()) {
      event_piles_not_correctly_executed.insert(time);
    }
    if (!pile_info.instigating_event) {
      event_piles.erase(event_pile_iter);
      return true;
    }
    return false;
  }
};

} //end namespace time_steward_system


