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
class siphash_id {
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
  
  // These constructors combine all the entropy from all the arguments.
  siphash_id(uint64_t i){
    uint64_t in[] = { i };
    data_[0] = siphash24((void*)in, sizeof(in), siphash_key);
    ++in[0];
    data_[1] = siphash24((void*)in, sizeof(in), siphash_key);
  }
  siphash_id(siphash_id a, uint64_t i){
    uint64_t in[] = { a.data_[0], a.data_[1], i };
    data_[0] = siphash24((void*)in, sizeof(in), siphash_key);
    ++in[0];
    data_[1] = siphash24((void*)in, sizeof(in), siphash_key);
  }
  siphash_id(siphash_id a, siphash_id b){
    uint64_t in[] = { a.data_[0], a.data_[1], b.data_[0], b.data_[1] };
    data_[0] = siphash24((void*)in, sizeof(in), siphash_key);
    ++in[0];
    data_[1] = siphash24((void*)in, sizeof(in), siphash_key);
  }
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
  bool operator==(siphash_id const& o)const {
    return (data_[0] == o.data_[0]) && (data_[1] == o.data_[1]);
  }
  bool operator!=(siphash_id const& o)const {
    return (data_[0] != o.data_[0]) || (data_[1] != o.data_[1]);
  }
private:
  constexpr siphash_id(data_type data):data_(data){}
  data_type data_;
  friend class std::hash<siphash_id>;
};

} //end namespace time_steward_system
namespace std {
  template<>
  class hash<time_steward_system::siphash_id> {
    public:
    size_t operator()(time_steward_system::siphash_id const& i)const {
      // Since the siphash is already a hash, just return some of it.
      return i.data_[0];
    }
  };
}
namespace time_steward_system {

template<typename IntType>
struct combine_hash_with_integer {
  siphash_id operator()(siphash_id const& hash, IntType i)const noexcept {
    return siphash_id(hash, static_cast<uint64_t>(i));
  }
};

//typedef siphash_id time_id;

template<typename TimeType = int64_t, TimeType Never = std::numeric_limits<TimeType>::min(), TimeType MinTime = Never+1, class CombineHashWithTimeTypeFuncType = combine_hash_with_integer<TimeType>>
struct time_type_info {
  static_assert(Never < MinTime, "Never must be less than MinTime");
  typedef TimeType time_type;
  static const TimeType never = Never;
  static const TimeType min_time = MinTime;
  typedef CombineHashWithTimeTypeFuncType combine_hash_with_time_type_func_type;
};

typedef siphash_id entity_id;
const entity_id global_object_id = siphash_id::least();
#if 0
template<typename EntitySubclass> struct entity_id;
namespace impl {
template<typename T> entity_id<T> create_entity_id(siphash_id id);
template<typename T> siphash_id get_untyped_id(entity_id<T> id);
}

template<typename EntitySubclass>
struct entity_id {
  entity_id() {}

  // entity_id<E1> can be converted to entity_id<E2> iff E1* can be converted to E2*
  template<typename DerivedEntity> entity_id(entity_id<DerivedEntity> id, EntitySubclass* = static_cast<DerivedEntity*>(nullptr))
    : id_(id.id_) {}
  
  template<class T, class U>
  friend inline entity_id<T> reinterpret_entity_id(entity_id<U> const& r);

  explicit operator bool()const {
    return id_ != siphash_id::null();
  }
private:
  siphash_id id_;

  template<typename T> friend struct entity_id;
  template<typename T> friend entity_id<T> impl::create_entity_id(siphash_id id);
  template<typename T> friend siphash_id impl::get_untyped_id(entity_id<T> id);
};
template<class T, class U>
inline entity_id<T> reinterpret_entity_id(entity_id<U> const& r) {
  return impl::create_entity_id<T>(r.id_);;
}
template<class T, class U>
bool operator==(entity_id<T> a, entity_id<U> b) { return impl::get_untyped_id(a) == impl::get_untyped_id(b); }
template<class T, class U>
bool operator!=(entity_id<T> a, entity_id<U> b) { return impl::get_untyped_id(a) != impl::get_untyped_id(b); }

namespace impl {
template<typename T>
entity_id<T> create_entity_id(siphash_id id) {
  entity_id<T> e;
  e.id_ = id;
  return e;
}
template<typename T>
siphash_id get_untyped_id(entity_id<T> id) {
  return id.id_;
}
}
template<typename T>
entity_id<T> global_object_id() { return impl::create_entity_id<T>(siphash_id::least()); }


template<typename TimeSteward>
class time_steward_accessor;

template<typename EntitySubclass>
class entity_ref {
public:
  // implicit conversion from non-const to const and derived to base
  template<typename DerivedEntity> entity_ref(entity_ref<DerivedEntity> er)
    : e_(er.e_), id_(er.id_) {}

  entity_id<EntitySubclass> id()const {
    return id_;
  }
  EntitySubclass* operator->()const {
    return e_;
  }
  EntitySubclass& operator*()const {
    return *e_;
  }
  explicit operator bool()const {
    return e_;
  }
  
  friend bool operator==(entity_ref a, entity_ref b) { return a.e_ == b.e_; }
  friend bool operator!=(entity_ref a, entity_ref b) { return a.e_ != b.e_; }

private:
  template<typename> friend class time_steward_accessor;
  template<typename> friend class entity_ref;
  template<class T, class U>
  friend entity_ref<T> dynamic_pointer_cast(entity_ref<U> const& r);

  entity_ref(EntitySubclass* e, entity_id<EntitySubclass> id) : e_(e), id_(id) {}

  EntitySubclass* e_;
  entity_id<EntitySubclass> id_;
};
template<class T, class U>
inline entity_ref<T> dynamic_pointer_cast(entity_ref<U> const& r) {
  return entity_ref<T>(dynamic_cast<T*>(r.e_), reinterpret_entity_id<T>(r.id()));
}

} //end namespace time_steward_system
namespace std {
  template<typename EntitySubclass>
  class hash<time_steward_system::entity_id<EntitySubclass>> {
    public:
    size_t operator()(time_steward_system::entity_id<EntitySubclass> const& i)const {
      return std::hash<time_steward_system::siphash_id>()(time_steward_system::impl::get_untyped_id(i));
    }
  };
}
namespace time_steward_system {

template<typename EntitySubclass>
using entity_const_ref = entity_ref<const EntitySubclass>;

class standard_entity_base {
public:
  virtual standard_entity_base* clone()const = 0;
  virtual ~standard_entity_base(){}
};
class null_entity : public standard_entity_base {
public:
  null_entity* clone()const override{ return new null_entity(*this); }
};
#endif

namespace fields_list_impl {
  class fields_list_nature_base { typedef fields_list_nature confirm; };
  class empty_fields_list : public fields_list_nature {};
  class nonempty_fields_list : public fields_list_nature {};
  class generic_input : public fields_list_nature {};
  
  struct two { char c[2]; };
  template<class T> std::enable_if<std::base_of<fields_list_nature_base, T::fields_list_nature>, char> test();
  template<class T> two test();
  template<class T, bool is_specified> struct get_fields_list_nature_impl { typedef T::fields_list_nature type; }
  template<class T> struct get_fields_list_nature_impl<false> { typedef generic_input type; }
  template<class T> struct get_fields_list_nature { typedef get_fields_list_nature_impl<T, sizeof(test<T>())==1>::type type; }

  template<typename ...Input> class fields_list;
  template<typename HeadNature, typename ...Input> class fields_list_contents;
  template<> class fields_list_contents<empty_fields_list> {
    typedef empty_fields_list fields_list_nature;
    static const size_t size = 0;
    template<typename T> static constexpr size_t idx_of() { static_assert(false, "No field of that type exists"); }
  }
  
  template<typename Head, typename ...Tail>
  class fields_list_contents<generic_input, Head, Tail...> {
    typedef nonempty_fields_list fields_list_nature;
    typedef Head head;
    typedef fields_list<Tail...> tail;
    static const size_t idx = tail::size;
    static const size_t size = tail::size + 1;
    template<typename T> static constexpr size_t idx_of() { return tail::idx_of<T>(); }
    template<> static constexpr size_t idx_of<head>() { return idx; }
  }
  
  template<typename Head, typename ...Tail>
  class fields_list_contents<nonempty_fields_list, Head, Tail...> {
    typedef nonempty_fields_list fields_list_nature;
    typedef Head Head::head;
    typedef fields_list<Head::tail, Tail...> tail;
  }
  
  template<>
  class fields_list<> : public fields_list_contents<empty_fields_list> {};
  template<typename Head, typename ...Tail>
  class fields_list<Head, Tail...> : public fields_list_contents<get_fields_list_nature<Head>, Head, Tail...> {};

  template<class FieldsList, template<typename> class repeated>
  class foreach_field {
    typedef repeated<FieldsList::head> head_type;
    head_type head;
    foreach_field<FieldsList::tail, repeated> tail;
    template<typename T> inline repeated<T>& get() { return tail.get<T>(); }
    template<typename T> inline repeated<T> const& get()const { return tail.get<T>(); }
    template<> inline repeated<head_type>& get<FieldsList::head>() { return head; }
    template<> inline repeated<head_type> const& get<FieldsList::head>()const { return head; }
  }
  template<template<typename> class repeated>
  class foreach_field<fields_list<>, repeated> {
    template<typename T> inline repeated<T>& get() { static_assert(false, "No field of that type exists"); }
    template<typename T> inline repeated<T> const& get()const { static_assert(false, "No field of that type exists"); }
  }
  
  template<class FieldsList, class repeated>
  class foreach_field_array {
    std::array<repeated, FieldsList::size> data;
    template<typename T> inline repeated& get() { return data[FieldsList::idx_of<T>()]; }
    template<typename T> inline repeated const& get()const { return data[FieldsList::idx_of<T>()]; }
  }
  
  template<class FieldsList, typename FuncType, template<typename> FuncType func_template>
  class function_array : foreach_field_array<FuncType*> {
    function_array() { add_fptr<FieldsList>(); }
    template<class SubList> inline void add_fptr() {
      get<SubList::head>() = &func_template<SubList::head>;
      add_fptr<SubList::tail>();
    }
    template<> inline void add_fptr<fields_list<>>() {}
  }
  
  /* e.g.
  template<typename T>
  void transfer(foreach_field<list, r1>& a, foreach_field<list, r2>& b) {
    a.get<T>.insert(std::make_pair(time, b.get<T>));
  }
  function_array<list, decltype<transfer<list::head>>, transfer> farray;
  farray.get<field_example>()(a, t);
  */
}
template<typename ...Input> using fields_list = fields_list_impl::fields_list<Input...>;

template<typename TimeSteward>
class time_steward_accessor {
  typedef typename TimeSteward::entity_fields entity_fields;
  typedef typename TimeSteward::entity entity;
  typedef typename TimeSteward::extended_time extended_time;
  typedef typename TimeSteward::time_type time_type;
  typedef typename TimeSteward::action action;
  typedef typename TimeSteward::consequential_event consequential_event;
  typedef typename TimeSteward::trigger trigger;
  typedef typename TimeSteward::trigger_id trigger_id;
  typedef uint32_t invalidation_counter_t;
  
  struct trigger_info {
    trigger_info(std::shared_ptr<const trigger> t, bool in_queue):t(t),in_queue(in_queue){}
    std::shared_ptr<const trigger> t;
    bool in_queue;
  };
  
  template<typename Field>
  using field_info = Field;
  /*struct field_info {
    Field data;
  };*/
  struct field_metadata {
    persistent_map<trigger_id, trigger_info> triggers;
    bool accessed_preexisting_state;
    bool ever_modified;
    bool all_triggers_queued;
    void accessed(time_steward_accessor const& acc, entity_field_id fid) {
      if (acc.event_whenfunc_access_tracker) { acc.event_whenfunc_access_tracker->insert(fid); }
      // If we already overwrote a field,
      // then this is not accessing the field's state *before* this time.
      if (!ever_modified && !accessed_preexisting_state) {
        accessed_preexisting_state = true;
        acc.entity_fields_preexisting_state_accessed.insert(fid);
        triggers = ts_.get_provisional_triggers_active_before(fid);
        return true;
      }
      return false;
    }
    void modified(time_steward_accessor& acc, entity_field_id fid) {
      if (!ever_modified) {
        ever_modified = true;
        acc.entity_fields_modified.insert(fid);
      }
      if (!all_triggers_queued) {
        for (std::pair<trigger_id, trigger_info> const& i : triggers) {
          if (!i.in_queue) {
            acc.triggers_queue.emplace(fid, i.first);
            i.in_queue = true;
          }
        }
        all_triggers_queued = true;
      }
    }
  };
    
  template<typename FieldsList>
  struct entity_info {
    entity_id id;
    fields_list_impl::foreach_field<entity_fields, field_info> fields;
    fields_list_impl::foreach_field_array<field_metadata> metadata;
  };
  
  struct entity_ref {
  public:
    entity_id id()const { return data->id; }
  private:
    entity_info* data;
    friend class time_steward_accessor;
  };
  
  template<typename Field>
  void copy_change_from_accessor(accessor::entity_info& src, entity_throughout_time_info& dst) {
    
  }
  
  template<typename Field>
  inline Field& get_impl(entity_ref e)const {
    const size_t idx = entity_fields::idx_of<Field>();
    Field& result = e.data->fields.get<Field>();
    const entity_field_id fid(e.id, idx);
    if (e.data->metadata[idx].accessed(*this, fid)) {
      result = ts_.get_provisional_entity_field_before<Field>(e.id());
    }
    return result;
  }
  
  struct queued_trigger_info {
    entity_field_id fid;
    trigger_id tid;
  };
    
  struct upcoming_event_info {
    std::shared_ptr<const consequential_event> e;
    extended_time time;
    std::unordered_set<entity_id> dependencies;
  };
  
public:
  inline entity_ref get(entity_id id)const {
    entity_info& e = entities[id];
    e.id = id;
    return &e;
  }
  
  template<typename Field> inline Field const& get    (entity_ref e)const { return get_impl(e); }
  template<typename Field> inline Field      & get_mut(entity_ref e) {
    Field& result = get_impl(e);
    const size_t idx = entity_fields::idx_of<Field>();
    e.data->metadata[idx].modified(*this, entity_field_id(e.id(), idx));
    return result;
  }
  template<typename Field>
  inline Field& set(entity_ref e, Field new_contents) {
    const size_t idx = entity_fields::idx_of<Field>();
    e.data->metadata[idx].modified(*this, entity_field_id(e.id(), idx));
    Field& result = e.data->fields.get<Field>();
    result = new_contents;
    return result;
  }
  
  void create_event(std::unique_ptr<const consequential_event> e) {
    fresh_new_upcoming_events.push_back(std::move(e));
  };
  void create_trigger(entity_field_id fid, trigger_id tid, std::shared_ptr<const trigger> t) {
    entity_ref e = get(fid.eid)
    field_metadata& metadata = e.data->metadata[fid.fid];
    trigger_info& inf = metadata.triggers[tid];
    inf.t = t;
    if (!inf.in_queue) {
      triggers_queue.emplace(fid, tid);
    }
    inf.in_queue = true;
  };
  void delete_trigger(entity_field_id fid, trigger_id tid) {
    entity_ref e = get(fid.eid)
    field_metadata& metadata = e.data->metadata[fid.fid];
    metadata.triggers.erase(tid);
  };
  time_type now()const {
    return time_.base_time;
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
  size_t ids_created_;
  mutable std::unordered_map<entity_id, entity_info> entities;
  std::queue<queued_trigger_info> triggers_queue;
  std::vector<std::unique_ptr<consequential_event>> fresh_new_upcoming_events;
  std::vector<upcoming_event_info> stale_new_upcoming_events;
  std::unordered_set<entity_field_id> *event_whenfunc_access_tracker;
  std::unordered_set<entity_field_id> entity_fields_modified;
  std::unordered_set<entity_field_id> entity_fields_preexisting_state_accessed;

  siphash_id create_id() { return siphash_id(time_.id(), ids_created_++); }
    
  time_steward_accessor(TimeSteward const* ts, extended_time const& time):ts_(ts),time_(time),ids_created_(0){}
  void process_event(event const* e) {
    (*e)(this);
    while (!triggers_queue.empty()) {
      const queued_trigger_info t = triggers_queue.front();
      triggers_queue.pop();
      entity_info& info = *entities.find(t.fid.eid)
      field_metadata& metadata = info.metadata[t.fid.fid];
      auto trigger_iter = metadata.triggers.find(t.tid);
      if (trigger_iter != metadata.triggers.end()) {
        trigger_iter->second.in_queue = false;
        metadata.all_triggers_queued = false;
        (*trigger_iter->second.t)(this, &info);
      }
    }
    
    for (auto new_upcoming_event : fresh_new_upcoming_events) {
      upcoming_event_info inf;
      event_whenfunc_access_tracker = &inf.dependencies;
      const time_type t = new_upcoming_event->when(this);
      event_whenfunc_access_tracker = nullptr;
      
      if (t != never) {
        assert(t >= time_.base_time); // TODO an exception
        inf.e.swap(new_upcoming_event);
        inf.time = (t == time_.base_time) ? extended_time(time_, create_id()) : extended_time(t, create_id());
        stale_new_upcoming_events.push_back(inf);
      }
    }
    fresh_new_upcoming_events.clear();
  }
  friend TimeSteward;
};

namespace impl {
template<typename TimeTypeInfo>
struct extended_time {
  typedef typename TimeTypeInfo::time_type time_type;

  time_type base_time;
  // We construct extended_times with new unique ids;
  // some of them extend old extended_times, so tiebreakers.front() is not unique among extended_times.
  // However, tiebreakers.back() is unique among extended_times.
  // TODO: a structure with better copy asymptotics than std::vector and better constant speed for short values.
  std::vector<siphash_id> tiebreakers;
  siphash_id id()const { return tiebreakers.back(); }
  struct first_t{}; struct last_t{}; static const first_t first; static const last_t last;
  extended_time(time_type base_time, first_t):base_time(base_time),tiebreakers(siphash_id::least()){}
  extended_time(time_type base_time,  last_t):base_time(base_time),tiebreakers(siphash_id::greatest()){}
  extended_time(time_type base_time, siphash_id tiebreaker):base_time(base_time),tiebreakers(id){}
  extended_time(extended_time base_exttime, siphash_id further_tiebreaker)
      :
      base_time(base_exttime.base_time),
      tiebreakers(base_exttime.tiebreakers)
  {
    tiebreakers.push_back(further_tiebreaker);
  }
  
  bool operator==(extended_time const& o)const { return id() == o.id(); }
  bool operator!=(extended_time const& o)const { return id() != o.id(); }
  bool operator<(extended_time const& o)const {
    if (base_time < o.base_time) { return true; }
    if (base_time > o.base_time) { return false; }
    for (size_t i = 0; ; ++i) {
      if (i == o.tiebreakers.size()) { return false; }
      else if (i == tiebreakers.size()) { return true; }
      else {
        if (tiebreakers[i] < o.tiebreakers[i]) { return true; }
        if (tiebreakers[i] > o.tiebreakers[i]) { return false; }
      }
    }
  }
  bool operator>(extended_time const& o)const { return o < *this; }
  bool operator>=(extended_time const& o)const { return !(*this < o); }
  bool operator<=(extended_time const& o)const { return !(o < *this); }
};
}

} //end namespace time_steward_system
namespace std {
  template<typename TimeTypeInfo>
  class hash<time_steward_system::impl::extended_time<TimeTypeInfo>> {
    public:
    size_t operator()(time_steward_system::impl::extended_time<TimeTypeInfo> const& t)const {
      return std::hash<time_steward_system::siphash_id>()(t.id());
    }
  };
}
namespace time_steward_system {

template<class TimeTypeInfo = time_type_info<>, class FieldsList/*, class Event = */>
class time_steward {
public:
  typedef time_steward_accessor<time_steward> accessor;
  typedef typename TimeTypeInfo::time_type time_type;
  static const time_type min_time = TimeTypeInfo::min_time;
  static const time_type never = TimeTypeInfo::never;
  typedef EntityData entity;
  // TODO can the action/event/trigger stuff be template parameters rather than always being virtual classes to be subclassed?
  class action {
  public:
    virtual void operator()(accessor* accessor)const = 0;
    virtual ~action(){}
  };
  typedef action event;
  typedef action trigger;
  class consequential_event : virtual public event {
  public:
    virtual void operator()(accessor* accessor)const = 0;
    virtual time_type when(accessor const* accessor)const = 0;
  };
private:
  typedef typename TimeTypeInfo::combine_hash_with_time_type_func_type combine_hash_with_time_type_func_type;
  typedef impl::extended_time<TimeTypeInfo> extended_time;
  typedef extended_time event_pile_id;
  typedef std::set<extended_time> time_set;
  
  struct event_pile_info {
    event_pile_info(std::shared_ptr<event> instigating_event) :
      instigating_event (instigating_event),
      num_instigating_event_creation_dependencies_cut_off (0),
      has_been_executed (false)
      {}
    event_pile_info(std::shared_ptr<event> instigating_event,
                    extended_time instigating_event_creation_time,
                    std::unordered_set<entity_field_id> instigating_event_creation_dependencies) :
      instigating_event (instigating_event),
      instigating_event_creation_time (instigating_event_creation_time),
      instigating_event_creation_dependencies (instigating_event_creation_dependencies),
      num_instigating_event_creation_dependencies_cut_off (0),
      has_been_executed (false)
      {}
    extended_time instigating_event_creation_time;
    std::unordered_set<entity_field_id> instigating_event_creation_dependencies;
    uint32_t num_instigating_event_creation_dependencies_cut_off;
    std::shared_ptr<event> instigating_event;
    std::unordered_set<entity_field_id> entity_fields_pile_accessed;
    std::unordered_set<entity_field_id> entity_fields_pile_modified;
    std::unordered_set<event_pile_id> instigating_events_pile_created;
    bool has_been_executed;
    bool should_be_executed()const { return instigating_event && (num_instigating_event_creation_dependencies_cut_off == 0); }
  };
  
  
  template<typename Field>
  using field_throughout_time = std::map<extended_time, Field>;
  template<typename Field>
  struct field_metadata_throughout_time {
    std::map<extended_time, persistent_map<trigger_id, trigger_info>> triggers_changes;
    std::unordered_set<event_pile_id> event_piles_which_accessed_this;
    std::set<event_pile_id> event_piles_whose_instigating_event_creation_accessed_this;
  };
  struct entity_throughout_time_info {
    fields_list_impl::foreach_field<FieldsList, field_throughout_time> fields;
    fields_list_impl::foreach_field_array<FieldsList, field_metadata_throughout_time> metadata;
  };
  
  void cut_off_events(bool cut, entity_throughout_time_info& dst, field_id id, extended_time const& cutoff_time, extended_time const& events_before_this_are_already_cut_off) {
    for (auto ev = dst.metadata[id].event_piles_whose_instigating_event_creation_accessed_this.upper_bound(time);
         ev != dst.metadata[id].event_piles_whose_instigating_event_creation_accessed_this.end(); ++ev) {
      event_pile_info& pile_info = *event_piles.find(*ev);
      if ((pile_info.instigating_event_creation_time < cutoff_time) &&
        (pile_info.instigating_event_creation_time >= events_before_this_are_already_cut_off)) {
        if (cut) {
          if (pile_info.num_instigating_event_creation_dependencies_cut_off == 0) {
            if (pile_info.has_been_executed) { event_piles_not_correctly_executed.insert(*ev); }
            else                              { event_piles_not_correctly_executed.erase (*ev); }
          }
          ++pile_info.num_instigating_event_creation_dependencies_cut_off;
        }
        else {
          --pile_info.num_instigating_event_creation_dependencies_cut_off;
          if (pile_info.should_be_executed()) {
            event_piles_not_correctly_executed.insert(*ev);
          }
        }
      }
    }
  }
  template<typename Field>
  void copy_field_change_from_accessor(extended_time const& time, decltype(accessor::entity_info::fields)& src, entity_throughout_time_info& dst) {
    const auto p = dst.fields.get<Field>().insert(std::make_pair(time, src.get<Field>()));
    assert(p.second);
    cut_off_events(true, dst, FieldsList::idx_of<Field>(), time, (p.first == dst.fields.get<Field>().begin()) ? min_time : boost::prior(p.first)->first);
  }
    
  
  void invalidate_event_piles_that_accessed_entity_when_it_was_defined_by_the_change_at(entity_throughout_time_info& e, extended_time const& change_time) {
    const auto next_change = e.changes.upper_bound(change_time);
    // upper_bound(), because references at the same extended_time as a event_pile
    // are that event_pile's own references, and refer to the previous state of the entity.
    const auto invalidated_event_piles_begin = e.event_piles_which_accessed_this.upper_bound(change_time);
    const auto invalidated_event_piles_end = (next_change == e.changes.end()) ? e.event_piles_which_accessed_this.end() : e.event_piles_which_accessed_this.upper_bound(next_change->first);
    // TODO: better to do this with something like std::move?
    for (auto i = invalidated_event_piles_begin; i != invalidated_event_piles_end; ) {
      assert(*i != change_time);
      event_piles_not_correctly_executed.insert(*i);
      e.event_piles_which_accessed_this.erase(i++);
    }
  }
  template<typename Field>
  void undo_field_change(extended_time const& time, entity_throughout_time_info& dst) {
    const auto change_iter = dst.fields.get<Field>().find(time);
    assert (change_iter != dst.fields.get<Field>().end());
    const auto change_end_iter = boost::next(change_iter);
    {
      const auto end_iter = (change_end_iter == dst.fields.get<Field>().end()) ?
        dst.metadata.event_piles_whose_instigating_event_creation_accessed_this.end() : 
        dst.metadata.event_piles_whose_instigating_event_creation_accessed_this.upper_bound(change_end_iter->first);
      for (auto i = dst.metadata.event_piles_whose_instigating_event_creation_accessed_this.upper_bound(time); i != end_iter; ) {
        erase_instigating_event(*(i++));
      }
    }
    {
      const auto begin_iter = dst.metadata.event_piles_which_accessed_this.upper_bound(time);
      const auto end_iter = (change_end_iter == dst.fields.get<Field>().end()) ?
        dst.metadata.event_piles_which_accessed_this.end() : 
        dst.metadata.event_piles_which_accessed_this.upper_bound(change_end_iter->first);
      // TODO is this std::move any faster? I feel that it is less-easy-to-understand code.
      //std::move(begin_iter, end_iter, std::inserter(event_piles_not_correctly_executed, event_piles_not_correctly_executed.end()));
      for (auto i = begin_iter; i != end_iter; ) {
        e.event_piles_which_accessed_this.erase(i++);
        event_piles_not_correctly_executed.insert(*i);
      }
    }
    cut_off_events(false, dst, FieldsList::idx_of<Field>(), time, (change_iter == dst.fields.get<Field>().begin()) ? min_time : boost::prior(change_iter)->first);
    dst.fields.get<Field>().erase(change_iter);
  }
  
  static const fields_list_impl::function_array<
    FieldsList,
    decltype(copy_field_change_from_accessor<FieldsList::head>),
    copy_field_change_from_accessor>
    copy_field_change_from_accessor_funcs;
  static const fields_list_impl::function_array<
    FieldsList,
    decltype(undo_field_change<FieldsList::head>),
    undo_field_change>
    undo_field_change_funcs;
  
  // The events map doesn't need to be ordered (even though it has a meaningful order)
  // because it's just for looking up events whose times we know.
  typedef std::unordered_map<event_pile_id, event_pile_info> event_piles_map;
  typedef std::unordered_map<entity_id, entity_throughout_time_info> entities_map;

  friend class time_steward_accessor<time_steward>;
public:
  time_steward() {}
  
  void insert_fiat_event(time_type time, uint64_t distinguisher, std::unique_ptr<event> e) {
    // This function is one of exactly two places where
    // (persistent) extended_times are created.
    // Uniqueness justification:
    // TODO
    const extended_time t(time, combine_hash_with_time_type_func_type()(siphash_id(distinguisher), time));
    // TODO throw an exception if the user inserts two events at the same time with the same distinguisher
    insert_instigating_event(t, event_pile_info(e));
  }
  void erase_fiat_event(time_type time, uint64_t distinguisher) {
    const extended_time t(time, siphash_id(distinguisher));
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
    extended_time et(time, extended_time::last);
    update_through_time(et);
    return std::unique_ptr<accessor>(new accessor(this, et));
  }
  
  // Some functions for external client code to regulate how much processing the time_steward
  //   does in advance of when it's needed.
  void update_until_time(time_type const& time) {
    // Note: For collision safety, these extended_times must not persist.
    update_until_time(extended_time(time, extended_time::first));
  }
  void update_through_time(time_type const& time) {
    // Note: For collision safety, these extended_times must not persist.
    update_through_time(extended_time(time, extended_time::last));
  }
private:
  entities_map entities;
  event_piles_map event_piles;
  time_set event_piles_not_correctly_executed;
  
  void update_through_time(extended_time const& time) {
    while (!is_updated_through(time)) execute_event_pile(*event_piles_not_correctly_executed.begin());
  }
  void update_until_time  (extended_time const& time) {
    while (!is_updated_until  (time)) execute_event_pile(*event_piles_not_correctly_executed.begin());
  }
  bool is_updated_until  (extended_time const& time)const {
    return event_piles_not_correctly_executed.empty() || *event_piles_not_correctly_executed.begin() >= time;
  }
  bool is_updated_through(extended_time const& time)const {
    return event_piles_not_correctly_executed.empty() || *event_piles_not_correctly_executed.begin() >  time;
  }
  /*entity const* get_actual_entity_data_before(entity_id<entity> const& id, extended_time const& time) {
    update_until_time(time);
    return get_provisional_entity_data_before(id, time);
  }*/
  template<typename Field>
  Field const& get_provisional_entity_field_before(entity_id id, extended_time const& time)const {
    const auto entity_stream_iter = entities.find(id);
    if (entity_stream_iter == entities.end()) {
      return Field();
    }
    
    entity_throughout_time_info const& e = entity_stream_iter->second;
    const auto next_change_iter = e.get<Field>().changes.upper_bound(time);
    if (next_change_iter == e.changes.begin()) {
      return Field();
    }
    return &*boost::prior(next_change_iter)->second;
  }
  
  void insert_instigating_event(extended_time const& time, event_pile_info& e) {
    const auto p1 = event_piles.insert(std::make_pair(time, e));
    assert(p1.second);
    const auto p2 = event_piles_not_correctly_executed.insert(time);
    assert(p2.second);
  }
  void erase_instigating_event(extended_time const& time) {
    const auto event_pile_iter = event_piles.find(time);
    if (event_pile_iter != event_piles.end()) {
      event_piles_not_correctly_executed.erase(time);
      event_pile_info& pile_info = event_pile_iter->second;
      pile_info.instigating_event = nullptr;
      for (entity_field_id id : pile_info.instigating_event_creation_dependencies) {
        const auto entity_stream_iter = entities.find(id.e);
        assert (entity_stream_iter != entities.end());
        field_metadata_throughout_time& metadata = entity_stream_iter->second.metadata[id.f];
        metadata.event_piles_whose_instigating_event_creation_accessed_this.erase(time);
      }
      if (pile_info.has_been_executed) {
        event_piles_not_correctly_executed.insert(time);
      }
      else {
        event_piles_not_correctly_executed.erase(time);
        event_piles.erase(event_iter);
      }
    }
  }
  // One of the essential, nonintuitive strengths of this system is that
  // execute_event_pile() is a safe operation.
  // If you execute a change that won't actually happen, it will just end
  // up getting invalidated later.
  void execute_event_pile(extended_time const& time) {
    const auto event_pile_iter = event_piles.find(time);
    assert(event_pile_iter != event_piles.end());
    event_pile_info& pile_info = event_pile_iter->second;
    
    // If the event already has been carried out in a different way,
    // we need to undo that before proceeding.
    const bool event_pile_deleted_for_being_out_of_date = unexecute_event_pile(time);
    if (event_pile_deleted_for_being_out_of_date) return;
    
    if (pile_info.should_be_executed()) {
      // Let's be very explicit about how long the accessor, which is a bit of a hack, is allowed to exist.
      accessor a(this, time);
      a.process_event(pile_info.instigating_event.get());
      
      for (entity_field_id const& id : a.entity_fields_preexisting_state_accessed) {
        pile_info.entity_fields_pile_accessed.insert(id);
        const auto p = entities[id.e].metadata[id.f].event_piles_which_accessed_this.insert(time);
        assert(p.second);
      }
      for (entity_field_id const& id : a.entity_fields_modified) {
        pile_info.entity_fields_pile_modified.insert(id);
        copy_change_from_accessor_funcs[id.f](time, a.entities.find(id.e)->fields, *entities.find(id.e));
      }
      
      for (accessor::upcoming_event_info ev : a.stale_new_upcoming_events) {
        insert_instigating_event(ev.time, event_pile_info(ev.e, time, ev.dependencies));
        pile_info.instigating_events_pile_created.insert(ev.time);
        for (entity_field_id id : ev.dependencies) {
          entities[id.e].metadata[id.f].event_piles_whose_instigating_event_creation_accessed_this.insert(ev.time);
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
  bool unexecute_event_pile(extended_time const& time) {
    const auto event_pile_iter = event_piles.find(time);
    assert (event_pile_iter != event_piles.end());
    event_pile_info& pile_info = event_pile_iter->second;
    if (!pile_info.has_been_executed) return;
    
    for (extended_time const& t : pile_info.instigating_events_pile_created) {
      erase_instigating_event(t);
    }
    for (entity_field_id const& id : pile_info.entity_fields_pile_accessed) {
      const auto i = entities.find(id.e); assert(i != entities.end());
      const auto j = i->metadata[id.f].event_piles_which_accessed_this.erase(time);
      assert(j);
    }
    for (entity_field_id const& id : pile_info.entity_fields_pile_modified) {
      undo_field_change_funcs[id.f](time, *entities.find(id.e));
    }
    pile_info.entity_fields_pile_accessed.clear();
    pile_info.entity_fields_pile_modified.clear();
    pile_info.instigating_events_pile_created.clear();
    
    pile_info.has_been_executed = false;
    if (pile_info.should_be_executed()) {
      event_piles_not_correctly_executed.insert(time);
    }
    if (!pile_info.instigating_event) {
      event_piles.erase(event_iter);
      return true;
    }
    return false;
  }
};

} //end namespace time_steward_system


