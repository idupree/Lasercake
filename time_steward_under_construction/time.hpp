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

#ifndef LASERCAKE_TIME_HPP__
#define LASERCAKE_TIME_HPP__

#include "siphash_id.hpp"
#include "manual_orderer.hpp"

typedef siphash_id trigger_id;

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

static const uint32_t not_a_trigger = std::numeric_limits<uint32_t>::max();
template<class TimeTypeInfo>
struct extended_time_metadata {
  struct dec_ref_callback {
    void operator()(extended_time_metadata& t, size_t old_ref_count)const {
      if (old_ref_count == 2 && t.parent && t.id != siphash_id::greatest() && t.base_time != TimeTypeInfo::never) {
        // i.e. no base time roots, sentinels, or temporary lookup-in-children hacks
        // TODO handle base time roots too
        t.delete_from_parent();
      }
    }
  };
  typedef typename manual_orderer<extended_time_metadata, dec_ref_callback>::entry_ref extended_time;
  typedef typename TimeTypeInfo::time_type time_type;
  
  void delete_from_parent() {
    //auto e = parent->children->erase(t);
    //assert (e);
    // hack
    auto parent_ = parent;
    parent = extended_time(); // prevent dec_ref_callback recursion
    const extended_time twin = extended_time::construct(time_type(TimeTypeInfo::never), parent_, id, trigger_iteration);
    const auto i = parent_->children->lower_bound(twin);
    assert (i != parent_->children->end());
    if (**i == *this) {
      parent_->children->erase(i);
      if (parent_->children->size() == 1) {
        parent_->children = nullptr;
      }
    }
  }
  
  extended_time_metadata(time_type base_time):
    base_time(base_time),parent(),id(),trigger_iteration(not_a_trigger){}
  extended_time_metadata(time_type base_time, extended_time parent, siphash_id id):
    base_time(base_time),parent(parent),id(id),trigger_iteration(not_a_trigger){}
  extended_time_metadata(time_type base_time, extended_time parent, siphash_id id, uint32_t trigger_iteration):
    base_time(base_time),parent(parent),id(id),trigger_iteration(trigger_iteration){}
  
  struct sort_sibling_extended_times_by_id_and_trigger_iteration {
    bool operator()(extended_time const& a, extended_time const& b)const {
      if (a->trigger_iteration != b->trigger_iteration) { return a->trigger_iteration < b->trigger_iteration; }
      return a->id < b->id;
    }
  };
  
  // TODO: can we optimize for size in common cases?
  time_type base_time;
  extended_time parent;
  siphash_id id; // unused for base_time_roots
  uint32_t trigger_iteration; // constant not_a_trigger for all non-triggers
  typedef std::set<extended_time, sort_sibling_extended_times_by_id_and_trigger_iteration> children_set;
  std::unique_ptr<children_set> children; // usually empty
  
  bool is_trigger()const {
    return trigger_iteration != not_a_trigger;
  }
  bool operator==(extended_time_metadata const& o)const {
    return id == o.id;
  }
};

template<class TimeTypeInfo>
struct extended_time_maker {
public:
  typedef typename TimeTypeInfo::time_type time_type;
  static const time_type min_time = TimeTypeInfo::min_time;
  static const time_type max_base_time = TimeTypeInfo::max_time;
  static const time_type never = TimeTypeInfo::never;
  
  typedef typename extended_time_metadata<TimeTypeInfo>::extended_time extended_time;
  static extended_time event_time(time_type t, siphash_id id) {
    return make_extended_time_impl(get_base_time_root(t), id);
  }
  static extended_time event_time(extended_time const& t, siphash_id id) {
    return make_extended_time_impl(t->is_trigger() ? t->parent : t, id);
  }
  static extended_time trigger_call_time(trigger_id tid, extended_time const& when_triggered) {
    const extended_time t = when_triggered->is_trigger() ? when_triggered->parent : when_triggered;
    uint32_t trigger_iteration = (when_triggered->trigger_iteration == not_a_trigger) ? 0 : when_triggered->trigger_iteration;
    
    siphash_id id = siphash_id::combining(t->id, tid, trigger_iteration);
    if ((when_triggered->trigger_iteration == trigger_iteration) && (when_triggered->id >= id)) {
      ++trigger_iteration;
      id = siphash_id::combining(t->id, tid, trigger_iteration);
    }
    
    const extended_time result = make_extended_time_impl(t, id, trigger_iteration);
    assert (result > when_triggered);
    assert (result->parent);
    return result;
  }
  static extended_time max_time() {
    return get_base_time_root(max_base_time);
  }
  static extended_time after_all_calls_triggered_at(extended_time const& when_triggered) {
    extended_time const& t = when_triggered->is_trigger() ? when_triggered->parent : when_triggered;
    return require_sentinel_child(t);
  }
  static extended_time base_time_end(time_type t) {
    extended_time root = get_base_time_root(t);
    return require_sentinel_child(root);
  }
private:
  struct sort_extended_times_by_base_time {
    bool operator()(extended_time const& a, extended_time const& b)const {
      return a->base_time < b->base_time;
    }
  };
  static manual_orderer<extended_time_metadata<TimeTypeInfo>, typename extended_time_metadata<TimeTypeInfo>::dec_ref_callback>& all_extended_times() {
    static manual_orderer<extended_time_metadata<TimeTypeInfo>, typename extended_time_metadata<TimeTypeInfo>::dec_ref_callback> a; return a;
  }
  static std::set<extended_time, sort_extended_times_by_base_time>& base_time_roots() {
    static std::set<extended_time, sort_extended_times_by_base_time> a; return a;
  }
  static extended_time require_sentinel_child(extended_time const& t) {
    if (!t->children) {
      // create a sentinel with no children at the end
      t->children.reset(new typename extended_time_metadata<TimeTypeInfo>::children_set());
      const extended_time sentinel = extended_time::construct(time_type(never), t, siphash_id::greatest());
      all_extended_times().put_after(sentinel, t);
      t->children->insert(sentinel);
      assert (sentinel->parent == t);
      sentinel->base_time = t->base_time;
      return sentinel;
    }
    return *boost::prior(t->children->end());
  }
  
  static extended_time get_base_time_root(time_type t) {
    if (base_time_roots().empty()) {
      // create a sentinel with no children at the end
      // time_type() is a hack - without it, compiler tries to pass max_base_time as reference and gets undefined reference
      const extended_time sentinel = extended_time::construct(time_type(max_base_time));
      all_extended_times().put_only(sentinel);
      base_time_roots().insert(sentinel);
    }
    const extended_time result = extended_time::construct(t);
    const auto i = base_time_roots().lower_bound(result);
    assert(i != base_time_roots().end());
    if ((*i)->base_time == t) { return *i; }
    all_extended_times().put_before(result, *i);
    base_time_roots().emplace_hint(i, result);
    return result;
  }
  template <class... Args>
  static extended_time make_extended_time_impl(extended_time const& parent, Args&&... args) {
    require_sentinel_child(parent);
    const extended_time result = extended_time::construct(time_type(never), parent, std::forward<Args>(args)...);
    const auto i = parent->children->lower_bound(result);
    assert(i != parent->children->end());
    if (**i == *result) { return *i; }
    all_extended_times().put_before(result, *i);
    parent->children->emplace_hint(i, result);
    assert (result->parent);
    result->base_time = parent->base_time;
    return result;
  }
};

#endif
