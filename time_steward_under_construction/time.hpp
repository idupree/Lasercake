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
  typedef typename manual_orderer<extended_time_metadata>::manually_orderable extended_time;
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
    base_time(base_time),parent(),id(),trigger_iteration(not_a_trigger),claims(0){}
  extended_time_metadata(time_type base_time, extended_time parent, siphash_id id):
    base_time(base_time),parent(parent),id(id),trigger_iteration(not_a_trigger),claims(0){}
  extended_time_metadata(time_type base_time, extended_time parent, siphash_id id, uint32_t trigger_iteration):
    base_time(base_time),parent(parent),id(id),trigger_iteration(trigger_iteration),claims(0){}
  
  struct sort_sibling_extended_times_by_id_and_trigger_iteration {
    bool operator()(extended_time a, extended_time b)const {
      if (a->trigger_iteration != b->trigger_iteration) { return a->trigger_iteration < b->trigger_iteration; }
      return a->id < b->id;
    }
  };
  
  // TODO: can we optimize for size in common cases?
  time_type base_time;
  extended_time parent;
  siphash_id id; // unused for base_time_roots
  uint32_t trigger_iteration; // constant not_a_trigger for all non-triggers
  uint32_t claims;
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
  static extended_time make_event_time(time_type t, siphash_id id) {
    return make_extended_time_impl(get_base_time_root(t), id);
  }
  static extended_time make_event_time(extended_time t, siphash_id id) {
    return make_extended_time_impl(t->is_trigger() ? t->parent : t, id);
  }
  static extended_time make_trigger_call_time(trigger_id tid, extended_time when_triggered) {
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
  static extended_time get_first_after(time_type t) {
    const extended_time hack = extended_time::construct(t);
    const auto i = base_time_roots().upper_bound(hack);
    if (i == base_time_roots().end()) {
      return max_time();
    }
    return *i;
  }

  static bool is_after_all_calls_triggered_at(extended_time t, extended_time when_triggered) {
    extended_time trigger_parent = when_triggered->is_trigger() ? when_triggered->parent : when_triggered;
    if (trigger_parent->children) {
      extended_time sentinel = *boost::prior(trigger_parent->children->end());
      return t > sentinel;
    }
    else {
      return t > trigger_parent;
    }
  }
  
  static void claim(extended_time t) {
    ++t->claims;
  }
  static void release(extended_time t) {
    // TODO delete unused base_time_roots (also base_time_roots are a hack, can we remove them?)
    caller_error_if(t->claims == 0, "released an extended_time too many times");
    --t->claims;
    if (t->claims == 0 && t->parent) {
      t->delete_from_parent();
      all_extended_times().erase(t);
      t.destroy();
    }
  }
  class scoped {
  public:
    scoped(extended_time data):data(data){}
    ~scoped(){ release(data); }
    extended_time get()const { return data; }
  private:
    extended_time data;
  };
private:
  struct sort_extended_times_by_base_time {
    bool operator()(extended_time a, extended_time b)const {
      return a->base_time < b->base_time;
    }
  };
  static manual_orderer<extended_time_metadata<TimeTypeInfo>>& all_extended_times() {
    static manual_orderer<extended_time_metadata<TimeTypeInfo>> a; return a;
  }
  static std::set<extended_time, sort_extended_times_by_base_time>& base_time_roots() {
    static std::set<extended_time, sort_extended_times_by_base_time> a; return a;
  }
  static extended_time require_sentinel_child(extended_time t) {
    if (!t->children) {
      // create a sentinel with no children at the end
      t->children.reset(new typename extended_time_metadata<TimeTypeInfo>::children_set());
      const extended_time sentinel = extended_time::construct(time_type(never), t, siphash_id::greatest());
      all_extended_times().insert_adjacent(sentinel, t, true);
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
      all_extended_times().insert_only(sentinel);
      base_time_roots().insert(sentinel);
    }
    const extended_time result = extended_time::construct(t);
    const auto i = base_time_roots().lower_bound(result);
    assert(i != base_time_roots().end());
    if ((*i)->base_time == t) { return *i; }
    all_extended_times().insert_adjacent(result, *i, false);
    base_time_roots().emplace_hint(i, result);
    return result;
  }
  template <class... Args>
  static extended_time make_extended_time_impl(extended_time parent, Args&&... args) {
    require_sentinel_child(parent);
    extended_time result = extended_time::construct(time_type(never), parent, std::forward<Args>(args)...);
    const auto i = parent->children->lower_bound(result);
    assert(i != parent->children->end());
    if (**i == *result) {
      result = *i;
    }
    else {
      all_extended_times().insert_adjacent(result, *i, false);
      parent->children->emplace_hint(i, result);
      assert (result->parent);
      result->base_time = parent->base_time;
    }
    ++result->claims;
    return result;
  }
};

#endif
