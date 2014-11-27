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


#include <assert.h>
#include <iostream>


#include "../bounded_int_calculus.hpp"
#include "../time_steward.hpp"
typedef ptrdiff_t num_coordinates_type;

typedef time_steward_system::time_steward<time_steward_system::fields_list<int>> hack_time_steward;
typedef hack_time_steward::time_type time_type;
const time_type never = hack_time_steward::never;
const time_type min_time = hack_time_steward::min_time;
typedef time_steward_system::entity_id entity_id;

typedef int64_t space_coordinate;
typedef int64_t tile_coordinate;
const num_coordinates_type num_dimensions = 2;
using time_steward_system::optional;
using time_steward_system::none;
typedef bounded_int_calculus::polynomial_with_origin<time_type, space_coordinate, 2> poly;
typedef bounded_int_calculus::polynomial_with_origin<time_type, space_coordinate, 3> poly3;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, space_coordinate> fd_vector;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, poly> poly_fd_vector;

const int view_rad = 10;
const int tile_size_shift = 20;
const space_coordinate tile_size = 1 << tile_size_shift;
const time_type second_time = 1 << 10;
constexpr inline tile_coordinate space_to_tile_min(space_coordinate c) { return c >> tile_size_shift; }
constexpr inline tile_coordinate space_to_tile_max(space_coordinate c) { return (c + 1) >> tile_size_shift; }
constexpr inline space_coordinate tile_to_space_min(tile_coordinate c) { return c << tile_size_shift; }
constexpr inline space_coordinate tile_to_space_max(tile_coordinate c) { return tile_to_space_min(c + 1); }
class shot_trajectory {};
class shot_tile {};
class player_center_reference_tile {};
class player_next_shot_time {};
class tile_shots {};
class cave_state {};
struct player_shape {
  player_shape(poly_fd_vector center, space_coordinate radius):center(center),radius(radius){}
  
  poly_fd_vector center;
  space_coordinate radius;
};
enum wall_state {
  UNDETERMINED,
  WALL,
  EMPTY,
};
struct wall_state_traits {
  static inline wall_state null() { return UNDETERMINED; }
};
struct cave_state_traits {
  static inline space_coordinate null() { return 0; }
};



using time_steward_system::field;
typedef time_steward_system::fields_list<
  player_shape,
  field<player_center_reference_tile, fd_vector>,
  field<player_next_shot_time, time_type, hack_time_steward::time_field_traits>,
  field<shot_trajectory, optional<poly_fd_vector>>,
  field<shot_tile, fd_vector>,
  field<tile_shots, persistent_siphash_id_set>,
  field<wall_state, wall_state, wall_state_traits>,
  field<cave_state, space_coordinate, cave_state_traits>
> fields;
typedef time_steward_system::time_steward<fields> time_steward;
typedef time_steward::accessor accessor;
typedef time_steward::event event;
typedef time_steward::trigger trigger;
typedef accessor::entity_ref entity_ref;


template<class Poly>
time_type when_nonpos(time_type start, Poly p) {
  if (p(start) <= 0) { return start; }
  auto i = p.sign_interval_boundaries_upper_bound(start);
  while ((i != p.sign_interval_boundaries_end()) && (sign(p(*i)) == 1)) { ++i; }
  return (i != p.sign_interval_boundaries_end()) ? *i : never;
}
template<class Poly>
time_type when_nonneg(time_type start, Poly p) {
  if (p(start) >= 0) { return start; }
  auto i = p.sign_interval_boundaries_upper_bound(start);
  while ((i != p.sign_interval_boundaries_end()) && (sign(p(*i)) == -1)) { ++i; }
  return (i != p.sign_interval_boundaries_end()) ? *i : never;
}

void anticipate_shot_moving(time_steward::accessor* accessor, entity_ref e, fd_vector tile);
entity_ref tile_entity(time_steward::accessor* accessor, fd_vector tile, bool req_wall_state = true, bool req_cave_state = false);
class shot_enters_new_tile : public event {
public:
  shot_enters_new_tile(entity_id id, fd_vector tile) : id(id),tile(tile) {}
  entity_id id;
  fd_vector tile;

  void operator()(time_steward::accessor* accessor)const override {
    auto e = accessor->get(id);
    auto old_te = tile_entity(accessor, accessor->get<shot_tile>(e));
    auto new_te = tile_entity(accessor, tile);
    auto& s = accessor->get_mut<tile_shots>(old_te);
    s = s.erase(id);
    if (accessor->get<wall_state>(new_te) == WALL) {
      accessor->set<shot_trajectory>(e, none);
      accessor->set<wall_state>(new_te, EMPTY);
    }
    else {
      anticipate_shot_moving(accessor, e, tile);
      accessor->set<shot_tile>(e, tile);
      auto& s = accessor->get_mut<tile_shots>(new_te);
      s = s.insert(id);
    }
  }
};

entity_id tile_entity_id(fd_vector tile) {
  return siphash_id::combining('t','i','l','e',tile(0),tile(1));
}

const int64_t max_cave_radius_in_tiles = 10;

void require_wall_state(time_steward::accessor* accessor, entity_ref e, fd_vector tile);
void require_cave_state(time_steward::accessor* accessor, entity_ref e, fd_vector tile);
entity_ref tile_entity(time_steward::accessor* accessor, fd_vector tile, bool req_wall_state, bool req_cave_state) {
  auto e = accessor->get(tile_entity_id(tile));
  if (req_wall_state) { require_wall_state(accessor, e, tile); }
  if (req_cave_state) { require_cave_state(accessor, e, tile); }
  return e;
}
void require_wall_state(time_steward::accessor* accessor, entity_ref e, fd_vector tile) {
  auto w = accessor->get<wall_state>(e);
  if (w == UNDETERMINED) {
    for (tile_coordinate x = tile(0) - max_cave_radius_in_tiles; x <= tile(0) + max_cave_radius_in_tiles; ++x) {
      for (tile_coordinate y = tile(1) - max_cave_radius_in_tiles; y <= tile(1) + max_cave_radius_in_tiles; ++y) {
        auto ce = tile_entity(accessor, fd_vector(x,y), false, true);
        auto c = accessor->get<cave_state>(ce);
        if ((c > 0) && ((x - tile(0))*(x - tile(0)) + (y - tile(1))*(y - tile(1)))*tile_size*tile_size <= c*c) {
          accessor->set<wall_state>(e, EMPTY);
    //std::cerr << "b";
          return;
        }
      }
    }
    accessor->set<wall_state>(e, WALL);
    //std::cerr << "a";
  }
}
void require_cave_state(time_steward::accessor* accessor, entity_ref e, fd_vector tile) {
  auto c = accessor->get<cave_state>(e);
  if (c == 0) {
    bool any_cave_here = ((e.id().data()[0] & 255) == 0) || ((tile(0) == 0) && (tile(1) == 0));
    if (any_cave_here) {
    //std::cerr << "q";
      accessor->set<cave_state>(e, tile_size*2 + (e.id().data()[1] & ((tile_size * 8)-1)));
    }
    else {
    //std::cerr << "n";
      accessor->set<cave_state>(e, -1);
    }
  }
}
void require_view_area(time_steward::accessor* accessor, fd_vector tile) {
  for (tile_coordinate tx = tile(0)-view_rad; tx <= tile(0)+view_rad; ++tx) {
    for (tile_coordinate ty = tile(1)-view_rad; ty <= tile(1)+view_rad; ++ty) {
    //std::cerr << "r";
      require_wall_state(accessor, accessor->get(tile_entity_id(fd_vector(tx,ty))), fd_vector(tx,ty));
    }
  }
}

void anticipate_shot_moving(time_steward::accessor* accessor, entity_ref e, fd_vector tile) {
  auto trajectory = *accessor->get<shot_trajectory>(e);
  time_type best_time = never;
  fd_vector best;
  for (num_coordinates_type dim = 0; dim < num_dimensions; ++dim) {
    auto t = trajectory(dim);
    time_type when = never;
    fd_vector where = tile;
    if (t.get_term(accessor->now(), 1) > 0) {
      ++where[dim];
      when = when_nonneg(accessor->now(), t - tile_to_space_max(tile(dim)));
      assert (when != never);
    }
    if (t.get_term(accessor->now(), 1) < 0) {
      --where[dim];
      when = when_nonpos(accessor->now(), t - tile_to_space_min(tile(dim)));
      assert (when != never);
    }
    if (when != never) {
      if ((best_time == never) || (best_time > when) || ((best_time >= when) && (accessor->random_bits(1)))) {
        best_time = when;
        best = where;
      }
    }
  }
  assert (best_time != never);
  accessor->anticipate_event(best_time, std::shared_ptr<event>(new shot_enters_new_tile(e.id(), best)));
}

class player_shoots : public event {
public:
  player_shoots(entity_id id, fd_vector v) : id(id),v(v) {}
  entity_id id;
  fd_vector v;

  void operator()(time_steward::accessor* accessor)const override {
    auto player = accessor->get(id);
    player_shape const& p = *accessor->get<player_shape>(player);
    auto shot = accessor->create_entity();
    accessor->set<shot_trajectory>(shot, poly_fd_vector(
          poly(accessor->now(), poly::without_origin_t(p.center(0)(accessor->now()),v(0))),
          poly(accessor->now(), poly::without_origin_t(p.center(1)(accessor->now()),v(1)))));
    fd_vector tile = accessor->get<player_center_reference_tile>(player);
    auto te = tile_entity(accessor, tile);
    auto& s = accessor->get_mut<tile_shots>(te);
    s = s.insert(shot.id());
    accessor->set<shot_tile>(shot, tile);
    anticipate_shot_moving(accessor, shot, tile);
    accessor->set<player_next_shot_time>(player, accessor->now() + (second_time/4));
  }
};

class player_accelerates : public event {
public:
  player_accelerates(entity_id id, fd_vector v) : id(id),v(v) {}
  entity_id id;
  fd_vector v;

  void operator()(time_steward::accessor* accessor)const override {
    auto player = accessor->get(id);
    player_shape& p = *accessor->get_mut<player_shape>(player);
    p.center.set_term(accessor->now(), 1, p.center.get_term<space_coordinate>(accessor->now(), 1) + v);
  }
};

// const  freq = 1;
// {
//   poly_fd_vector new_acceleration_func_magradrad(0,0);
//   const poly_fd_vector collision_vec;
//   const fd_vector collision_dir_magrad = collision_vec.get_term<space_coordinate>(accessor->now(), 0) magnitude radius;
//   const poly x = collision_vec.dot(collision_dir_magrad) - radius*radius;
//   const poly v = x.derivative();
//   const poly directed_acc_func = -(v*freq + x*freq*freq);
//   new_acceleration_func_magradrad += collision_dir_magrad * directed_acc_func;
//   
//   
//   
//   c0.center.set_term(accessor->now(), 2, new_acceleration_func_magradrad.get_term<space_coordinate>(accessor->now(), 0) / (radius*radius));
// }

class player_strikes_tile : public event {
public:
  player_strikes_tile(entity_id id, fd_vector tile) : id(id),tile(tile) {}
  entity_id id;
  fd_vector tile;

  void operator()(time_steward::accessor* accessor)const override {
    //std::cerr << "player_strikes_tile\n";
    auto e = accessor->get(id);
    player_shape& p = *accessor->get_mut<player_shape>(e);
    auto c = p.center.get_term<space_coordinate>(accessor->now(), 0);
    auto r = p.radius;
    bool oob0 = (c(0) + r < tile_to_space_min(tile(0))) || (c(0) - r > tile_to_space_max(tile(0)));
    bool oob1 = (c(1) + r < tile_to_space_min(tile(1))) || (c(1) - r > tile_to_space_max(tile(1)));
    if (oob0 && !oob1) {
      p.center[0].set_term(accessor->now(), 1, 0);
    }
    else if (oob1 && !oob0) {
      p.center[1].set_term(accessor->now(), 1, 0);
    }
    else {
      space_coordinate mid_x = (tile_to_space_min(tile(0)) + tile_to_space_max(tile(0))) / 2;
      space_coordinate corner_x = (c(0) < mid_x) ? tile_to_space_min(tile(0)) : tile_to_space_max(tile(0));
      space_coordinate mid_y = (tile_to_space_min(tile(1)) + tile_to_space_max(tile(1))) / 2;
      space_coordinate corner_y = (c(1) < mid_y) ? tile_to_space_min(tile(1)) : tile_to_space_max(tile(1));
      fd_vector d = c - fd_vector(corner_x, corner_y);
      auto v = p.center.get_term<space_coordinate>(accessor->now(), 1);
      auto dmagsq = d(0)*d(0)+d(1)*d(1);
      auto removed = d * v.dot(d);
      removed[0] = divide(removed[0], dmagsq, rounding_strategy<round_up, negative_mirrors_positive>());
      removed[1] = divide(removed[1], dmagsq, rounding_strategy<round_up, negative_mirrors_positive>());
      p.center.set_term(accessor->now(), 1, v-removed);
    }
  }
};

class player_could_hit_walls : public trigger {
public:
  player_could_hit_walls(entity_id id) : id(id) {}
  entity_id id;

  void operator()(time_steward::accessor* accessor)const override {
    //std::cerr << "player_could_hit_walls\n";
    auto e = accessor->get(id);
    fd_vector ct = accessor->get<player_center_reference_tile>(e);
    player_shape const& p = *accessor->get<player_shape>(e);
    
    for (tile_coordinate tx = ct(0)-2; tx <= ct(0)+2; ++tx) {
      for (tile_coordinate ty = ct(1)-2; ty <= ct(1)+2; ++ty) {
        auto tile = fd_vector(tx,ty);
        auto te = tile_entity(accessor, tile);
        if (accessor->get<wall_state>(te) == WALL) {
          for (space_coordinate corner_x = tile_to_space_min(tx); corner_x <= tile_to_space_max(tx); corner_x += tile_size) {
            for (space_coordinate corner_y = tile_to_space_min(ty); corner_y <= tile_to_space_max(ty); corner_y += tile_size) {
              poly3 distish = (p.center(0)-corner_x)*(p.center(0)-corner_x) + (p.center(1)-corner_y)*(p.center(1)-corner_y) - p.radius*p.radius;
              assert(distish.get_term(accessor->now(), 0) > 0);
              time_type when = when_nonpos(accessor->now(), distish);
              if (when != never) {
                assert (when > accessor->now());
                accessor->anticipate_event(when-1, std::shared_ptr<event>(new player_strikes_tile(e.id(), tile)));
              }
            }
          }
          for (num_coordinates_type dim0 = 0; dim0 < num_dimensions; ++dim0) {
            num_coordinates_type dim1 = !dim0;
            tile_coordinate t0 = dim0 ? ty : tx;
            tile_coordinate t1 = dim0 ? tx : ty;
            time_type when = never;
            if (p.center(dim0).get_term(accessor->now(), 1) > 0) {
              when = when_nonneg(accessor->now(), p.center(dim0)+p.radius - tile_to_space_min(t0));
            }
            if (p.center(dim0).get_term(accessor->now(), 1) < 0) {
              when = when_nonpos(accessor->now(), p.center(dim0)-p.radius - tile_to_space_max(t0));
            }
            if (when != never) {
              if (when > accessor->now()) {
                space_coordinate s = p.center(dim1)(when);
                if (s >= tile_to_space_min(t1) && s <= tile_to_space_max(t1)) {
                  accessor->anticipate_event(when-1, std::shared_ptr<event>(new player_strikes_tile(e.id(), tile)));
                }
              }
            }
          }
        }
      }
    }
  }
};

class player_enters_new_tile : public event {
public:
  player_enters_new_tile(entity_id id, fd_vector tile) : id(id),tile(tile) {}
  entity_id id;
  fd_vector tile;

  void operator()(time_steward::accessor* accessor)const override {
    //std::cerr << "player_enters_new_tile\n";
    accessor->set<player_center_reference_tile>(accessor->get(id), tile);
    require_view_area(accessor, tile);
  }
};

class player_moves_around : public trigger {
public:
  player_moves_around(entity_id id) : id(id) {}
  entity_id id;

  void operator()(time_steward::accessor* accessor)const override {
    //std::cerr << "player_moves_around\n";
    auto e = accessor->get(id);
    fd_vector ct = accessor->get<player_center_reference_tile>(e);
    player_shape const& p = *accessor->get<player_shape>(e);
    
    for (num_coordinates_type dim = 0; dim < num_dimensions; ++dim) {
      fd_vector where = ct;
      auto t = p.center(dim);
      //std::cerr << accessor->now() << ": " << tile_to_space_min(ct(dim)) << ", " << t << ", " << tile_to_space_max(ct(dim)) << "\n";
      if (t.get_term(accessor->now(), 1) > 0) {
        ++where[dim];
        accessor->anticipate_event(when_nonneg(accessor->now(), t - tile_to_space_max(ct(dim))), std::shared_ptr<event>(new player_enters_new_tile(e.id(), where)));
        //std::cerr << "c " << when_nonneg(accessor->now(), t - tile_to_space_max(ct(dim))) << "\n";
      }
      if (t.get_term(accessor->now(), 1) < 0) {
        --where[dim];
        accessor->anticipate_event(when_nonpos(accessor->now(), t - tile_to_space_min(ct(dim))), std::shared_ptr<event>(new player_enters_new_tile(e.id(), where)));
        //std::cerr << "d " << when_nonpos(accessor->now(), t - tile_to_space_min(ct(dim))) << "\n";
      }
    }
  }
};

class initialize_world : public event {
public:
  void operator()(time_steward::accessor* accessor)const override {
    entity_ref player = accessor->get(time_steward_system::global_object_id);
    accessor->set<player_center_reference_tile>(player, fd_vector(0,0));
    require_view_area(accessor, fd_vector(0,0));
    accessor->set<player_shape>(player, player_shape(
        poly_fd_vector(
          poly(0, poly::without_origin_t(0,0)),
          poly(0, poly::without_origin_t(0,0))),
        tile_size/3
    ));
    accessor->set_trigger(player.id(), std::shared_ptr<trigger>(new player_could_hit_walls(player.id())));
    accessor->set_trigger(siphash_id::combining(player.id(), 2), std::shared_ptr<trigger>(new player_moves_around(player.id())));
  }
};

class history_tree {
public:
  struct fiat_event {
    fiat_event(uint64_t distinguisher, std::shared_ptr<event> e):distinguisher(distinguisher),e(e){}
    uint64_t distinguisher;
    std::shared_ptr<event> e;
  };
  struct node {
    //node* parent;
    node(time_type start_time):start_time(start_time){}
    time_type start_time;
    std::multimap<time_type, fiat_event> events;
    std::multimap<time_type, node> children;
  };
  time_steward steward;
  node history_root;
  typedef std::vector<node*> history;
  history current_history;
  time_type current_time;
  size_t place_in_current_history;
  
  history_tree():
    history_root(-1),
    current_time(-1),
    place_in_current_history(0)
  {
    current_history.push_back(&history_root);
    validate();
  }
  
  struct spatial_representation_entry {
    spatial_representation_entry(double height, history h):height(height),h(h){}
    double height;
    history h;
  };
  struct spatial_representation_column {
    time_type time;
    std::vector<spatial_representation_entry> entries;
  };
  struct spatial_representation {
    std::vector<spatial_representation_column> columns;
  };
  struct queue_entry {
    time_type time;
    bool starting;
    history h;
    queue_entry(history h, bool starting):starting(starting),h(h){
      if (starting) {
        time = h.back()->start_time;
      }
      else if (h.back()->events.empty()) {
        time = 1;
      }
      else {
        time = boost::prior(h.back()->events.end())->first;
        if (time < 1) { time = 1; }
      }
    }
    bool operator<(queue_entry const& other)const { return time > other.time; }
  };
  spatial_representation spatialize() {
    spatial_representation result;
    std::priority_queue<queue_entry> q; //node*, std::vector<node*>, sort_nodes_desc
    std::vector<history> current;
    history bh;
    bh.push_back(&history_root);
    q.push(queue_entry(bh, true));
      
    while (!q.empty()) {
      queue_entry t = q.top(); q.pop();
      
      spatial_representation_column col;
      col.time = t.time;
      
      if (t.starting) {
        if (current.empty()) {
          current.push_back(t.h);
          col.entries.push_back(spatial_representation_entry(0.5, t.h));
        }
        else {
          double inc = double(1) / double(current.size() + 1);
          double height = 0;
          for (size_t i = 0; i < current.size(); ++i) {
            if (current[i].back() == t.h[t.h.size()-2]) {
              current.insert(current.begin()+i+1, t.h);
            }
            if (current[i].back() != t.h.back()) {
              height += inc;
            }
            col.entries.push_back(spatial_representation_entry(height, current[i]));
          }
        }
        
        for (auto& p : t.h.back()->children) {
          history h = t.h;
          h.push_back(&p.second);
          q.push(queue_entry(h, true));
        }
        q.push(queue_entry(t.h, false));
      }
      else {
        double inc = double(1) / double(current.size() + 1);
        double height = 0;
        for (size_t i = 0; i < current.size(); ++i) {
          height += inc;
          col.entries.push_back(spatial_representation_entry(height, current[i]));
        }
        for (size_t i = 0; i < current.size(); ++i) {
          if (current[i].back() == t.h.back()) {
            current.erase(current.begin()+i);
          }
        }
      }
      result.columns.push_back(col);
    }
    assert (current.empty());
    return result;
  }
  
  std::unique_ptr<accessor> accessor_after(time_type time) {
    if (time > current_time) {
      set_time(time);
    }
    return steward.accessor_after(time);
  }
  
  void insert_fiat_event(time_type time, uint64_t distinguisher, std::shared_ptr<event> e) {
    if (time-1 < current_time) {
      set_time(time-1);
    }
    node* cur_node = current_history[place_in_current_history];
    if ((place_in_current_history+1 >= current_history.size()) && (cur_node->events.empty() || boost::prior(cur_node->events.end())->first <= time)) {
      cur_node->events.insert(std::make_pair(time, fiat_event(distinguisher, e)));
    }
    else {
      while (place_in_current_history+1 < current_history.size()) { current_history.pop_back(); }
      auto p = cur_node->children.insert(std::make_pair(time, node(time)));
      node* new_node = &p->second;
      current_history.push_back(new_node);
      new_node->events.insert(std::make_pair(time, fiat_event(distinguisher, e)));
    }
  }
  
  void validate_node(node const& n) {
    for (auto const& p : n.events) {
      assert (p.first >= n.start_time);
    }
    for (auto const& p : n.children) {
      assert (p.first > n.start_time);
      assert (p.second.start_time == p.first);
      validate_node(p.second);
    }
  }
  void validate() {
    assert (place_in_current_history < current_history.size());
    assert (current_history[place_in_current_history]->start_time <= current_time);
    assert ((place_in_current_history+1 >= current_history.size()) || current_history[place_in_current_history+1]->start_time > current_time);
    node* n = &history_root;
    for (size_t i = 1; i < current_history.size(); ++i) {
      assert(n->children.find(current_history[i]->start_time) != n->children.end());
      n = current_history[i];
    }
    validate_node(history_root);
  }
  
  void set_time(time_type new_time) {
    validate();
    if (new_time > current_time) {
      while (true) {
        node* cur_node = current_history[place_in_current_history];
        node* next_node = (place_in_current_history+1 < current_history.size()) ? current_history[place_in_current_history+1] : nullptr;
        for (auto i = cur_node->events.upper_bound(current_time);
            i != cur_node->events.end() && (i->first <= new_time) &&
            ((!next_node) || (i->first <= next_node->start_time)); ++i) {
          steward.insert_fiat_event(i->first, i->second.distinguisher, i->second.e);
        }
        if (next_node && new_time >= next_node->start_time) {
          ++place_in_current_history;
        }
        else {
          break;
        }
      }
    }
    if (new_time < current_time) {
      while (true) {
        node* cur_node = current_history[place_in_current_history];
        node* next_node = (place_in_current_history+1 < current_history.size()) ? current_history[place_in_current_history+1] : nullptr;
        for (auto i = cur_node->events.upper_bound(new_time);
            i != cur_node->events.end() && (i->first <= current_time) &&
            ((!next_node) || (i->first <= next_node->start_time)); ++i) {
          steward.erase_fiat_event(i->first, i->second.distinguisher);
        }
        if (place_in_current_history > 0 && new_time < cur_node->start_time) {
          --place_in_current_history;
        }
        else {
          break;
        }
      }
    }
    current_time = new_time;
    validate();
  }
  void set_history(history const& new_history) {
    validate();
    size_t first_difference = 0;
    while (current_history.size() < first_difference ||
               new_history.size() < first_difference ||
               current_history[first_difference+1] == new_history[first_difference+1]) {
      ++first_difference;
    }
    if (current_history.size() < first_difference &&
            new_history.size() < first_difference) { return; }
    const time_type s0 = (current_history.size() < first_difference) ? min_time : current_history[first_difference]->start_time;
    const time_type s1 = (    new_history.size() < first_difference) ? min_time :     new_history[first_difference]->start_time;
    const time_type branch_point = std::max(s0, s1);
    if (branch_point-1 < current_time) {
      set_time(branch_point-1);
    }
    /*while (current_history->start_time >= branch_point) {
      current_history.pop_back();
    }
    
    size_t next_difference = first_difference;
    while (new_history.size() < next_difference) {
      current_history.push_back(new_history[next_difference++]);
    }*/
    current_history = new_history;
    validate();
  }
  
  /*history_tree_node* inc_sibling(history_tree_node* old_history, bool direction)const {
    if (!old_history->parent) { return nullptr; }
    auto i = old_history->parent->children.find(old_history);
    assert (i != old_history->parent->children.end());
    if (direction) {
      ++i;
      if (i == old_history->parent->children.end()) { return nullptr; }
      return *i;
    }
    else {
      if (i == old_history->parent->children.begin()) { return nullptr; }
      --i;
      return *i;
    }
  }
  history inc_history(history_tree_node* old_history, time_type switch_time, bool direction)const {
    history_tree_node* new_history = old_history;
    while (new_history->start_time >= switch_time) {
      new_history = new_history->parent;
    }
    while (true) {
      history_tree_node* next_history = inc_sibling(new_history, direction);
      if (next_history) { return next_history; }
      new_history = new_history->parent;
      if (!new_history) { return old_history; }
    }
  }*/
};

time_type shot_tail_delay = second_time/10;
template<class DrawFuncsType>
void draw_green_caves(history_tree& w, time_type time, DrawFuncsType& draw) {
  std::unique_ptr<time_steward::accessor> accessor = w.accessor_after(time);
  auto player = accessor->get(time_steward_system::global_object_id);
  fd_vector ct = accessor->get<player_center_reference_tile>(player);
  player_shape const& p = *accessor->get<player_shape>(player);
  
  const space_coordinate cx = p.center(0)(accessor->now());
  const space_coordinate cy = p.center(1)(accessor->now());
  draw.circle(0, 0, p.radius);
  
  for (tile_coordinate tx = ct(0)-view_rad; tx <= ct(0)+view_rad; ++tx) {
    for (tile_coordinate ty = ct(1)-view_rad; ty <= ct(1)+view_rad; ++ty) {
      auto tile = fd_vector(tx,ty);
      auto te = accessor->get(tile_entity_id(tile));
      
      if (accessor->get<wall_state>(te) == WALL) {
        draw.rect(
          tile_to_space_min(tx) - cx,
          tile_to_space_min(ty) - cy,
          tile_to_space_max(tx) - cx,
          tile_to_space_max(ty) - cy
        );
      }
      for (entity_id shot : accessor->get<tile_shots>(te)) {
        assert(accessor->get<shot_tile>(accessor->get(shot)) == tile);
        auto trajectory = *accessor->get<shot_trajectory>(accessor->get(shot));
        draw.segment(
          trajectory(0)(accessor->now()) - cx,
          trajectory(1)(accessor->now()) - cy,
          trajectory(0)(accessor->now() - shot_tail_delay) - cx,
          trajectory(1)(accessor->now() - shot_tail_delay) - cy
        );
      }
    }
  }
  
  history_tree::spatial_representation sprep = w.spatialize();
  const double max = sprep.columns.back().time;
  const double vqqq = tile_size*view_rad;
  size_t place_in_current_history = 0;
  for (size_t i = 1; i < sprep.columns.size(); ++i) {
    auto& cur = sprep.columns[i];
    auto& prev = sprep.columns[i-1];
    for (auto e : cur.entries) {
      for (auto f : prev.entries) {
        if (f.h.back() == e.h.back()) {
          double x0 = (double(prev.time * vqqq * 2) / max) - vqqq;
          double y0 = f.height * vqqq * 0.3 + vqqq * 0.7;
          double x1 = (double(cur.time * vqqq * 2) / max) - vqqq;
          double y1 = e.height * vqqq * 0.3 + vqqq * 0.7;
          draw.segment(x0, y0, x1, y1);
          if ((place_in_current_history+1 < w.current_history.size()) && w.current_history[place_in_current_history+1]->start_time == prev.time) {
            ++place_in_current_history;
          }
          if (w.current_history[place_in_current_history] == e.h.back()) {
            draw.segment(x1, y1, x0, y0);
            if (prev.time <= time && time < cur.time) {
              draw.circle(
                x0 + (x1-x0)*double(time-prev.time)/double(cur.time-prev.time),
                y0 + (y1-y0)*double(time-prev.time)/double(cur.time-prev.time),
                vqqq * 0.02);
            }
          }
        }
      }
    }
  }
}
