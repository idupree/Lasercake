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
#include "../history_tree.hpp"
typedef ptrdiff_t num_coordinates_type;

typedef time_steward_system::time_steward<time_steward_system::fields_list<int>> hack_time_steward;
typedef hack_time_steward::time_type time_type;
const time_type never = hack_time_steward::never;
const time_type min_time = hack_time_steward::min_time;
const time_type max_time = hack_time_steward::max_time;
typedef time_steward_system::entity_id entity_id;

typedef int64_t space_coordinate;
typedef int64_t tile_coordinate;
const num_coordinates_type num_dimensions = 2;
using time_steward_system::optional;
using time_steward_system::none;
typedef bounded_int_calculus::polynomial_with_origin<time_type, space_coordinate, 2> poly;
typedef bounded_int_calculus::polynomial_with_origin<time_type, space_coordinate, 3> poly3;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, space_coordinate> fd_vector;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, double> double_vector;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, poly> poly_fd_vector;

const int view_rad = 10;
const int tile_size_shift = 20;
const space_coordinate tile_size = 1 << tile_size_shift;
const time_type second_time = 1 << 10;
constexpr inline tile_coordinate space_to_tile_min(space_coordinate c) { return c >> tile_size_shift; }
constexpr inline tile_coordinate space_to_tile_max(space_coordinate c) { return (c + 1) >> tile_size_shift; }
constexpr inline space_coordinate tile_to_space_min(tile_coordinate c) { return c * tile_size; }
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

typedef history_tree<time_steward, second_time> gc_history_tree;


struct view_rect {
  fd_vector screen_min;
  fd_vector screen_size;
  double_vector to_screen(double_vector loc) {
    double_vector r(double(screen_min(0)), double(screen_min(1)));
    r[0] += screen_size(0) * loc(0);
    r[1] += screen_size(1) * loc(1);
    return r;
  }
  double_vector from_screen(fd_vector loc) {
    fd_vector r0 = loc - screen_min;
    double_vector r(double(r0(0)), double(r0(1)));
    r[0] /= screen_size(0);
    r[1] /= screen_size(1);
    return r;
  }
};
const double vqqq = tile_size*view_rad*2;

struct draw_green_caves_metadata {
  view_rect main_view;
  view_rect hist_view;
  num_coordinates_type hist_time_dim;
  time_type focus_time;
  gc_history_tree* w;
  double_vector main_to_screen(space_coordinate x, space_coordinate y) {
    double_vector r;
    r[0] = 0.5 + double(x) / vqqq;
    r[1] = 0.5 + double(y) / vqqq;
    return main_view.to_screen(r);
  }
  std::pair<time_type, gc_history_tree::history> hist_from_screen(fd_vector loc) {
    double_vector r = hist_view.from_screen(loc);
    if ((r(0) < 0) || (r(0) > 1) || (r(1) < 0) || (r(1) > 1)) {
      return std::pair<time_type, gc_history_tree::history>(never, gc_history_tree::history());
    }
    return w->hist_from_draw_coords(r(hist_time_dim), r(!hist_time_dim), focus_time);
  }
  draw_green_caves_metadata() {}
  draw_green_caves_metadata(fd_vector screen_size) {
    hist_time_dim = (screen_size(1) > screen_size(0)) ? 0 : 1;
    main_view.screen_min = fd_vector(0,0);
    hist_view.screen_size[hist_time_dim] = screen_size(hist_time_dim);
    hist_view.screen_size[!hist_time_dim] = screen_size(hist_time_dim) / 3;
    hist_view.screen_min[hist_time_dim] = 0;
    int64_t remaining = screen_size(!hist_time_dim) - hist_view.screen_size(!hist_time_dim);
    hist_view.screen_min[!hist_time_dim] = remaining;
    int64_t mss = std::min(remaining, screen_size(hist_time_dim));
    main_view.screen_size = fd_vector(mss,mss);
  }
};

time_type shot_tail_delay = second_time/10;
template<class DrawFuncsType>
struct hist_line_func {
  draw_green_caves_metadata* m;
  DrawFuncsType* draw;
  void operator()(double x0, double y0, double x1, double y1, bool in_current_history) {
    double_vector v0;
    double_vector v1;
    v0[m->hist_time_dim] = x0;
    v0[!m->hist_time_dim] = y0;
    v1[m->hist_time_dim] = x1;
    v1[!m->hist_time_dim] = y1;
    v0 = m->hist_view.to_screen(v0);
    v1 = m->hist_view.to_screen(v1);
    draw->segment(v0(0), v0(1), v1(0), v1(1), in_current_history ? 4.0 : 1.5);
  }
};
template<class DrawFuncsType>
draw_green_caves_metadata draw_green_caves(fd_vector screen_size, gc_history_tree& w, time_type time, DrawFuncsType& draw) {
  draw_green_caves_metadata metadata(screen_size);
  metadata.focus_time = time;
  metadata.w = &w;
  std::unique_ptr<time_steward::accessor> accessor = w.accessor_after(time);
  auto player = accessor->get(time_steward_system::global_object_id);
  fd_vector ct = accessor->get<player_center_reference_tile>(player);
  player_shape const& p = *accessor->get<player_shape>(player);
  
  const space_coordinate cx = p.center(0)(accessor->now());
  const space_coordinate cy = p.center(1)(accessor->now());
  double_vector pv = metadata.main_to_screen(0, 0);
  draw.circle(pv(0), pv(1), double(p.radius) * double(metadata.main_view.screen_size(0)) / vqqq);
  
  for (tile_coordinate tx = ct(0)-view_rad; tx <= ct(0)+view_rad; ++tx) {
    for (tile_coordinate ty = ct(1)-view_rad; ty <= ct(1)+view_rad; ++ty) {
      auto tile = fd_vector(tx,ty);
      auto te = accessor->get(tile_entity_id(tile));
      
      if (accessor->get<wall_state>(te) == WALL) {
        double_vector v0 = metadata.main_to_screen(tile_to_space_min(tx) - cx, tile_to_space_min(ty) - cy);
        double_vector v1 = metadata.main_to_screen(tile_to_space_max(tx) - cx, tile_to_space_max(ty) - cy);
        draw.rect(v0(0), v0(1), v1(0), v1(1));
      }
      for (entity_id shot : accessor->get<tile_shots>(te)) {
        assert(accessor->get<shot_tile>(accessor->get(shot)) == tile);
        auto trajectory = *accessor->get<shot_trajectory>(accessor->get(shot));
        double_vector v0 = metadata.main_to_screen(trajectory(0)(accessor->now()) - cx, trajectory(1)(accessor->now()) - cy);
        double_vector v1 = metadata.main_to_screen(trajectory(0)(accessor->now() - shot_tail_delay) - cx, trajectory(1)(accessor->now() - shot_tail_delay) - cy);
        draw.segment(v0(0), v0(1), v1(0), v1(1), 1.5);
      }
    }
  }
  
  hist_line_func<DrawFuncsType> line;
  line.m = &metadata;
  line.draw = &draw;
  w.draw_tree(line, time);
  
  double_vector current_hist_pos;
  current_hist_pos[metadata.hist_time_dim] = w.time_coord(time, time);
  current_hist_pos[!metadata.hist_time_dim] = w.height_coord(w.current_history, time, time);
  current_hist_pos = metadata.hist_view.to_screen(current_hist_pos);
  draw.circle(current_hist_pos(0), current_hist_pos(1), 8);
              
  return metadata;
}

const int64_t acc_updates_per_second = 50;
const space_coordinate acc = tile_size*30/(second_time*acc_updates_per_second);

struct green_caves_ui_backend {
  fd_vector screen_size;
  draw_green_caves_metadata last_metadata;

  gc_history_tree hist;
  time_type current_time = 0;
  int64_t last_milliseconds = 0;
  bool shooting = false;
  int mouse_x;
  int mouse_y;
  bool left = false;
  bool right = false;
  bool up = false;
  bool down = false;
  green_caves_ui_backend() {
    hist.insert_fiat_event(0, 0, std::shared_ptr<event>(new initialize_world()));
  }
  void update_to_real_time(int64_t milliseconds) {
    const int64_t dur = milliseconds-last_milliseconds;
    last_milliseconds = milliseconds;
    if (dur > 0 && dur < 400) {
      const time_type new_time = current_time + second_time * dur / 1000;
      for (int64_t i = current_time*acc_updates_per_second/second_time; i < new_time*acc_updates_per_second/second_time; ++i) {
        const time_type ut = i * second_time / acc_updates_per_second;
        if (   up && ! down) { hist.insert_fiat_event(ut, 1, std::shared_ptr<event>(new player_accelerates(time_steward_system::global_object_id, fd_vector(0, acc)))); }
        if ( down && !   up) { hist.insert_fiat_event(ut, 2, std::shared_ptr<event>(new player_accelerates(time_steward_system::global_object_id, fd_vector(0, -acc)))); }
        if ( left && !right) { hist.insert_fiat_event(ut, 3, std::shared_ptr<event>(new player_accelerates(time_steward_system::global_object_id, fd_vector(-acc, 0)))); }
        if (right && ! left) { hist.insert_fiat_event(ut, 4, std::shared_ptr<event>(new player_accelerates(time_steward_system::global_object_id, fd_vector(acc, 0)))); }
      }
      while (shooting) {
        std::unique_ptr<time_steward::accessor> accessor = hist.accessor_after(new_time);
        time_type st = accessor->get<player_next_shot_time>(accessor->get(time_steward_system::global_object_id));
        if (st < current_time) { st = current_time; }
        if (st <= new_time) {
          double_vector v0 = last_metadata.main_view.from_screen(fd_vector(mouse_x, screen_size(1)-mouse_y)) - double_vector(0.5, 0.5);
          //std::cerr << v0;
          fd_vector v(space_coordinate(v0(0)*100000), space_coordinate(v0(1)*100000));
          space_coordinate mag = isqrt((v(0) * v(0)) + (v(1) * v(1)));
          if (mag > 0) {
            v[0] = divide(v[0] * tile_size*5, second_time*mag, rounding_strategy<round_up, negative_mirrors_positive>());
            v[1] = divide(v[1] * tile_size*5, second_time*mag, rounding_strategy<round_up, negative_mirrors_positive>());
            hist.insert_fiat_event(st, 5, std::shared_ptr<event>(new player_shoots(time_steward_system::global_object_id, v)));
          }
        }
        else {
          break;
        }
      }
      hist.expand_to_time(new_time);
      current_time = new_time;
    }
  }
  void mouse_down(int x, int y) {
    mouse_x = x;
    mouse_y = y;
    auto h = last_metadata.hist_from_screen(fd_vector(mouse_x,screen_size(1)-mouse_y));
    if (!h.second.empty()) {
      bool dif = true;
      for (size_t i = 1; i < hist.current_history.size(); ++i) {
        if (hist.current_history[i]->start_time > h.first) {
          if (hist.current_history[i-1] == h.second.back()) {
            dif = false;
          }
          break;
        }
      }
      if (dif) {
        hist.set_history(h.second);
      }
      current_time = h.first;
      shooting = false;
    }
    else {
      shooting = true;
    }
  }
  void mouse_up(int x, int y) {
    mouse_x = x;
    mouse_y = y;
    shooting = false;
  }
  void mouse_moves(int x, int y) {
    mouse_x = x;
    mouse_y = y;
  }
  template<class DrawFuncsType>
  void draw(DrawFuncsType& draw) {
    last_metadata = draw_green_caves(screen_size, hist, current_time, draw);
  }
};

