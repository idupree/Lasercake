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



#include <vector>
#include <boost/optional.hpp>

#include <memory>
#include <vector>
#include <assert.h>

#include <iostream>

#include <cmath>
#include <limits>
#include <limits.h>


#include "../bounded_int_calculus.hpp"
#include "../bbox_collision_detector.hpp"

typedef int64_t space_coordinate;
const num_coordinates_type num_dimensions = 2;
using time_steward_system::optional;
using time_steward_system::none;
typedef time_steward_system::time_steward<time_steward_system::fields_list<int>> hack_time_steward;
typedef hack_time_steward::time_type time_type;
const time_type never = hack_time_steward::never;
typedef bounded_int_calculus::polynomial_with_origin<time_type, space_coordinate> poly;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, space_coordinate> fd_vector;
typedef bounded_int_calculus::finite_dimensional_vector<num_dimensions, poly> poly_fd_vector;
typedef time_steward_system::entity_id entity_id;

uint64_t space_to_bbox(space_coordinate const& input) {
  return uint64_t(int64_t(input)) + 0x0101010101010101ULL;
}
space_coordinate bbox_to_space(uint64_t input) {
  return space_coordinate(int64_t(input - 0x0101010101010101ULL));
}

const int64_t position_units_per_zolttt_per_square_gloppp = 1LL << 6; // plenty to avoid rounding error
const int64_t time_units_per_gloppp = 1LL << 6;
const int64_t position_units_per_zolttt_per_gloppp = position_units_per_zolttt_per_square_gloppp * time_units_per_gloppp;
const int64_t position_units_per_zolttt = position_units_per_zolttt_per_gloppp * time_units_per_gloppp;
const int64_t arena_width = position_units_per_zolttt << 10;
const int64_t view_width = arena_width >> 3;

struct circle_shape {
  circle_shape(poly_fd_vector center, poly radius):center(center),radius(radius){}
  
  poly_fd_vector center;
  poly radius;
};

struct circle_overlaps {
  persistent_set<entity_id> overlaps;
};

poly distish(circle_shape const& c1, circle_shape const& c2) {
  const auto d = c1.center - c2.center;
  const auto d_mag_sq = d.dot(d);
  const auto radsum = c1.radius + c2.radius;
  const auto r_sq = (radsum*radsum);
  return d_mag_sq - r_sq;
}

struct global_data {
  global_data(entity_id bbcd_id):bbcd_id(bbcd_id){}
  entity_id bbcd_id;
};

typedef bbox_collision_detector_system<64, num_dimensions> bbcd_system;
typedef time_steward_system::fields_list<global_data, circle_shape, circle_overlaps, bbcd_system::fields> fields;
typedef time_steward_system::time_steward<fields> time_steward;
typedef time_steward::accessor accessor;
typedef time_steward::event event;
typedef accessor::entity_ref entity_ref;

  
void update_acceleration(accessor* accessor, entity_ref e) {
  fd_vector new_acceleration(0,0);
  circle_shape& c0 = *accessor->get_mut<circle_shape>(e);
  for (entity_id other_id : accessor->get<circle_overlaps>(e)->overlaps) {
    auto other = accessor->get(other_id);
    circle_shape const& c1 = *accessor->get<circle_shape>(other);
    const fd_vector diff = (c0.center - c1.center).get_term<space_coordinate>(accessor->now(), 0);
    const space_coordinate magsq = diff.dot(diff);
    const space_coordinate radsum = c0.radius(accessor->now()) + c1.radius(accessor->now());
    if (magsq - (radsum*radsum) < 0) { // not necessarily <0, if there are two interacts at the same time
      // simple hack!
      const space_coordinate mag = isqrt(magsq);
      if (mag != 0) {
        new_acceleration += fd_vector(
          diff(0)*position_units_per_zolttt_per_square_gloppp*1000/mag,
          diff(1)*position_units_per_zolttt_per_square_gloppp*1000/mag);
      }
    }
  }
  c0.center.set_term(accessor->now(), 2, new_acceleration);
}
  

class circles_interact : public event {
public:
  circles_interact(entity_id id0, entity_id id1) : id0(id0),id1(id1) {}
  entity_id id0, id1;

  void operator()(time_steward::accessor* accessor)const override {
    auto e0 = accessor->get(id0);
    auto e1 = accessor->get(id1);
    auto c0 = accessor->get_mut<circle_shape>(e0);
    auto c1 = accessor->get_mut<circle_shape>(e1);
    const space_coordinate d = distish(*c0,*c1)(accessor->now());
    
    const fd_vector diff = (c0->center - c1->center).get_term<space_coordinate>(accessor->now(), 0);
    const space_coordinate magsq = diff.dot(diff);
    const space_coordinate radsum = c0->radius(accessor->now()) + c1->radius(accessor->now());
    assert ((d < 0) == (magsq - (radsum*radsum) < 0));
    
    persistent_set<entity_id>& c0_overlaps = accessor->get_mut<circle_overlaps>(e0)->overlaps;
    persistent_set<entity_id>& c1_overlaps = accessor->get_mut<circle_overlaps>(e1)->overlaps;
    if (d < 0) {
      c0_overlaps.insert(id1);
      c1_overlaps.insert(id0);
    }
    else {
      c0_overlaps.erase(id1);
      c1_overlaps.erase(id0);
    }
    update_acceleration(accessor, e0);
    update_acceleration(accessor, e1);
  }
};

class bbcd_funcs {
public:
  typedef bbcd_system::bounding_box bounding_box;
  typedef bbcd_system::coordinate_array coordinate_array;
  bounding_box bbox(accessor const* accessor, entity_ref e) {
    auto ci = accessor->get<circle_shape>(e);
    coordinate_array min;
    coordinate_array max;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      space_coordinate c = ci->center(i)(accessor->now());
      space_coordinate r = ci->radius(accessor->now());
      min[i] = space_to_bbox(c-r);
      max[i] = space_to_bbox(c+r);
    }
    return bounding_box::min_and_max(min, max);
  }
  
  time_type escape_time(accessor const* accessor, entity_ref e, bounding_box const& bbox) {
    auto c = accessor->get<circle_shape>(e);
    time_type result = never;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      const poly min_poly = c->center(i) - c->radius - bbox_to_space(bbox.min(i));
      for (auto j = min_poly.sign_interval_boundaries_upper_bound(accessor->now()); j != min_poly.sign_interval_boundaries_end(); ++j) {
        if ((min_poly(*j) < 0) && ((result == never) || (*j < result))) {
          result = *j;
        }
      }
      const poly max_poly = c->center(i) + c->radius - bbox_to_space(bbox.max(i));
      for (auto j = max_poly.sign_interval_boundaries_upper_bound(accessor->now()); j != max_poly.sign_interval_boundaries_end(); ++j) {
        if ((max_poly(*j) > 0) && ((result == never) || (*j < result))) {
          result = *j;
        }
      }
    }
    return result;
  }
  
  void anticipate_interactions(accessor* accessor, entity_ref e0, entity_ref e1) {
    auto c0 = accessor->get<circle_shape>(e0);
    auto c1 = accessor->get<circle_shape>(e1);
    
    const poly d = distish(*c0, *c1);
    time_type when = never;
    const bool overlapping = d(0) < 0;
    for (auto i = d.sign_interval_boundaries_upper_bound(accessor->now()); i != d.sign_interval_boundaries_end(); ++i) {
      if ((d(*i) < 0) != overlapping) {
        when = *i;
        break;
      }
    }
    if (overlapping) {
      // TODO: this isn't idempotent when called at different times
      when = std::max(when, accessor->now() + (time_units_per_gloppp>>6));
    }
    if (when != never) {
      accessor->anticipate_event(when, std::shared_ptr<event>(new circles_interact(e0.id(), e1.id())));
    }
  }
};

typedef typename bbcd_system::operations<time_steward, bbcd_funcs> bbcd_operations;


struct circles_overlapping_bbox_filter {
  circles_overlapping_bbox_filter(time_steward::accessor const* accessor, bbcd_system::bounding_box bbox):accessor(accessor),bbox(bbox){}
  time_steward::accessor const* accessor;
  bbcd_system::bounding_box bbox;
  
  bool min_cost(bbcd_system::bounding_box bbox2) { return bbox.overlaps(bbox2); }
  bool cost(entity_id id) {
    auto e = accessor->get<circle_shape>(accessor->get(id));
    auto c = e->center.get_term<space_coordinate>(accessor->now(), 0);
    auto r = e->radius(accessor->now());
    space_coordinate square_sum = 0;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      space_coordinate min = bbox_to_space(bbox.min()[i]);
      space_coordinate max = bbox_to_space(bbox.max()[i]);
      if (c[i] < min) { square_sum += (c[i]-min)*(c[i]-min); }
      if (c[i] > max) { square_sum += (c[i]-max)*(c[i]-max); }
    }
    return r*r >= square_sum;
  }
};

class initialize_world : public time_steward::event {
public:
  void operator()(time_steward::accessor* accessor)const override {
    auto bbcd = bbcd_operations::create_bbox_collision_detector(accessor);
    auto g = accessor->set<global_data>(accessor->get(time_steward_system::global_object_id), global_data(bbcd.id()));
    srand(0);
    for (int i = 0; i < 150; ++i) {
      auto e = accessor->create_entity();
      std::vector<space_coordinate> xv; xv.push_back(-(arena_width/2)+(rand()%1024)*(arena_width/1024));
      xv.push_back((-10000+(rand()%20000))*arena_width/(100000*time_units_per_gloppp)); // "Might cross the whole arena in only ten gloppps"
      std::vector<space_coordinate> yv; yv.push_back(-(arena_width/2)+(rand()%1024)*(arena_width/1024));
      yv.push_back((-10000+(rand()%20000))*arena_width/(100000*time_units_per_gloppp)); // "Might cross the whole arena in only ten gloppps"
      accessor->set<circle_shape>(e, circle_shape(
        poly_fd_vector(poly(0, bounded_int_calculus::polynomial<time_type,space_coordinate>(xv)), poly(0, bounded_int_calculus::polynomial<time_type,space_coordinate>(yv))),
        poly(0, bounded_int_calculus::polynomial<time_type,space_coordinate>((50+(rand()%100))*arena_width/5000)) // 1-3% the width of the arena
      ));
      accessor->set<circle_overlaps>(e, circle_overlaps());
      bbcd_operations::insert(accessor, bbcd.id(), e);
    }
  }
};
