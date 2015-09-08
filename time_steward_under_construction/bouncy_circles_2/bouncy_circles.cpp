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


#include <boost/optional.hpp>

#include <assert.h>

#include <iostream>

#include "../bounded_int_calculus.hpp"

struct circle_shape {
poly_vector center;
space_coordinate radius;
};

MAKE_FIELD (circle_shape);
MAKE_FIELD (circle_tile, FD_vector);
MAKE_FIELD (circles_here, persistent_ID_map);

MAKE_EVENT (circle_switches_tile, 1) {
FD_vector old_tile = GET (entity, circle_tile);
FD_vector new_tile = space_to_tile (GET_JUST (entity, circle_shape).center.get_term <space_coordinate> (NOW, 0));
assert (new_tile != old_tile);
SET (entity, circle_tile, new_tile);
auto old_tile_entity = tile_entity (accessor, old_tile);
auto new_tile_entity = tile_entity (accessor, new_tile);
SET (old_tile_entity, circles_here, GET (old_tile_entity, circles_here).erase (entity .id ()));
SET (new_tile_entity, circles_here, GET (new_tile_entity, circles_here).insert (entity .id ()));
};

MAKE_TRIGGER (circle_may_switch_tile,1) {
circle_shape const & circle = GET_JUST (entity, circle_shape);
FD_vector current_tile = space_to_tile (circle.center.get_term <space_coordinate> (NOW, 0));
assert (GET (entity, circle_tile)== current_tile);
time_type when = never;
for (num_coordinates_type dimension = 0; dimension <num_dimensions;++ dimension) {
when = earliest(when, next_time_positive (NOW, circle.center (dimension) - tile_to_space_max (current_tile (dimension)));
when = earliest(when, next_time_negative (NOW, circle.center (dimension) - tile_to_space_min (current_tile (dimension)));
}
accessor-> anticipate_event (when, new circle_switches_tile (entity .id ()));
};

MAKE_EVENT (circles_interact, 2) {

};

MAKE_TRIGGER (circle_may_interact,1) {
FD_vector tile = GET (entity, circle_tile);
for (tile_coordinate coordinate_0 = tile (0)-1; coordinate_0 <= tile (0) +1; ++ coordinate_0) {
for (tile_coordinate coordinate_1 = tile (1)-1; coordinate_1 <= tile (1) +1; ++ coordinate_1) {
for (entity_id id1 : GET (tile_entity (accessor, FD_vector (coordinate_0, coordinate_1)), circles_here) {
CAPTURE (id1, entity_1);
auto shape_0 = GET (entity, circle_shape);
auto shape_1 = GET (entity_1, circle_shape);
auto center_difference = shape_1.center - shape_2.center;
auto center_distance_squared = center_difference (0).dot (center_difference (0)) + center_difference (1).dot(center_difference (1));
space_coordinate radius_sum = shape_1.radius - shape_0.radius;
space_coordinate current_distance_rounded_down = isqrt (center_distance_squared (NOW))-radius_sum;
space_coordinate current_overlap_category = divide (current_distance_rounded_down, overlap_increment_size, rounding_strategy <round_down, negative_continuous_with_positive> ());
if (current_overlap_category >0) {current_overlap_category = 0;}
time_type when = next_time_negative (NOW, center_distance_squared - square (radius_sum + current_overlap_category*overlap_increment_size));
if (current_overlap_category < 0) {
when = earliest (when, next_time_nonnegative (NOW, center_distance_squared - square (radius_sum + (current_overlap_category +1)*overlap_increment_size)));
}
accessor-> anticipate_event (when, new circles_interact (entity .id (), entity_1 .id ());
}
}
} 
}

