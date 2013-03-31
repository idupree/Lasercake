/*

    Copyright Eli Dupree and Isaac Dupree, 2011, 2012, 2013

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

#define USE_BOUNDS_CHECKED_INTS 1

#include <GL/glew.h>
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "SDL.h"
#include "GL/gl.h"
#include "GL/glu.h"

#include <iostream>
#include <cmath>

#include <array>
#include <vector>
#include <queue>

#include <boost/shared_ptr.hpp>

#include "../units.hpp"
#include "../gl_data_format.hpp"
#include "../utils.hpp"
#include "../data_structures/misc_structures.hpp"


#ifdef LASERCAKE_HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>
#include <boost/chrono/thread_clock.hpp>

#include <gmpxx.h>

#include "../utils.hpp"
//#include "../data_structures/geometry.hpp"

#if 0
namespace std {
template<>
struct numeric_limits<mpz_class>{
  static const bool is_specialized = true;
  static const int digits = INT_MAX;
  static const bool is_integer = true;
  static const bool is_signed = true;
  static const bool radix = 2;
};
template<typename T1, typename T2>
struct numeric_limits<__gmp_expr<T1, T2>>{
  static const bool is_specialized = true;
  static const int digits = INT_MAX;
  static const bool is_integer = true;
  static const bool is_signed = true;
  static const bool radix = 2;
};
}
#endif
/*
void fffff() {
  //const acceleration1d gravity_acceleration_magnitude = mpz(9806650) * (micro*meters) / (seconds*seconds) * identity(distance_units / (micro*meters));
const auto aaaaa = mpz(9806650) * (micro*meters);
const auto bbbbb = aaaaa / (seconds * seconds);
const auto iiiii = identity(distance_units / (micro*meters));
const auto ccccc = bbbbb * iiiii;
}*/

// prevents mpz_class's clever optimizations,
// so that physical_quantity can handle it better.
// TODO do better.
class mpz {
  mpz_class data_;
public:
  mpz(){}
  mpz(mpz_class d):data_(d){}
  template<typename T> mpz(T t) : data_(t) {}

  explicit operator bool()const { return data_ != 0; }
  explicit operator float()const { return data_.get_d(); }
  explicit operator double()const { return data_.get_d(); }

friend inline mpz operator*(mpz a, mpz b) { return mpz(a.data_ * b.data_); }
friend inline mpz operator/(mpz a, mpz b) { return mpz(a.data_ / b.data_); }
friend inline mpz operator%(mpz a, mpz b) { return mpz(a.data_ % b.data_); }
friend inline mpz operator+(mpz a) { return mpz(+a.data_); }
friend inline mpz operator-(mpz a) { return mpz(-a.data_); }
friend inline mpz abs(mpz a) { return abs(a.data_); }
friend inline mpz operator+(mpz a, mpz b) { return mpz(a.data_ + b.data_); }
friend inline mpz operator-(mpz a, mpz b) { return mpz(a.data_ - b.data_); }

friend inline mpz operator&(mpz a, mpz b) { return mpz(a.data_ & b.data_); }
friend inline mpz operator|(mpz a, mpz b) { return mpz(a.data_ | b.data_); }
friend inline mpz operator^(mpz a, mpz b) { return mpz(a.data_ ^ b.data_); }
friend inline mpz operator~(mpz a) { return mpz(~a); }

friend inline bool operator==(mpz a, mpz b) { return a.data_ == b.data_; }
friend inline bool operator!=(mpz a, mpz b) { return a.data_ != b.data_; }
friend inline bool operator<(mpz a, mpz b) { return a.data_ < b.data_; }
friend inline bool operator>(mpz a, mpz b) { return a.data_ > b.data_; }
friend inline bool operator>=(mpz a, mpz b) { return a.data_ >= b.data_; }
friend inline bool operator<=(mpz a, mpz b) { return a.data_ <= b.data_; }

friend inline mpz operator<<(mpz a, uint32_t shift)
{ return mpz(a.data_ << shift); }
friend inline mpz operator>>(mpz a, uint32_t shift)
{ return mpz(a.data_ >> shift); }

friend inline mpz& operator+=(mpz& a, mpz b) { a = a + b; return a; }
friend inline mpz& operator-=(mpz& a, mpz b) { a = a - b; return a; }
friend inline mpz& operator*=(mpz& a, mpz b) { a = a * b; return a; }
friend inline mpz& operator&=(mpz& a, mpz b) { a = a & b; return a; }
friend inline mpz& operator|=(mpz& a, mpz b) { a = a | b; return a; }
friend inline mpz& operator^=(mpz& a, mpz b) { a = a ^ b; return a; }
friend inline mpz& operator<<=(mpz& a, uint32_t shift) { a = a << shift; return a; }
friend inline mpz& operator>>=(mpz& a, uint32_t shift) { a = a >> shift; return a; }

friend inline mpz& operator++(mpz& a) { ++a.data_; return a; }
friend inline mpz operator++(mpz& a, int) { mpz old = a; ++a.data_; return old; }
friend inline mpz& operator--(mpz& a) { --a.data_; return a; }
friend inline mpz operator--(mpz& a, int) { mpz old = a; ++a.data_; return old; }

//Imagine taking the sqrt of a signed 8 bit value.  Can it fit into 4 bits?
// sqrt(127) is 11.something, which is less than 15 (unsigned 4bit max) but
// greater than 7 (signed 4bit max).  So we conservatively leave enough space here.
friend inline mpz isqrt(mpz a)
{ return mpz(sqrt(a.data_)); }
//friend inline int32_t ilog2(mpz a) { return ; }

friend std::ostream& operator<<(std::ostream& os, mpz a) {
  return os << a.data_;
}
friend float get_primitive_float(mpz a) {
  return a.data_.get_d();
}
};
namespace std {
template<> struct numeric_limits<mpz> {
  static const bool is_specialized = true;
  static const bool is_integer = true;
  static const bool is_signed = true;
  static const int digits = INT_MAX;
  static const int radix = 2;
};
}
namespace boost {
template<> struct make_signed<mpz> { typedef mpz type; };
// HACK: used because rounding_strategy negative_is_forbidden efficiency hack.
template<> struct make_unsigned<mpz> { typedef mpz type; };
}


static_assert(get_units<mpz>::is_nonunit_scalar, "erejrqrq");

namespace YO {
using boost::shared_ptr;
using boost::dynamic_pointer_cast;


struct gl_polygon {
  std::vector<gl_data_format::vertex_with_color> vertices_;
};
// we want to sort threevertices together: triangles: not individual vertices
struct gl_triangle {
  array<gl_data_format::vertex_with_color, 3> vertices_;
};
typedef std::vector<gl_triangle> gl_triangles;
glm::vec3 v_to_gv(gl_data_format::vertex v) {
  return glm::vec3(v.x, v.y, v.z);
}
float gl_triangle_distance_order(glm::vec3 from, gl_triangle const& triangle) {
  return glm::distance(from, v_to_gv(triangle.vertices_[0].v))
       + glm::distance(from, v_to_gv(triangle.vertices_[1].v))
       + glm::distance(from, v_to_gv(triangle.vertices_[2].v));
}
void sort_gl_triangles_near_to_far(glm::vec3 view_center, gl_triangles& ts) {
  std::sort(ts.begin(), ts.end(), [view_center](gl_triangle const& t1, gl_triangle const& t2) {
    return gl_triangle_distance_order(view_center, t1) < gl_triangle_distance_order(view_center, t2);
  });
}
void sort_gl_triangles_far_to_near(glm::vec3 view_center, gl_triangles& ts) {
  std::sort(ts.begin(), ts.end(), [view_center](gl_triangle const& t1, gl_triangle const& t2) {
    return gl_triangle_distance_order(view_center, t2) < gl_triangle_distance_order(view_center, t1);
  });
}



///duplicate code
glm::vec3 average_direction(glm::vec3 d1, glm::vec3 d2) {
  return glm::normalize(
      glm::normalize(d1) + glm::normalize(d2)
    );
}
gl_data_format::color ceethrough(gl_data_format::color c) {
  c.a = 0;
  return c;
}
gl_data_format::vertex_with_color ceethrough(gl_data_format::vertex_with_color vc) {
  vc.c.a = 0;
  return vc;
}
// makes it gradually translucent towards the centre.
void push_wireframe_polygon(
      gl_triangles& coll, GLfloat width,
      gl_polygon polygon) {
  const auto vs = polygon.vertices_;
  const int num_vertices = vs.size();
  const int last = num_vertices - 1;
  caller_correct_if(num_vertices >= 3, "that's not a polygon");

  auto vcenters(vs);
  // todo don't round to float before subtraction wait its already float
  vcenters[0].v = v_to_gv(vcenters[0].v) + width * average_direction(v_to_gv(vs[last].v) - v_to_gv(vs[0].v), v_to_gv(vs[1].v) - v_to_gv(vs[0].v));
  for(int i = 1; i < last; ++i) {
    vcenters[i].v = v_to_gv(vcenters[i].v) + width * average_direction(v_to_gv(vs[i-1].v) - v_to_gv(vs[i].v), v_to_gv(vs[i+1].v) - v_to_gv(vs[i].v));
  }
  vcenters[last].v = v_to_gv(vcenters[last].v) + width * average_direction(v_to_gv(vs[last-1].v) - v_to_gv(vs[last].v), v_to_gv(vs[0].v) - v_to_gv(vs[last].v));
  for(int i = 0; i != num_vertices; ++i) {
    vcenters[i].c.a = 0;
  }
  for(int i = 0; i < last; ++i) {
    coll.push_back(gl_triangle{{{ vs[i], vcenters[i], vs[i+1] }}});
    coll.push_back(gl_triangle{{{ vcenters[i], vs[i+1], vcenters[i+1] }}});
  }
  coll.push_back(gl_triangle{{{ vs[last], vcenters[last], vs[0] }}});
  coll.push_back(gl_triangle{{{ vcenters[last], vs[0], vcenters[0] }}});
}



//typedef seconds<non_normalized_rational<64>> time_type;
//TODO make this size better for long durations or something
//long distances
//128
typedef typename units_prod<pico_t, seconds_t>::type time_units_t;
constexpr time_units_t time_units = time_units_t();
typedef physical_quantity<mpz, time_units_t> time_type;

typedef typename units_prod<nano_t, meters_t>::type distance_units_t;
constexpr distance_units_t distance_units = distance_units_t();
typedef physical_quantity<mpz, distance_units_t> distance;
typedef physical_quantity<mpz, typename units_prod<distance_units_t, dim::second<(-1)> >::type> velocity1d;
typedef physical_quantity<mpz, typename units_prod<distance_units_t, dim::second<(-2)> >::type> acceleration1d;

// Standard (Earth-equivalent) gravity: precisely 9.80665 m/s2
const acceleration1d gravity_acceleration_magnitude = mpz(9806650) * (micro*meters) / (seconds*seconds) * identity(distance_units / (micro*meters));
const vector3<acceleration1d> gravity_acceleration(0, 0, -gravity_acceleration_magnitude);

typedef uint32_t region_idx_type;
typedef uint32_t face_idx_type;
typedef uint32_t vertex_idx_type;
typedef uint64_t revision_number_type;
//constexpr region_idx_type no_region_idx = std::numeric_limits<region_idx_type>::max();
//constexpr face_idx_type no_face_idx = std::numeric_limits<region_idx_type>::max();
//constexpr vertex_idx_type no_vertex_idx = std::numeric_limits<vertex_idx_type>::max();

const int triangle_sides = 3;

// A vertex is an imaginary thing that moves around
// independent of the substances it is helping to represent.
struct vertex {
  time_type base_time_;
  vector3<distance> vertex_position_;
  vector3<velocity1d> vertex_velocity_;
  vector3<acceleration1d> vertex_acceleration_;
};
struct edge {
};
// All faces are triangles.
struct face {
  face():base_time_(0),revision_number_(0),ABC(0),D(0),D_velocity(0),D_acceleration(0){}
  time_type base_time_;
  revision_number_type revision_number_;
  vector3<mpz> ABC;
  distance D;
  velocity1d D_velocity;
  acceleration1d D_acceleration;
  std::vector<face_idx_type> neighboring_faces_;
  std::vector<region_idx_type> neighboring_regions_;
  face updated_to_time(time_type t)const {
    assert(t >= base_time_);
    face result(*this);
    const time_type relative_time = t - base_time_;
    const auto relative_timem = relative_time/identity(units_factor<1, 1000000000>());
    //TODO deal with times somehow more correctly
    result.D          += D_velocity    *relative_timem                / identity(units_factor<1,    1000>())
                       + D_acceleration*relative_timem*relative_timem / identity(units_factor<1, 1000000>())/2;
    result.D_velocity += D_acceleration*relative_timem                / identity(units_factor<1,    1000>());
    result.base_time_ = t;
    ++result.revision_number_;
    return result;
  }
};

mpz scalar_triple_product(vector3<mpz> v1, vector3<mpz> v2, vector3<mpz> v3) {
  return v1(X)*v2(Y)*v3(Z) - v1(X)*v3(Y)*v2(Z) + v2(X)*v3(Y)*v1(Z) - v2(X)*v1(Y)*v3(Z) + v3(X)*v1(Y)*v2(Z) - v3(X)*v2(Y)*v1(Z);
}
vector3<distance> approx_loc_of_triple_intersection_of_up_to_date_faces(face const& f1, face const& f2, face const& f3) {
  mpz p = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  assert(p != 0);
  return vector3<distance>(
    f1.D * (f2.ABC(1)*f3.ABC(2) - f3.ABC(1)*f2.ABC(2)) + f2.D * (f3.ABC(1)*f1.ABC(2) - f1.ABC(1)*f3.ABC(2)) + f3.D * (f1.ABC(1)*f2.ABC(2) - f2.ABC(1)*f1.ABC(2)),
    f1.D * (f2.ABC(2)*f3.ABC(0) - f3.ABC(2)*f2.ABC(0)) + f2.D * (f3.ABC(2)*f1.ABC(0) - f1.ABC(2)*f3.ABC(0)) + f3.D * (f1.ABC(2)*f2.ABC(0) - f2.ABC(2)*f1.ABC(0)),
    f1.D * (f2.ABC(0)*f3.ABC(1) - f3.ABC(0)*f2.ABC(1)) + f2.D * (f3.ABC(0)*f1.ABC(1) - f1.ABC(0)*f3.ABC(1)) + f3.D * (f1.ABC(0)*f2.ABC(1) - f2.ABC(0)*f1.ABC(1))
                 ) / p;
}

// TODO : What about the cases where two of the planes are parallel? That's legit (handle them separately? they're easier)
faux_optional<time_type> how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(face const& f1, face const& f2, face const& f3, face const& f4) {
  // When the 4x4 determinant is 0.
  // That's
  // + D1 * scalar_triple_product(f2.ABC, f3.ABC, f4.ABC)
  // - D2 * scalar_triple_product(f3.ABC, f4.ABC, f1.ABC)
  // + D3 * scalar_triple_product(f4.ABC, f1.ABC, f2.ABC)
  // - D4 * scalar_triple_product(f1.ABC, f2.ABC, f3.ABC)
  // But Dn is Dn + dDn*t + .5ddDn*t^2
  // so
  const mpz t1 = scalar_triple_product(f2.ABC, f3.ABC, f4.ABC);
  const mpz t2 = scalar_triple_product(f3.ABC, f4.ABC, f1.ABC);
  const mpz t3 = scalar_triple_product(f4.ABC, f1.ABC, f2.ABC);
  const mpz t4 = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  // in at^2 + bt + c = 0
  acceleration1d a_times_2 = f1.D_acceleration*t1 - f2.D_acceleration*t2 + f3.D_acceleration*t3 - f4.D_acceleration*t4;
  velocity1d     b         = f1.D_velocity    *t1 - f2.D_velocity    *t2 + f3.D_velocity    *t3 - f4.D_velocity    *t4;
  distance       c         = f1.D             *t1 - f2.D             *t2 + f3.D             *t3 - f4.D             *t4;

#if 0
if(!(f4.ABC(X) != 0 || f4.ABC(Y) != 0 || f4.ABC(Z) <= 0 || f3.ABC(Z) == 0)){

  std::cerr << f1.ABC << "\n";
  std::cerr << f2.ABC << "\n";
  std::cerr << f3.ABC << "\n";
  std::cerr << f4.ABC << "\n\n";
  std::cerr << t1 << "\n";
  std::cerr << t2 << "\n";
  std::cerr << t3 << "\n";
  std::cerr << t4 << "\n\n";
  std::cerr << a_times_2 << "\n";
  std::cerr << b << "\n";
  std::cerr << c << "\n\n";
}
#endif
  
  // (-b +/- sqrt(b^2 - 2(a_times_2)c)) / a_times_2
  const auto discriminant = b*b - a_times_2*2*c;
  if (discriminant < 0) {
    return boost::none;
  }
  else {
    const rounding_strategy<round_down, negative_is_forbidden> strat;
    if (a_times_2 == 0) {
      if (b == 0) {
        // Eww. A constant. They're either ALWAYS intersecting or NEVER.
        if (c == 0) {
          // Uhh... apparently sometimes c=0 can occur when they're never intersecting, too.
          // I don't really know. Returning none here is techincally a false-negative sometimes.
          // But probably only in border cases.
          return boost::none;
        }
        else return boost::none;
      }
      if (b < 0) {
        b = -b;
        c = -c;
      }
      if (c > 0) return boost::none;
      // hack - should just be
      // return divide(-c * identity(time_units / seconds), b, strat);
      // but overflow stuff
      // also see below
      const time_type zero = divide(-c * mpz(identity(time_units / seconds)*seconds/time_units), b, strat)*time_units/seconds;
      //std::cerr << a_times_2 << ", " << b << ", " << c << ": " << zero << ", " << get(zero,time_units)*get(zero,time_units)*get(a_times_2,typename units_prod<distance_units_t, dim::second<(-2)> >::type()) / 2 + get(zero,time_units)*get(b,typename units_prod<distance_units_t, dim::second<(-1)> >::type())*mpz(1e12) + get(c,distance_units)*mpz(1e24) << "\n";
      return zero;
    }
    else {
      if (a_times_2 < 0) {
        a_times_2 = -a_times_2;
        b = -b;
        c = -c;
      }
      // we want the earlier time, but not if it's negative
      const velocity1d sqrt_disc = isqrt(discriminant);
      const velocity1d  lesser_numerator = -b - sqrt_disc;
      if ( lesser_numerator >= 0)  {
        const time_type zero = divide(
            lesser_numerator * mpz(identity(time_units / seconds)*seconds/time_units),
            a_times_2,
            strat)*time_units/seconds;
        //std::cerr << a_times_2 << ", " << b << ", " << c << ": " << zero << ", " << get(zero,time_units)*get(zero,time_units)*get(a_times_2,typename units_prod<distance_units_t, dim::second<(-2)> >::type()) / 2 + get(zero,time_units)*get(b,typename units_prod<distance_units_t, dim::second<(-1)> >::type())*mpz(1e12) + get(c,distance_units)*mpz(1e24) << "\n";
        return zero;
      }
      const velocity1d greater_numerator = -b + sqrt_disc;
      if (greater_numerator >= 0) {
        const time_type zero = divide(
            greater_numerator * mpz(identity(time_units / seconds)*seconds/time_units),
            a_times_2,
            strat)*time_units/seconds;
        //std::cerr << a_times_2 << ", " << b << ", " << c << ": " << zero << ", " << get(zero,time_units)*get(zero,time_units)*get(a_times_2,typename units_prod<distance_units_t, dim::second<(-2)> >::type()) / 2 + get(zero,time_units)*get(b,typename units_prod<distance_units_t, dim::second<(-1)> >::type())*mpz(1e12) + get(c,distance_units)*mpz(1e24) << "\n";
        return zero;
      }
      return boost::none;
    }
  }
}
faux_optional<time_type> when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(face const& f1, face const& f2, face const& f3, face const& f4) {
  if (faux_optional<time_type> result = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1,f2,f3,f4)) {
    return *result + f1.base_time_;
  }
  else return boost::none;
}

struct region_vertex {
  vertex_idx_type shared_vertex_data_;
  //other data
};
struct region {
  // Regions are arbitrary, non-convex polyhedra which can have holes.
  //std::vector<vertex_idx_type, tetrahedron_sides> vertices_;
  std::vector<face_idx_type> faces_;
};

struct event {
  virtual ~event(){}
  time_type when_event_occurs_;
  event(time_type t):when_event_occurs_(t){}
};
struct collision : public event {
  face_idx_type f1_;
  face_idx_type f2_;
  face_idx_type f3_;
  face_idx_type f4_;
  revision_number_type r1_;
  revision_number_type r2_;
  revision_number_type r3_;
  revision_number_type r4_;
  
  collision(time_type t,
            face_idx_type f1, revision_number_type r1,
            face_idx_type f2, revision_number_type r2,
            face_idx_type f3, revision_number_type r3,
            face_idx_type f4, revision_number_type r4):
  event(t),f1_(f1),f2_(f2),f3_(f3),f4_(f4),r1_(r1),r2_(r2),r3_(r3),r4_(r4){}
};
struct vertex_face_collision : public collision {
  face_idx_type vertex_face_1()const { return f1_; }
  face_idx_type vertex_face_2()const { return f2_; }
  face_idx_type vertex_face_3()const { return f3_; }
  face_idx_type struck_face  ()const { return f4_; }
  vertex_face_collision(time_type t,
                        face_idx_type vf1, revision_number_type vf1r,
                        face_idx_type vf2, revision_number_type vf2r,
                        face_idx_type vf3, revision_number_type vf3r,
                        face_idx_type sf , revision_number_type sfr):
    collision(t, vf1, vf1r, vf2, vf2r, vf3, vf3r, sf, sfr){}
};
struct edge_edge_collision : public collision {
  face_idx_type edge_1_face_1()const { return f1_; }
  face_idx_type edge_1_face_2()const { return f2_; }
  face_idx_type edge_2_face_1()const { return f3_; }
  face_idx_type edge_2_face_2()const { return f4_; }
  edge_edge_collision(time_type t,
                      face_idx_type e11, revision_number_type e11r,
                      face_idx_type e12, revision_number_type e12r,
                      face_idx_type e21, revision_number_type e21r,
                      face_idx_type e22, revision_number_type e22r):
    collision(t, e11, e11r, e12, e12r, e21, e21r, e22, e22r){}
};

struct event_ptr_compare {
  bool operator()(shared_ptr<event> const& a, shared_ptr<event> const& b)const { return a->when_event_occurs_ > b->when_event_occurs_; }
};

// Maybe the structures should be required to follow some kind of right-hand-rule
// order? Then asserts could check that it is still correct, and OpenGL would
// be more easily able to draw directional surfaces.
class grand_structure_of_lasercake {
  std::vector<region> regions_;
  std::vector<face> faces_;
  std::priority_queue<shared_ptr<event>, std::vector<shared_ptr<event>>, event_ptr_compare> next_events_;

  bool bounded_edges_cross__hack(time_type t, face const& e11_old, face const& e12_old, size_t neighbor_idx_in_e11, face const& e21_old, face const& e22_old, size_t neighbor_idx_in_e21)const {
    const face e11 = e11_old.updated_to_time(t);
    const face e12 = e12_old.updated_to_time(t);
    const face e21 = e21_old.updated_to_time(t);
    const face e22 = e22_old.updated_to_time(t);
    const vector3<distance> approx_cross_location = approx_loc_of_triple_intersection_of_up_to_date_faces(e11, e12, e21); // should be about the same location for any three
    for (int i = 0; i < 2; ++i) {
      face const& e1 = i ? e21 : e11;
      face const& e2 = i ? e22 : e12;
      size_t const& neighbor_idx = i ? neighbor_idx_in_e21 : neighbor_idx_in_e11;
      const size_t next_neighbor_idx = (neighbor_idx + 1) % e1.neighboring_faces_.size();
      const size_t prev_neighbor_idx = (neighbor_idx+e1.neighboring_faces_.size() - 1) % e1.neighboring_faces_.size();

      const vector3<distance> end1 = approx_loc_of_triple_intersection_of_up_to_date_faces(
        e1, e2, faces_[e1.neighboring_faces_[prev_neighbor_idx]].updated_to_time(t));
      const vector3<distance> end2 = approx_loc_of_triple_intersection_of_up_to_date_faces(
        e1, e2, faces_[e1.neighboring_faces_[next_neighbor_idx]].updated_to_time(t));

      // This edge is okay if the intersection point is between the two ends of the edge.
      // We would normally have to check only one dimension. Do it for all because sometimes
      // the edge is parallel to an axis.
      if ((((end1.x-approx_cross_location.x) * (end2.x-approx_cross_location.x)) < 0)
        || (((end1.y-approx_cross_location.y) * (end2.y-approx_cross_location.y)) < 0)
        || (((end1.z-approx_cross_location.z) * (end2.z-approx_cross_location.z)) < 0)) {
        // ...
      }
      else {
        return false;
      }
    }
    return true;
  }

  bool vertex_is_in_bounded_face__hack(time_type t, face const& v1_old, face const& v2_old, face const& v3_old, face const& f_old)const {
    const face v1 = v1_old.updated_to_time(t);
    const face v2 = v2_old.updated_to_time(t);
    const face v3 = v3_old.updated_to_time(t);
    const face  f =  f_old.updated_to_time(t);
    
    // we prefer to eliminate the face that's most foreshortened, i.e. the one for which the normal's value is biggest
    which_dimension_type dim1 = ((abs(f.ABC(X)) > abs(f.ABC(Y))) && (abs(f.ABC(X)) > abs(f.ABC(Z)))) ? Y : X;
    which_dimension_type dim2 = ((abs(f.ABC(Z)) > abs(f.ABC(Y))) && (abs(f.ABC(Z)) > abs(f.ABC(X)))) ? Y : Z;

    const vector3<distance> v = approx_loc_of_triple_intersection_of_up_to_date_faces(v1, v2, v3);

    int crosses = 0;
    
    for (size_t i = 0; i < f.neighboring_faces_.size(); ++i) {
      const face_idx_type neighbor_id_1 = f.neighboring_faces_[i];
      const face_idx_type neighbor_id_2 = f.neighboring_faces_[(i+1)%f.neighboring_faces_.size()];
      const face_idx_type neighbor_id_3 = f.neighboring_faces_[(i+2)%f.neighboring_faces_.size()];
      const face present_neighbor_1 = faces_[neighbor_id_1].updated_to_time(t);
      const face present_neighbor_2 = faces_[neighbor_id_2].updated_to_time(t);
      const face present_neighbor_3 = faces_[neighbor_id_3].updated_to_time(t);
      const vector3<distance> pv1 = approx_loc_of_triple_intersection_of_up_to_date_faces(f, present_neighbor_1, present_neighbor_2);
      const vector3<distance> pv2 = approx_loc_of_triple_intersection_of_up_to_date_faces(f, present_neighbor_2, present_neighbor_3);
      if ((pv1(dim1) > v(dim1)) != (pv2(dim1) > v(dim1))) {
        const vector3<distance> pv1_to_v = (v - pv1);
        const vector3<distance> pv1_to_pv2 = (pv2 - pv1);
        if (pv1_to_pv2(dim2) == 0) {
          if (pv1_to_v(dim2) <= 0) {
            ++crosses;
          }
        }
        else if (pv1_to_v(dim2) == 0) {
          if (pv1_to_pv2(dim2) >= 0) {
            ++crosses;
          }
        }
        else {
          const mpz sign_correction = (pv1_to_v(dim1)*pv1_to_pv2(dim1) > 0) ? 1 : -1;
          if (pv1_to_v(dim2)*pv1_to_pv2(dim1)*sign_correction >= pv1_to_pv2(dim2)*pv1_to_v(dim1)*sign_correction) {
            ++crosses;
          }
        }
      }
    }
    return (crosses % 2);
  }
  
  void insert_events_involving(face_idx_type fi) {
    // Two kinds of events: Vertex-face collisions and edge-edge collisions.
    // From one face, that's three kinds (f being part of the vertex and f being the face.)
    // Also there's an arguable three extra categories depending on whether the face
    // neighbors zero, one, two, or three of the faces establishing the vertex.
    // TODO: try to de-duplicate this "iterate through the present vertices" system which is used elsewhere in the code.
    face const& f = faces_[fi];
    
    for (size_t i = 0; i < f.neighboring_faces_.size(); ++i) {
      const face_idx_type neighbor_id_1 = f.neighboring_faces_[i];
      const face_idx_type neighbor_id_2 = f.neighboring_faces_[(i+1)%f.neighboring_faces_.size()];
      face const& old_neighbor_1 = faces_[neighbor_id_1];
      face const& old_neighbor_2 = faces_[neighbor_id_2];
      const face present_neighbor_1 = old_neighbor_1.updated_to_time(f.base_time_);
      const face present_neighbor_2 = old_neighbor_2.updated_to_time(f.base_time_);

      for (region_idx_type ri : f.neighboring_regions_) {
        region const& r = regions_[ri];
        for (face_idx_type fi2 : r.faces_) {
          // A vertex doesn't collide with one of its own faces
          if ((fi2 != fi) && (fi2 != neighbor_id_1)) {
            face const& old_other_face = faces_[fi2];
            const face present_other_face = old_other_face.updated_to_time(f.base_time_);
            if (fi2 != neighbor_id_2) {
              if (faux_optional<time_type> collision_time = when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(
                f, present_neighbor_1, present_neighbor_2, present_other_face)) {
                if (vertex_is_in_bounded_face__hack(*collision_time, f, present_neighbor_1, present_neighbor_2, present_other_face)) {
                  next_events_.push(shared_ptr<event>(new vertex_face_collision(
                      *collision_time,
                                 fi,              f.revision_number_,
                      neighbor_id_1, old_neighbor_1.revision_number_,
                      neighbor_id_2, old_neighbor_2.revision_number_,
                                fi2, old_other_face.revision_number_)));
                }
              }
            }
            
            // Also catch edge-edge collisions...
            for (size_t j = 0; j < present_other_face.neighboring_faces_.size(); ++j) {
              const face_idx_type other_neighbor_id = present_other_face.neighboring_faces_[j];
              if ((other_neighbor_id != fi) && (other_neighbor_id != neighbor_id_1)) {
                // only consider each edge once. (without this if, each one would be considered twice...)
                if (other_neighbor_id > fi2) {
                  face const& old_other_neighbor = faces_[other_neighbor_id];
                  const face present_other_neighbor = old_other_neighbor.updated_to_time(f.base_time_);
                  if (faux_optional<time_type> collision_time = when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(
                    f, present_neighbor_1, present_other_face, present_other_neighbor)) {
                    if (bounded_edges_cross__hack(*collision_time, f, present_neighbor_1, i, present_other_face, present_other_neighbor, j)) {
                      next_events_.push(shared_ptr<event>(new edge_edge_collision(
                          *collision_time,
                                         fi,                  f.revision_number_,
                              neighbor_id_1,     old_neighbor_1.revision_number_,
                                        fi2,     old_other_face.revision_number_,
                          other_neighbor_id, old_other_neighbor.revision_number_)));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    for (region_idx_type ri : f.neighboring_regions_) {
      region const& r = regions_[ri];
      for (face_idx_type fi2 : r.faces_) {
        if (fi2 != fi) {
          face const& old_other_face = faces_[fi2];
          const face present_other_face = old_other_face.updated_to_time(f.base_time_);
          for (size_t i = 0; i < present_other_face.neighboring_faces_.size(); ++i) {
            const face_idx_type neighbor_id_1 = present_other_face.neighboring_faces_[i];
            const face_idx_type neighbor_id_2 = present_other_face.neighboring_faces_[(i+1)%present_other_face.neighboring_faces_.size()];
            // a vertex doesn't collide with one of its own faces
            if ((neighbor_id_1 != fi) && (neighbor_id_2 != fi)) {
              // only consider each vertex once. (without this if, each one would be considered three times...)
              if ((neighbor_id_1 > fi2) && (neighbor_id_2 > fi2)) {
                face const& old_neighbor_1 = faces_[neighbor_id_1];
                face const& old_neighbor_2 = faces_[neighbor_id_2];
                const face present_neighbor_1 = old_neighbor_1.updated_to_time(f.base_time_);
                const face present_neighbor_2 = old_neighbor_2.updated_to_time(f.base_time_);
                if (faux_optional<time_type> collision_time = when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(
                  present_other_face, present_neighbor_1, present_neighbor_2, f)) {
                  if (vertex_is_in_bounded_face__hack(*collision_time, present_other_face, present_neighbor_1, present_neighbor_2, f)) {
                    next_events_.push(shared_ptr<event>(new vertex_face_collision(
                        *collision_time,
                                  fi2, old_other_face.revision_number_,
                        neighbor_id_1, old_neighbor_1.revision_number_,
                        neighbor_id_2, old_neighbor_2.revision_number_,
                                   fi,              f.revision_number_)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  void hack_insert_rock(vector3<distance> loc) {
    face_idx_type first_face_idx = faces_.size();
    region_idx_type region_idx = regions_.size();

    static const vector3<mpz> diffs[4] = {
      vector3<mpz>(0, 0, 10000),
      vector3<mpz>(0, 8300, -4500),
      vector3<mpz>(6000, -4000, -5000),
      vector3<mpz>(-6200, -4200, -4300)
    };

    regions_.push_back(region());
    region& air = regions_[0];
    region& r = regions_.back();
    for (int i = 0; i < 4; ++i) {
      faces_.push_back(face());
      face& f = faces_.back();
      f.base_time_ = 0;
      f.ABC = diffs[i] + vector3<mpz>(first_face_idx+1, first_face_idx+1, first_face_idx+1);
      f.D = loc.dot<mpz>(f.ABC) + f.ABC.magnitude_using<mpz>()*100*(centi*meters)*identity(distance_units/(centi*meters));
      f.D_velocity = 0;
      f.D_acceleration = f.ABC.dot<mpz>(gravity_acceleration);
      f.neighboring_regions_.push_back(0);
      f.neighboring_regions_.push_back(region_idx);
      r.faces_.push_back(first_face_idx + i);
      air.faces_.push_back(first_face_idx + i);
      //r.vertices_.push_back(first_vertex_idx + i);
      for (int j = 0; j < 4; ++j) {
        if (j != i) f.neighboring_faces_.push_back(first_face_idx + j);
      }
    }
  }
  void hack_insert_fixed_cube(vector3<distance> loc) {
    face_idx_type first_face_idx = faces_.size();
    region_idx_type region_idx = regions_.size();
    regions_.push_back(region());
    region& air = regions_[0];
    region& r = regions_.back();
    for (int i = 0; i < 6; ++i) {
      faces_.push_back(face());
      face& f = faces_.back();
      f.base_time_ = 0;
      f.ABC = vector3<mpz>(cardinal_direction_vectors[i]);
      f.D = loc.dot<mpz>(f.ABC) + f.ABC.magnitude_using<mpz>()*1000*(centi*meters)*identity(distance_units/(centi*meters));
      f.D_velocity = 0;
      f.D_acceleration = 0;
      f.neighboring_regions_.push_back(0);
      f.neighboring_regions_.push_back(region_idx);
      r.faces_.push_back(first_face_idx + i);
      air.faces_.push_back(first_face_idx + i);
      //r.vertices_.push_back(first_vertex_idx + i);
      for (int j = 0; j < 6; ++j) {
        // Hack - relying on the order of the cardinal directions
        if ((j % 3) != (i % 3)) f.neighboring_faces_.push_back(first_face_idx + j);
      }
    }
  }
  
public:
  grand_structure_of_lasercake() {
    regions_.push_back(region()); // the air
    hack_insert_rock(vector3<mpz>(2, 7, 11)*meters*identity(distance_units/meters));
    hack_insert_rock(vector3<mpz>(2, 14, 11)*meters*identity(distance_units/meters));
    hack_insert_rock(vector3<mpz>(15, 14, 21)*meters*identity(distance_units/meters));
    hack_insert_fixed_cube(vector3<mpz>(8, 10, -11)*meters*identity(distance_units/meters));

    // TODO don't duplicate events here.
    for (face_idx_type fi = 0; fi < faces_.size(); ++fi) {
      insert_events_involving(fi);
    }
    
    this->debug_check_consistency();
  }
  void debug_check_consistency()const {
    // region internal consistency
    for(region const& r : regions_) {
      // no two vertices of a region have the same identity
      /*for(size_t j = 0; j < r.vertices_.size(); ++j) {
        for(size_t k = j + 1; k != r.vertices_.size(); ++k) {
          assert(r.vertices_[j] != r.vertices_[k]);
        }
      }*/
      // no two faces of a region have the same identity
      for(size_t j = 0; j < r.faces_.size(); ++j) {
        for(size_t k = j + 1; k != r.faces_.size(); ++k) {
          assert(r.faces_[j] != r.faces_[k]);
        }
      }
    }
    // reference in-bound-ness
    for(region const& r : regions_) {
      //for(vertex_idx_type vi : r.vertices_) { assert(vi < vertices_.size()); }
      for(face_idx_type fi : r.faces_) { assert(fi < faces_.size()); }
    }
    for(face const& f : faces_) {
      for(vertex_idx_type fi : f.neighboring_faces_) { assert(fi < faces_.size()); }
      for(region_idx_type ri : f.neighboring_regions_) { assert(ri < regions_.size()); }
    }
    // region-face data consistency
    for(size_t i = 0; i != regions_.size(); ++i) {
      region const& r = regions_[i];
      for(face_idx_type fi : r.faces_) {
        face const& f = faces_[fi];
        bool found_reference = false;
        for (region_idx_type ri : f.neighboring_regions_) {
          if (ri == i) {
            found_reference = true;
          }
        }
        assert(found_reference);
      }
    }
    for(size_t i = 0; i != faces_.size(); ++i) {
      face const& f = faces_[i];
      for(region_idx_type ri : f.neighboring_regions_) {
        region const& r = regions_[ri];
        bool found_reference = false;
        for (face_idx_type fi : r.faces_) {
          if (fi == i) {
            found_reference = true;
          }
        }
        assert(found_reference);
      }
    }
    // TODO: check stuff about edges as well.
  }
  // TODO thing-ness e.g. robots
  //void player_input_becomes(time_type when, );
  //void insert_event(time_type when, );
  gl_triangles display(time_type when, vector3<distance> where) {
    //assert(when < next_event_time)
    gl_triangles triangles;
    for(face const& f : faces_) {
      gl_polygon polygon;
      const face present_face = f.updated_to_time(when);
      for(size_t i = 0; i < f.neighboring_faces_.size(); ++i) {
        const size_t j = (i+1)%f.neighboring_faces_.size();
        const face present_neighbor_1 = faces_[f.neighboring_faces_[i]].updated_to_time(when);
        const face present_neighbor_2 = faces_[f.neighboring_faces_[j]].updated_to_time(when);
        const vector3<distance> loc = approx_loc_of_triple_intersection_of_up_to_date_faces(present_face, present_neighbor_1, present_neighbor_2);
        polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          get_primitive_float(loc.x/distance_units),
          get_primitive_float(loc.y/distance_units),
          get_primitive_float(loc.z/distance_units),
          gl_data_format::color(0xffff0080)));
        //std::cerr << vertices[i] << '\n';
      }
      push_wireframe_polygon(triangles, 1e9, polygon);
    }
    sort_gl_triangles_far_to_near(
      glm::vec3(where.x/distance_units, where.y/distance_units, where.z/distance_units),
      triangles);
    return triangles;
  }

  time_type time_of_next_event()const {
    if (next_events_.empty()) return 0;
    return next_events_.top()->when_event_occurs_;
  }
private:
  
  void do_next_event() {
    const shared_ptr<event> ev = next_events_.top();
    next_events_.pop();
    if (const shared_ptr<collision> c = dynamic_pointer_cast<collision>(ev)) {
      if (
        (faces_[c->f1_].revision_number_ == c->r1_) &&
        (faces_[c->f2_].revision_number_ == c->r2_) &&
        (faces_[c->f3_].revision_number_ == c->r3_) &&
        (faces_[c->f4_].revision_number_ == c->r4_)
      ) {
        if (const shared_ptr<vertex_face_collision> vfc = dynamic_pointer_cast<vertex_face_collision>(c)) {
          // The vertex stays existent; N new vertices also appear where it is,
          // where N is the number of edges connected to that vertex.
          // This requires us to split the faces incident to this vertex into
          // triangles (N more faces).
          // The face (which is a triangle) disappears and is replaced with
          // 3+N faces (triangles), plus N faces inside the ring of new vertices.
        }
        if (const shared_ptr<edge_edge_collision> eec = dynamic_pointer_cast<edge_edge_collision>(c)) {
          // The two edges do not share a vertex, because we make sure to
          // generate only a vertex--face collision in that case.
        }
      }
    }
  }
  
};



void do_gl(grand_structure_of_lasercake& simulated_world, uint64_t frame, time_type when) {
  //std::cerr<<"hi.\n";
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  //gluPerspective(90, 1, 1, 1000);
  //gluLookAt(20,20,20,0,0,0,0,0,1);
  //gluLookAt(0,0,0,1,1,1,0,0,1);
  const double wid = 100e9;
  gluPerspective(80, 1, 1e9, 2.0*wid);
  const vector3<double> viewcenter(0+wid*std::cos(double(frame) / 200), 0+wid*std::sin(double(frame) / 200), 0);
  gluLookAt(viewcenter.x, viewcenter.y, viewcenter.z,
            0,0,0, 0,0,1);

  gl_triangles data = simulated_world.display(when, vector3<distance>(viewcenter*distance_units));
  if(const size_t count = data.size()*3) {
    GLuint triangles_VBO_name;
    glGenBuffersARB(1, &triangles_VBO_name);
    glBindBufferARB(GL_ARRAY_BUFFER, triangles_VBO_name);
    glBufferDataARB(GL_ARRAY_BUFFER, count*sizeof(gl_data_format::vertex_with_color), &data[0], GL_STREAM_DRAW);
    glInterleavedArrays(GL_C4UB_V3F, 0, (void*)(0));
    glDrawArrays(GL_TRIANGLES, 0, count);
    glDeleteBuffersARB(1, &triangles_VBO_name);
  }
}











namespace chrono = boost::chrono;

typedef int64_t microseconds_t;

microseconds_t get_this_thread_microseconds() {
#if defined(BOOST_CHRONO_HAS_THREAD_CLOCK)
  return chrono::duration_cast<chrono::microseconds>(chrono::thread_clock::now().time_since_epoch()).count();
#else
  return 0;
#endif
}

static SDL_Surface *gScreen;

static void initAttributes ()
{
    // Setup attributes we want for the OpenGL context
    
    int value;
    
    // Don't set color bit sizes (SDL_GL_RED_SIZE, etc)
    //    Mac OS X will always use 8-8-8-8 ARGB for 32-bit screens and
    //    5-5-5 RGB for 16-bit screens
    
    // Request a 16-bit depth buffer (without this, there is no depth buffer)
    value = 16;
    SDL_GL_SetAttribute (SDL_GL_DEPTH_SIZE, value);
    
    
    // Request double-buffered OpenGL
    //     The fact that windows are double-buffered on Mac OS X has no effect
    //     on OpenGL double buffering.
    value = 1;
    SDL_GL_SetAttribute (SDL_GL_DOUBLEBUFFER, value);
}

static void printAttributes ()
{
    // Print out attributes of the context we created
    int nAttr;
    int i;
    
    int  attr[] = { SDL_GL_RED_SIZE, SDL_GL_BLUE_SIZE, SDL_GL_GREEN_SIZE,
                    SDL_GL_ALPHA_SIZE, SDL_GL_BUFFER_SIZE, SDL_GL_DEPTH_SIZE };
                    
    const char *desc[] = { "Red size: %d bits\n", "Blue size: %d bits\n", "Green size: %d bits\n",
                     "Alpha size: %d bits\n", "Color buffer size: %d bits\n", 
                     "Depth bufer size: %d bits\n" };

    nAttr = sizeof(attr) / sizeof(int);
    
    for (i = 0; i < nAttr; i++) {
    
        int value;
        SDL_GL_GetAttribute ((SDL_GLattr)attr[i], &value);
        printf (desc[i], value);
    } 
}

static void createSurface (int fullscreen)
{
    Uint32 flags = 0;
    
    flags = SDL_OPENGL;
    if (fullscreen)
        flags |= SDL_FULLSCREEN;
    
    // Create window
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    gScreen = SDL_SetVideoMode (640, 640, 32, flags);
    if (gScreen == NULL) {
		
        fprintf (stderr, "Couldn't set 640x640 OpenGL video mode: %s\n",
                 SDL_GetError());
		SDL_Quit();
		exit(2);
	}
}


static void mainLoop (std::string /*scenario*/)
{
  SDL_Event event;
  int done = 0;
  int p_mode = 0;
  bool moving = false;
large_fast_noncrypto_rng rng(time(NULL));

  
  
int frame = 0;
time_type when = 0;

  grand_structure_of_lasercake simulated_world;
  
  while ( !done ) {

    /* Check for events */
    while ( SDL_PollEvent (&event) ) {
      switch (event.type) {
        case SDL_MOUSEMOTION:
          break;
          
        case SDL_MOUSEBUTTONDOWN:
          break;
          
        case SDL_KEYDOWN:
          if(event.key.keysym.sym == SDLK_p) ++p_mode;
          if(event.key.keysym.sym == SDLK_m) moving = !moving;
          if(event.key.keysym.sym == SDLK_e) when = simulated_world.time_of_next_event();
          //if(event.key.keysym.sym == SDLK_r) ++view_dist;
          //if(event.key.keysym.sym == SDLK_f) --view_dist;
          if(event.key.keysym.sym != SDLK_ESCAPE)break;
          
        case SDL_QUIT:
          done = 1;
          break;
          
        default:
          break;
      }
    }
    if(p_mode == 1)continue;
    if(p_mode > 1)--p_mode;
    __attribute__((unused)) int before_drawing = SDL_GetTicks();

    
    __attribute__((unused)) int before_GL = SDL_GetTicks();
    
    do_gl(simulated_world, frame, when);
    glFinish();	
    SDL_GL_SwapBuffers();
   
    __attribute__((unused)) int before_processing = SDL_GetTicks();
    
    //doing stuff code here
    
    ++frame;
    if (moving) when += int64_t(1000000000LL)*pico*seconds;
    
    
    __attribute__((unused)) int after = SDL_GetTicks();
    //std::cerr << (after - before_processing) << ", " << (before_GL - before_drawing) << ", " << (before_processing - before_GL) << "\n";

//    SDL_Delay(50);
  }
}

}
using namespace YO;

int main(int argc, char *argv[])
{
	// Init SDL video subsystem
	if ( SDL_Init (SDL_INIT_VIDEO | SDL_INIT_NOPARACHUTE) < 0 ) {
		
        fprintf(stderr, "Couldn't initialize SDL: %s\n",
			SDL_GetError());
		exit(1);
	}

    // Set GL context attributes
    initAttributes ();
    
    // Create GL context
    createSurface (0);
    
    // Get GL context attributes
    printAttributes ();
    
    const GLenum glew_init_err = glewInit();
    if(glew_init_err != GLEW_OK) {
      throw "glew failed";
    }
    
    // Draw, get events...
    if (argc < 2) {
      std::cerr << "You didn't give an argument saying which scenario to use! Using default value...\n";
      mainLoop ("default");
    }
    else mainLoop (argv[1]);
    
    // Cleanup
	SDL_Quit();
	
    return 0;
}
