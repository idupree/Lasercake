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
typedef typename units_prod<milli_t, seconds_t>::type time_units_t;
constexpr time_units_t time_units = time_units_t();
typedef physical_quantity<mpz, time_units_t> time_type;

typedef typename units_prod<nano_t, meters_t>::type distance_units_t;
typedef typename units_prod<distance_units_t, typename units_pow<time_units_t, -1>::type>::type velocity_units_t;
typedef typename units_prod<velocity_units_t, typename units_pow<time_units_t, -1>::type>::type acceleration_units_t;
typedef typename units_prod<acceleration_units_t, dim::ratio<2>>::type acceleration_coefficient_units_t;
constexpr distance_units_t distance_units = distance_units_t();
constexpr velocity_units_t velocity_units = velocity_units_t();
constexpr acceleration_units_t acceleration_units = acceleration_units_t();
constexpr acceleration_coefficient_units_t acceleration_coefficient_units = acceleration_coefficient_units_t();
typedef physical_quantity<mpz, distance_units_t> distance;
typedef physical_quantity<mpz, velocity_units_t> velocity1d;
typedef physical_quantity<mpz, acceleration_units_t> acceleration1d;
typedef physical_quantity<mpz, acceleration_coefficient_units_t> acceleration_coefficient;

// Standard (Earth-equivalent) gravity: precisely 9.80665 m/s2
const acceleration1d gravity_acceleration_magnitude = mpz(9806650) * (micro*meters) / (seconds*seconds) * identity(distance_units / (micro*meters)) / (identity(time_units * time_units / seconds / seconds));
const vector3<acceleration1d> gravity_acceleration(0, 0, -gravity_acceleration_magnitude);

typedef uint32_t region_idx_type;
typedef uint32_t face_idx_type;
typedef uint32_t vertex_idx_type;
typedef uint64_t revision_number_type;


struct face {
  face():base_time_(0),revision_number_(0),ABC(0),D(0),D_velocity(0),D_acc_coeff(0){}
  time_type base_time_;
  revision_number_type revision_number_;
  vector3<mpz> ABC;
  distance D;
  velocity1d D_velocity;
  acceleration_coefficient D_acc_coeff;
  std::vector<face_idx_type> neighboring_faces_;
  std::vector<region_idx_type> neighboring_regions_;
  face updated_to_time(time_type t)const {
    assert(t >= base_time_);
    face result(*this);
    const time_type relative_time = t - base_time_;
    //const auto relative_timem = relative_time/identity(units_factor<1, 1000000000>());
    //TODO deal with times somehow more correctly
    assert((identity(acceleration_units / acceleration_coefficient_units) / 2) == (1 * acceleration_units / acceleration_coefficient_units));
    result.D          += D_velocity * relative_time
                       + D_acc_coeff * relative_time * relative_time * (identity(acceleration_units / acceleration_coefficient_units) / 2);
    result.D_velocity += D_acc_coeff * relative_time * identity(acceleration_units / acceleration_coefficient_units);
    result.base_time_ = t;
    ++result.revision_number_;
    return result;
  }
};

struct silly_rational_loc {
  vector3<distance> nums;
  mpz shared_denom;
  bool operator<(silly_rational_loc const& o)const {
    mpz negative_factor = (is_negative(shared_denom) == is_negative(o.shared_denom)) ? 1 : -1;
    for (which_dimension_type dim = 0; dim < 3; ++dim) {
      const auto prod = negative_factor * (nums(X) * o.shared_denom - o.nums(X) * shared_denom);
      if (prod < 0) return true;
      if (prod > 0) return false;
    }
    return false;
  }
  bool operator==(silly_rational_loc const& o)const {
    for (which_dimension_type dim = 0; dim < 3; ++dim) {
      const auto prod = (nums(X) * o.shared_denom - o.nums(X) * shared_denom);
      if (prod != 0) return false;
    }
    return true;
  }
  bool operator!=(silly_rational_loc const& o)const {
    return !(o == *this);
  }
};

mpz scalar_triple_product(vector3<mpz> v1, vector3<mpz> v2, vector3<mpz> v3) {
  return v1(X)*v2(Y)*v3(Z) - v1(X)*v3(Y)*v2(Z) + v2(X)*v3(Y)*v1(Z) - v2(X)*v1(Y)*v3(Z) + v3(X)*v1(Y)*v2(Z) - v3(X)*v2(Y)*v1(Z);
}
silly_rational_loc exact_loc_of_triple_intersection_of_up_to_date_faces(face const& f1, face const& f2, face const& f3) {
  mpz p = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  assert(p != 0);
  return silly_rational_loc{ vector3<distance>(
    f1.D * (f2.ABC(1)*f3.ABC(2) - f3.ABC(1)*f2.ABC(2)) + f2.D * (f3.ABC(1)*f1.ABC(2) - f1.ABC(1)*f3.ABC(2)) + f3.D * (f1.ABC(1)*f2.ABC(2) - f2.ABC(1)*f1.ABC(2)),
    f1.D * (f2.ABC(2)*f3.ABC(0) - f3.ABC(2)*f2.ABC(0)) + f2.D * (f3.ABC(2)*f1.ABC(0) - f1.ABC(2)*f3.ABC(0)) + f3.D * (f1.ABC(2)*f2.ABC(0) - f2.ABC(2)*f1.ABC(0)),
    f1.D * (f2.ABC(0)*f3.ABC(1) - f3.ABC(0)*f2.ABC(1)) + f2.D * (f3.ABC(0)*f1.ABC(1) - f1.ABC(0)*f3.ABC(1)) + f3.D * (f1.ABC(0)*f2.ABC(1) - f2.ABC(0)*f1.ABC(1))
                 ), p};
}
vector3<distance> approx_loc_of_triple_intersection_of_up_to_date_faces(face const& f1, face const& f2, face const& f3) {
  silly_rational_loc l = exact_loc_of_triple_intersection_of_up_to_date_faces(f1, f2, f3);
  return l.nums / l.shared_denom;
}
// Note: This function will need to change a LOT more than the above function if the change stuff gets more complicated.
// Perhaps it should just be approximated experimentally (dv/dt for very small t).
vector3<velocity1d> approx_velocity_of_triple_intersection_of_up_to_date_faces(face const& f1, face const& f2, face const& f3) {
  mpz p = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  assert(p != 0);
  return vector3<velocity1d>(
    f1.D_velocity * (f2.ABC(1)*f3.ABC(2) - f3.ABC(1)*f2.ABC(2)) + f2.D_velocity * (f3.ABC(1)*f1.ABC(2) - f1.ABC(1)*f3.ABC(2)) + f3.D_velocity * (f1.ABC(1)*f2.ABC(2) - f2.ABC(1)*f1.ABC(2)),
    f1.D_velocity * (f2.ABC(2)*f3.ABC(0) - f3.ABC(2)*f2.ABC(0)) + f2.D_velocity * (f3.ABC(2)*f1.ABC(0) - f1.ABC(2)*f3.ABC(0)) + f3.D_velocity * (f1.ABC(2)*f2.ABC(0) - f2.ABC(2)*f1.ABC(0)),
    f1.D_velocity * (f2.ABC(0)*f3.ABC(1) - f3.ABC(0)*f2.ABC(1)) + f2.D_velocity * (f3.ABC(0)*f1.ABC(1) - f1.ABC(0)*f3.ABC(1)) + f3.D_velocity * (f1.ABC(0)*f2.ABC(1) - f2.ABC(0)*f1.ABC(1))
                 ) / p;
}

// TODO : What about the cases where two of the planes are parallel? That's legit (handle them separately? they're easier)
std::vector<time_type> how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(face const& f1, face const& f2, face const& f3, face const& f4/*, bool hack_recurse_test = true*/) {
#if 0
  if (hack_recurse_test) {
    const faux_optional<time_type> a = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1, f2, f3, f4, false);
    const faux_optional<time_type> b = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1, f2, f4, f3, false);
    const faux_optional<time_type> c = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1, f3, f2, f4, false);
    const faux_optional<time_type> d = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1, f3, f4, f2, false);
    const faux_optional<time_type> e = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1, f4, f3, f2, false);
    const faux_optional<time_type> f = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1, f4, f2, f3, false);
    const faux_optional<time_type> g = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f2, f1, f3, f4, false);
    const faux_optional<time_type> h = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f2, f1, f4, f3, false);
    const faux_optional<time_type> i = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f2, f3, f1, f4, false);
    const faux_optional<time_type> j = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f2, f3, f4, f1, false);
    const faux_optional<time_type> k = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f2, f4, f3, f1, false);
    const faux_optional<time_type> l = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f2, f4, f1, f3, false);
    const faux_optional<time_type> m = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f3, f2, f1, f4, false);
    const faux_optional<time_type> n = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f3, f2, f4, f1, false);
    const faux_optional<time_type> o = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f3, f1, f2, f4, false);
    const faux_optional<time_type> p = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f3, f1, f4, f2, false);
    const faux_optional<time_type> q = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f3, f4, f1, f2, false);
    const faux_optional<time_type> r = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f3, f4, f2, f1, false);
    const faux_optional<time_type> s = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f4, f2, f3, f1, false);
    const faux_optional<time_type> t = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f4, f2, f1, f3, false);
    const faux_optional<time_type> u = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f4, f3, f2, f1, false);
    const faux_optional<time_type> v = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f4, f3, f1, f2, false);
    const faux_optional<time_type> w = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f4, f1, f3, f2, false);
    const faux_optional<time_type> x = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f4, f1, f2, f3, false);
    if (a) {
      assert(*a == *b);
      assert(*a == *c);
      assert(*a == *d);
      assert(*a == *e);
      assert(*a == *f);
      assert(*a == *g);
      assert(*a == *h);
      assert(*a == *i);
      assert(*a == *j);
      assert(*a == *k);
      assert(*a == *l);
      assert(*a == *m);
      assert(*a == *n);
      assert(*a == *o);
      assert(*a == *p);
      assert(*a == *q);
      assert(*a == *r);
      assert(*a == *s);
      assert(*a == *t);
      assert(*a == *u);
      assert(*a == *v);
      assert(*a == *w);
      assert(*a == *x);
    }
    else {
      assert(!b);
      assert(!c);
      assert(!d);
      assert(!e);
      assert(!f);
      assert(!g);
      assert(!h);
      assert(!i);
      assert(!j);
      assert(!k);
      assert(!l);
      assert(!m);
      assert(!n);
      assert(!o);
      assert(!p);
      assert(!q);
      assert(!r);
      assert(!s);
      assert(!t);
      assert(!u);
      assert(!v);
      assert(!w);
      assert(!x);
    }
  }
#endif
  // When the 4x4 determinant is 0.
  // That's
  // + D1 * scalar_triple_product(f2.ABC, f3.ABC, f4.ABC)
  // - D2 * scalar_triple_product(f3.ABC, f4.ABC, f1.ABC)
  // + D3 * scalar_triple_product(f4.ABC, f1.ABC, f2.ABC)
  // - D4 * scalar_triple_product(f1.ABC, f2.ABC, f3.ABC)
  // But Dn is Dn + dDn*t + .5ddDn*t^2
  // so
  std::vector<time_type> result;
  const mpz t1 = scalar_triple_product(f2.ABC, f3.ABC, f4.ABC);
  const mpz t2 = scalar_triple_product(f3.ABC, f4.ABC, f1.ABC);
  const mpz t3 = scalar_triple_product(f4.ABC, f1.ABC, f2.ABC);
  const mpz t4 = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  if ((t1 == 0) || (t2 == 0) || (t3 == 0) || (t4 == 0)) {
    //std::cerr << "Warning: Linearly dependent normals passed to how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point()\n";
    return result;
  }
  // in at^2 + bt + c = 0
  acceleration_coefficient a = (f1.D_acc_coeff*t1 - f2.D_acc_coeff*t2 + f3.D_acc_coeff*t3 - f4.D_acc_coeff*t4);
  velocity1d               b =  f1.D_velocity *t1 - f2.D_velocity *t2 + f3.D_velocity *t3 - f4.D_velocity *t4;
  distance                 c =  f1.D          *t1 - f2.D          *t2 + f3.D          *t3 - f4.D          *t4;

  acceleration1d a_times_2 = a * identity(acceleration_units / acceleration_coefficient_units);
  
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
    return result;
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
          return result;
        }
        else return result;
      }
      if (b < 0) {
        b = -b;
        c = -c;
      }
      if (c > 0) return result;
      const time_type zero = divide(-c, b, strat);
      //std::cerr << a_times_2 << ", " << b << ", " << c << ": " << zero << ", " << get(zero,time_units)*get(zero,time_units)*get(a_times_2,typename units_prod<distance_units_t, dim::second<(-2)> >::type()) / 2 + get(zero,time_units)*get(b,typename units_prod<distance_units_t, dim::second<(-1)> >::type())*mpz(1e12) + get(c,distance_units)*mpz(1e24) << "\n";
      result.push_back(zero);
      return result;
    }
    else {
      if (a_times_2 < 0) {
        a_times_2 = -a_times_2;
        b = -b;
        c = -c;
      }
      // we want the earlier time, but not if it's negative
      const velocity1d sqrt_disc = isqrt(discriminant);
      const velocity1d  lesser_numerator = -b - sqrt_disc - ((sqrt_disc*sqrt_disc == discriminant) ? 0 : 1)*distance_units/time_units;
      if ( lesser_numerator >= 0)  {
        const time_type zero = divide(lesser_numerator, a_times_2, strat);
        //std::cerr << a_times_2 << ", " << b << ", " << c << ": " << zero << ", " << get(zero,time_units)*get(zero,time_units)*get(a_times_2,typename units_prod<distance_units_t, dim::second<(-2)> >::type()) / 2 + get(zero,time_units)*get(b,typename units_prod<distance_units_t, dim::second<(-1)> >::type())*mpz(1e12) + get(c,distance_units)*mpz(1e24) << "\n";
        result.push_back(zero);
      }
      const velocity1d greater_numerator = -b + sqrt_disc;
      if (greater_numerator >= 0) {
        const time_type zero = divide(greater_numerator, a_times_2, strat);
        //std::cerr << a_times_2 << ", " << b << ", " << c << ": " << zero << ", " << get(zero,time_units)*get(zero,time_units)*get(a_times_2,typename units_prod<distance_units_t, dim::second<(-2)> >::type()) / 2 + get(zero,time_units)*get(b,typename units_prod<distance_units_t, dim::second<(-1)> >::type())*mpz(1e12) + get(c,distance_units)*mpz(1e24) << "\n";
        result.push_back(zero);
      }
      return result;
    }
  }
}
std::vector<time_type> when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(face const& f1, face const& f2, face const& f3, face const& f4) {
  assert(f1.base_time_ == f2.base_time_);
  assert(f1.base_time_ == f3.base_time_);
  assert(f1.base_time_ == f4.base_time_);
  std::vector<time_type> result = how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(f1,f2,f3,f4);
  for (time_type& r : result) {
    r += f1.base_time_;
  }
  return result;
}

enum region_contents {
  AIR,
  ROCK
};

struct region {
  // Regions are arbitrary, non-convex polyhedra which can have holes.
  //std::vector<vertex_idx_type, tetrahedron_sides> vertices_;
  std::vector<face_idx_type> faces_;
  region_contents contents;
  region(region_contents contents):contents(contents){}
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

struct face_triple {
  // standardized to be in increasing order
  face_idx_type f1_;
  face_idx_type f2_;
  face_idx_type f3_;
  bool operator==(face_triple const& o)const{return (f1_ == o.f1_) && (f2_ == o.f2_) && (f3_ == o.f3_);}
  bool operator!=(face_triple const& o)const{return (f1_ != o.f1_) || (f2_ != o.f2_) || (f3_ != o.f3_);}
};

} // namespace YO 
namespace HASH_NAMESPACE {
  template<> struct hash<YO::face_triple> {
    inline size_t operator()(YO::face_triple const& t) const {
      size_t seed = 0;
      boost::hash_combine(seed, t.f1_);
      boost::hash_combine(seed, t.f2_);
      boost::hash_combine(seed, t.f3_);
      return seed;
    }
  };
}
namespace YO {

template<typename JoinedType>
struct segment_collection {
  // TODO: Evaluate data structure choice.
  // This structure will usually contain 3 segments at most, and very rarely more than, let's say, eight.
  std::multimap<JoinedType, JoinedType> contents_;
  void insert(std::pair<JoinedType, JoinedType> const& o) {
    contents_.insert(o);
    contents_.insert(std::make_pair(o.second, o.first));
  }
  void erase(std::pair<JoinedType, JoinedType> const& o) {
    {
      auto r = contents_.equal_range(o.first);
      for (auto i = r.first; i != r.second; ++i) { if (i->second == o.second) { contents_.erase(i); break; } }
    }
    {
      auto r = contents_.equal_range(o.second);
      for (auto i = r.first; i != r.second; ++i) { if (i->second == o.first) { contents_.erase(i); break; } }
    }
  }
  bool empty() {
    return contents_.empty();
  }
  JoinedType const* arbitrary_other_end(JoinedType const& t)const {
    auto result = contents_.find(t);
    if (result == contents_.end()) return nullptr;
    else return &result->second;
  }
  JoinedType arbitrary_triple()const {
    return contents_.begin()->first;
  }
};

// Maybe the structures should be required to follow some kind of right-hand-rule
// order? Then asserts could check that it is still correct, and OpenGL would
// be more easily able to draw directional surfaces.
class grand_structure_of_lasercake {
  std::vector<region> regions_;
  std::vector<face> faces_;
  time_type present_time_;
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

      // Hack - sometimes the line can read as right in some dimensions and wrong in others.
      // This mostly happens if a face is very nearly parallel to one of the lines.
      // If it does, eliminate it rather than including it, because it might be wrong and isn't super important if it's right.
      if ((((end1.x-approx_cross_location.x) * (end2.x-approx_cross_location.x)) > 0)
        || (((end1.y-approx_cross_location.y) * (end2.y-approx_cross_location.y)) > 0)
        || (((end1.z-approx_cross_location.z) * (end2.z-approx_cross_location.z)) > 0)) {
        return false;
      }
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
          const mpz sign_correction = (pv1_to_v(dim1) > 0) ? 1 : -1;
          if (pv1_to_pv2(dim2)*pv1_to_v(dim1)*sign_correction >= pv1_to_v(dim2)*pv1_to_pv2(dim1)*sign_correction) {
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
              std::vector<time_type> collision_times = when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(
                  f, present_neighbor_1, present_neighbor_2, present_other_face);
              for (auto collision_time : collision_times) {
                next_events_.push(shared_ptr<event>(new vertex_face_collision(
                    collision_time,
                               fi,              f.revision_number_,
                    neighbor_id_1, old_neighbor_1.revision_number_,
                    neighbor_id_2, old_neighbor_2.revision_number_,
                              fi2, old_other_face.revision_number_)));
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
                  std::vector<time_type> collision_times = when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(
                      f, present_neighbor_1, present_other_face, present_other_neighbor);
                  for (auto collision_time : collision_times) {
                    next_events_.push(shared_ptr<event>(new edge_edge_collision(
                        collision_time,
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
                std::vector<time_type> collision_times = when_will_planes_of_up_to_date_faces_be_coincident_at_a_point(
                    present_other_face, present_neighbor_1, present_neighbor_2, f);
                for (auto collision_time : collision_times) {
                  next_events_.push(shared_ptr<event>(new vertex_face_collision(
                      collision_time,
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
  
  void hack_insert_rock(vector3<distance> loc) {
    face_idx_type first_face_idx = faces_.size();
    region_idx_type region_idx = regions_.size();

    static const vector3<mpz> diffs[4] = {
      vector3<mpz>(0, 0, 10000),
      vector3<mpz>(0, 8300, -4500),
      vector3<mpz>(6000, -4000, -5000),
      vector3<mpz>(-6200, -4200, -4300)
    };

    regions_.push_back(region(ROCK));
    region& air = regions_[0];
    region& r = regions_.back();
    for (int i = 0; i < 4; ++i) {
      faces_.push_back(face());
      face& f = faces_.back();
      f.base_time_ = 0;
      f.ABC = diffs[i] + vector3<mpz>(first_face_idx+1, first_face_idx+1, first_face_idx+1);
      f.D = loc.dot<mpz>(f.ABC) + f.ABC.magnitude_using<mpz>()*100*(centi*meters)*identity(distance_units/(centi*meters));
      f.D_velocity = 0;
      f.D_acc_coeff = divide(f.ABC.dot<mpz>(gravity_acceleration), identity(acceleration_units / acceleration_coefficient_units), rounding_strategy<round_to_nearest_with_ties_rounding_to_even>());
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
    regions_.push_back(region(ROCK));
    region& air = regions_[0];
    region& r = regions_.back();
    for (int i = 0; i < 6; ++i) {
      faces_.push_back(face());
      face& f = faces_.back();
      f.base_time_ = 0;
      f.ABC = vector3<mpz>(cardinal_direction_vectors[i]);
      f.D = loc.dot<mpz>(f.ABC) + f.ABC.magnitude_using<mpz>()*1000*(centi*meters)*identity(distance_units/(centi*meters));
      f.D_velocity = 0;
      f.D_acc_coeff = 0;
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
    regions_.push_back(region(AIR));
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
  void update_to_time(time_type when, bool do_events = true) {
    assert(when >= present_time_);
    while(!next_events_.empty() && next_events_.top()->when_event_occurs_ < when) {
      if (do_events) do_next_event();
      else next_events_.pop();
    }
    present_time_ = when;
  }

  
  std::vector<std::vector<silly_rational_loc>> find_self_overlaps(region const& r)const {
    segment_collection<silly_rational_loc> face_overlaps;

    
    //for (each pair (f1, f2) of nonadjacent faces of r) {
    for (size_t f1ii =        0; f1ii < r.faces_.size(); ++f1ii) {
      const face_idx_type f1i = r.faces_[f1ii];
      const face f1 = faces_[f1i].updated_to_time(present_time_);
    for (size_t f2ii = f1ii + 1; f2ii < r.faces_.size(); ++f2ii) {
      const face_idx_type f2i = r.faces_[f2ii];
      bool is_neighboring = false;
      for (auto q : f1.neighboring_faces_) { if (q == f2i) { is_neighboring = true; break; } }
      if (is_neighboring) { continue; }
      const face f2 = faces_[f2i].updated_to_time(present_time_);
      
      
      std::priority_queue<silly_rational_loc> f1_transition_points;
      std::priority_queue<silly_rational_loc> f2_transition_points;

      for (int i = 0; i < 2; ++i) {
        face const& fA = i ? f2 : f1;
        face const& fB = i ? f1 : f2;
        std::priority_queue<silly_rational_loc>& fA_transition_points = i ? f2_transition_points : f1_transition_points;

        
        //for (each triple (fAn1, fAn2, fAn3) of consecutive faces adjacent to fA, constituting a pair of adjacent vertices (v1, v2) of fA) {
        for (size_t fAn1ii = 0; fAn1ii < fA.neighboring_faces_.size(); ++fAn1ii) {
          const size_t fAn2ii = (fAn1ii + 1) % fA.neighboring_faces_.size();
          const size_t fAn3ii = (fAn1ii + 2) % fA.neighboring_faces_.size();
          const face_idx_type fAn1i = fA.neighboring_faces_[fAn1ii];
          const face_idx_type fAn2i = fA.neighboring_faces_[fAn2ii];
          const face_idx_type fAn3i = fA.neighboring_faces_[fAn3ii];
          const face fAn1 = faces_[fAn1i].updated_to_time(present_time_);
          const face fAn2 = faces_[fAn2i].updated_to_time(present_time_);
          const face fAn3 = faces_[fAn3i].updated_to_time(present_time_);
          const silly_rational_loc v1 = exact_loc_of_triple_intersection_of_up_to_date_faces(fA, fAn1, fAn2);
          const silly_rational_loc v2 = exact_loc_of_triple_intersection_of_up_to_date_faces(fA, fAn2, fAn3);
          
          //if (v1 and v2 are on opposite sides of fB) {
          const mpz d1neg = is_negative(v1.shared_denom) ? -1 : 1;
          const mpz d2neg = is_negative(v2.shared_denom) ? -1 : 1;
          if ((v1.nums.dot<mpz>(fB.ABC) * d1neg < fB.D * v1.shared_denom * d1neg) != (v2.nums.dot<mpz>(fB.ABC) * d2neg < fB.D * v2.shared_denom * d2neg)) {

            fA_transition_points.push(v1);
            fA_transition_points.push(v2);
            // this comment is wrong - "record the triple intersection (fA, fB, fAn2) in fA_transition_points;"

          //}
          }
          
        //}
        }
      
      }

      
      //Compute the intersection of the line-subsets described in f1_transition_points, f2_transition_points; insert them into face_overlaps;
      bool f1_on = false;
      bool f2_on = false;
      faux_optional<silly_rational_loc> segment_forming;
      assert(!(f1_transition_points.size() % 2));
      assert(!(f2_transition_points.size() % 2));
      while ((!f1_transition_points.empty()) && (!f2_transition_points.empty())) {
        bool really_transitioning;
        silly_rational_loc l;
        if (f1_transition_points.top() < f2_transition_points.top()) {
          l = f2_transition_points.top();
          f2_on = !f2_on;
          really_transitioning = f1_on;
          f2_transition_points.pop();
        }
        else {
          l = f1_transition_points.top();
          f1_on = !f1_on;
          really_transitioning = f2_on;
          f1_transition_points.pop();
        }
        if (really_transitioning) {
          if (segment_forming) {
            assert((f1_on || f2_on) && (!(f1_on && f2_on)));
            face_overlaps.insert(std::make_pair(l, *segment_forming));
            segment_forming = boost::none;
          }
          else {
            assert(f1_on && f2_on);
            segment_forming = l;
          }
        }
      }
      assert(!(f1_on && f2_on));
      assert(!segment_forming);
    
    //}
    }}
    
    std::vector<std::vector<silly_rational_loc>> result;
    while (!face_overlaps.empty()) {
      std::vector<silly_rational_loc> loop;
      loop.push_back(face_overlaps.arbitrary_triple());
      do {
        silly_rational_loc const* nextp = face_overlaps.arbitrary_other_end(loop.back());
        assert(nextp);
        silly_rational_loc next = *nextp;
        face_overlaps.erase(std::make_pair(loop.back(), next));
        loop.push_back(next);
      }
      while (loop.back() != loop.front());
    }
    return result;
  }

  void display_self_overlaps(region const& r, gl_triangles& triangles, gl_data_format::color c, float width) {
    auto foo = find_self_overlaps(r);
    for (auto const& bar : foo) {
      gl_polygon polygon;
      for (auto const& baz : bar) {
        const vector3<distance> loc = baz.nums / baz.shared_denom;
        polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          get_primitive_float(loc.x/distance_units),
          get_primitive_float(loc.y/distance_units),
          get_primitive_float(loc.z/distance_units),
          c));
      }
      push_wireframe_polygon(triangles, width, polygon);
    }
  }
  
  void display_face(face const& f, gl_triangles& triangles, gl_data_format::color c, float width) {
    gl_polygon polygon;
    const face present_face = f.updated_to_time(present_time_);
    for(size_t i = 0; i < f.neighboring_faces_.size(); ++i) {
      const size_t j = (i+1)%f.neighboring_faces_.size();
      const face present_neighbor_1 = faces_[f.neighboring_faces_[i]].updated_to_time(present_time_);
      const face present_neighbor_2 = faces_[f.neighboring_faces_[j]].updated_to_time(present_time_);
      const vector3<distance> loc = approx_loc_of_triple_intersection_of_up_to_date_faces(present_face, present_neighbor_1, present_neighbor_2);
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
        get_primitive_float(loc.x/distance_units),
        get_primitive_float(loc.y/distance_units),
        get_primitive_float(loc.z/distance_units),
        c));
      //std::cerr << vertices[i] << '\n';
    }
    push_wireframe_polygon(triangles, width, polygon);
  }
  
  gl_triangles display(vector3<distance> where, shared_ptr<event> current_event) {
    //assert(when < next_event_time)
    gl_triangles triangles;
    for (face const& f : faces_) {
      display_face(f, triangles, gl_data_format::color(0xffff0080), 1e9);
    }
    for (region const& r : regions_) {
      display_self_overlaps(r, triangles, gl_data_format::color(0xffffff80), 3e9);
    }
    if (current_event) {
      if (const shared_ptr<collision> c = dynamic_pointer_cast<collision>(current_event)) {
        display_face(faces_[c->f1_], triangles, gl_data_format::color(0xff000080), 3e9);
        display_face(faces_[c->f2_], triangles, gl_data_format::color(0xff008080), 3e9);
        display_face(faces_[c->f3_], triangles, gl_data_format::color(0x8000ff80), 3e9);
        display_face(faces_[c->f4_], triangles, gl_data_format::color(0x0000ff80), 3e9);
        const face f1 = faces_[c->f1_].updated_to_time(present_time_);
        const face f2 = faces_[c->f2_].updated_to_time(present_time_);
        const face f3 = faces_[c->f3_].updated_to_time(present_time_);
        const face f4 = faces_[c->f4_].updated_to_time(present_time_);
        const vector3<mpz> v1 = approx_loc_of_triple_intersection_of_up_to_date_faces(f2,f3,f4)/distance_units;
        const vector3<mpz> v2 = approx_loc_of_triple_intersection_of_up_to_date_faces(f1,f3,f4)/distance_units;
        const vector3<mpz> v3 = approx_loc_of_triple_intersection_of_up_to_date_faces(f1,f2,f4)/distance_units;
        const vector3<mpz> v4 = approx_loc_of_triple_intersection_of_up_to_date_faces(f1,f2,f3)/distance_units;
        //std::cerr << v4 << "\n";
        /// haaaaaack
        glBegin(GL_POINTS);
          glVertex3f((float)v1.x, (float)v1.y, (float)v1.z);
          glVertex3f((float)v2.x, (float)v2.y, (float)v2.z);
          glVertex3f((float)v3.x, (float)v3.y, (float)v3.z);
          glVertex3f((float)v4.x, (float)v4.y, (float)v4.z);
        glEnd();
      }
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
  shared_ptr<event> advance_to_next_real_event() {
    while(!next_events_.empty()) {
      const shared_ptr<event> e = next_events_.top();
      if (do_next_event()) return e;
    }
    return shared_ptr<event>();
  }
private:

  // Returns true if an event happened, false if the event was invalid for one reason or another.
  // I can't think of a real reason for it to do that, it's just that way for debugging purposes.
  bool do_next_event() {
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
          face& vf1 = faces_[vfc->vertex_face_1()];
          face& vf2 = faces_[vfc->vertex_face_2()];
          face& vf3 = faces_[vfc->vertex_face_3()];
          face& sf  = faces_[vfc->struck_face()];
          if (vertex_is_in_bounded_face__hack(vfc->when_event_occurs_, vf1, vf2, vf3, sf)) {
            vector3<mpz> normal = sf.updated_to_time(vfc->when_event_occurs_).ABC;
            vector3<velocity1d> problem_velocity = -normal * 10000000 * distance_units / time_units / normal.magnitude_using<mpz>(); //normal*(approx_velocity_of_triple_intersection_of_up_to_date_faces(vf1.updated_to_time(vfc->when_event_occurs_), vf2.updated_to_time(vfc->when_event_occurs_), vf3.updated_to_time(vfc->when_event_occurs_)).dot<mpz>(normal) - sf.D_velocity) / normal.dot<mpz>(normal);
            assert(problem_velocity != 0);
            // Hack - this standardizes "normals point outwards from rock"
            if (problem_velocity.dot<mpz>(normal) < 0) {
              //std::cerr << c->r1_ << "," << c->r2_ << "," << c->r2_ << "," << c->r2_ << "\n";
              //std::cerr << "vfc.\n";
              // Haaaaaaack
              region_idx_type ri1 = 0;
              region_idx_type ri2 = 0;
              for (auto ri : vf1.neighboring_regions_) {
                if (regions_[ri].contents == ROCK) ri1 = ri;
              }
              for (auto ri : sf.neighboring_regions_) {
                if (regions_[ri].contents == ROCK) ri2 = ri;
              }
              region& r1 = regions_[ri1];
              region& r2 = regions_[ri2];
              assert(r1.contents = ROCK); assert(r2.contents = ROCK);
              for (face_idx_type fi : r1.faces_) {
                faces_[fi] = faces_[fi].updated_to_time(vfc->when_event_occurs_);
              }
              for (face_idx_type fi : r2.faces_) {
                faces_[fi] = faces_[fi].updated_to_time(vfc->when_event_occurs_);
              }
              for (face_idx_type fi : r1.faces_) {
                faces_[fi].D_velocity -= problem_velocity.dot<mpz>(faces_[fi].ABC) * 2 / 3;
              }
              for (face_idx_type fi : r2.faces_) {
                faces_[fi].D_velocity += problem_velocity.dot<mpz>(faces_[fi].ABC) * 2 / 3;
              }
              // TODO : have the recomputation be automated somehow
              for (face_idx_type fi : r1.faces_) {
                insert_events_involving(fi);
              }
              for (face_idx_type fi : r2.faces_) {
                insert_events_involving(fi);
              }
              return true;
            }
          }
        }
        if (const shared_ptr<edge_edge_collision> eec = dynamic_pointer_cast<edge_edge_collision>(c)) {
          face& e11 = faces_[eec->edge_1_face_1()];
          face& e12 = faces_[eec->edge_1_face_2()];
          face& e21 = faces_[eec->edge_2_face_1()];
          face& e22 = faces_[eec->edge_2_face_2()];
          size_t n1 = 0;
          size_t n2 = 0;
          for (size_t i = 0; i < e11.neighboring_faces_.size(); ++i) {
            if (e11.neighboring_faces_[i] == eec->edge_1_face_2()) {
              n1 = i;
              break;
            }
          }
          for (size_t i = 0; i < e21.neighboring_faces_.size(); ++i) {
            if (e21.neighboring_faces_[i] == eec->edge_2_face_2()) {
              n2 = i;
              break;
            }
          }
          if (bounded_edges_cross__hack(eec->when_event_occurs_, e11, e12, n1, e21, e22, n2)) {
            face p11 = e11.updated_to_time(eec->when_event_occurs_);
            face p12 = e12.updated_to_time(eec->when_event_occurs_);
            face p21 = e21.updated_to_time(eec->when_event_occurs_);
            face p22 = e22.updated_to_time(eec->when_event_occurs_);
#if 0
            vector3<mpz> normal = sf.updated_to_time(vfc->when_event_occurs_).ABC;
            vector3<velocity1d> problem_velocity = normal*(approx_velocity_of_triple_intersection_of_up_to_date_faces(vf1.updated_to_time(vfc->when_event_occurs_), vf2.updated_to_time(vfc->when_event_occurs_), vf3.updated_to_time(vfc->when_event_occurs_)).dot<mpz>(normal) - sf.D_velocity) / normal.dot<mpz>(normal);
            if () {
              //std::cerr << "eec.\n";
              // Haaaaaaack
              region_idx_type ri1 = 0;
              region_idx_type ri2 = 0;
              for (auto ri : e11.neighboring_regions_) {
                if (regions_[ri].contents == ROCK) ri1 = ri;
              }
              for (auto ri : e21.neighboring_regions_) {
                if (regions_[ri].contents == ROCK) ri2 = ri;
              }
              region& r1 = regions_[ri1];
              region& r2 = regions_[ri2];
              assert(r1.contents = ROCK); assert(r2.contents = ROCK);
              //vector3<mpz> normal = ;
              // Hack, this is COMPLETELY wrong
              //vector3<velocity1d> vel_thing = normal * (mpz(1000000000)*distance_units/seconds) / normal.magnitude_using<mpz>();
              for (face_idx_type fi : r1.faces_) {
                faces_[fi] = faces_[fi].updated_to_time(eec->when_event_occurs_);
                faces_[fi].D_velocity *= -1; //= faces_[fi].ABC.dot<mpz>(vel_thing);
              }
              for (face_idx_type fi : r2.faces_) {
                faces_[fi] = faces_[fi].updated_to_time(eec->when_event_occurs_);
                faces_[fi].D_velocity *= -1; // = -faces_[fi].ABC.dot<mpz>(vel_thing);
              }
              // TODO : have the recomputation be automated somehow
              for (face_idx_type fi : r1.faces_) {
                insert_events_involving(fi);
              }
              for (face_idx_type fi : r2.faces_) {
                insert_events_involving(fi);
              }
              return true;
            }
#endif
          }
        }
      }
    }
    return false;
  }
  
};



void do_gl(grand_structure_of_lasercake& simulated_world, uint64_t frame, shared_ptr<event> current_event) {
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

  gl_triangles data = simulated_world.display(vector3<distance>(viewcenter*distance_units), current_event);
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
  SDL_Event sdle;
  int done = 0;
  int p_mode = 0;
  bool moving = false;
large_fast_noncrypto_rng rng(time(NULL));

  
  
int frame = 0;
time_type when = 0;
bool do_events = true;
shared_ptr<event> current_event;

  grand_structure_of_lasercake simulated_world;
  
  while ( !done ) {

    /* Check for events */
    while ( SDL_PollEvent (&sdle) ) {
      switch (sdle.type) {
        case SDL_MOUSEMOTION:
          break;
          
        case SDL_MOUSEBUTTONDOWN:
          break;
          
        case SDL_KEYDOWN:
          if(sdle.key.keysym.sym == SDLK_p) ++p_mode;
          if(sdle.key.keysym.sym == SDLK_m) moving = !moving;
          if(sdle.key.keysym.sym == SDLK_e) when = simulated_world.time_of_next_event();
          if(sdle.key.keysym.sym == SDLK_r) {
            current_event = simulated_world.advance_to_next_real_event();
            if (current_event) when = current_event->when_event_occurs_;
          }
          if(sdle.key.keysym.sym == SDLK_q) do_events = false;
          //if(sdle.key.keysym.sym == SDLK_r) ++view_dist;
          //if(sdle.key.keysym.sym == SDLK_f) --view_dist;
          if(sdle.key.keysym.sym != SDLK_ESCAPE)break;
          
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

    simulated_world.update_to_time(when, do_events);
    do_gl(simulated_world, frame, current_event);
    glFinish();	
    SDL_GL_SwapBuffers();
   
    __attribute__((unused)) int before_processing = SDL_GetTicks();
    
    //doing stuff code here
    
    ++frame;
    if (moving) when += int64_t(10LL)*milli*seconds;
    //std::cerr << when << "\n";
    
    
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
