
// TODO:
// This entire file is just bits of other files hacked together.
// It only exists for historical reasons. It should be destroyed.



#include <array>
using std::array;
#include <utility>
using std::declval;


#ifndef ATTRIBUTE_NORETURN
    // from http://www.boost.org/doc/libs/1_48_0/boost/exception/detail/attribute_noreturn.hpp
#if defined(_MSC_VER)
#define ATTRIBUTE_NORETURN __declspec(noreturn)
#elif defined(__GNUC__)
#define ATTRIBUTE_NORETURN __attribute__((noreturn))
#else
#define ATTRIBUTE_NORETURN
#endif
#endif

#ifndef LASERCAKE_CONFIG_HPP__
// It's not polite for library functions to assert() because the library's users
// misused a correct library; use these for that case.
inline ATTRIBUTE_NORETURN void caller_error(const char* error) {
  // If exceptions prove worse for debugging than asserts/segfaults,
  // feel free to comment this out and use asserts/segfaults/breakpoints.
  //boost::throw_exception(std::logic_error(error));
  throw std::logic_error(error);
}
// You must provide an explanatory string so that the user of the library
// will know what *they* did wrong, and not have to interpret an assert() expression
// to find out.
inline void caller_error_if(bool cond, const char* error) {
  if(cond) {
    caller_error(error);
  }
}
inline void caller_correct_if(bool cond, const char* error) {
  if(!cond) {
    caller_error(error);
  }
}
#endif

// This has to be a non-function because ret should only
// be evaluated if cond is true.
#define constexpr_require_and_return(cond, str, ret) ((cond) ? (ret) : throw std::logic_error((str)))
// This has to be a non-function so that the compiler can deduce its
// return type as anything for use in ?: results. It has to be
// non-parenthesized because Clang, at least, will cast it to void
// (rather than its magic whatever-the-other-thing-in-the-?:-is type)
// if parenthesized.
#define constexpr_caller_error(str)  throw std::logic_error((str))

// If you need to pass around a dimension for use as an index to
// a vector3, you can optionally use these names.
enum {
  X = 0, Y = 1, Z = 2, vector3_num_dimensions = 3
};
typedef int which_dimension_type;

template<typename ScalarType> class vector3 {
public:
  ScalarType x, y, z;
  constexpr vector3():x(0),y(0),z(0){}
  // implicit conversion from literal 0, so you can write 0 for any zero vector.
  constexpr vector3(decltype(nullptr)):x(0),y(0),z(0){}
  constexpr vector3(ScalarType x, ScalarType y, ScalarType z):x(x),y(y),z(z){}
  template<typename OtherType> explicit constexpr vector3(vector3<OtherType> const& other):
    x(other.x),y(other.y),z(other.z){}

  // implicit conversions to/from array:
  BOOST_FORCEINLINE operator array<ScalarType, 3>()const {
    array<ScalarType, 3> result = {{ x, y, z }};
    return result;
  }
  BOOST_FORCEINLINE vector3(array<ScalarType, 3> const& arr)
    : x(arr[0]), y(arr[1]), z(arr[2]) {}
  
  BOOST_FORCEINLINE ScalarType& operator[](which_dimension_type index) {
    if(index == X){return x;} if(index == Y){return y;} if(index == Z){return z;}
      constexpr_caller_error("Trying to index a vector3 with an out-of-bounds index!");
  }
  BOOST_FORCEINLINE constexpr ScalarType operator[](which_dimension_type index)const {
    return (index == X) ? (x) : (index == Y) ? (y) : (index == Z) ? (z) :
      constexpr_caller_error("Trying to index a vector3 with an out-of-bounds index!");
  }
  BOOST_FORCEINLINE constexpr ScalarType operator()(which_dimension_type index)const { return (*this)[index]; }
  BOOST_FORCEINLINE void set(which_dimension_type index, ScalarType value) { (*this)[index] = value; }

  template<typename OtherType> constexpr auto operator+(vector3<OtherType> const& other)const
  -> vector3<decltype(x + other.x)> {
    return vector3<decltype(x + other.x)>(x + other.x, y + other.y, z + other.z);
  }
  template<typename OtherType> vector3& operator+=(vector3<OtherType> const& other) {
    x += other.x; y += other.y; z += other.z; return *this;
  }
  template<typename OtherType> constexpr auto operator-(vector3<OtherType> const& other)const
  -> vector3<decltype(x - other.x)> {
    return vector3<decltype(x - other.x)>(x - other.x, y - other.y, z - other.z);
  }
  template<typename OtherType> vector3& operator-=(vector3<OtherType> const& other) {
    x -= other.x; y -= other.y; z -= other.z; return *this;
  }
  template<typename OtherType> constexpr auto operator*(OtherType const& other)const
  -> vector3<decltype(x * other)> {
    return vector3<decltype(x * other)>(x * other, y * other, z * other);
  }
  vector3& operator*=(ScalarType other) {
    x *= other; y *= other; z *= other; return *this;
  }
  template<typename OtherType, typename RoundingStrategy>
  friend inline constexpr auto divide(vector3 const& v, OtherType const& other, RoundingStrategy strat)
  -> vector3<decltype(declval<ScalarType>() / other)>{
    return vector3<decltype(v.x / other)>(
      divide(v.x, other, strat),
      divide(v.y, other, strat),
      divide(v.z, other, strat));
  }
  // Default to rounding towards zero. (TODO: is it wise to have any default here?
  // It's not like we use division much. But we don't want to use shifting
  // without considering that shifting rounds down towards negative infinity, too.)
  //typedef rounding_strategy<round_down, negative_mirrors_positive> default_rounding_strategy;

  // In C++11 integer division rounds towards zero,
  // which is often what we want for vectors; IEEE754 floating point division,
  // by default, rounds to nearest and to even for ties.
  template<typename OtherType> constexpr auto operator/(OtherType const& other)const
  -> vector3<decltype(x / other)> {
    return vector3<decltype(x / other)>(x/other, y/other, z/other);
  }
  vector3& operator/=(ScalarType other) {
    x /= other; y /= other; z /= other;
    return *this;
  }
  // Multiplying two vectors is usually a type-error mistake, so
  // you have to say you're doing it in words:
  template<typename OtherType> constexpr auto multiply_piecewise_by(vector3<OtherType> const& other)const
  -> vector3<decltype(x * other.x)> {
    return vector3<decltype(x * other.x)>(x * other.x, y * other.y, z * other.z);
  }
  template<typename OtherType, typename RoundingStrategy>
  constexpr auto divide_piecewise_by(vector3<OtherType> const& other, RoundingStrategy strat)const
  -> vector3<decltype(x / other.x)> {
    return vector3<decltype(x / other.x)>(
      divide(x, other.x, strat),
      divide(y, other.y, strat),
      divide(z, other.z, strat));
  }
  // Careful, shift operators on builtin types (ScalarType?) are only
  // defined for shift >= 0 && shift < bits_in_type
  constexpr vector3 operator<<(int shift)const {
    return vector3(x << shift, y << shift, z << shift);
  }
  vector3& operator<<=(int shift) {
    x <<= shift; y <<= shift; z <<= shift; return *this;
  }
  constexpr vector3 operator>>(int shift)const {
    return vector3(x >> shift, y >> shift, z >> shift);
  }
  vector3& operator>>=(int shift) {
    x >>= shift; y >>= shift; z >>= shift; return *this;
  }
  constexpr vector3 operator^(ScalarType other)const {
    return vector3(x ^ other, y ^ other, z ^ other);
  }
  constexpr vector3 operator|(ScalarType other)const {
    return vector3(x | other, y | other, z | other);
  }
  constexpr vector3 operator&(ScalarType other)const {
    return vector3(x & other, y & other, z & other);
  }
  constexpr vector3 operator~()const {
    return vector3(~x, ~y, ~z);
  }

  constexpr vector3 operator+()const { return *this; } // unary plus
  constexpr vector3 operator-()const { // unary minus
    return vector3(-x, -y, -z);
  }

  constexpr bool operator==(vector3 const& other)const {return x == other.x && y == other.y && z == other.z; }
  constexpr bool operator!=(vector3 const& other)const {return x != other.x || y != other.y || z != other.z; }
  
  // Do not try to use this if either vector has an unsigned ScalarType.
  // It might work in some situations, but why would you ever do that anyway?
  //
  // You are required to specify an output representation type, because of the
  // risk of overflow. Make sure to choose one that can fit the squares of the
  // numbers you're dealing with.
#if 0
  template<typename OutputRepr, typename OtherType>
  constexpr auto dot(vector3<OtherType> const& other)const
  -> decltype(numeric_representation_cast<OutputRepr>(x) * numeric_representation_cast<OutputRepr>(other.x)) {
    return
      numeric_representation_cast<OutputRepr>(x) * numeric_representation_cast<OutputRepr>(other.x) +
      numeric_representation_cast<OutputRepr>(y) * numeric_representation_cast<OutputRepr>(other.y) +
      numeric_representation_cast<OutputRepr>(z) * numeric_representation_cast<OutputRepr>(other.z);
  }

  typedef lint64_t int64_type_to_use_with_dot;
  ScalarType magnitude_within_32_bits()const {
    return ScalarType(i64sqrt(dot<int64_type_to_use_with_dot>(*this)));
  }
  
  // Choose these the way you'd choose dot's output type (see the comment above)
  // we had trouble making these templates, so now they just always use int64_t
  constexpr bool magnitude_within_32_bits_is_less_than(ScalarType amount)const {
    return dot<int64_type_to_use_with_dot>(*this) <
          numeric_representation_cast<int64_type_to_use_with_dot>(amount)
        * numeric_representation_cast<int64_type_to_use_with_dot>(amount);
  }
  constexpr bool magnitude_within_32_bits_is_greater_than(ScalarType amount)const {
    return dot<int64_type_to_use_with_dot>(*this) >
          numeric_representation_cast<int64_type_to_use_with_dot>(amount)
        * numeric_representation_cast<int64_type_to_use_with_dot>(amount);
  }
#endif
  constexpr bool operator<(vector3 const& other)const {
    return (x < other.x) || ((x == other.x) && ((y < other.y) || ((y == other.y) && (z < other.z))));
  }

  friend inline std::ostream& operator<<(std::ostream& os, vector3 const& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
  }
#if 0
  friend inline size_t hash_value(vector3 const& v) {
      size_t seed = 0;
      boost::hash_combine(seed, v.x);
      boost::hash_combine(seed, v.y);
      boost::hash_combine(seed, v.z);
      return seed;
  }
#endif
};

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

#ifndef LASERCAKE_GL_DATA_FORMAT_HPP__
#define LASERCAKE_GL_DATA_FORMAT_HPP__

// This header does not use GL types (GLfloat, GLubyte) in order to
// be agnostic between GLEW and Qt's opinions of OpenGL headers.
// gl_data_preparation.cpp static_asserts that they're the types we
// expect.
// (This might not be necessary currently, but doesn't hurt anything
// (that I know of).)

#include <string>
#include <vector>
//#include "cxx11/array.hpp"
#include <array>
using std::array;
#include <ostream>
#include <stdlib.h> //malloc
#include <string.h> //memcpy

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/type_ptr.hpp>

//#include "world_constants.hpp"


namespace gl_data_format {

typedef float header_GLfloat;
typedef unsigned char header_GLubyte;

// These are to be passed as part of arrays to OpenGL.
// Thus their data format/layout makes a difference.
// TODO maybe make vector3<GLfloat> and vertex the same thing?
// Is vector3<GLfloat>'s layout required to be the same as that of
// this struct? I believe so.  typedef vector3<GLfloat> vertex; ?
// TODO use glm types?
struct vertex {
  vertex() {}

  vertex(header_GLfloat x, header_GLfloat y, header_GLfloat z) : x(x), y(y), z(z) {}
  /* implicit conversion */ vertex(vector3<header_GLfloat> const& v) : x(v.x), y(v.y), z(v.z) {}
  /* implicit conversion */ vertex(glm::vec3 const& v) : x(v.x), y(v.y), z(v.z) {}
  /* implicit conversion */ vertex(glm::dvec3 const& v) : x(v.x), y(v.y), z(v.z) {}

  header_GLfloat x, y, z;
};
static_assert(sizeof(vertex) == 3*sizeof(header_GLfloat), "OpenGL needs this data layout.");
inline std::ostream& operator<<(std::ostream& os, vertex const& v) {
  // TODO does the float precision we output here matter?
  return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

struct color {
  color() {}

  // Use hex RGBA values as familiar from e.g. CSS.
  // e.g. 0x00ff0077 is pure green at 50% opacity.
  explicit color(uint32_t rgba)
    : r(header_GLubyte(rgba >> 24)), g(header_GLubyte(rgba >> 16)),
      b(header_GLubyte(rgba >>  8)), a(header_GLubyte(rgba)) {}

  color(header_GLubyte r, header_GLubyte g, header_GLubyte b, header_GLubyte a)
    : r(r), g(g), b(b), a(a) {}

  // random Web forums thought a factor of 255.0 was correct, but is it exactly the right conversion?
  color(header_GLfloat r, header_GLfloat g, header_GLfloat b, header_GLfloat a)
   : r(header_GLubyte(r*255)), g(header_GLubyte(g*255)), b(header_GLubyte(b*255)), a(header_GLubyte(a*255)) {}

  header_GLubyte r, g, b, a;
};
static_assert(sizeof(color) == 4*sizeof(header_GLubyte), "OpenGL needs this data layout.");
inline std::ostream& operator<<(std::ostream& os, color const& c) {
  // TODO does the float precision we output here matter?
  const char oldfill = os.fill('0');
  const std::ios::fmtflags oldflags = os.setf(std::ios::hex, std::ios::basefield);
  os << "rgba(";
  os.width(2); os << int(c.r);
  os.width(2); os << int(c.g);
  os.width(2); os << int(c.b);
  os.width(2); os << int(c.a);
  os << ')';
  os.flags(oldflags);
  os.fill(oldfill);
  return os;
}

struct vertex_with_color {
  vertex_with_color() {}

  vertex_with_color(vertex v, color c) : c(c), v(v) {}
  vertex_with_color(header_GLfloat x, header_GLfloat y, header_GLfloat z, color c)
    : c(c), v(x,y,z) {}
  color c;
  vertex v;
};
static_assert(sizeof(vertex_with_color) == 16, "OpenGL needs this data layout.");
inline std::ostream& operator<<(std::ostream& os, vertex_with_color const& vc) {
  // TODO does the float precision we output here matter?
  return os << vc.c << '~' << vc.v;
}

struct gl_call_data {
  typedef uint32_t size_type;
  size_type count;
  size_type alloced;
  vertex_with_color* vertices;
  //vertex_with_color* vertices_end;
  static const size_type default_size = 5;
  static const size_type expand_multiplier = 4;
  gl_call_data()
    : count(0),
      alloced(default_size),
      vertices((vertex_with_color*)malloc(default_size*sizeof(vertex_with_color))) {
    assert(vertices);
  }
  gl_call_data(gl_call_data&& other)
    : count(other.count),
      alloced(other.alloced),
      vertices(other.vertices) {
    other.vertices = nullptr;
  }
  gl_call_data(gl_call_data const& other)
    : count(other.count),
      alloced(other.alloced),
      vertices((vertex_with_color*)malloc(other.alloced*sizeof(vertex_with_color))) {
    assert(vertices);
    memcpy(vertices, other.vertices, count*sizeof(vertex_with_color));
  }
  gl_call_data& operator=(gl_call_data const& other) {
    free(vertices);
    count = other.count;
    alloced = other.alloced;
    vertices = ((vertex_with_color*)malloc(other.alloced*sizeof(vertex_with_color)));
    assert(vertices);
    memcpy(vertices, other.vertices, count*sizeof(vertex_with_color));
    return *this;
  }
  gl_call_data& operator=(gl_call_data&& other) {
    free(vertices);
    count = other.count;
    alloced = other.alloced;
    vertices = other.vertices;
    other.vertices = nullptr;
    return *this;
  }
  ~gl_call_data() { free(vertices); }
  void push_vertex(vertex_with_color const& v) {
    if(count == alloced) do_realloc(alloced * expand_multiplier);
    vertices[count] = v;
    ++count;
  }
  void push_vertex(vertex const& v, color const& c) {
    if(count == alloced) do_realloc(alloced * expand_multiplier);
    vertices[count].c = c;
    vertices[count].v = v;
    ++count;
  }
  void reserve_new_slots(size_type num_slots) {
    const size_type new_count = count + num_slots;
    if(new_count > alloced) do_realloc(new_count * expand_multiplier);
    count = new_count;
  }
  size_type size() const { return count; }
  void do_realloc(size_type new_size) {
    assert(new_size > alloced);
    vertex_with_color* new_vertices = (vertex_with_color*)malloc(new_size * sizeof(vertex_with_color));
    assert(new_vertices);
    memcpy(new_vertices, vertices, count*sizeof(vertex_with_color));
    free(vertices);
    vertices = new_vertices;
    alloced = new_size;
  }
};

struct gl_collection {
  // Points and lines are for debugging, because they don't change size at distance,
  // and TODO can be represented reasonably by triangles.
  // Triangles are legit.
  // Quads are TODO OpenGL prefers you to use triangles (esp. OpenGL ES) and
  // they'll be converted at some point.
  gl_call_data points;
  gl_call_data lines;
  gl_call_data triangles;
  gl_call_data quads;
};

//The gl_collection:s with higher indices here are intended to be
//further away and rendered first (therefore covered up most
//by everything else that's closer).
typedef std::vector<gl_collection> gl_collectionplex;

struct heads_up_display_text {
  // text may contain newlines, and will also
  // be soft-wrapped to fit on the screen.
  std::string text;
  color c;
  std::string font_name;
  int point_size;
  int horizontal_margin_in_pixels;
  int vertical_margin_in_pixels;
};

struct gl_all_data {
  gl_collectionplex stuff_to_draw_as_gl_collections_by_distance;
  color tint_everything_with_this_color;
  heads_up_display_text hud_text;
  vector3<header_GLfloat> facing;
  vector3<header_GLfloat> facing_up;
};



const int32_t fovy_degrees = 80;
//const non_normalized_rational<int32_t> pretend_aspect_ratio_value_when_culling(2,1);
//const distance near_clipping_plane = tile_width / 10;
//const distance far_clipping_plane = tile_width * 300;
// TODO: we can, with some more work,
// use non-floating-point matrices for several things.
// Which would allow main-simulation code to filter based on
// view frustums, for example.
inline glm::mat4 make_projection_matrix(float aspect_ratio) {
  return glm::perspective(
    float(fovy_degrees),
    float(aspect_ratio),
    float(1), //get_primitive_float(near_clipping_plane/fine_distance_units),
    float(3000) //get_primitive_float(far_clipping_plane/fine_distance_units)
  );
}
const vector3<float> view_from(0);
inline glm::mat4 make_view_matrix(vector3<float> view_towards, vector3<float> up) {
  return glm::lookAt(
    glm::vec3(view_from.x, view_from.y, view_from.z),
    glm::vec3(view_towards.x, view_towards.y, view_towards.z),
    glm::vec3(up.x, up.y, up.z)
  );
}

// TODO glm has glm::detail::tmat4x4<> etc, needed to templatize
// this, if I want to templatize it.
struct frustum {
  enum direction {
    // Windows headers define NEAR and FAR
    // and I don't want to tempt any dragons by #undef'ing
    // them, so just add an underscore to those names here.
    LEFT, RIGHT, BOTTOM, TOP, NEAR_, FAR_
  };
  array<glm::vec4, 6> half_spaces;
};
inline frustum make_frustum_from_matrix(glm::mat4 m) {
  frustum result;
  result.half_spaces[frustum::LEFT]   = glm::normalize(glm::row(m, 3) + glm::row(m, 0));
  result.half_spaces[frustum::RIGHT]  = glm::normalize(glm::row(m, 3) - glm::row(m, 0));
  result.half_spaces[frustum::BOTTOM] = glm::normalize(glm::row(m, 3) + glm::row(m, 1));
  result.half_spaces[frustum::TOP]    = glm::normalize(glm::row(m, 3) - glm::row(m, 1));
  result.half_spaces[frustum::NEAR_]   = glm::normalize(glm::row(m, 3) + glm::row(m, 2));
  result.half_spaces[frustum::FAR_]    = glm::normalize(glm::row(m, 3) - glm::row(m, 2));
  return result;
}

} // end namespace gl_data_format

#endif
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

#if 0
#ifdef LASERCAKE_HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>
#include <boost/chrono/thread_clock.hpp>
#endif


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



#if 0
namespace chrono = boost::chrono;

typedef int64_t microseconds_t;

microseconds_t get_this_thread_microseconds() {
#if defined(BOOST_CHRONO_HAS_THREAD_CLOCK)
  return chrono::duration_cast<chrono::microseconds>(chrono::thread_clock::now().time_since_epoch()).count();
#else
  return 0;
#endif
}
#endif

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
#if 0
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
#endif
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

