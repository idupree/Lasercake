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




#include "../utils.hpp"
//#include "../data_structures/geometry.hpp"


namespace YO {
using boost::shared_ptr;


// we want to sort threevertices together: triangles: not individual vertices
struct gl_triangle {
  array<gl_data_format::vertex_with_color, 3> vertices;
};
typedef std::vector<gl_triangle> gl_triangles;
glm::vec3 v_to_gv(gl_data_format::vertex v) {
  return glm::vec3(v.x, v.y, v.z);
}
float gl_triangle_distance_order(glm::vec3 from, gl_triangle const& triangle) {
  return glm::distance(from, v_to_gv(triangle.vertices[0].v))
       + glm::distance(from, v_to_gv(triangle.vertices[1].v))
       + glm::distance(from, v_to_gv(triangle.vertices[2].v));
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
void push_wireframe_triangle(
      gl_triangles& coll, GLfloat width,
      gl_triangle triangle) {
  const auto vs = triangle.vertices;
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
typedef physical_quantity<lint64_t, time_units_t> time_type;

typedef typename units_prod<nano_t, meters_t>::type distance_units_t;
constexpr distance_units_t distance_units = distance_units_t();
typedef physical_quantity<lint64_t, distance_units_t> distance;
typedef physical_quantity<lint64_t, typename units_prod<distance_units_t, dim::second<(-1)> >::type> velocity1d;
typedef physical_quantity<lint64_t, typename units_prod<distance_units_t, dim::second<(-2)> >::type> acceleration1d;

// Standard (Earth-equivalent) gravity: precisely 9.80665 m/s2
constexpr acceleration1d gravity_acceleration_magnitude = 9806650LL * (micro*meters) / (seconds*seconds) * identity(distance_units / (micro*meters));
constexpr vector3<acceleration1d> gravity_acceleration(0, 0, -gravity_acceleration_magnitude);

typedef uint32_t region_idx_type;
typedef uint32_t face_idx_type;
typedef uint32_t vertex_idx_type;
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
  time_type base_time_;
  vector3<lint64_t> ABC;
  distance D;
  velocity1d D_velocity;
  acceleration1d D_acceleration;
  std::vector<face_idx_type> neighboring_faces_;
  std::vector<region_idx_type> neighboring_regions_;
  face updated_to_time(time_type t)const {
    face result(*this);
    const time_type relative_time = t - base_time_;
    const auto relative_timem = relative_time/identity(units_factor<1, 1000000000>());
    //TODO deal with times somehow more correctly
    result.D          += D_velocity    *relative_timem                / identity(units_factor<1,    1000>())
                       + D_acceleration*relative_timem*relative_timem / identity(units_factor<1, 1000000>())/2;
    result.D_velocity += D_acceleration*relative_timem                / identity(units_factor<1,    1000>());
    return result;
  }
};

lint64_t scalar_triple_product(vector3<lint64_t> v1, vector3<lint64_t> v2, vector3<lint64_t> v3) {
  return v1(X)*v2(Y)*v3(Z) - v1(X)*v3(Y)*v2(Z) + v2(X)*v3(Y)*v1(Z) - v2(X)*v1(Y)*v3(Z) + v3(X)*v1(Y)*v2(Z) - v3(X)*v2(Y)*v1(Z);
}
vector3<distance> approx_loc_of_triple_intersection_of_up_to_date_faces(face const& f1, face const& f2, face const& f3) {
  lint64_t p = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  assert(p != 0);
  return vector3<distance>(
    f1.D * (f2.ABC(1)*f3.ABC(2) - f3.ABC(1)*f2.ABC(2)) + f2.D * (f3.ABC(1)*f1.ABC(2) - f1.ABC(1)*f3.ABC(2)) + f3.D * (f1.ABC(1)*f2.ABC(2) - f2.ABC(1)*f1.ABC(2)),
    f1.D * (f2.ABC(2)*f3.ABC(0) - f3.ABC(2)*f2.ABC(0)) + f2.D * (f3.ABC(2)*f1.ABC(0) - f1.ABC(2)*f3.ABC(0)) + f3.D * (f1.ABC(2)*f2.ABC(0) - f2.ABC(2)*f1.ABC(0)),
    f1.D * (f2.ABC(0)*f3.ABC(1) - f3.ABC(0)*f2.ABC(1)) + f2.D * (f3.ABC(0)*f1.ABC(1) - f1.ABC(0)*f3.ABC(1)) + f3.D * (f1.ABC(0)*f2.ABC(1) - f2.ABC(0)*f1.ABC(1))
                 ) / p;
}

faux_optional<time_type> how_long_from_now_will_planes_of_up_to_date_faces_be_coincident_at_a_point(face const& f1, face const& f2, face const& f3, face const& f4) {
  // When the 4x4 determinant is 0.
  // That's
  // + D1 * scalar_triple_product(f2.ABC, f3.ABC, f4.ABC)
  // - D2 * scalar_triple_product(f3.ABC, f4.ABC, f1.ABC)
  // + D3 * scalar_triple_product(f4.ABC, f1.ABC, f2.ABC)
  // - D4 * scalar_triple_product(f1.ABC, f2.ABC, f3.ABC)
  // But Dn is Dn + dDn*t + .5ddDn*t^2
  // so
  const lint64_t t1 = scalar_triple_product(f2.ABC, f3.ABC, f4.ABC);
  const lint64_t t2 = scalar_triple_product(f3.ABC, f4.ABC, f1.ABC);
  const lint64_t t3 = scalar_triple_product(f4.ABC, f1.ABC, f2.ABC);
  const lint64_t t4 = scalar_triple_product(f1.ABC, f2.ABC, f3.ABC);
  // in at^2 + bt + c = 0
  acceleration1d a_times_2 = f1.D_acceleration*t1 + f2.D_acceleration*t2 + f3.D_acceleration*t3 + f4.D_acceleration*t4;
  velocity1d     b         = f1.D_velocity    *t1 + f2.D_velocity    *t2 + f3.D_velocity    *t3 + f4.D_velocity    *t4;
  distance       c         = f1.D             *t1 + f2.D             *t2 + f3.D             *t3 + f4.D             *t4;
  
  // (-b +/- sqrt(b^2 - 2(a_times_2)c)) / a_times_2
  const auto discriminant = b*b - a_times_2*2*c;
  if (discriminant < 0) {
    return boost::none;
  }
  else {
    const rounding_strategy<round_down, negative_is_forbidden> strat;
    if (a_times_2 == 0) {
      if (b == 0) {
        // Eww. A constant. They're either ALWAYS intersecting (c == 0) or NEVER (c != 0).
        if (c == 0) {
          assert(false); // I don't think we can handle this case.
          return time_type(0);
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
      return divide(-c * lint64_t(identity(time_units / seconds)*seconds/time_units), b, strat)*time_units/seconds;
    }
    else {
      if (a_times_2 < 0) {
        a_times_2 = -a_times_2;
        b = -b;
        c = -c;
      }
      // we want the earlier time, but not if it's negative
      const velocity1d sqrt_disc = i64sqrt(discriminant);
      const velocity1d  lesser_numerator = -b - sqrt_disc;
      if ( lesser_numerator >= 0) return divide(lesser_numerator * lint64_t(identity(time_units / seconds)*seconds/time_units), a_times_2, strat)*time_units/seconds;
      const velocity1d greater_numerator = -b + sqrt_disc;
      if (greater_numerator >= 0) return divide(greater_numerator * lint64_t(identity(time_units / seconds)*seconds/time_units), a_times_2, strat)*time_units/seconds;
      return boost::none;
    }
  }
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
};
struct collision : public event {
};
struct vertex_face_collision : public collision {
  //This includes polyhedra that involute (hopefully)
};
struct edge_edge_collision : public collision {
};


// Maybe the structures should be required to follow some kind of right-hand-rule
// order? Then asserts could check that it is still correct, and OpenGL would
// be more easily able to draw directional surfaces.
class grand_structure_of_lasercake {
  std::vector<region> regions_;
  std::vector<face> faces_;

  void hack_insert_rock(vector3<distance> loc) {
    face_idx_type first_face_idx = faces_.size();
    region_idx_type region_idx = regions_.size();

    static const vector3<lint64_t> diffs[4] = {
      vector3<lint64_t>(0, 0, 100),
      vector3<lint64_t>(0, 83, -45),
      vector3<lint64_t>(60, -40, -50),
      vector3<lint64_t>(-62, -42, -43)
    };

    regions_.push_back(region());
    region& air = regions_[0];
    region& r = regions_.back();
    for (int i = 0; i < 4; ++i) {
      faces_.push_back(face());
      face& f = faces_.back();
      f.base_time_ = 0;
      f.ABC = diffs[i];
      f.D = loc.dot<lint64_t>(f.ABC) + f.ABC.magnitude_within_32_bits()*100*(centi*meters)*identity(distance_units/(centi*meters));
      f.D_velocity = 0;
      f.D_acceleration = f.ABC.dot<lint64_t>(gravity_acceleration);
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
  
public:
  grand_structure_of_lasercake() {
    regions_.push_back(region()); // the air
    hack_insert_rock(vector3<lint64_t>(2, 7, 11)*meters*identity(distance_units/meters));
    hack_insert_rock(vector3<lint64_t>(2, 14, 11)*meters*identity(distance_units/meters));
    hack_insert_rock(vector3<lint64_t>(15, 14, 21)*meters*identity(distance_units/meters));
    
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
  std::priority_queue<shared_ptr<event>> next_events_;
  //void player_input_becomes(time_type when, );
  //void insert_event(time_type when, );
  gl_triangles/*triangles*/ display(time_type when, vector3<distance> where) {
    //assert(when < next_event_time)
    gl_triangles triangles;
    for(face const& f : faces_) {
      gl_triangle triangle;
      const face present_face = f.updated_to_time(when);
      for(size_t i = 0; i < f.neighboring_faces_.size(); ++i) {
        size_t j = (i+1)%f.neighboring_faces_.size();
        const face present_neighbor_1 = faces_[f.neighboring_faces_[i]].updated_to_time(when);
        const face present_neighbor_2 = faces_[f.neighboring_faces_[j]].updated_to_time(when);
        vector3<distance> loc = approx_loc_of_triple_intersection_of_up_to_date_faces(present_face, present_neighbor_1, present_neighbor_2);
        triangle.vertices[i] = gl_data_format::vertex_with_color(
          get_primitive_int(loc.x/distance_units),
          get_primitive_int(loc.y/distance_units),
          get_primitive_int(loc.z/distance_units),
          gl_data_format::color(0xffff0080));
        //std::cerr << vertices[i] << '\n';
      }
      push_wireframe_triangle(triangles, 0.5e9, triangle);
    }
    sort_gl_triangles_far_to_near(
      glm::vec3(where.x/distance_units, where.y/distance_units, where.z/distance_units),
      triangles);
    return triangles;
  }
private:
  /*
  void do_next_event() {
    const shared_ptr<event> ev = next_events.top();
    next_events.pop();
    if(const shared_ptr<vertex_face_collision> vfc = dynamic_pointer_cast<vertex_face_collision>(ev)) {
      // The vertex stays existent; N new vertices also appear where it is,
      // where N is the number of edges connected to that vertex.
      // This requires us to split the faces incident to this vertex into
      // triangles (N more faces).
      // The face (which is a triangle) disappears and is replaced with
      // 3+N faces (triangles), plus N faces inside the ring of new vertices.
    }
    if(const shared_ptr<edge_edge_collision> eec = dynamic_pointer_cast<edge_edge_collision>(ev)) {
      // The two edges do not share a vertex, because we make sure to
      // generate only a vertex--face collision in that case.
    }
  }*/
  
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
          /*
          if(event.key.keysym.sym == SDLK_z) insert_objects ++;
          if(event.key.keysym.sym == SDLK_x) insert_objects += 50;
          if(event.key.keysym.sym == SDLK_c) insert_objects += 2500;
          if(event.key.keysym.sym == SDLK_a) do_2d_test_scenario(root);
          if(event.key.keysym.sym == SDLK_s) do_3d_test_scenario(root3d);
          //if(event.key.keysym.sym == SDLK_z) draw_poly = !draw_poly;*/
          //if(event.key.keysym.sym == SDLK_x) draw_normals = !draw_normals;
          //if(event.key.keysym.sym == SDLK_c) draw_endp = !draw_endp;
          //if(event.key.keysym.sym == SDLK_v) draw_coll_stuff = !draw_coll_stuff;
          //if(event.key.keysym.sym == SDLK_b) use_foo1 = !use_foo1;
          //if(event.key.keysym.sym == SDLK_q) ++velocity[X];
          //if(event.key.keysym.sym == SDLK_a) --velocity[X];
          //if(event.key.keysym.sym == SDLK_w) ++velocity[Y];
          //if(event.key.keysym.sym == SDLK_s) --velocity[Y];
          //if(event.key.keysym.sym == SDLK_e) ++velocity[Z];
          //if(event.key.keysym.sym == SDLK_d) --velocity[Z];
          //if(event.key.keysym.sym == SDLK_r) obstacle.translate(vector3<geometry_int_type>(1,0,0));
          //if(event.key.keysym.sym == SDLK_f) obstacle.translate(vector3<geometry_int_type>(-1,0,0));
          //if(event.key.keysym.sym == SDLK_t) obstacle.translate(vector3<geometry_int_type>(0,1,0));
          //if(event.key.keysym.sym == SDLK_g) obstacle.translate(vector3<geometry_int_type>(0,-1,0));
          //if(event.key.keysym.sym == SDLK_y) obstacle.translate(vector3<geometry_int_type>(0,0,1));
          //if(event.key.keysym.sym == SDLK_h) obstacle.translate(vector3<geometry_int_type>(0,0,-1));
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
    if (moving) when += 1000000000LL*pico*seconds;
    
    
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
