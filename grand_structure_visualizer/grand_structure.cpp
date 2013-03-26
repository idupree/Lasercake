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
void sort_gl_triangles_by_distance_to(glm::vec3 view_center, gl_triangles& ts) {
  std::sort(ts.begin(), ts.end(), [view_center](gl_triangle const& t1, gl_triangle const& t2) {
    return gl_triangle_distance_order(view_center, t1) < gl_triangle_distance_order(view_center, t2);
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


typedef luint32_t region_idx_type;
typedef luint32_t vertex_idx_type;
constexpr region_idx_type no_region_idx = std::numeric_limits<region_idx_type>::max();
constexpr vertex_idx_type no_vertex_idx = std::numeric_limits<vertex_idx_type>::max();

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
};


struct region_vertex {
  vertex_idx_type shared_vertex_data_;
  //other data
};
const int tetrahedron_sides = 4;
struct region {
  // Regions are tetrahedrons.
  array<region_vertex, tetrahedron_sides> vertices_;
  // Each face has the array index of the vertex it does not contain.
  // Regions touching the edge of the world have a TODO DEFINE WHICH
  // idx value (0? region_idx_type(-1)? the same region?)
  array<region_idx_type, tetrahedron_sides> neighbors_;
};

struct event {
  virtual ~event(){}
  time_type when_event_occurs_;
};
struct collision : public event {
};
struct vertex_face_collision : public collision {
  //This includes polyhedrons that involute (hopefully)
};
//not with tetrahedron-mesh!
//struct edge_edge_collision : public collision {
//};


// What happens when two neighboring regions overlap?
// They can only become overlapping, if the initial state is consistent, by
// one of them becoming flat for a moment.
// Also, maybe region::vertices_ should be required to follow a right-hand-rule
// order? Then asserts could check that it is still correct, and OpenGL would
// be more easily able to draw directional surfaces.
class grand_structure_of_lasercake {
  std::vector<region> regions_;
  std::vector<vertex> vertices_;
public:
  grand_structure_of_lasercake() {
    vertices_.push_back(vertex{0*time_units,
      vector3<distance>(2*distance_units, 7*distance_units, 11*distance_units), 0, 0});
    vertices_.push_back(vertex{0*time_units,
      vector3<distance>(2*distance_units, 70*distance_units, 11*distance_units), 0, 0});
    vertices_.push_back(vertex{0*time_units,
      vector3<distance>(15*distance_units, 70*distance_units, 51*distance_units), 0, 0});
    vertices_.push_back(vertex{0*time_units,
      vector3<distance>(90*distance_units, 7*distance_units, 51*distance_units), 0, 0});
    vertices_.push_back(vertex{0*time_units,
      vector3<distance>(300*distance_units, 200*distance_units, 100*distance_units), 0, 0});
    regions_.push_back(region{
      {{region_vertex{0}, region_vertex{1}, region_vertex{2}, region_vertex{3}}},
      {{region_idx_type{1}, no_region_idx, no_region_idx, no_region_idx}}});
    regions_.push_back(region{
      {{region_vertex{4}, region_vertex{1}, region_vertex{2}, region_vertex{3}}},
      {{region_idx_type{0}, no_region_idx, no_region_idx, no_region_idx}}});
    
    this->debug_check_consistency();
  }
  void debug_check_consistency()const {
    // region internal consistency
    for(size_t i = 0; i != regions_.size(); ++i) {
      region const& r = regions_[i];
      for(int j = 0; j != tetrahedron_sides; ++j) {
        // a region has all four vertices
        assert(r.vertices_[j].shared_vertex_data_ != no_vertex_idx);
        for(int k = 0; k != tetrahedron_sides; ++k) {
          if(j != k) {
            // no two vertices or neighbors of a region have the same identity
            assert(r.vertices_[j].shared_vertex_data_ != r.vertices_[k].shared_vertex_data_);
            if(r.neighbors_[j] != no_region_idx && r.neighbors_[k] != no_region_idx) {
              assert(r.neighbors_[j] != r.neighbors_[k]);
            }
          }
        }
      }
    }
    // region reference in-bound-ness
    for(size_t i = 0; i != regions_.size(); ++i) {
      region const& r = regions_[i];
      for(int j = 0; j != tetrahedron_sides; ++j) {
        assert(r.vertices_[j].shared_vertex_data_ == no_vertex_idx || r.vertices_[j].shared_vertex_data_ < vertices_.size());
        assert(r.neighbors_[j] == no_region_idx || r.neighbors_[j] < regions_.size());
      }
    }
    // region neighbor data consistency
    for(size_t i = 0; i != regions_.size(); ++i) {
      region const& r = regions_[i];
      for(int l = 0; l != tetrahedron_sides; ++l) {
        const region_idx_type ni = r.neighbors_[l];
        if(ni != no_region_idx) {
          assert(ni < regions_.size());
          // neighboring tetrahedra link to each other
          region const& n = regions_[ni];
          int reciprocal_links = 0;
          int reciprocal_link;
          for(int j = 0; j != tetrahedron_sides; ++j) {
            if(n.neighbors_[j] == i) {
              reciprocal_links += 1;
              reciprocal_link = j;
            }
          }
          assert(reciprocal_links == 1);

          // neighboring tetrahedra share three vertices
          int shared_vertices = 0;
          for(int j = 0; j != tetrahedron_sides; ++j) {
            for(int k = 0; k != tetrahedron_sides; ++k) {
              shared_vertices += (r.vertices_[j].shared_vertex_data_ == n.vertices_[k].shared_vertex_data_);
            }
          }
          assert(shared_vertices == 3);
          
          // ordering of neighbor links and of vertices are consistent
          // with each other
          for(int j = 0; j != tetrahedron_sides; ++j) {
            assert(r.vertices_[l].shared_vertex_data_ != n.vertices_[j].shared_vertex_data_);
            assert(r.vertices_[j].shared_vertex_data_ != n.vertices_[reciprocal_link].shared_vertex_data_);
          }
        }
      }
    }
  }
  // TODO thing-ness e.g. robots
  std::priority_queue<shared_ptr<event>> next_events_;
  //void player_input_becomes(time_type when, );
  //void insert_event(time_type when, );
  gl_triangles/*triangles*/ display(time_type when, vector3<distance> where) {
    //assert(when < next_event_time)
    gl_triangles triangles;
    for(region const& r : regions_) {
      array<gl_data_format::vertex_with_color, tetrahedron_sides> vertices;
      for(int i = 0; i < tetrahedron_sides; ++i) {
        region_vertex const& rv = r.vertices_[i];
        vertex const& v = vertices_[rv.shared_vertex_data_];
        const time_type relative_time = when - v.base_time_;
        const auto relative_timem = relative_time/identity(units_factor<1, 1000000000>());
        //TODO deal with times somehow more correctly
        const vector3<distance> loc =
              v.vertex_position_
            + v.vertex_velocity_*relative_timem/identity(units_factor<1, 1000>())
            + v.vertex_acceleration_*relative_timem*relative_timem/identity(units_factor<1, 1000000>())/2
            ;
        vertices[i] = gl_data_format::vertex_with_color(
          loc.x/distance_units, loc.y/distance_units, loc.z/distance_units,
          gl_data_format::color(0xffff0080));
        //std::cerr << vertices[i] << '\n';
      }
      for(int i = 0; i < tetrahedron_sides; ++i) {
        // draw edges (for debugging)
        for(int j = i+1; j < tetrahedron_sides; ++j) {
        }
        // draw sides
        gl_triangle triangle;
        int k = 0;
        //std::cerr << "\n\n";
        for(int j = 0; j < tetrahedron_sides; ++j) {
          if(i != j) {
            triangle.vertices[k++] = vertices[j];
            //std::cerr << triangle.vertices[k-1] << "  ";
          }
          //std::cerr << " ?"<<i << " ?"<<j << "\n";
        }
        //std::cerr << "\n";
        //triangles.push_back(triangle);
        push_wireframe_triangle(triangles, 10, triangle);
      }
    }
    sort_gl_triangles_by_distance_to(
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



void do_gl(grand_structure_of_lasercake& simulated_world, uint64_t frame) {
  //std::cerr<<"hi.\n";
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  //gluPerspective(90, 1, 1, 1000);
  //gluLookAt(20,20,20,0,0,0,0,0,1);
  //gluLookAt(0,0,0,1,1,1,0,0,1);
  const double wid = 500;
  gluPerspective(80, 1, 1, 2*wid);
  const vector3<double> viewcenter(0+wid*std::cos(double(frame) / 200), 0+wid*std::sin(double(frame) / 200), 0);
  gluLookAt(viewcenter.x, viewcenter.y, viewcenter.z,
            0,0,0, 0,0,1);

  gl_triangles data = simulated_world.display(7*pico*seconds, vector3<distance>(viewcenter*distance_units));
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
large_fast_noncrypto_rng rng(time(NULL));

  
  
int frame = 0;

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
          /*
          if(event.key.keysym.sym == SDLK_p) ++p_mode;
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
    
    do_gl(simulated_world, frame);
    glFinish();	
    SDL_GL_SwapBuffers();
   
    __attribute__((unused)) int before_processing = SDL_GetTicks();
    
    //doing stuff code here
    
    
    
	++frame;
    
    
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
