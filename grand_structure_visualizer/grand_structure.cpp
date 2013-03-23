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


#ifdef LASERCAKE_HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>
#include <boost/chrono/thread_clock.hpp>


#include "../utils.hpp"
//#include "../data_structures/geometry.hpp"

//typedef seconds<non_normalized_rational<64>> time_type;
//TODO make this size better for long durations or something
//long distances
//128
typedef picoseconds<lint64_t> time_type;
typedef nanometers<lint64_t> distance;
typedef nanometers/second<lint64_t> velocity1d;
typedef nanometers/second/second<lint64_t> acceleration1d;

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
struct region {
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
struct edge_edge_collision : public collision {
};

class grand_structure_of_lasercake {
  std::priority_queue<shared_ptr<event>> next_events;
private:
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
  }
};















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
    gScreen = SDL_SetVideoMode (640, 640, 0, flags);
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

    //drawing code here
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    gluPerspective(80, 1, 100000, 10000000);
    gluLookAt(0+2500000*std::cos(double(frame) / 200),0+2500000*std::sin(double(frame) / 200),1500000,0,0,0,0,0,1);

    
    __attribute__((unused)) int before_GL = SDL_GetTicks();
    
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
    
    // Init GL state
	gluPerspective(90, 1, 10, 1000);
	gluLookAt(20,20,20,0,0,0,0,0,1);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
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
