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



#define GLM_FORCE_RADIANS

#include <iostream>
#include <cstddef>

#include "bouncy_circles.cpp"
#include "../testcase_draw.cpp"
#include <GL/glut.h>

static const double cube_size = 1;
static const time_type view_duration = time_units_per_gloppp * 10;

void draw_ztree_node(gl_triangles& triangles, time_steward::accessor const* accessor, entity_id id, gl_data_format::color color) {
  if (id) {
    auto const& node = accessor->get<bbcd_system::ztree_node>(accessor->get(id));
    bbcd_system::bounding_box bbox = node->here.get_bbox();
    if (bbox.size_minus_one()[0]<arena_width){
      gl_polygon polygon;
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
        (cube_size*accessor->now()/view_duration)-cube_size/2,
        bbox_to_space(bbox.min()[0])*cube_size/view_width,
        bbox_to_space(bbox.min()[1])*cube_size/view_width,
        color));
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
        (cube_size*accessor->now()/view_duration)-cube_size/2,
        bbox_to_space(bbox.min()[0])*cube_size/view_width,
        bbox_to_space(bbox.max()[1])*cube_size/view_width,
        color));
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
        (cube_size*accessor->now()/view_duration)-cube_size/2,
        bbox_to_space(bbox.max()[0])*cube_size/view_width,
        bbox_to_space(bbox.max()[1])*cube_size/view_width,
        color));
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
        (cube_size*accessor->now()/view_duration)-cube_size/2,
        bbox_to_space(bbox.max()[0])*cube_size/view_width,
        bbox_to_space(bbox.min()[1])*cube_size/view_width,
        color));
      push_wireframe_polygon(triangles, bbox.size_minus_one()[0]*cube_size/view_width/30, polygon);
    }
    for (entity_id i : node->children) {
      draw_ztree_node(triangles, accessor, i, color);
    }
  }
}

void draw_time(gl_triangles& triangles, time_steward& w, time_type time, gl_data_format::color color) {
  bbcd_system::coordinate_array min;
  bbcd_system::coordinate_array max;
  for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
    min[i] = space_to_bbox(-view_width/2);
    max[i] = space_to_bbox(view_width/2);
  }
  const bbcd_system::bounding_box view_box = bbcd_system::bounding_box::min_and_max(min,max);
  
  std::unique_ptr<time_steward::accessor> accessor = w.accessor_after(time);
  auto bbcd_id = accessor->get<global_data>(accessor->get(time_steward_system::global_object_id))->bbcd_id;
  std::unordered_set<entity_id> visible_circles = bbcd_operations::filter(accessor.get(), bbcd_id, circles_overlapping_bbox_filter(accessor.get(), view_box));
  for (entity_id id : visible_circles) {
    auto e = accessor->get(id);
    auto c = accessor->get<circle_shape>(e);
    gl_polygon polygon;
    double x = (double)c->center(0)(accessor->now())*cube_size/view_width;
    double y = (double)c->center(1)(accessor->now())*cube_size/view_width;
    double rad = (double)c->radius(accessor->now())*cube_size/view_width;
    std::cerr<<x<<","<<y<<","<<rad<<"\n";
    for (double theta = 0; theta < 6.283; theta += 0.3) {
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
        (cube_size*time/view_duration)-cube_size/2,
        x + std::cos(theta)*rad,
        y + std::sin(theta)*rad,
        accessor->get<circle_overlaps>(e)->empty() ? color : gl_data_format::color(0x00ffffb0)));
    }
    push_wireframe_polygon(triangles, rad/5, polygon);
  }
  draw_ztree_node(triangles, accessor.get(), accessor->get<bbcd_system::bbox_collision_detector_root_node>(accessor->get(bbcd_id))->root_node_id_, color);
}

gl_triangles display(vector3<double> const& where, time_steward& w, time_type focus_time) {
  gl_triangles triangles;
  {
    gl_polygon polygon;
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          -cube_size/2, -cube_size/2, -cube_size/2,
        gl_data_format::color(0x00800080)));
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          cube_size/2, -cube_size/2, -cube_size/2,
        gl_data_format::color(0xff800080)));
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          cube_size/2, cube_size/2, -cube_size/2,
        gl_data_format::color(0xff80ff80)));
      polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          -cube_size/2, cube_size/2, -cube_size/2,
        gl_data_format::color(0x0080ff80)));
    push_wireframe_polygon(triangles, cube_size/10, polygon);
  }
  
  draw_time(triangles, w, focus_time, gl_data_format::color(0x0000ffb0));
  for (time_type t = 0; t < view_duration; t += time_units_per_gloppp>>2) {
    //draw_time(triangles, w, t, gl_data_format::color(0xff0000b0));
  }
  
#if 0
  //for (time_type i = 1; i < 2500; i += 10) {
    std::vector<entity_id> const& list = w.get_entity_data_after(w.global_object_id(),time)->as<asteroid_ids_list>()->data;
    for (size_t j = 0; j < list.size(); ++j) {
      entity_id const& e = list[j];
      asteroid const& a = *w.get_entity_data_after(e,time)->as<asteroid>();
      gl_polygon polygon;
      double x = (double)a.center(0)(time);
      double y = (double)a.center(1)(time);
      double rad = (double)a.radius(time);
      for (double theta = 0; theta < 6.283; theta += 0.1) {
        polygon.vertices_.push_back(gl_data_format::vertex_with_color(
          0,//double(i*SCALE/5 - 25*SCALE),
          x + std::cos(theta)*rad,
          y + std::sin(theta)*rad,
          gl_data_format::color(j ? 0xff0000b0 : 0x0000ffb0)));
      }
      push_wireframe_polygon(triangles, rad/5, polygon);
    }
  //}
#endif
  sort_gl_triangles_far_to_near(
    glm::vec3(where.x, where.y, where.z),
    triangles);
  return triangles;
}

#if 0
std::string read_file(const char* filename) {
  // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
  // " Turns out that this method, while following STL idioms well, is actually surprisingly inefficient! Don't do this with large files. "
  std::ifstream t(filename);
  std::string str((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());
  return str;
}
#endif

void shader_debug_info(GLuint obj) {
  int infologLength = 0;
  int charsWritten  = 0;
  char* infoLog;

  glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &infologLength);

  if(infologLength > 0) {
    infoLog = (char *)malloc(infologLength);
    glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
    printf("%s\n",infoLog);
    free(infoLog);
  }
}

void program_debug_info(GLuint obj) {
  int infologLength = 0;
  int charsWritten  = 0;
  char* infoLog;

  glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infologLength);

  if(infologLength > 0) {
    infoLog = (char *)malloc(infologLength);
    glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
    printf("%s\n",infoLog);
    free(infoLog);
  }
}
GLuint shader_program_name;
void add_shaders() {
  GLuint v,f,f2,p;

  v = glCreateShader(GL_VERTEX_SHADER);
  f = glCreateShader(GL_FRAGMENT_SHADER);
 // f2 = glCreateShader(GL_FRAGMENT_SHADER);

  //std::string vs = read_file("horrible.vert");
  //std::string fs = read_file("horrible.frag");
  const char* vs =
#ifdef EMSCRIPTEN
"precision mediump float;\n"
#endif
"attribute vec3 position;\n"
"attribute vec4 color;\n"
"varying vec4 colorV;\n"
"void main() {\n"
"  colorV = color;\n"
"  gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1.0);\n"
"}\n"
;
  const char* fs =
#ifdef EMSCRIPTEN
"precision mediump float;\n"
#endif
"varying vec4 colorV;\n"
"void main() {\n"
"  gl_FragColor = colorV;\n"
"}\n"
;

  glShaderSource(v, 1, &vs, NULL);
  glShaderSource(f, 1, &fs, NULL);

  glCompileShader(v);
  glCompileShader(f);

  shader_debug_info(v);
  shader_debug_info(f);

  p = glCreateProgram();
  glAttachShader(p,v);
  glAttachShader(p,f);

  glLinkProgram(p);
  program_debug_info(p);

  glUseProgram(p);

  glEnableVertexAttribArray(0);
  shader_program_name = p;
}

void do_gl(time_steward& w, time_type time, double wid, double height_angle, double rot) {
  static bool first_time = true;//hack
  if(first_time) {
    add_shaders();
    first_time = false;
  }
  //std::cerr<<"hi.\n";
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(0.1, 1.0, 0.5, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  //gluPerspective(90, 1, 1, 1000);
  //gluLookAt(20,20,20,0,0,0,0,0,1);
  //gluLookAt(0,0,0,1,1,1,0,0,1);
  gluPerspective(80, 1, wid/1000, wid*2);
  double horiz_wid = wid * std::cos(height_angle);
  double vert_wid = wid * std::sin(height_angle);
  const vector3<double> viewcenter(0+horiz_wid*std::cos(rot), 0+horiz_wid*std::sin(rot), vert_wid);
  gluLookAt(viewcenter.x, viewcenter.y, viewcenter.z,
            0,0,0, 0,0,1);

  gl_triangles data = display(viewcenter, w, time);
  if(const size_t count = data.size()*3) {
    GLuint triangles_VBO_name;
    glGenBuffers(1, &triangles_VBO_name);
    glBindBuffer(GL_ARRAY_BUFFER, triangles_VBO_name);
    glBufferData(GL_ARRAY_BUFFER, count*sizeof(gl_data_format::vertex_with_color), &data[0], GL_STREAM_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, count*sizeof(vertex_with_color), data.vertices, GL_STREAM_DRAW);
    auto position_attrib_location = glGetAttribLocation(shader_program_name, "position");
    auto color_attrib_location = glGetAttribLocation(shader_program_name, "color");
        if(GLEW_VERSION_2_0) {
          std::cerr << "GLEW_VERSION_2_0\n";
#define BUFFER_OFFSET(i) ((void*)(i))
          glVertexAttribPointer(color_attrib_location, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(gl_data_format::vertex_with_color), BUFFER_OFFSET(0+offsetof(gl_data_format::vertex_with_color, c)));
          glVertexAttribPointer(position_attrib_location, 3, GL_FLOAT, GL_FALSE, sizeof(gl_data_format::vertex_with_color), BUFFER_OFFSET(0+offsetof(gl_data_format::vertex_with_color, v)));
          glEnableVertexAttribArray(color_attrib_location);
          glEnableVertexAttribArray(position_attrib_location);
          glBindBuffer(GL_ARRAY_BUFFER, 0);
        }else {
          glInterleavedArrays(GL_C4UB_V3F, 0, (void*)(0));
        }
        glDrawArrays(GL_TRIANGLES, 0, count);
        //glDrawElements for indexed mode
        if(GLEW_VERSION_2_0) {
          glDisableVertexAttribArray(position_attrib_location);
          glDisableVertexAttribArray(color_attrib_location);
        }
    //glVertexAttribPointer(, 3, GL_FLOAT, GL_FALSE, sizeof(gl_data_format::vertex_with_color), (void*)4);
    //glEnableVertexAttribArray(glGetAttribLocation(shader_program_name, "position"));
    //      glBindBuffer(GL_ARRAY_BUFFER, 0);
    //glInterleavedArrays(GL_C4UB_V3F, 0, (void*)(0));
    //glDrawArrays(GL_TRIANGLES, 0, count);
    glDeleteBuffers(1, &triangles_VBO_name);
    //glDisableVertexAttribArray(glGetAttribLocation(shader_program_name, "position"));
  }
}


time_steward steward;
double view_length = 2*cube_size;
double height_angle = 0;
double rot = 0;
time_type min_time = 0;
time_type gtime = min_time;
time_type max_time = 2*time_units_per_gloppp;

static void Key(unsigned char key, int x, int y) {
  switch (key) {
  case 27:
    exit(0); break;/*
  case 'q':
          if(sdle.key.keysym.sym == SDLK_q) { time -= 100; if (time < 0) time = 0; }
  case 'w':
          if(sdle.key.keysym.sym == SDLK_w) time += 100;*/
  case 'a':
    view_length *= 0.9; break;
  case 's':
    view_length /= 0.9; break;
  case 'e':
    height_angle += 0.2; break;
  case 'd':
    height_angle -= 0.2; break;
  }
  glutPostRedisplay();
}

static void Idle(void) {
  rot = rot + 0.01;
  gtime = gtime + (time_units_per_gloppp>> 6);
  if (gtime>max_time) gtime=min_time;
  glutPostRedisplay();
}

static void Draw(void) {
  do_gl(steward, gtime, view_length, height_angle, rot);
  glutSwapBuffers();
}

int main(int argc, char **argv)
{
  bounded_int_calculus::test();
  steward.insert_fiat_event(0, 0, std::shared_ptr<event>(new initialize_world()));
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutCreateWindow("bouncy circles");
  const GLenum glew_init_err = glewInit();
  if(glew_init_err != GLEW_OK) { throw "glew failed"; }
  glutKeyboardFunc(Key);
  glutDisplayFunc(Draw);
  glutIdleFunc(Idle);
  glutReshapeWindow(800,800);
  glutMainLoop();
  return 0;
}
