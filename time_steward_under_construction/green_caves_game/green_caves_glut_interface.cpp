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
#include "boost/date_time/posix_time/posix_time.hpp"

#include "green_caves.cpp"
#include "../testcase_draw.cpp"
#include <GL/glut.h>

green_caves_ui_backend backend;

gl_data_format::color green_color(0x00ff00ff);

struct draw_funcs {
  draw_funcs(gl_triangles& triangles) : triangles(triangles){}
  gl_triangles& triangles;
  void circle(double cx, double cy, double r) {
    gl_data_format::vertex_with_color first;
    gl_data_format::vertex_with_color prev;
    gl_data_format::vertex_with_color center(cx, cy, 0, green_color);
    for (double theta = 0; ; theta += 0.3) {
      gl_data_format::vertex_with_color cur(
        cx + std::cos(theta)*r,
        cy + std::sin(theta)*r,
        0,
        green_color);
      if (theta == 0) {
        first = cur;
      }
      else {
        if (theta >= 6.283) {
          triangles.push_back(gl_triangle{{{ center, prev, first }}});
          break;
        }
        else {
          triangles.push_back(gl_triangle{{{ center, prev, cur }}});
        }
      }
      prev = cur;
    }
  }
  void rect(double x0, double y0, double x1, double y1) {
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0, y0, 0, green_color),
      gl_data_format::vertex_with_color(x1, y0, 0, green_color),
      gl_data_format::vertex_with_color(x1, y1, 0, green_color) }}});
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0, y0, 0, green_color),
      gl_data_format::vertex_with_color(x0, y1, 0, green_color),
      gl_data_format::vertex_with_color(x1, y1, 0, green_color) }}});
  }
  void segment(double x0, double y0, double x1, double y1, double width) {
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0+width, y0, 0, green_color),
      gl_data_format::vertex_with_color(x0-width, y0, 0, green_color),
      gl_data_format::vertex_with_color(x1, y1, 0, green_color) }}});
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0, y0+width, 0, green_color),
      gl_data_format::vertex_with_color(x0, y0-width, 0, green_color),
      gl_data_format::vertex_with_color(x1, y1, 0, green_color) }}});
  }
};

gl_triangles display() {
  gl_triangles triangles;
  draw_funcs draw(triangles);
  backend.draw(draw);
  
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

void do_gl() {
  static bool first_time = true;//hack
  if(first_time) {
    add_shaders();
    first_time = false;
  }
  //std::cerr<<"hi.\n";
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  gl_triangles data = display();
  if(const size_t count = data.size()*3) {
    GLuint triangles_VBO_name;
    glGenBuffers(1, &triangles_VBO_name);
    glBindBuffer(GL_ARRAY_BUFFER, triangles_VBO_name);
    glBufferData(GL_ARRAY_BUFFER, count*sizeof(gl_data_format::vertex_with_color), &data[0], GL_STREAM_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, count*sizeof(vertex_with_color), data.vertices, GL_STREAM_DRAW);
    auto position_attrib_location = glGetAttribLocation(shader_program_name, "position");
    auto color_attrib_location = glGetAttribLocation(shader_program_name, "color");
        if(GLEW_VERSION_2_0) {
          //std::cerr << "GLEW_VERSION_2_0\n";
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


boost::posix_time::ptime ptmin(boost::posix_time::min_date_time);
int64_t ms() { return (boost::posix_time::microsec_clock::universal_time() - ptmin).total_milliseconds(); }
void update() {
  backend.update_to_real_time(ms());
}

static void keydown(unsigned char key, int /*x*/, int /*y*/) {
  switch (key) {
    case 27: exit(0); break;
    case 'w': backend.set_key(ms(), UP, true); break;
    case 'a': backend.set_key(ms(), LEFT, true); break;
    case 's': backend.set_key(ms(), DOWN, true); break;
    case 'd': backend.set_key(ms(), RIGHT, true); break;
  }
  if (key >= '0' && key <= '9') { backend.set_time_rate(key-'0'); }
}
static void keyup(unsigned char key, int /*x*/, int /*y*/) {
  switch (key) {
    case 27: exit(0); break;
    case 'w': backend.set_key(ms(), UP, false); break;
    case 'a': backend.set_key(ms(), LEFT, false); break;
    case 's': backend.set_key(ms(), DOWN, false); break;
    case 'd': backend.set_key(ms(), RIGHT, false); break;
  }
}

static void mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      backend.mouse_down(ms(), x, y);
    }
    else {
      backend.mouse_up(ms(), x, y);
    }
  }
}
static void mouse2(int x, int y) {
  backend.mouse_moves(ms(), x, y);
}
void reshape(int width, int height) {
  update();
  backend.screen_size[0] = width;
  backend.screen_size[1] = height;
  glViewport(0,0,width,height);
  glLoadIdentity();
  gluOrtho2D(0, GLdouble(width), 0, GLdouble(height));
}

static void Idle(void) {
  update();
  glutPostRedisplay();
}

static void Draw(void) {
  update();
  do_gl();
  glutSwapBuffers();
}

int main(int argc, char **argv)
{
  srand(0);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutCreateWindow("green caves");
  const GLenum glew_init_err = glewInit();
  if(glew_init_err != GLEW_OK) { throw "glew failed"; }
  glutKeyboardFunc(keydown);
  glutKeyboardUpFunc(keyup);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse2);
  glutReshapeFunc(reshape);
  glutPassiveMotionFunc(mouse2);
  glutDisplayFunc(Draw);
  glutIdleFunc(Idle);
  glutReshapeWindow(800,800);
  glutMainLoop();
  return 0;
}
