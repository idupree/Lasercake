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
//#define GLM_PRECISION_MEDIUMP_FLOAT

#include <iostream>
#include <cstddef>

#include "green_caves.cpp"
#include "../testcase_draw.cpp"

gl_data_format::color green_color(0x00ff00ff);
gl_data_format::color black_color(0x000000ff);

struct gl_draw_funcs {
  gl_draw_funcs(gl_triangles& triangles) : triangles(triangles){}
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
  void rect(double x0, double y0, double x1, double y1, bool color) {
    gl_data_format::color c = color ? green_color : black_color;
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0, y0, 0, c),
      gl_data_format::vertex_with_color(x1, y0, 0, c),
      gl_data_format::vertex_with_color(x1, y1, 0, c) }}});
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0, y0, 0, c),
      gl_data_format::vertex_with_color(x0, y1, 0, c),
      gl_data_format::vertex_with_color(x1, y1, 0, c) }}});
  }
  void segment(double x0, double y0, double x1, double y1, double width) {
    double m = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
    double dx = (y1-y0)*width/m/2;
    double dy = -(x1-x0)*width/m/2;
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x0+dx, y0+dy, 0, green_color),
      gl_data_format::vertex_with_color(x0-dx, y0-dy, 0, green_color),
      gl_data_format::vertex_with_color(x1+dx, y1+dy, 0, green_color) }}});
    triangles.push_back(gl_triangle{{{
      gl_data_format::vertex_with_color(x1+dx, y1+dy, 0, green_color),
      gl_data_format::vertex_with_color(x1-dx, y1-dy, 0, green_color),
      gl_data_format::vertex_with_color(x0-dx, y0-dy, 0, green_color) }}});
  }
};

gl_triangles display(green_caves_ui_backend& backend) {
  gl_triangles triangles;
  gl_draw_funcs draw(triangles);
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
"uniform mat4 transformation;\n"
"attribute vec3 position;\n"
"attribute vec4 color;\n"
"varying vec4 colorV;\n"
"void main() {\n"
"  colorV = color;\n"
// This is the deprecated way to do it that doesn't work at all in OpenGL ES 2.0 / WebGL:
//"  gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1.0);\n"
"  gl_Position = transformation * vec4(position, 1.0);\n"
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

void do_gl(green_caves_ui_backend& backend) {
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

  gl_triangles data = display(backend);
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

void gl_reshape(int width, int height) {
  glViewport(0,0,width,height);
  const GLint uniform_id = glGetUniformLocation(shader_program_name, "transformation");
  const glm::mat4 transformation = glm::ortho(0.0f, GLfloat(width), 0.0f, GLfloat(height));
  glUniformMatrix4fv(uniform_id, 1, GL_FALSE, glm::value_ptr(transformation));
// This is the deprecated way to do it that doesn't work at all in OpenGL ES 2.0 / WebGL:
//  glLoadIdentity();
//  gluOrtho2D(0, GLdouble(width), 0, GLdouble(height));
}
