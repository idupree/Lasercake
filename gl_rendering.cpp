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

#include <cstddef>
#include <vector>
#include "cxx11/array.hpp"

#include <glm/gtc/type_ptr.hpp>

#include "gl.hpp"
#include "gl_rendering.hpp"
#include "gl_data_format.hpp"

using namespace gl_data_format;

#define BUFFER_OFFSET(i) ((void*)(i))
const GLuint INVALID_BUFFER_ID = 0;

// by_distance_VBO_* are indexed by distance*DISTANCE_IDX_FACTOR + this enum:
enum by_distance_idx_adjustment_enum {
  QUADS_IDX = 0, TRIANGLES_IDX, LINES_IDX, POINTS_IDX,
  DISTANCE_IDX_FACTOR
};

typedef array<vertex_with_color, 6> rect_type;

struct gl_renderer::state_t_ {
  GLuint rect_VBO_name;
  std::vector<GLuint> by_distance_VBO_names;
  std::vector<size_t> by_distance_VBO_sizes;
  //GLuint vertex_shader_name;
  //GLuint fragment_shader_name;
  GLuint shader_program_name;
  GLint matrix_uniform_location;
  GLint position_attrib_location;
  GLint color_attrib_location;
};

gl_renderer::gl_renderer() {}
gl_renderer::~gl_renderer() {}

GLuint compile_shader(GLenum shader_type, const char* program) {
  const GLuint shader_name = glCreateShader(shader_type);
  glShaderSource(shader_name, 1, &program, nullptr);
  glCompileShader(shader_name);
  GLint good;
  glGetShaderiv(shader_name, GL_COMPILE_STATUS, &good);
  if(!good) {
    GLint log_length;
    glGetShaderiv(shader_name, GL_INFO_LOG_LENGTH, &log_length);
    std::vector<char> log(log_length);
    glGetShaderInfoLog(shader_name, log_length, nullptr, &log[0]);
    LOG << "Shader compile error: " << &log[0];
    glDeleteShader(shader_name);
    return 0;
  }
  return shader_name;
}

template<typename CollectionOfGLuint>
GLuint link_shaders(CollectionOfGLuint shader_names)
{
  for(GLuint shader_name : shader_names) {
    if(shader_name == 0) {
      return 0;
    }
  }
  const GLuint program_name = glCreateProgram();
  for(GLuint shader_name : shader_names) {
    glAttachShader(program_name, shader_name);
  }
  glLinkProgram(program_name);
  GLint good;
  glGetProgramiv(program_name, GL_LINK_STATUS, &good);
  if(!good) {
    GLint log_length;
    glGetProgramiv(program_name, GL_INFO_LOG_LENGTH, &log_length);
    std::vector<char> log(log_length);
    glGetProgramInfoLog(program_name, log_length, nullptr, &log[0]);
    LOG << "Shader link error: " << &log[0];
    glDeleteProgram(program_name);
    return 0;
  }
  return program_name;
}

GLuint make_shaders(const char* vertex_program, const char* fragment_program) {
  const GLuint vertex_shader_name = compile_shader(GL_VERTEX_SHADER, vertex_program);
  const GLuint fragment_shader_name = compile_shader(GL_FRAGMENT_SHADER, fragment_program);
  const array<GLuint, 2> shader_names = {{ vertex_shader_name, fragment_shader_name }};
  const GLuint program_name = link_shaders(shader_names);
  for(GLuint shader_name : shader_names) {
    glDeleteShader(shader_name);
  }
  return program_name;
}


void gl_renderer::output_gl_data_to_OpenGL(
    abstract_gl_data const& abstract_gl_data,
    viewport_dimension viewport_width,
    viewport_dimension viewport_height,
    LasercakeGLWidget* gl_widget,
    // volatile: Ensure that every load is in fact done,
    // even though they look redundant.
    atomic::atomic_bool const volatile& interrupt
) {
  // This code is intended to have no operations in it that can possibly
  // throw exceptions, to simplify dealing with OpenGL context state.
  // When allocating, e.g. via std::vector, wrap it in a try/catch
  // and do something sensible if there's an exception.

  gl_all_data const& gl_data = abstract_gl_data.data();

  if(interrupt.load(atomic::memory_order_relaxed)) {return;};

  if(!state_) {
    try {
      state_.reset(new state_t_);
    }
    catch(std::bad_alloc const&) {
      return;
    }

    if(!init_gl_library())
    {
      // give up
      state_.reset();
      return;
    }
    if(interrupt.load(atomic::memory_order_relaxed)) {return;};

    // webgl uses shader #version 100
    const char* vertex_shader_source =
    "#version 100\n"
    "attribute vec3 position;\n"
    "attribute vec4 color;\n"
    "varying vec4 colorV;\n"
    "uniform mat4 matrix;\n"
    "void main() {\n"
    "  colorV = color;\n"
#if 1
    "  gl_Position = matrix * vec4(position, 1.0);\n"
#else /* useless, but entertaining to see the potential of webgl */
    "  vec4 hm = matrix * vec4(position, 1.0);\n"
    "  float d = distance(hm, vec4(1000.0, 500.0, 0.0, 1.0));\n"
    "  float f = d*d*0.00001 - 10000.0;\n"
    "  gl_Position = hm + vec4(f, f, 0.0, 0.0);\n"
#endif
    "}\n"
    ;
    const char* fragment_shader_source =
    "#version 100\n"
    "precision mediump float;\n"
    "varying vec4 colorV;\n"
    "void main() {\n"
    "  gl_FragColor = colorV;\n"/*
    //"  gl_FragColor = vec4 ( 0.0, 1.0, 0.0, 1.0 );\n"
    //"  gl_FragColor = colorV*colorV;\n"
    "  vec4 col = colorV*colorV;\n"
    "  col.xyz *= (gl_FragCoord.x*0.003);\n"
    "  gl_FragColor = col;\n"*/
    "}\n"
    ;
    state_->shader_program_name = make_shaders(vertex_shader_source, fragment_shader_source);
    assert(state_->shader_program_name);
    if(!state_->shader_program_name) {
      // give up
      state_.reset();
      return;
    }
    state_->matrix_uniform_location = glGetUniformLocation(state_->shader_program_name, "matrix");
    assert(state_->matrix_uniform_location != -1);
    state_->position_attrib_location = glGetAttribLocation(state_->shader_program_name, "position");
    assert(state_->position_attrib_location != -1);
    state_->color_attrib_location = glGetAttribLocation(state_->shader_program_name, "color");
    assert(state_->color_attrib_location != -1);

    glGenBuffers(1, &state_->rect_VBO_name);
    glBindBuffer(GL_ARRAY_BUFFER, state_->rect_VBO_name);
    glBufferData(GL_ARRAY_BUFFER, sizeof(rect_type), nullptr, GL_STREAM_DRAW);

    if(interrupt.load(atomic::memory_order_relaxed)) {return;};
  }

  glViewport(0, 0, viewport_width, viewport_height);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  // Depth func LEQUAL not LESS.  We try to draw objects in back-to-front
  // order, so rounding error means the front (later) one should win.
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#if MODERN_GL_PLEASE
  glUseProgram(state_->shader_program_name);

  const glm::mat4 ortho_matrix = glm::ortho<float>(0, 1, 0, 1);
  //glUniformMatrix4fv(state_->matrix_uniform_location, 1, false, glm::value_ptr(ortho_matrix));

  const glm::mat4 projection_matrix =
    make_projection_matrix(float(viewport_width) / float(viewport_height));
  const glm::mat4 view_matrix =
    make_view_matrix(gl_data.facing, gl_data.facing_up);
  const glm::mat4 matrix = projection_matrix * view_matrix;
  glUniformMatrix4fv(state_->matrix_uniform_location, 1, false, glm::value_ptr(matrix));

#else
// The OpenGL ES headers don't let this compile, and GLEW doesn't
// work out of the box with Emscripten, so it's ifdef'd instead of
// testing GLEW_VERSION_2_0 at runtime.

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  const glm::mat4 projection_matrix =
    make_projection_matrix(float(viewport_width) / float(viewport_height));
  glLoadMatrixf(glm::value_ptr(projection_matrix));

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  const glm::mat4 view_matrix =
    make_view_matrix(gl_data.facing, gl_data.facing_up);
  glLoadMatrixf(glm::value_ptr(view_matrix));
#endif

#if !MODERN_GL_PLEASE
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
#endif

  if(interrupt.load(atomic::memory_order_relaxed)) {return;};
size_t vertices = 0;
  if(gl_data.stuff_to_draw_as_gl_collections_by_distance.size() > state_->by_distance_VBO_names.size()) {
    const size_t new_buffers_base = state_->by_distance_VBO_names.size()*DISTANCE_IDX_FACTOR;
    const size_t num_new_buffers = gl_data.stuff_to_draw_as_gl_collections_by_distance.size()*DISTANCE_IDX_FACTOR - new_buffers_base;
    try {
      state_->by_distance_VBO_names.resize(new_buffers_base + num_new_buffers);
      state_->by_distance_VBO_sizes.resize(new_buffers_base + num_new_buffers);
    }
    catch(std::bad_alloc const&) {
      return;
    }
    glGenBuffers(num_new_buffers, &state_->by_distance_VBO_names[new_buffers_base]);
    for(size_t i = 0; i != num_new_buffers; ++i) {
      if(interrupt.load(atomic::memory_order_relaxed)) {return;};
      state_->by_distance_VBO_sizes[new_buffers_base + i] = 0;
    }
  }
//TODO close->far nearby opaque things with z buffer
// to help out GL throw out pixels etc
  for(size_t dist_plus_one = gl_data.stuff_to_draw_as_gl_collections_by_distance.size(); dist_plus_one != 0; --dist_plus_one) {
    if(interrupt.load(atomic::memory_order_relaxed)) {return;};
    const size_t dist = dist_plus_one - 1;
    gl_collection const& coll = gl_data.stuff_to_draw_as_gl_collections_by_distance[dist];
    struct polygon_type {
      GLenum gl_type;
      by_distance_idx_adjustment_enum our_idx_adj;
      gl_call_data gl_collection::* gl_data_container_ptr_to_member;
    };
#if MODERN_GL_NO_QUADS
    caller_error_if(coll.quads.size() > 0, "quads don't exist in OpenGL ES; use two triangles");
#endif
    const array<polygon_type, DISTANCE_IDX_FACTOR> types = {{
#if !MODERN_GL_NO_QUADS
      { GL_QUADS, QUADS_IDX, &gl_collection::quads },
#endif
      { GL_TRIANGLES, TRIANGLES_IDX, &gl_collection::triangles },
      { GL_LINES, LINES_IDX, &gl_collection::lines },
      { GL_POINTS, POINTS_IDX, &gl_collection::points },
    }};
    for (polygon_type type : types) {
      if(interrupt.load(atomic::memory_order_relaxed)) {return;};
      gl_call_data const& data = coll.*(type.gl_data_container_ptr_to_member);
      if(const size_t count = data.size()) {
        vertices += count;
        const size_t buf_name_idx = dist*DISTANCE_IDX_FACTOR + type.our_idx_adj;
        glBindBuffer(GL_ARRAY_BUFFER, state_->by_distance_VBO_names[buf_name_idx]);
        if(state_->by_distance_VBO_sizes[buf_name_idx] < count) {
          glBufferData(GL_ARRAY_BUFFER, count*sizeof(vertex_with_color), data.vertices, GL_STREAM_DRAW);
          state_->by_distance_VBO_sizes[buf_name_idx] = count;
        }
        else {
          glBufferSubData(GL_ARRAY_BUFFER, 0, count*sizeof(vertex_with_color), data.vertices);
        }
        #if MODERN_GL_PLEASE
          glVertexAttribPointer(state_->position_attrib_location, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_color), BUFFER_OFFSET(0+offsetof(vertex_with_color, v)));
          glEnableVertexAttribArray(state_->position_attrib_location);
          glVertexAttribPointer(state_->color_attrib_location, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(vertex_with_color), BUFFER_OFFSET(0+offsetof(vertex_with_color, c)));
          glEnableVertexAttribArray(state_->color_attrib_location);
          glBindBuffer(GL_ARRAY_BUFFER, 0);
        #else
          glInterleavedArrays(GL_C4UB_V3F, 0, BUFFER_OFFSET(0));
        #endif
        glDrawArrays(type.gl_type, 0, count);
        //glDrawElements for indexed mode
        #if MODERN_GL_PLEASE
          glDisableVertexAttribArray(state_->position_attrib_location);
          glDisableVertexAttribArray(state_->color_attrib_location);
        #endif
      }
    }
  }
  if(interrupt.load(atomic::memory_order_relaxed)) {return;};
LOG << "\n!!!!!!!!! " << vertices << " VERTICES!!!!!!!!\n";
  // Is there a simpler way to tint the whole screen a color?
  const color tint = gl_data.tint_everything_with_this_color;
  const rect_type rect = {{
    vertex_with_color(0, 0, 0, tint),
    vertex_with_color(0, 1, 0, tint),
    vertex_with_color(1, 1, 0, tint),
    vertex_with_color(1, 1, 0, tint),
    vertex_with_color(1, 0, 0, tint),
    vertex_with_color(0, 0, 0, tint)
  }};
  // This screen-tinting code may not even be working (neither
  // with nor without MODERN_GL_PLEASE). No visual change shows up.
  // (Even on the rare occasion when it is supposed to.)
  glDisable(GL_DEPTH_TEST);
  glDepthFunc(GL_ALWAYS);
#if MODERN_GL_PLEASE
  glUniformMatrix4fv(state_->matrix_uniform_location, 1, false, glm::value_ptr(ortho_matrix));
#else
  glLoadIdentity();
  glOrtho(0, 1, 0, 1, -1, 1);
#endif
  glBindBuffer(GL_ARRAY_BUFFER, state_->rect_VBO_name);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(rect), &rect[0]);
#if !MODERN_GL_PLEASE
  glInterleavedArrays(GL_C4UB_V3F, 0, BUFFER_OFFSET(0));
#endif
  glDrawArrays(GL_TRIANGLES, 0, 6);
  glBindBuffer(GL_ARRAY_BUFFER, INVALID_BUFFER_ID);

  if(interrupt.load(atomic::memory_order_relaxed)) {return;};

  if(gl_widget) {
    render_2d_text_overlay_(abstract_gl_data, viewport_width, viewport_height, *gl_widget);
  }
}

void gl_renderer::fini() {
  if(state_) {
#if MODERN_GL_PLEASE
    glDeleteProgram(state_->shader_program_name);
#endif
    glDeleteBuffers(1, &state_->rect_VBO_name);
    glDeleteBuffers(state_->by_distance_VBO_names.size(),
                       &state_->by_distance_VBO_names[0]);
    // by_distance_VBO_sizes does not contain OpenGL-owned data
    state_.reset();
  }
}

