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

#ifndef LASERCAKE_GL_HPP__
#define LASERCAKE_GL_HPP__

// for portable header inclusion...

#if LASERCAKE_USE_QT

#include <GL/glew.h>

inline bool init_gl_library() {
  return glewInit() == GLEW_OK;
}

#else

#include <GLES2/gl2.h>

inline bool init_gl_library() {
  return true;
}

#endif

#if EXTREMELY_OBSOLETE_GL_BUFFER_FUNC_NAMES
// Using the *ARB versions is needed for OS X 10.6 macbook2,1,
// but slows down Nouveau because it probably gets confused between
// old functions names and new GLSL code.
#define glGenBuffers glGenBuffersARB
#define glBindBuffer glBindBufferARB
#define glBufferData glBufferDataARB
#define glBufferSubData glBufferSubDataARB
#define glDeleteBuffers glDeleteBuffersARB
#endif

#define MODERN_GL_PLEASE 1
#define MODERN_GL_NO_QUADS 1

#endif

