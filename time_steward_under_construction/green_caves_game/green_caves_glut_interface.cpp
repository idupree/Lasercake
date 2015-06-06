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


#include <iostream>
#include <cstddef>
#include "boost/date_time/posix_time/posix_time.hpp"

#include "green_caves_gl.hpp"
#include <GL/glut.h>

green_caves_ui_backend backend;

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
  gl_reshape(width, height);
}

static void Idle(void) {
  update();
  glutPostRedisplay();
}

static void Draw(void) {
  update();
  do_gl(backend);
  glutSwapBuffers();
}

int main(int argc, char **argv)
{
  for(int i = 1; i < argc; ++i) {
    if(strcmp(argv[i], "-b") == 0) {
      assert(false);
    }
  }
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
