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


#include "green_caves.cpp"

green_caves_ui_backend backend;

extern "C" {

extern void draw_circle(double cx, double cy, double r);
extern void draw_rect(double x0, double y0, double x1, double y1);
extern void draw_segment(double x0, double y0, double x1, double y1, double width);

} // extern "C"

struct draw_funcs {
  void circle(double cx, double cy, double r) {
    draw_circle(cx, cy, r);
  }
  void rect(double x0, double y0, double x1, double y1) {
    draw_rect(x0, y0, x1, y1);
  }
  void segment(double x0, double y0, double x1, double y1, double width) {
    draw_segment(x0, y0, x1, y1, width);
  }
};

extern "C" {

// when you add something here, you need to also add to the EXPORTED_FUNCTIONS
// in the build script

void draw() { draw_funcs draw; backend.draw(draw); }
void update_to_real_time(double milliseconds) { backend.update_to_real_time(int64_t(milliseconds)); }
void mouse_down(int x, int y) { backend.mouse_down(x, y); }
void mouse_up(int x, int y) { backend.mouse_up(x, y); }
void mouse_moves(int x, int y) { backend.mouse_moves(x, y); }
void set_left(bool b) { backend.left = b; }
void set_right(bool b) { backend.right = b; }
void set_up(bool b) { backend.up = b; }
void set_down(bool b) { backend.down = b; }

} // extern "C"
