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

(function(){
  "use strict";
  mergeInto(LibraryManager.library, {
    draw_rect: function(x,y,w,h) {
      green_caves.canvas_context.fillStyle = "#00ff00";
      green_caves.canvas_context.fillRect(x,y,w,h);
    },
    draw_circle: function(x,y,r) {
      green_caves.canvas_context.fillStyle = "#00ff00";
      green_caves.canvas_context.arc(x, y, r, 0, 2 * Math.PI);
      green_caves.canvas_context.fill();
    },
    draw_segment: function(x0,y0,x1,y1,w) {
      green_caves.canvas_context.strokeStyle = "#00ff00";
      green_caves.canvas_context.lineWidth = w;
      green_caves.canvas_context.beginPath();
      green_caves.canvas_context.moveTo(x0, y0);
      green_caves.canvas_context.lineTo(x1, y1);
      green_caves.canvas_context.stroke();
    },
  });
})();
