#!/bin/sh
em++ -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-function -ferror-limit=5 -fdiagnostics-show-option -std=c++1y -I../../bundled_libs/boostbcp/ \
-s EXPORTED_FUNCTIONS="['_draw', '_update_to_real_time', '_mouse_down', '_mouse_up', '_mouse_moves', '_set_left', 'set_right', 'set_up', 'set_down']" \
--js-library green_caves_lib.js -o green_caves.js green_caves_emscripten_interface.cpp ../csiphash.c "$@"
