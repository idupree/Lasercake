#!/bin/sh
em++ -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-function -ferror-limit=5 -fdiagnostics-show-option -std=c++1y -I../../bundled_libs/boostbcp/ --js-library green_caves_lib.js -o green_caves.js green_caves_emscripten_interface.cpp ../csiphash.c "$@"
