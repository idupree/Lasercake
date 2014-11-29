#!/bin/sh
emcc --js-library green_caves_lib.js -o green_caves.js green_caves_emscripten_interface.cpp ../csiphash.c "$@"
