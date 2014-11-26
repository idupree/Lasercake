#!/bin/sh
g++ -fdiagnostics-color=auto -o sim -ggdb -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-function -Wno-unused-local-typedefs -fmax-errors=5 -fdiagnostics-show-option -std=c++1y -I/usr/include/SDL/ -lSDL -lGL -lGLU -lGLEW -lglut -lgmpxx -lgmp ../csiphash.c green_caves_glut_interface.cpp "$@"
