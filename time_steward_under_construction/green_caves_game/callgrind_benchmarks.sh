#!/bin/sh
set -eu
mkdir -p zzzuncommitted/callgrindz
valgrind --tool=callgrind --dump-instr=yes --dump-line=yes --branch-sim=yes --cache-sim=yes --callgrind-out-file='zzzuncommitted/callgrindz/callgrind.out.%p' ./benchmarks "$@"
