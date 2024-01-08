#!/bin/bash

# Compile the program
make gomp-only

# Optional: Add any initial test commands here
./gomp-only i1.coord coStudent1.dat foStudent1.dat noStudent1.dat
./gomp-only i2.coord coStudent2.dat foStudent2.dat noStudent2.dat
./gomp-only i3.coord coStudent3.dat foStudent3.dat noStudent3.dat