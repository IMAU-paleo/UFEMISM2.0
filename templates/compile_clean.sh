#!/bin/bash

# Go to src/ make, and come back
cd src
make clean
make all
cd ..

# Update program
rm -f UFEMISM_program
mv src/UFEMISM_program .
