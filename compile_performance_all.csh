#! /bin/csh -f

# Go to src/, make, and come back
cd src
make all DO_ASSERTIONS=no DO_RESOURCE_TRACKING=no
cd ..

# Update program
rm -f UFEMISM_program
mv src/UFEMISM_program .
