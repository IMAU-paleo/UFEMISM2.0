#! /bin/csh -f

# Safety
if ($#argv != 2) goto usage
if ($1 != 'dev' && $1 != 'perf') goto usage
if ($2 != 'changed' && $2 != 'clean') goto usage

echo ""

#  Confirm user's compilation choices
if ($1 == 'dev') then
  echo "dev: compiling UFEMISM, developers' version"
else if ($1 == 'perf') then
  echo "perf: compiling UFEMISM, performance version"
endif

if ($2 == 'changed') then
  echo "changed: (re)compiling changed modules only"
else if ($2 == 'clean') then
  echo "clean: recompiling all modules"
endif
echo ""

# Go to src/, make, and come back
cd src
if ($2 == 'clean') make clean
if ($1 == 'dev') then
  make all   DO_ASSERTIONS=yes   DO_RESOURCE_TRACKING=yes   DO_INCLUDE_COMPILER_CHECKS=yes
else if ($1 == 'perf') then
  make all   DO_ASSERTIONS=no    DO_RESOURCE_TRACKING=no    DO_INCLUDE_COMPILER_CHECKS=no
endif
cd ..

# Copy compiled program
if ($1 == 'dev') then

  rm -f UFEMISM_program_dev
  mv src/UFEMISM_program UFEMISM_program_dev
  rm -f UFEMISM_program
  cp UFEMISM_program_dev UFEMISM_program

else if ($1 == 'perf') then

  rm -f UFEMISM_program_perf
  mv src/UFEMISM_program UFEMISM_program_perf
  rm -f UFEMISM_program
  cp UFEMISM_program_perf UFEMISM_program

endif

exit 1

usage:

echo ""
echo "Usage: ./compile_UFEMISM.csh  [VERSION]  [SELECTION]"
echo ""
echo "  [VERSION]: dev, perf"
echo ""
echo "    dev : developer's version; extra compiler flags (see Makefile_include, COMPILER_FLAGS_CHECK),"
echo "          plus run-time assertions (i.e. DO_ASSERTIONS = true) and resource tracking (i.e."
echo "          DO_RESOURCE_TRACKING = true). Useful for tracking coding errors, but slow."
echo ""
echo "    perf: performance version; all of the above disabled. Coding errors will much more often"
echo "          simply result in segmentation faults, but runs much faster."
echo ""
echo "  [SELECTION]: changed, clean"
echo ""
echo "    changed: (re)compile only changed modules. Works 99% of the time, fails when you change the"
echo "             definition of a derived type without recompiling all of the modules thst use that"
echo "             type. In that case, better to just do a clean compilation."
echo ""
echo "    clean: recompile all modules. Always works, but slower."
echo ""