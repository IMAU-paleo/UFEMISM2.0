#! /bin/csh -f

# Compile all the UFEMISM code.
#
# Usage: ./compile_UFEMISM.csh  [VERSION]  [SELECTION]
#
#   [VERSION]: dev, perf
#
#     dev : developer's version; extra compiler flags (see Makefile_include, COMPILER_FLAGS_CHECK),
#           plus run-time assertions (i.e. DO_ASSERTIONS = true) and resource tracking (i.e.
#           DO_RESOURCE_TRACKING = true). Useful for tracking coding errors, but slow.
#
#     perf: performance version; all of the above disabled. Coding errors will much more often
#           simply result in segmentation faults, but runs much faster.
#
#   [SELECTION]: changed, clean
#
#     changed: (re)compile only changed modules. Works 99% of the time, fails when you change the
#              definition of a derived type without recompiling all of the modules thst use that
#              type. In that case, better to just do a clean compilation.
#
#     clean: recompile all modules. Always works, but slower.
#

# Safety
if ($#argv != 2) goto usage

echo ""

#  Confirm user's compilation choices
set version   = $argv[1]
if ($version == 'dev') then
  echo "dev: compiling UFEMISM, developers' version"
else if ($version == 'perf') then
  echo "perf: compiling UFEMISM, performance version"
else
  goto usage
endif

set selection = $argv[2]
if ($selection == 'changed') then
  echo "changed: (re)compiling changed modules only"
else if ($selection == 'clean') then
  echo "clean: recompiling all modules"
else
  goto usage
endif

echo ""

# Go to src/, make, and come back
cd src
if ($selection == 'clean') make clean
if ($version == 'dev') then
  make all   DO_ASSERTIONS=yes   DO_RESOURCE_TRACKING=yes   DO_INCLUDE_COMPILER_CHECKS=yes
else if ($version == 'perf') then
  make all   DO_ASSERTIONS=no    DO_RESOURCE_TRACKING=no    DO_INCLUDE_COMPILER_CHECKS=no
endif
cd ..

# Copy compiled program
if ($version == 'dev') then

  rm -f UFEMISM_program_dev
  mv src/UFEMISM_program UFEMISM_program_dev
  rm -f UFEMISM_program
  cp UFEMISM_program_dev UFEMISM_program

else if ($version == 'perf') then

  rm -f UFEMISM_program_perf
  mv src/UFEMISM_program UFEMISM_program_perf
  rm -f UFEMISM_program
  cp UFEMISM_program_perf UFEMISM_program

endif

exit 0

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

exit 1