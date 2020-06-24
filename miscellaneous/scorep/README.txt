To instrument the code with scorep.
Use the profiling-environment docker image.
Start from a pristine build.

One difficulty commes from the fact CMake prohibits to change the compiler after the CMake step
and it also prohibits that the compiler variable value includes any flags
Cf. http://scorepci.pages.jsc.fz-juelich.de/scorep-pipelines/docs/scorep-6.0/html/scorepwrapper.html

So we have to use the score-p compiler wrapper with SCOREP_WRAPPER set to off
during the configuration stage but set to on (or unset) during the compilation stage.

For the time being this can be achieved in two ways:
- either the OLD_CMAKE_INSTALL option with cmake/with-scorep cmake ...
- either use two steps:
  1/ cmake/with-scorep pip install -e .
  2/ locate the build directory (usually build/temp-xxxxxx) and from there do : make clean and make install
