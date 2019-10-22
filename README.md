[![Codacy Badge](https://api.codacy.com/project/badge/Grade/dc251db65f7a4c06ae07380544ea08fc)](https://www.codacy.com/manual/hrosiak/catima?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hrosiak/catima&amp;utm_campaign=Badge_Grade)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/hrosiak/catima.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/hrosiak/catima/context:cpp)
[![Documentation Status](https://readthedocs.org/projects/catima/badge/?version=latest)](https://catima.readthedocs.io/en/latest/?badge=latest)

CATima
=======
C++ library for caclulalaton of energy loss, range, angular scattering and time of flight of the particle passing through matter.
The library is based on physics used in the ATIMA code,however its not 100% copy of ATIMA physics.
 see CREDITS for more details.

The WebAtima UI to this library can be found here:
  * https://web-docs.gsi.de/~aprochaz/webatima (only inside GSI)
  * https://isotopea.com/webatima

Installation
------------
CMake is used to build the library. For default build use:

```
> mkdir build
> cd build
> cmake ../
> make
```

after the compilation the libraries and headers must be either installed system-wide by make install or PATH and LD_LIBRARY_PATH must be adjusted to point to headers and library files.
The default install path can be change, ie: cmake -DCMAKE_INSTALL_PREFIX=/opt/catima

the option to system-wide installation is to adjust library path and include paths.
This can be done sourcing the init.sh file, which is generated in the build directory:
```
> source init.sh
```

cmake options
-------------
compile options, enable or disable with cmake:
> cmake ../ -D[OPTION]

available options:
  * BUILD_SHARED_LIBS - if ON shared library is build, otherwise static
  * PYTHON_MODULE - enable/disable building of the python bindigs, pybind11 is required to build the catima python module, default OFF
  * APPS - build command line app, default ON
  * TESTS - build tests, default OFF
  * EXAMPLES - build examples, default OFF
  * GENERATE_DATA - makes program to re-generate precalculated tables (ie precalculated LS coefficients), default:OFF
  * THIN_TARGET_APPROXIMATION - compile the library with thin target approximation, default: ON
  * GSL_INTEGRATION - use GSL integration functions, otherwise use built-in integrator, default: OFF
  * GLOBAL - compile with GLOBAL code (source not included at the moment, needs to be manually added to __global__ directory, default:OFF)
  * STORE_SPLINES - store splines in cache, if disabled datapoints are stored and splines are recreated, default ON

ie:
> cmake -DPYTHON_MODULE=ON ../
