[![Codacy Badge](https://api.codacy.com/project/badge/Grade/dc251db65f7a4c06ae07380544ea08fc)]()
[![Build](https://github.com/agoose77/catima/actions/workflows/build.yml/badge.svg)](https://github.com/agoose77/catima/actions/workflows/build.yml)
[![Documentation Status](https://readthedocs.org/projects/catima/badge/?version=latest)](https://catima.readthedocs.io/en/latest/?badge=latest)

CATima
=======
C++ library for calculation of energy loss, range, angular scattering and time of flight of the particle passing through matter.
The library is based on physics used in the ATIMA code,however its not a 100% copy of the ATIMA physics.
 see CREDITS for more details.

The WebAtima UI to this library can be found here:
  * https://web-docs.gsi.de/~aprochaz/webatima (only inside GSI)
  * https://isotopea.com/webatima

Installation
------------
CMake is used to build the library. For a default build, first prepare the build location

```bash
mkdir build
cd build
```

## System-wide install
Simply configure, build, and install the library
```bash
cmake ../
make
make install
```

## Local install
First, configure, build, and install the library into a prefix; replacing `<PATH-TO-INSTALL>` with the install location.
```bash
cmake ../ -DCMAKE_INSTALL_PREFIX=<PATH-TO-INSTALL>
make
make install
```

Then, update the `PATH` and `LD_LIBRARY_PATH` variables:
```bash
PATH="<PATH-TO-INSTALL>/bin:$PATH"
LD_LIBRARY_PATH="<PATH-TO-INSTALL>/lib:$LD_LIBRARY_PATH"
```


Python Module
-------------
Python module can be installed directly from PyPI using pip:
```
pip install pycatima
```

Or installed from the local Git checkout with
```
pip install ./pymodule
```


cmake options
-------------
compile options, enable or disable with cmake:
```bash
cmake ../ -D[OPTION]
```

available options:
  * BUILD_SHARED_LIBS - if ON shared library is build, otherwise static
  * APPS - build command line app, default ON
  * TESTS - build tests, default OFF
  * EXAMPLES - build examples, default OFF
  * DOCS - prepare doxygen documentation (after cmake, __make docs__ needs to be executed)
  * GENERATE_DATA - makes program to re-generate precalculated tables (ie precalculated LS coefficients), default:OFF
  * THIN_TARGET_APPROXIMATION - compile the library with thin target approximation, default: ON
  * GSL_INTEGRATION - use GSL integration functions, otherwise use built-in integrator, default: OFF
  * GLOBAL - compile with GLOBAL code (source not included at the moment, needs to be manually added to __global__ directory, default:OFF)
  * STORE_SPLINES - store splines in cache, if disabled datapoints are stored and splines are recreated, default ON

ie:
```bash
cmake -DEXAMPLES=ON ../
```

