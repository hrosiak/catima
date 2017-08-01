CATima
=======
C++ library for caclulalaton of energy loss, range, angular scattering and time of flight of the particle passing through matter.
The library is based on physics used in the ATIMA code,however its not 100% copy of ATIMA physics.
 see CREDITS for more details.


Installation
------------
CMake is used to build the library. For default build use:

```
> mkdir build
> cd build
> cmake ../
> make
```

cmake options
-------------
compile options, enable or disable with cmake:
> cmake ../ -D[OPTION]

available options:
  * CATIMA_PYTHON - enable/disable building of the python bindigs, cython and numpy are required to build the catima python module, default OFF
  * TESTS - build tests
  * EXAMPLES - build examples
  * DOCS - prepare doxygen documentation (after cmake, __make docs__ needs to be executed)
  * GENERATE_DATA - makes program to re-generate precalculated tables (ie precalculated LS coefficients), default:OFF
  * THIN_TARGET_APPROXIMATION - compile the library with thin target approximation, default: ON

ie:
> cmake -DCATIMA_PYTHON=ON -DEXAMPLES=ON ../


after the compilation the libraries and headers must be either installed system-wide by make install or PATH and LD_LIBRARY_PATH must be adjusted to point to headers and library files.
The default install path can be change, ie: cmake -DCMAKE_INSTALL_PREFIX=/opt/catima

the option to system-wide installation is to adjust library path and include paths.
This can be done sourcing the init.sh file, which is generated in the build directory:
> source init.sh
