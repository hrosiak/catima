# BSD 3-Clause License; see https://github.com/scikit-hep/awkward-1.0/blob/main/LICENSE

cmake_minimum_required(VERSION 3.15...3.24)

# Project must be near the top
project(
  ${SKBUILD_PROJECT_NAME}
  LANGUAGES CXX
  VERSION ${SKBUILD_PROJECT_VERSION})

message(STATUS "CMake version ${CMAKE_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Defaults for properties in this directory (and below)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build a static library
set(BUILD_SHARED_LIBS OFF)
# TODO: remove this once scikit-build-core supports changing the source-dir
set(PYTHON_NO_INCLUDE_PYMODULE "")

# Hide all unused symbols
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# Add main library
set(APPS OFF)
add_subdirectory(../ catima EXCLUDE_FROM_ALL)

# Build extension
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_ext src/_ext.cpp)
target_link_libraries(_ext PRIVATE catima)

# Install python bindings
install(
  TARGETS _ext
  LIBRARY DESTINATION "${SKBUILD_PROJECT_NAME}"
  RUNTIME DESTINATION "${SKBUILD_PROJECT_NAME}"
  ARCHIVE DESTINATION "${SKBUILD_PROJECT_NAME}")
