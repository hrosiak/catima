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
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# Set RPATH
if(APPLE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH "@loader_path")
else()
  set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
endif()

# Add main library
set(APPS OFF)
add_subdirectory(../ catima EXCLUDE_FROM_ALL)

# Install python bindings
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_ext src/_ext.cpp)
target_link_libraries(_ext PRIVATE catima)

install(
  TARGETS _ext catima
  LIBRARY DESTINATION "${SKBUILD_PROJECT_NAME}/lib"
  RUNTIME DESTINATION "${SKBUILD_PROJECT_NAME}/lib"
  ARCHIVE DESTINATION "${SKBUILD_PROJECT_NAME}/lib")
