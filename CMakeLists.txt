cmake_minimum_required(VERSION 3.12.0)
project(catima)

############ options #############
option(BUILD_SHARED_LIBS "build as shared library" ON)
option(PYTHON_MODULE "compile the Catima python module(requires numpy and cython installed)" OFF)
option(TESTS "build tests" OFF)
option(EXAMPLES "build examples" OFF)
option(APPS "build catima applications" ON)
option(GLOBAL "build with global, sources are required" OFF)
option(REACTIONS "enable/disable nuclear reaction rate" ON)
option(STORE_SPLINES "store splines, if disables splines are always recreated" ON)
option(GSL_INTEGRATION "use GSL integration" OFF)
option(GSL_INTERPOLATION "use GSL inteRPOLATION" OFF)
option(THIN_TARGET_APPROXIMATION "thin target approximation" ON)
option(GENERATE_DATA "make data tables generator" OFF)
option(PYTHON_WHEEL "make python wheel" OFF)
######## build type ############
if(NOT CMAKE_BUILD_TYPE  STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    MESSAGE(STATUS "Build type Release")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wfatal-errors -Wno-unused-parameter -Wno-sign-compare")
    endif()
    MESSAGE(STATUS "Build type Debug")
endif()
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
################################

######### compiler flags ###########
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
MESSAGE(STATUS "install prefix: " ${CMAKE_INSTALL_PREFIX})

if(APPLE)
  set(RPATH_VARIABLE "DYLD_LIBRARY_PATH")
else()
  set(RPATH_VARIABLE "LD_LIBRARY_PATH")
endif()

############# Requirements ##################
if(GSL_INTEGRATION OR GSL_INTERPOLATION)
    find_package(GSL REQUIRED)
    MESSAGE(STATUS "GSL include dirs: " ${GSL_INCLUDE_DIRS})
    list(APPEND EXTRA_LIBS ${GSL_LIBRARIES} )
endif()

#find_package(nurex QUIET)
#if(nurex_FOUND)
#message(STATUS "nurex library found")
#set(NUREX ON)
#list(APPEND EXTRA_LIBS nurex::nurex)
#endif(nurex_FOUND)

find_package(fmt)
list(APPEND EXTRA_LIBS fmt::fmt)

configure_file(
                    "${CMAKE_CURRENT_SOURCE_DIR}/build_config.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/include/catima/build_config.h"
                    )

configure_file("${PROJECT_SOURCE_DIR}/init.sh.in"
               "${PROJECT_BINARY_DIR}/init.sh"
              )
############### main build ###########################

file(GLOB SOURCES *.cpp)
if(GLOBAL)
file(GLOB GLOBAL_SOURCES global/*.c)
LIST (APPEND SOURCES ${GLOBAL_SOURCES})
endif(GLOBAL)
file(GLOB HEADERS *.h libs/*.h)

add_library(catima ${SOURCES})
set_target_properties(catima PROPERTIES
		      POSITION_INDEPENDENT_CODE ON
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
                      )

target_link_libraries(catima PUBLIC ${EXTRA_LIBS})

target_include_directories(catima
                           PUBLIC  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
                                   $<BUILD_INTERFACE:${GSL_INCLUDE_DIRS}>
                          )

add_library(catima::catima ALIAS catima)

FILE(COPY ${HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include/catima)

# the compiler used for C++ files
MESSAGE( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )

######## for python module
find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
    message(STATUS "Python found: ${Python_EXECUTABLE}")
endif()
if(PYTHON_MODULE)
    if(NOT Python_FOUND)
        MESSAGE(SEND_ERROR "Python is required to build nurex python modules")
    endif(NOT Python_FOUND)
    find_package(pybind11 QUIET)
    if(NOT pybind11_FOUND)
	    message(INFO "pybind11 not found, trying to dowload")
	    include(FetchContent)
	    FetchContent_Declare(
		    pybind11
		    GIT_REPOSITORY https://github.com/pybind/pybind11.git
		    GIT_TAG v2.6.2
		    )
	    FetchContent_MakeAvailable(pybind11)
    endif(NOT pybind11_FOUND)
    
    #set(PYBIND11_CPP_STANDARD -std=c++14)
    pybind11_add_module(pycatima pymodule/pycatima.cpp)
    target_include_directories(pycatima PUBLIC
                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs>
                                    $<INSTALL_INTERFACE:include>)
    target_link_libraries(pycatima PRIVATE catima)
endif(PYTHON_MODULE )
if(PYTHON_WHEEL)
	execute_process(COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/pymodule/setup.py bdist_wheel)
endif(PYTHON_WHEEL)

########## Sub Directories ###########
if(EXAMPLES)
    file(GLOB EXAMPLES examples/*.cpp)
    FILE(COPY ${EXAMPLES} DESTINATION ${PROJECT_BINARY_DIR}/examples)
    FILE(COPY examples/makefile DESTINATION ${PROJECT_BINARY_DIR}/examples)
    #add_subdirectory("examples")
endif(EXAMPLES)
if(TESTS)
enable_testing()
add_subdirectory("tests")
endif(TESTS)

########## data generator ########
if(GENERATE_DATA)
add_executable(generate_ls_coeff utils/generator.cpp)
target_link_libraries(generate_ls_coeff catima)
#add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR/include/generated_LS_coeff.h}
#    COMMAND
#)
endif(GENERATE_DATA)

###### subdirectories ######
if(APPS)
add_subdirectory("bin")
endif(APPS)

####### install part #######
FILE(GLOB headers "*.h")
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(catimaConfigVersion.cmake VERSION 1.7 COMPATIBILITY AnyNewerVersion)

install (TARGETS catima
        EXPORT catimaConfig
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (FILES ${headers} DESTINATION include/catima)

install(EXPORT catimaConfig
        NAMESPACE catima::
        DESTINATION lib/cmake/catima
         )

export(TARGETS catima NAMESPACE catima:: FILE catimaConfig.cmake)
export(PACKAGE catima)

###### packaging #######
set(CPACK_PACKAGE_NAME "catima")
set(CPACK_PACKAGE_VENDOR "A. Prochazka")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "catima")
include(CPack)
