cmake_minimum_required(VERSION 3.2.0)
project(catima)

############ options #############
#option(THREADS "Use multi-threading" ON)
option(CATIMA_PYTHON "compile the Catima python module(requires numpy and cython installed)" OFF)
option(TESTS "build tests" OFF)
option(EXAMPLES "build examples" ON)
option(GSL_INTEGRATION "use GSL integration" ON)
option(GENERATE_DATA "make data tables generator" OFF)
option(THIN_TARGET_APPROXIMATION "thin target approximation" ON)
option(DOCS "build documentation (requires doxygen)" OFF)
option(GLOBAL "build with global, sources are required" OFF)
option(APPS "build catima applications" ON)

######## build type ############
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
################################

######### compiler flags ###########
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(warnings "-Wall -Wextra -Werror")
endif()

if(APPLE)
  set(RPATH_VARIABLE "DYLD_LIBRARY_PATH")
else()
  set(RPATH_VARIABLE "LD_LIBRARY_PATH")
endif()

MESSAGE(STATUS "install prefix: " ${CMAKE_INSTALL_PREFIX})

############# Requirements ##################
find_package(GSL REQUIRED)
MESSAGE(STATUS "GSL include dirs: " ${GSL_INCLUDE_DIRS})
#if(THREADS)
#    find_package(Threads REQUIRED)
#    set (EXTRA_LIBS ${EXTRA_LIBS} ${CMAKE_THREAD_LIBS_INIT})
#    set (USE_THREADS ON)
#    MESSAGE(STATUS "Nurex will use threads")
#endif(THREADS)

find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
    message("-- Python found: ${PYTHON_EXECUTABLE}")
endif()

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
file(GLOB HEADERS *.h)

add_library(catima SHARED ${SOURCES})
add_library(catima_static STATIC ${SOURCES})
set_target_properties(catima_static PROPERTIES OUTPUT_NAME catima POSITION_INDEPENDENT_CODE ON)
target_link_libraries(catima ${EXTRA_LIBS} ${GSL_LIBRARIES} )
target_link_libraries(catima_static ${EXTRA_LIBS} ${GSL_LIBRARIES})

target_include_directories(catima 
                           PUBLIC  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
                                   $<BUILD_INTERFACE:${GSL_INCLUDE_DIRS}>
                          )
target_include_directories(catima_static
                           PUBLIC  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
                                   $<BUILD_INTERFACE:${GSL_INCLUDE_DIRS}>
                          )
add_library(catima::catima ALIAS catima)
add_library(catima::catima_static ALIAS catima_static)

FILE(COPY ${HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include/catima)

# the compiler used for C++ files 
MESSAGE( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )
# the compiler flags for compiling C++ sources 
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS_RELEASE} )

######## for python module
if(CATIMA_PYTHON)
    if(NOT PYTHONINTERP_FOUND)
        MESSAGE(SEND_ERROR "Python is required to build nurex python modules") 
    endif(NOT PYTHONINTERP_FOUND)
    find_program(CYTHON_EXECUTABLE 
                 NAMES cython cython2 cython3 cython.bat
                 DOC "path to the cython executable"
                )
    if(NOT CYTHON_EXECUTABLE)
        MESSAGE(SEND_ERROR "Cython not found, it is required to build nurex python modules") 
    endif(NOT CYTHON_EXECUTABLE)
    MESSAGE(STATUS "Cython found: " ${CYTHON_EXECUTABLE})
    
     ### build libraries string
    foreach(entry ${EXTRA_LIBS} ${GSL_LIBRARIES} catima)
        LIST (APPEND EXTRA_PYTHON_LIBS \"${entry}\")
    endforeach(entry ${EXTRA_LIBS} ${GSL_LIBRARIES} catima)
    string (REPLACE ";" "," EXTRA_PYTHON_LIBS "${EXTRA_PYTHON_LIBS}")
#    if(THREADS)
#    set (CYTHON_DEFINES "-DUSE_THREADS=1")
#    endif(THREADS)
     ### insert libraries string and create setup.py
    FILE(COPY catimac.pxd catima.pyx DESTINATION ${PROJECT_BINARY_DIR})
    set(CATIMA_LIB ${CMAKE_SHARED_LIBRARY_PREFIX}catima${CMAKE_SHARED_LIBRARY_SUFFIX})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${PROJECT_BINARY_DIR}/setup.py)
    add_custom_target(target_python ALL DEPENDS catima)
    add_custom_command(TARGET target_python
                        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/setup.py  build_ext ${CYTHON_DEFINES} -i                         
                        )
endif(CATIMA_PYTHON )

########## Sub Directories ###########
if(EXAMPLES)
    file(GLOB EXAMPLES examples/*.cpp)
    FILE(COPY ${EXAMPLES} DESTINATION ${PROJECT_BINARY_DIR}/examples)
    FILE(COPY examples/makefile DESTINATION ${PROJECT_BINARY_DIR}/examples)
    #add_subdirectory("examples")
endif(EXAMPLES)
if(TESTS)
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

####### DOCS generation #########
if(DOCS)
find_package(Doxygen REQUIRED)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(docs 
                   COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile COMMENT "Generating API documentation with Doxygen" VERBATIM )
endif(DOXYGEN_FOUND)
endif(DOCS)

###### subdirectories ######
if(APPS)
add_subdirectory("bin")
endif(APPS)

####### install part #######
FILE(GLOB headers "*.h")
include(GNUInstallDirs)
install (TARGETS catima catima_static 
        EXPORT catimaConfig
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
        
install (FILES ${headers} DESTINATION include/catima)

install(EXPORT catimaConfig
        NAMESPACE catima::
        DESTINATION lib/cmake/catima
         )
         
export(TARGETS catima catima_static NAMESPACE catima:: FILE catimaConfig.cmake)
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
