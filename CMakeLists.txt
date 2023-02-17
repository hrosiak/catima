cmake_minimum_required(VERSION 3.14...3.25)
project(catima VERSION 1.0
               DESCRIPTION "C++ library for calculation of energy loss, range, angular scattering and time of flight of the particle passing through matter. "
               LANGUAGES CXX)

############ options #############
option(BUILD_SHARED_LIBS "build as shared library" ON)
option(TESTS "build tests" OFF)
option(EXAMPLES "build examples" OFF)
option(APPS "build catima applications" ON)
option(REACTIONS "enable/disable nuclear reaction rate" ON)
option(STORE_SPLINES "store splines, if disables splines are always recreated" ON)
option(GSL_INTEGRATION "use GSL integration" OFF)
option(GSL_INTERPOLATION "use GSL interpolation" OFF)
option(THIN_TARGET_APPROXIMATION "thin target approximation" ON)
option(ET_CALCULATED_INDEX "calculate energy table index, otherwise search" ON)
option(GENERATE_DATA "make data tables generator" OFF)
######## build type ############
if(NOT CMAKE_BUILD_TYPE  STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
    MESSAGE(STATUS "Build type Release")
else()
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wfatal-errors -Wno-unused-parameter -Wno-sign-compare")
    endif()
    MESSAGE(STATUS "Build type Debug")
endif()
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
################################

######### compiler flags ###########
set(CMAKE_CXX_EXTENSIONS OFF)
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


configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/build_config.in"
                "${CMAKE_CURRENT_BINARY_DIR}/generated/catima/build_config.h")

configure_file("${PROJECT_SOURCE_DIR}/init.sh.in"
               "${PROJECT_BINARY_DIR}/init.sh")
############### main build ###########################


add_library(catima src/abundance_database.cpp
                   src/calculations.cpp
                   src/catima.cpp
                   src/cwrapper.cpp
                   src/integrator.cpp
                   src/material_database.cpp
                   src/nucdata.cpp
                   src/reactions.cpp
                   src/spline.cpp
                   src/srim.cpp
                   src/storage.cpp
                   src/structures.cpp)
set_target_properties(catima PROPERTIES
		    POSITION_INDEPENDENT_CODE ON
            LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
            ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
                      )

target_link_libraries(catima PUBLIC ${EXTRA_LIBS})
target_compile_features(catima PRIVATE cxx_std_17)
target_include_directories(catima
                           PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
                                   $<BUILD_INTERFACE:${GSL_INCLUDE_DIRS}>
                          )

add_library(catima::catima ALIAS catima)

# the compiler used for C++ files
MESSAGE( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )


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
FILE(GLOB HEADERS "include/**/*.h")
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(catimaConfigVersion.cmake VERSION 1.7 COMPATIBILITY AnyNewerVersion)

install (TARGETS catima
        EXPORT catimaConfig
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (FILES ${HEADERS} DESTINATION include/catima)

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
