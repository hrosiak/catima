cmake_minimum_required(VERSION 3.14...3.25)

# TODO: remove this once scikit-build-core supports changing the source-dir
# If we're actually building the Python extension, via this file
if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND DEFINED SKBUILD AND NOT DEFINED PYTHON_NO_INCLUDE_PYMODULE)
    message(STATUS "A build of pymodule was initiated from the root directory. Loading pymodule")
    add_subdirectory(pymodule)
    return()
endif ()

set(LOCAL_PROJECT_VERSION "1.0.0")
set(LOCAL_PROJECT_VENDOR "A. Prochazka")
set(LOCAL_PROJECT_NAMESPACE "catima")
set(LOCAL_PROJECT_NAME "catima")
set(LOCAL_PROJECT_OUTPUT_NAME "catima")
set(LOCAL_PROJECT_DESCRIPTION "C++ library for calculation of energy loss, range, angular scattering and time of flight of the particle passing through matter. ")


project(${LOCAL_PROJECT_NAME}
        VERSION ${LOCAL_PROJECT_VERSION}
        DESCRIPTION ${LOCAL_PROJECT_DESCRIPTION}
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
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
    MESSAGE(STATUS "Build type Release")
else ()
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wfatal-errors -Wno-unused-parameter -Wno-sign-compare")
    endif ()
    MESSAGE(STATUS "Build type Debug")
endif ()
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
################################

######### compiler flags ###########
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
MESSAGE(STATUS "install prefix: " ${CMAKE_INSTALL_PREFIX})


############# Requirements ##################
if (GSL_INTEGRATION OR GSL_INTERPOLATION)
    find_package(GSL REQUIRED)
    MESSAGE(STATUS "GSL include dirs: " ${GSL_INCLUDE_DIRS})
    list(APPEND EXTRA_LIBS ${GSL_LIBRARIES})
endif ()


configure_file("${CMAKE_CURRENT_SOURCE_DIR}/build_config.in"
        "${CMAKE_CURRENT_BINARY_DIR}/generated/catima/build_config.h")

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
set_target_properties(catima PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(catima PUBLIC ${EXTRA_LIBS})
target_compile_features(catima PRIVATE cxx_std_17)
target_include_directories(catima
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
        $<BUILD_INTERFACE:${GSL_INCLUDE_DIRS}>
        )
# Exported alias
add_library(${LOCAL_PROJECT_NAMESPACE}::${LOCAL_PROJECT_NAME} ALIAS catima)

# the compiler used for C++ files
MESSAGE(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER})


########## Tests ###########
if (TESTS)
    enable_testing()
    add_subdirectory("tests")
endif (TESTS)

########## data generator ########
if (GENERATE_DATA)
    add_executable(generate_ls_coeff utils/generator.cpp)
    target_link_libraries(generate_ls_coeff catima)
    #add_custom_command(
    #    OUTPUT ${CMAKE_CURRENT_BINARY_DIR/include/generated_LS_coeff.h}
    #    COMMAND
    #)
endif (GENERATE_DATA)

###### subdirectories ######
if (APPS)
    add_subdirectory("bin")
    find_package(nlohmann_json QUIET)
    if (NOT nlohmann_json_FOUND)
        set(JSON_BuildTests OFF CACHE INTERNAL "")
        add_subdirectory(extern/json)
        find_package(nlohmann_json QUIET)
    endif (NOT nlohmann_json_FOUND)

endif (APPS)

####### install part #######
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${LOCAL_PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY AnyNewerVersion
)

include(GNUInstallDirs)
install(TARGETS catima
        EXPORT ${LOCAL_PROJECT_NAME}Targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if (EXAMPLES)
    install(DIRECTORY examples DESTINATION ${CMAKE_INSTALL_DATADIR}/catima)
endif (EXAMPLES)

export(TARGETS catima NAMESPACE ${LOCAL_PROJECT_NAMESPACE}:: FILE ${LOCAL_PROJECT_NAME}Targets.cmake)
export(PACKAGE ${LOCAL_PROJECT_NAME})

# Packaging
set(CPACK_PACKAGE_VENDOR ${LOCAL_PROJECT_VENDOR})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${LOCAL_PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
