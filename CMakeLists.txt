cmake_minimum_required(VERSION 3.14...3.25)

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
option(GLOBAL "build with global, sources are required" OFF)
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
################################

######### compiler flags ###########
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


############# Requirements ##################
if (GSL_INTEGRATION OR GSL_INTERPOLATION)
    find_package(GSL REQUIRED)
    list(APPEND EXTRA_LIBS GSL::gsl)
endif ()

if (GLOBAL)
    add_subdirectory("extern/global")
endif (GLOBAL)


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
               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

# Exported alias
add_library(${LOCAL_PROJECT_NAMESPACE}::${LOCAL_PROJECT_NAME} ALIAS catima)


########## Tests ###########
if (TESTS)
    enable_testing()
    add_subdirectory("tests")
endif (TESTS)

########## Data generator ########
if (GENERATE_DATA)
    add_executable(generate_ls_coeff utils/generator.cpp)
    target_link_libraries(generate_ls_coeff catima)
    #add_custom_command(
    #    OUTPUT ${CMAKE_CURRENT_BINARY_DIR/include/generated_LS_coeff.h}
    #    COMMAND
    #)
endif (GENERATE_DATA)

###### Apps ######
if (APPS)
    add_subdirectory("bin")
endif (APPS)

####### install part #######
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY AnyNewerVersion
)

include(GNUInstallDirs)

# Generate targets 
install(TARGETS catima
        EXPORT ${LOCAL_PROJECT_NAME}Targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install header files (TODO, this can be replaced with FILE_SET in cmake 3.23)
install(DIRECTORY include/catima "${CMAKE_CURRENT_BINARY_DIR}/generated/catima" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install examples
if (EXAMPLES)
    install(DIRECTORY examples DESTINATION ${CMAKE_INSTALL_DATADIR}/catima)
endif (EXAMPLES)

# Generate CMake file into the install tree (rather than export(TARGETS) which installs into build tree)
# Directly generate `Config` instead of `Targets`, as we have no additional configuration
install(EXPORT ${LOCAL_PROJECT_NAME}Targets 
        NAMESPACE ${LOCAL_PROJECT_NAMESPACE}:: 
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${LOCAL_PROJECT_NAME}
        FILE ${LOCAL_PROJECT_NAME}Config.cmake)
# Install version
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${LOCAL_PROJECT_NAME}
)
        
# Packaging
set(CPACK_PACKAGE_VENDOR ${LOCAL_PROJECT_VENDOR})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${LOCAL_PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
