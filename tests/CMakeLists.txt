enable_testing()
add_definitions(-Dlest_FEATURE_COLOURISE=1)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests)

set(CATIMA_TESTS test_calculations test_generated test_storage test_structures test_dedx_range test_abundances)

foreach(entry ${CATIMA_TESTS})
    add_executable(${entry} ${entry}.cpp)
    target_link_libraries(${entry} catima)
    add_test(${entry} ${PROJECT_BINARY_DIR}/tests/${entry})
endforeach(entry in ${CATIMA_TESTS})

if(NUREX)
    add_executable(test_reaction test_reaction.cpp)
    target_link_libraries(test_reaction nurex::nurex catima)
    add_test(test_reaction ${PROJECT_BINARY_DIR}/tests/test_reaction)
endif(NUREX)

set(C_TESTS test_c)
MESSAGE( STATUS "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER} )
foreach(entry ${C_TESTS})
    add_executable(${entry} ${entry}.c)
    target_link_libraries(${entry} catima)
endforeach(entry in ${C_TESTS})
