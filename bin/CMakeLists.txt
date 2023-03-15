include(FetchContent)

set(CATIMA_APPS catima_calculator)

# Don't install nlohmann_json for private libraries
set(JSON_Install OFF CACHE INTERNAL "")
FetchContent_Declare(nlohmann_json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(nlohmann_json)

foreach(entry ${CATIMA_APPS})
    add_executable(${entry} ${entry}.cpp)
    target_link_libraries(${entry} PRIVATE catima nlohmann_json::nlohmann_json)
endforeach(entry in ${CATIMA_APPS})

install (TARGETS ${CATIMA_APPS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
