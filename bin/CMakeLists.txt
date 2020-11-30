set(CATIMA_APPS catima_calculator)

find_package(nlohmann_json QUIET)

if(nlohmann_json_FOUND)
MESSAGE(STATUS "system nlohmann_json found ")
else(nlohmann_json_FOUND)
#MESSAGE(STATUS "nlohmann_json looking ")
include(ExternalProject)
ExternalProject_Add(
	json_hpp
	URL https://github.com/nlohmann/json/releases/download/v3.9.1/include.zip	
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	PREFIX "external"
	)
ExternalProject_Get_Property(json_hpp SOURCE_DIR)
endif(nlohmann_json_FOUND)

foreach(entry ${CATIMA_APPS})
    add_executable(${entry} ${entry}.cpp)
    target_link_libraries(${entry} catima )
	#target_link_libraries(${entry} PRIVATE nlohmann_json::nlohmann_json)
	target_include_directories(${entry} PRIVATE ${SOURCE_DIR}/include) 
endforeach(entry in ${CATIMA_APPS})

install (TARGETS ${CATIMA_APPS} RUNTIME DESTINATION bin)
