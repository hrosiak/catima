set(CATIMA_APPS catima_calculator)
include(ExternalProject)
ExternalProject_Add(
	json_hpp
	URL https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp
	DOWNLOAD_NO_EXTRACT 1
	DOWNLOAD_NO_PROGRESS 
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	PREFIX "external"
	)
ExternalProject_Get_Property(json_hpp download_dir)

foreach(entry ${CATIMA_APPS})
    add_executable(${entry} ${entry}.cpp)
    target_link_libraries(${entry} catima)
    target_include_directories(${entry} PRIVATE ${download_dir})
endforeach(entry in ${CATIMA_APPS})

install (TARGETS ${CATIMA_APPS} RUNTIME DESTINATION bin)
