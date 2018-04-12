set(CATIMA_APPS catima_calculator)

foreach(entry ${CATIMA_APPS})
    add_executable(${entry} ${entry}.cpp)
    target_link_libraries(${entry} catima)
endforeach(entry in ${CATIMA_APPS})

install (TARGETS ${CATIMA_APPS} RUNTIME DESTINATION bin)
