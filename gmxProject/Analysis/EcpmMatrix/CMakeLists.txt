
add_subdirectory(./src)

add_executable(${sourcefile} ${sourcefile})

target_link_libraries(${sourcefile} conp)

install(TARGETS ${sourcefile} DESTINATION bin)
