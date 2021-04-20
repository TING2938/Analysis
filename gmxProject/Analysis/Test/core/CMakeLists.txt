

add_executable(${sourcefile} ${sourcefile})

target_link_libraries(${sourcefile} gtest_main)

install(TARGETS ${sourcefile} DESTINATION bin)

