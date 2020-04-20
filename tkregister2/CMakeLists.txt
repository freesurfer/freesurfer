project(tkregister2)

include_directories(${FS_INCLUDE_DIRS})

add_executable(tkregister2_cmdl tkregister2.cpp)
target_compile_definitions(tkregister2_cmdl PRIVATE -DNO_GUI)
target_link_libraries(tkregister2_cmdl utils)
install(TARGETS tkregister2_cmdl DESTINATION bin)
symlink(tkregister2_cmdl ${CMAKE_INSTALL_PREFIX}/bin/tkregister2)
