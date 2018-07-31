project(tkregister2_cmdl)

add_definitions(-DNO_GUI)

include_directories(${FS_INCLUDE_DIRS})

add_executable(tkregister2_cmdl ../tkregister2.c)
target_link_libraries(tkregister2_cmdl utils)
install(TARGETS tkregister2_cmdl DESTINATION tktools)
