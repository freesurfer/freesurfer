project(mris_reposition_surface)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_reposition_surface mris_reposition_surface.cpp)
add_help(mris_reposition_surface mris_reposition_surface.help.xml)
target_link_libraries(mris_reposition_surface utils)
install(TARGETS mris_reposition_surface DESTINATION bin) 
