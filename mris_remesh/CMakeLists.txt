project(mris_remesh)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mris_remesh mris_remesh.cpp remesher.cpp)
add_help(mris_remesh mris_remesh.help.xml)
target_link_libraries(mris_remesh utils)

install(TARGETS mris_remesh DESTINATION bin)
