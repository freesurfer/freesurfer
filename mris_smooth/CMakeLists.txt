project(mris_smooth)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_smooth mris_smooth.c)
add_help(mris_smooth mris_smooth.help.xml)
target_link_libraries(mris_smooth utils)

install(TARGETS mris_smooth DESTINATION bin)
