project(mris_make_template)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_make_template mris_make_template.c)
target_link_libraries(mris_make_template utils)

install(TARGETS mris_make_template DESTINATION bin)
