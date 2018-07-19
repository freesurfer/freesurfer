project(mris_add_template)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_add_template mris_add_template.c)
target_link_libraries(mris_add_template utils)

install(TARGETS mris_add_template DESTINATION bin)
