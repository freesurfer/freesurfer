project(mris_talairach)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_talairach mris_talairach.c)
target_link_libraries(mris_talairach utils)

install(TARGETS mris_talairach DESTINATION bin)
