project(mris_shrinkwrap)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_shrinkwrap mris_shrinkwrap.c)
target_link_libraries(mris_shrinkwrap utils)
install(TARGETS mris_shrinkwrap DESTINATION bin)

add_executable(mris_AA_shrinkwrap mris_AA_shrinkwrap.c)
target_link_libraries(mris_AA_shrinkwrap utils)
install(TARGETS mris_AA_shrinkwrap DESTINATION bin)
