project(mri_topologycorrection)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_topologycorrection mri_topologycorrection.c)
target_link_libraries(mri_topologycorrection utils)

install(TARGETS mri_topologycorrection DESTINATION bin)
