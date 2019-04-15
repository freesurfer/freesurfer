project(mri_topologycorrection)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_topologycorrection mri_topologycorrection.cpp)
target_link_libraries(mri_topologycorrection utils)

install(TARGETS mri_topologycorrection DESTINATION bin)
