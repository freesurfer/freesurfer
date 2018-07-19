project(mri_gtmseg)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_gtmseg mri_gtmseg.c)
target_link_libraries(mri_gtmseg utils)

install(TARGETS mri_gtmseg DESTINATION bin)
