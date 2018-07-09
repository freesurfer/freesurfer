project(mri_probe_ima)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_probe_ima mri_probe_ima.c)
target_link_libraries(mri_probe_ima utils)

install(TARGETS mri_probe_ima DESTINATION bin)
