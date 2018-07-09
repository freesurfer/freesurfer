project(mri_volsynth)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_volsynth mri_volsynth.c)
target_link_libraries(mri_volsynth utils)

install(TARGETS mri_volsynth DESTINATION bin)
