project(mri_volsynth)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_volsynth mri_volsynth.cpp)
target_link_libraries(mri_volsynth utils)

install(TARGETS mri_volsynth DESTINATION bin)
