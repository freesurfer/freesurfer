project(mri_morphology)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_morphology mri_morphology.c)
target_link_libraries(mri_morphology utils)

install(TARGETS mri_morphology DESTINATION bin)
