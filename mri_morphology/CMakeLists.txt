project(mri_morphology)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_morphology mri_morphology.c)
target_link_libraries(mri_morphology utils)

install(TARGETS mri_morphology DESTINATION bin)
