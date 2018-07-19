project(mri_correct_segmentations)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_correct_segmentations mri_correct_segmentations.c)
target_link_libraries(mri_correct_segmentations utils)

install(TARGETS mri_correct_segmentations DESTINATION bin)
