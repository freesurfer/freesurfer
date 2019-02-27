project(mri_correct_segmentations)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_correct_segmentations mri_correct_segmentations.cpp)
target_link_libraries(mri_correct_segmentations utils)

install(TARGETS mri_correct_segmentations DESTINATION bin)
