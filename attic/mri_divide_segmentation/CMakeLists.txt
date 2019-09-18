project(mri_divide_segmentation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_divide_segmentation mri_divide_segmentation.cpp)
target_link_libraries(mri_divide_segmentation utils)

install(TARGETS mri_divide_segmentation DESTINATION bin)
