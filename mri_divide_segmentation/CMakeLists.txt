project(mri_divide_segmentation)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_divide_segmentation mri_divide_segmentation.c)
target_link_libraries(mri_divide_segmentation utils)

install(TARGETS mri_divide_segmentation DESTINATION bin)
