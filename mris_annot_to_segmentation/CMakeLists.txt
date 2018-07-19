project(mris_annot_to_segmentation)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_annot_to_segmentation mris_annot_to_segmentation.c)
target_link_libraries(mris_annot_to_segmentation utils)

install(TARGETS mris_annot_to_segmentation DESTINATION bin)
