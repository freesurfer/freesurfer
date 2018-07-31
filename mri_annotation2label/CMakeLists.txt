project(mri_annotation2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_annotation2label mri_annotation2label.c)
target_link_libraries(mri_annotation2label utils)

install(TARGETS mri_annotation2label DESTINATION bin)
