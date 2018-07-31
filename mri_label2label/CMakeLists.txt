project(mri_label2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label2label mri_label2label.c)
target_link_libraries(mri_label2label utils)

install(TARGETS mri_label2label DESTINATION bin)
