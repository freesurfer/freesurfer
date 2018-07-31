project(mri_otl)

include_directories(${FS_INCLUDE_DIRS})

add_executable(list_otl_labels list_otl_labels.c)
target_link_libraries(list_otl_labels utils)
install(TARGETS list_otl_labels DESTINATION bin)
