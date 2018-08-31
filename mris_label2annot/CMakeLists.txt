project(mris_label2annot)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_label2annot mris_label2annot.c)
target_link_libraries(mris_label2annot utils)

install(TARGETS mris_label2annot DESTINATION bin)
