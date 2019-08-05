project(mris_annot_diff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_annot_diff mris_annot_diff.cpp)
add_help(mris_annot_diff mris_annot_diff.help.xml)
target_link_libraries(mris_annot_diff utils)
install(TARGETS mris_annot_diff DESTINATION bin) 
