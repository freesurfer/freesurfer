project(mri_label_histo)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label_histo mri_label_histo.cpp)
target_link_libraries(mri_label_histo utils)

install(TARGETS mri_label_histo DESTINATION bin)
