project(mri_aseg_edit_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_aseg_edit_train mri_aseg_edit_train.cpp)
target_link_libraries(mri_aseg_edit_train utils svm)

install(TARGETS mri_aseg_edit_train DESTINATION bin)
