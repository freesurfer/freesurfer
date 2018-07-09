project(mri_aseg_edit_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_aseg_edit_train mri_aseg_edit_train.c)
target_link_libraries(mri_aseg_edit_train utils)

install(TARGETS mri_aseg_edit_train DESTINATION bin)
