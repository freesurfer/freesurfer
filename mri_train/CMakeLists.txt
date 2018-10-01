project(mri_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_train mri_train.c)
target_link_libraries(mri_train utils)

install(TARGETS mri_train DESTINATION bin)
