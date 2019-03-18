project(mri_ca_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ca_train mri_ca_train.cpp)
target_link_libraries(mri_ca_train utils)

install(TARGETS mri_ca_train DESTINATION bin)
