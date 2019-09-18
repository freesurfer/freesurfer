project(mri_gcab_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_gcab_train mri_gcab_train.cpp)
target_link_libraries(mri_gcab_train utils)

install(TARGETS mri_gcab_train DESTINATION bin)
