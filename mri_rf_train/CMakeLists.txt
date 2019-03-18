project(mri_rf_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_rf_train mri_rf_train.cpp)
target_link_libraries(mri_rf_train utils)

install(TARGETS mri_rf_train DESTINATION bin)
