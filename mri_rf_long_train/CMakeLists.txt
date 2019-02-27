project(mri_rf_long_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_rf_long_train mri_rf_long_train.cpp)
target_link_libraries(mri_rf_long_train utils)

install(TARGETS mri_rf_long_train DESTINATION bin)
