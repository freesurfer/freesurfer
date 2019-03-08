project(mri_rf_long_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_rf_long_label mri_rf_long_label.cpp)
target_link_libraries(mri_rf_long_label utils)

install(TARGETS mri_rf_long_label DESTINATION bin)
