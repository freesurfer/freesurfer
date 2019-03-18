project(mri_rf_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_rf_label mri_rf_label.cpp)
add_help(mri_rf_label mri_rf_label.help.xml)
target_link_libraries(mri_rf_label utils)

install(TARGETS mri_rf_label DESTINATION bin)
