project(mri_rf_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_rf_label mri_rf_label.c)
add_help(mri_rf_label mri_rf_label.help.xml)
target_link_libraries(mri_rf_label utils)

install(TARGETS mri_rf_label DESTINATION bin)
