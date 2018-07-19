project(mri_rf_long_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_rf_long_label mri_rf_long_label.c)
target_link_libraries(mri_rf_long_label utils)

install(TARGETS mri_rf_long_label DESTINATION bin)
