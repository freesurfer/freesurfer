project(mri_rf_long_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_rf_long_train mri_rf_long_train.c)
target_link_libraries(mri_rf_long_train utils)

install(TARGETS mri_rf_long_train DESTINATION bin)
