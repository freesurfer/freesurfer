project(mri_rf_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_rf_train mri_rf_train.c)
target_link_libraries(mri_rf_train utils)

install(TARGETS mri_rf_train DESTINATION bin)
