project(mri_ca_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ca_train mri_ca_train.c)
target_link_libraries(mri_ca_train utils)

install(TARGETS mri_ca_train DESTINATION bin)
