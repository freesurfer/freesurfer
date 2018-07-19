project(mri_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_train mri_train.c)
target_link_libraries(mri_train utils)

install(TARGETS mri_train DESTINATION bin)
