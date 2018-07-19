project(mri_gcab_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_gcab_train mri_gcab_train.c)
target_link_libraries(mri_gcab_train utils)

install(TARGETS mri_gcab_train DESTINATION bin)
