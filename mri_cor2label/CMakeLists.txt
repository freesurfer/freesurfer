project(mri_cor2label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_cor2label mri_cor2label.c)
target_link_libraries(mri_cor2label utils)

install(TARGETS mri_cor2label DESTINATION bin)
