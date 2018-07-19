project(mri_path2label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_path2label mri_path2label.c)
target_link_libraries(mri_path2label utils)

install(TARGETS mri_path2label DESTINATION bin)
