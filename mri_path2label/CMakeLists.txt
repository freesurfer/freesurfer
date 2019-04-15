project(mri_path2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_path2label mri_path2label.cpp)
target_link_libraries(mri_path2label utils)

install(TARGETS mri_path2label DESTINATION bin)
