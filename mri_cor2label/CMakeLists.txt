project(mri_cor2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_cor2label mri_cor2label.cpp)
target_link_libraries(mri_cor2label utils)
install(TARGETS mri_cor2label DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mri_cor2label DESTINATION bin RENAME mri_vol2label)
