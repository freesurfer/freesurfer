project(mri_label_accuracy)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label_accuracy mri_label_accuracy.cpp)
target_link_libraries(mri_label_accuracy utils)

install(TARGETS mri_label_accuracy DESTINATION bin)
