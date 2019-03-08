project(mri_classify)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_classify mri_classify.cpp)
target_link_libraries(mri_classify utils)

install(TARGETS mri_classify DESTINATION bin)
