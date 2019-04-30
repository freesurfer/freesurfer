project(mri_transform_to_COR)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_transform_to_COR mri_transform_to_COR.cpp)
target_link_libraries(mri_transform_to_COR utils)

install(TARGETS mri_transform_to_COR DESTINATION bin)
