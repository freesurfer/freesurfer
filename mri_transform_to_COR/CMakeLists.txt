project(mri_transform_to_COR)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_transform_to_COR mri_transform_to_COR.c)
target_link_libraries(mri_transform_to_COR utils)

install(TARGETS mri_transform_to_COR DESTINATION bin)
