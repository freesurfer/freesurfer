project(mri_extract_conditions)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_extract_conditions mri_extract_conditions.cpp)
target_link_libraries(mri_extract_conditions utils)

install(TARGETS mri_extract_conditions DESTINATION bin)
