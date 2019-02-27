project(mri_extract_fcd_features)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_extract_fcd_features mri_extract_fcd_features.cpp)
target_link_libraries(mri_extract_fcd_features utils)

install(TARGETS mri_extract_fcd_features DESTINATION bin)
