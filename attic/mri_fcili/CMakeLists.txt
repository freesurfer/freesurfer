project(mri_fcili)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fcili mri_fcili.cpp)
target_link_libraries(mri_fcili utils)

install(TARGETS mri_fcili DESTINATION bin)
