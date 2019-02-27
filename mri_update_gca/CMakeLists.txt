project(mri_update_gca)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_update_gca mri_update_gca.cpp)
target_link_libraries(mri_update_gca utils)

install(TARGETS mri_update_gca DESTINATION bin)
