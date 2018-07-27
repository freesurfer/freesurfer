project(mri_map_cpdat)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_map_cpdat mri_map_cpdat.cpp)
target_link_libraries(mri_map_cpdat utils)

install(TARGETS mri_map_cpdat DESTINATION bin)
