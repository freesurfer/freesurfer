project(mri_mcsim)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mcsim mri_mcsim.cpp)
target_link_libraries(mri_mcsim utils)

install(TARGETS mri_mcsim DESTINATION bin)

add_executable(mri_maps2csd mri_maps2csd.cpp)
target_link_libraries(mri_maps2csd utils)

install(TARGETS mri_maps2csd DESTINATION bin)
