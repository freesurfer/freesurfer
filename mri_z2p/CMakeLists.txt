project(mri_z2p)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_z2p mri_z2p.cpp)
target_link_libraries(mri_z2p utils)

install(TARGETS mri_z2p DESTINATION bin)
