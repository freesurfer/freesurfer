project(mri_z2p)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_z2p mri_z2p.c)
target_link_libraries(mri_z2p utils)

install(TARGETS mri_z2p DESTINATION bin)
