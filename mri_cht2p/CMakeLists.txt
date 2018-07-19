project(mri_cht2p)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_cht2p mri_cht2p.c)
target_link_libraries(mri_cht2p utils)

install(TARGETS mri_cht2p DESTINATION bin)
