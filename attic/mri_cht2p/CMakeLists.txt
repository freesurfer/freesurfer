project(mri_cht2p)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_cht2p mri_cht2p.cpp)
target_link_libraries(mri_cht2p utils)

install(TARGETS mri_cht2p DESTINATION bin)
