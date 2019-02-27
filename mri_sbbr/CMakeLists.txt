project(mri_sbbr)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_sbbr mri_sbbr.cpp)
target_link_libraries(mri_sbbr utils)

install(TARGETS mri_sbbr DESTINATION bin)
