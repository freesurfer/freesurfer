project(mri_gtmseg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_gtmseg mri_gtmseg.cpp)
target_link_libraries(mri_gtmseg utils)

install(TARGETS mri_gtmseg DESTINATION bin)
