project(mri_xvolavg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_xvolavg mri_xvolavg.cpp)
target_link_libraries(mri_xvolavg utils)

install(TARGETS mri_xvolavg DESTINATION bin)
