project(mri_wmfilter)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_wmfilter mri_wmfilter.cpp)
target_link_libraries(mri_wmfilter utils)

install(TARGETS mri_wmfilter DESTINATION bin)
