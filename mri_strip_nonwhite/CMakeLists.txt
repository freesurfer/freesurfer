project(mri_strip_nonwhite)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_strip_nonwhite mri_strip_nonwhite.cpp)
target_link_libraries(mri_strip_nonwhite utils)

install(TARGETS mri_strip_nonwhite DESTINATION bin)
