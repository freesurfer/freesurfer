project(mri_polv)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_polv mri_polv.c)
target_link_libraries(mri_polv utils)

install(TARGETS mri_polv DESTINATION bin)
