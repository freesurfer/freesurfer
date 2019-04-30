project(mri_auto_fill)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_auto_fill mri_auto_fill.cpp)
target_link_libraries(mri_auto_fill utils)

install(TARGETS mri_auto_fill DESTINATION bin)
