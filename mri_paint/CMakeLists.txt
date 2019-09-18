project(mri_paint)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_paint mri_paint.cpp)
target_link_libraries(mri_paint utils)

install(TARGETS mri_paint DESTINATION bin)
