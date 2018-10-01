project(mri_surfacemask)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surfacemask mri_surfacemask.cpp)
target_link_libraries(mri_surfacemask utils)

install(TARGETS mri_surfacemask DESTINATION bin)
