project(mri_synthesize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_synthesize mri_synthesize.cpp)
target_link_libraries(mri_synthesize utils)

install(TARGETS mri_synthesize DESTINATION bin)
