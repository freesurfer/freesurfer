project(mri_extract)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_extract mri_extract.cpp)
target_link_libraries(mri_extract utils)

install(TARGETS mri_extract DESTINATION bin)
