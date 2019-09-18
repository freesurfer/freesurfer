project(mri_vol2surf)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_vol2surf mri_vol2surf.cpp)
target_link_libraries(mri_vol2surf utils)

install(TARGETS mri_vol2surf DESTINATION bin)
