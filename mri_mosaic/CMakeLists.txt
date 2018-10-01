project(mri_mosaic)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mosaic mri_mosaic.c)
target_link_libraries(mri_mosaic utils)

install(TARGETS mri_mosaic DESTINATION bin)
