project(mri_aparc2wmseg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_aparc2wmseg mri_aparc2wmseg.cpp)
target_link_libraries(mri_aparc2wmseg utils)

install(TARGETS mri_aparc2wmseg DESTINATION bin)
