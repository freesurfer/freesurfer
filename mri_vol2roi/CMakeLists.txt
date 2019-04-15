project(mri_vol2roi)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_vol2roi mri_vol2roi.cpp)
target_link_libraries(mri_vol2roi utils)

install(TARGETS mri_vol2roi DESTINATION bin)
