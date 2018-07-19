project(mri_vol2roi)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_vol2roi mri_vol2roi.c)
target_link_libraries(mri_vol2roi utils)

install(TARGETS mri_vol2roi DESTINATION bin)
