project(mri_xcorr)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_xcorr mri_xcorr.c)
target_link_libraries(mri_xcorr utils)

install(TARGETS mri_xcorr DESTINATION bin)
