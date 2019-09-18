project(mri_xcorr)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_xcorr mri_xcorr.cpp)
target_link_libraries(mri_xcorr utils)

install(TARGETS mri_xcorr DESTINATION bin)
