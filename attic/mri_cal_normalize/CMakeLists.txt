project(mri_cal_normalize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_cal_normalize mri_cal_normalize.cpp)
target_link_libraries(mri_cal_normalize utils)
install(TARGETS mri_cal_normalize DESTINATION bin)
