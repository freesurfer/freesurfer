project(mri_cal_renormalize_gca)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_cal_renormalize_gca mri_cal_renormalize_gca.cpp)
add_help(mri_cal_renormalize_gca mri_cal_renormalize_gca.help.xml)
target_link_libraries(mri_cal_renormalize_gca utils)

install(TARGETS mri_cal_renormalize_gca DESTINATION bin)
