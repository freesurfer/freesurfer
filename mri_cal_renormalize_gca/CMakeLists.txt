project(mri_cal_renormalize_gca)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_cal_renormalize_gca mri_cal_renormalize_gca.c)
add_help(mri_cal_renormalize_gca mri_cal_renormalize_gca.help.xml)
target_link_libraries(mri_cal_renormalize_gca utils)

install(TARGETS mri_cal_renormalize_gca DESTINATION bin)
