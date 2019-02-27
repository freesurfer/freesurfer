project(mri_apply_inu_correction)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_apply_inu_correction mri_apply_inu_correction.cpp)
target_link_libraries(mri_apply_inu_correction utils)

install(TARGETS mri_apply_inu_correction DESTINATION bin)
