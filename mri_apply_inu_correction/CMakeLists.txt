project(mri_apply_inu_correction)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_apply_inu_correction mri_apply_inu_correction.c)
target_link_libraries(mri_apply_inu_correction utils)

install(TARGETS mri_apply_inu_correction DESTINATION bin)
