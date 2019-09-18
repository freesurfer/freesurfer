project(mri_estimate_tissue_parms)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_estimate_tissue_parms mri_estimate_tissue_parms.cpp)
target_link_libraries(mri_estimate_tissue_parms utils)

install(TARGETS mri_estimate_tissue_parms DESTINATION bin)
