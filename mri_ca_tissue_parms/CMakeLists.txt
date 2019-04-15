project(mri_ca_tissue_parms)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ca_tissue_parms mri_ca_tissue_parms.cpp)
target_link_libraries(mri_ca_tissue_parms utils)

install(TARGETS mri_ca_tissue_parms DESTINATION bin)
