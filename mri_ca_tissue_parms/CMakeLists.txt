project(mri_ca_tissue_parms)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ca_tissue_parms mri_ca_tissue_parms.c)
target_link_libraries(mri_ca_tissue_parms utils)

install(TARGETS mri_ca_tissue_parms DESTINATION bin)
