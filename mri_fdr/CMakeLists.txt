project(mri_fdr)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fdr mri_fdr.c)
target_link_libraries(mri_fdr utils)

install(TARGETS mri_fdr DESTINATION bin)
