project(mri_fdr)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fdr mri_fdr.cpp)
target_link_libraries(mri_fdr utils)

install(TARGETS mri_fdr DESTINATION bin)
