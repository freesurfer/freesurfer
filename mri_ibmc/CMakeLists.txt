project(mri_ibmc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ibmc mri_ibmc.c)
target_link_libraries(mri_ibmc utils)

install(TARGETS mri_ibmc DESTINATION bin)
