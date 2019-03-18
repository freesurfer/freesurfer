project(mri_gtmpvc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_gtmpvc mri_gtmpvc.cpp)
target_link_libraries(mri_gtmpvc utils)

install(TARGETS mri_gtmpvc DESTINATION bin)
