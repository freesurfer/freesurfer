project(mri_gtmpvc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_gtmpvc mri_gtmpvc.c)
target_link_libraries(mri_gtmpvc utils)

install(TARGETS mri_gtmpvc DESTINATION bin)
