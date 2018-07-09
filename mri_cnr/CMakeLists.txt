project(mri_cnr)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_cnr mri_cnr.c)
target_link_libraries(mri_cnr utils)

install(TARGETS mri_cnr DESTINATION bin)
