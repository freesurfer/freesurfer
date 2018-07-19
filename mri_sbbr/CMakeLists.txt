project(mri_sbbr)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_sbbr mri_sbbr.c)
target_link_libraries(mri_sbbr utils)

install(TARGETS mri_sbbr DESTINATION bin)
