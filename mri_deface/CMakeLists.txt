project(mri_deface)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_deface mri_deface.c)
target_link_libraries(mri_deface utils)

install(TARGETS mri_deface DESTINATION bin)
