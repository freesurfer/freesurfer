project(mri_voldiff)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_voldiff mri_voldiff.c)
target_link_libraries(mri_voldiff utils)

install(TARGETS mri_voldiff DESTINATION bin)
