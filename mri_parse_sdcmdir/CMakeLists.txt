project(mri_parse_sdcmdir)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_parse_sdcmdir mri_parse_sdcmdir.c)
target_link_libraries(mri_parse_sdcmdir utils)

install(TARGETS mri_parse_sdcmdir DESTINATION bin)
