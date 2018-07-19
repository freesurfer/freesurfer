project(mri_extract)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_extract mri_extract.c)
target_link_libraries(mri_extract utils)

install(TARGETS mri_extract DESTINATION bin)
