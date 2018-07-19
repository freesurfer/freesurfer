project(mri_vol2surf)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_vol2surf mri_vol2surf.c)
target_link_libraries(mri_vol2surf utils)

install(TARGETS mri_vol2surf DESTINATION bin)
