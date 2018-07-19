project(mri_mcsim)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_mcsim mri_mcsim.c)
target_link_libraries(mri_mcsim utils)

install(TARGETS mri_mcsim DESTINATION bin)
