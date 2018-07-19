project(mri_segreg)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_segreg mri_segreg.c)
target_link_libraries(mri_segreg utils)

install(TARGETS mri_segreg DESTINATION bin)
