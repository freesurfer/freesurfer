project(mri_wbc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_wbc mri_wbc.c)
target_link_libraries(mri_wbc utils)

install(TARGETS mri_wbc DESTINATION bin)
