project(mri_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_register mri_register.c)
target_link_libraries(mri_register utils)

install(TARGETS mri_register DESTINATION bin)
