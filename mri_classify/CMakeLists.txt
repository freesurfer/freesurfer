project(mri_classify)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_classify mri_classify.c)
target_link_libraries(mri_classify utils)

install(TARGETS mri_classify DESTINATION bin)
