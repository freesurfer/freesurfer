project(mri_histo_normalize)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_histo_normalize mri_histo_normalize.c)
target_link_libraries(mri_histo_normalize utils)

install(TARGETS mri_histo_normalize DESTINATION bin)
