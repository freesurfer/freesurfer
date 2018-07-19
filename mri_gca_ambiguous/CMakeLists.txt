project(mri_gca_ambiguous)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_gca_ambiguous mri_gca_ambiguous.c)
target_link_libraries(mri_gca_ambiguous utils)

install(TARGETS mri_gca_ambiguous DESTINATION bin)
