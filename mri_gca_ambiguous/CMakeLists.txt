project(mri_gca_ambiguous)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_gca_ambiguous mri_gca_ambiguous.cpp)
target_link_libraries(mri_gca_ambiguous utils)

install(TARGETS mri_gca_ambiguous DESTINATION bin)
