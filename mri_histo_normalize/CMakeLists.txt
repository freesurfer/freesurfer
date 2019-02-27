project(mri_histo_normalize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_histo_normalize mri_histo_normalize.cpp)
target_link_libraries(mri_histo_normalize utils)

install(TARGETS mri_histo_normalize DESTINATION bin)
