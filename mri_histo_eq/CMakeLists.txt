project(mri_histo_eq)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_histo_eq mri_histo_eq.cpp)
target_link_libraries(mri_histo_eq utils)

install(TARGETS mri_histo_eq DESTINATION bin)
