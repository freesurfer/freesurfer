project(mri_wbc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_wbc mri_wbc.cpp)
target_link_libraries(mri_wbc utils)

install(TARGETS mri_wbc DESTINATION bin)
