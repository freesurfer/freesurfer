project(mri_segreg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segreg mri_segreg.cpp)
target_link_libraries(mri_segreg utils)

install(TARGETS mri_segreg DESTINATION bin)
