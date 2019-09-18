project(mri_make_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_make_register mri_make_register.cpp)
target_link_libraries(mri_make_register utils)

install(TARGETS mri_make_register DESTINATION bin)
