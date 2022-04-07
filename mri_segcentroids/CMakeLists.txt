project(mri_segcentroids)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segcentroids mri_segcentroids.cpp)
add_help(mri_segcentroids mri_segcentroids.help.xml)
target_link_libraries(mri_segcentroids utils)

add_test_script(NAME mri_segcentroids_test SCRIPT test.sh)

install(TARGETS mri_segcentroids DESTINATION bin) 

