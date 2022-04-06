project(mri_surfcluster)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surfcluster mri_surfcluster.cpp)
target_link_libraries(mri_surfcluster utils)

add_test_script(NAME mri_surfcluster_test SCRIPT test.sh)

install(TARGETS mri_surfcluster DESTINATION bin)
