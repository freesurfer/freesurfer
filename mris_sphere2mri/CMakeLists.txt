project(mris_sphere2mri)

add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-unused-variable -Wno-unused-local-typedefs -Wno-reorder)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/resurf/Code)

add_executable(mris_sphere2mri mris_sphere2mri.cpp)
add_help(mris_sphere2mri mris_sphere2mri.help.xml mri_stopmask.help.xml)
target_link_libraries(mris_sphere2mri utils)
install(TARGETS mris_sphere2mri DESTINATION bin)

# add_test_script(NAME mris_sphere2mri_test SCRIPT test.sh)

