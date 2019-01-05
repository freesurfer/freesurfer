project(mris_curvature)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_curvature mris_curvature.c)
add_help(mris_curvature mris_curvature.help.xml)
target_link_libraries(mris_curvature utils)

add_test_script(NAME mris_curvature_test SCRIPT test_mris_curvature DEPENDS mris_curvature)
configure_file(testdata.tar.gz . COPYONLY)

install(TARGETS mris_curvature DESTINATION bin)
