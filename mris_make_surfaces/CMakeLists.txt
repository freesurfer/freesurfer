project(mris_make_surfaces)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_make_surfaces mris_make_surfaces.c)
add_help(mris_make_surfaces mris_make_surfaces.help.xml)
target_link_libraries(mris_make_surfaces utils)
install(TARGETS mris_make_surfaces DESTINATION bin)

add_test_script(NAME mris_make_surfaces_test SCRIPT test.sh DEPENDS mris_make_surfaces)

add_executable(mris_refine_surfaces mris_refine_surfaces.c)
target_link_libraries(mris_refine_surfaces utils)
install(TARGETS mris_refine_surfaces DESTINATION bin)

add_executable(mris_mef_surfaces mris_mef_surfaces.c)
target_link_libraries(mris_mef_surfaces utils)
install(TARGETS mris_mef_surfaces DESTINATION bin)

add_executable(mris_exvivo_surfaces mris_exvivo_surfaces.c)
target_link_libraries(mris_exvivo_surfaces utils)
install(TARGETS mris_exvivo_surfaces DESTINATION bin)
