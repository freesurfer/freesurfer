project(mri_normalize)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_normalize mri_normalize.c)
add_help(mri_normalize mri_normalize.help.xml)
target_link_libraries(mri_normalize utils)
install(TARGETS mri_normalize DESTINATION bin)

add_test_script(NAME mri_normalize_test SCRIPT test DEPENDS mri_normalize)

add_executable(mri_long_normalize mri_long_normalize.c)
add_help(mri_long_normalize mri_long_normalize.help.xml)
target_link_libraries(mri_long_normalize utils)
install(TARGETS mri_long_normalize DESTINATION bin)
