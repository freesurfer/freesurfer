project(mri_tessellate)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_tessellate mri_tessellate.c)
add_help(mri_tessellate mri_tessellate.help.xml)
target_link_libraries(mri_tessellate utils)
install(TARGETS mri_tessellate DESTINATION bin)

add_executable(ico_supersample ico_supersample.c)
target_link_libraries(ico_supersample utils)
install(TARGETS ico_supersample DESTINATION bin)
