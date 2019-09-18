project(mris_defects_pointset)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_defects_pointset mris_defects_pointset.cpp)
add_help(mris_defects_pointset mris_defects_pointset.help.xml)
target_link_libraries(mris_defects_pointset utils)
install(TARGETS mris_defects_pointset DESTINATION bin) 
