project(mris_watershed)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_watershed mris_watershed.cpp)
add_help(mris_watershed mris_watershed.help.xml)
target_link_libraries(mris_watershed utils)

install(TARGETS mris_watershed DESTINATION bin)
