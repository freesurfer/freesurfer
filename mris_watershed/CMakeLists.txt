project(mris_watershed)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_watershed mris_watershed.c)
add_help(mris_watershed mris_watershed.help.xml)
target_link_libraries(mris_watershed utils)

install(TARGETS mris_watershed DESTINATION bin)
