project(mris_ca_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_ca_label mris_ca_label.c)
add_help(mris_ca_label mris_ca_label.help.xml)
target_link_libraries(mris_ca_label utils)

install(TARGETS mris_ca_label DESTINATION bin)
