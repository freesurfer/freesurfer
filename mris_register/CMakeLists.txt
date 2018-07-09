project(mris_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_register mris_register.c)
add_help(mris_register mris_register.help.xml)
target_link_libraries(mris_register utils)

install(TARGETS mris_register DESTINATION bin)
