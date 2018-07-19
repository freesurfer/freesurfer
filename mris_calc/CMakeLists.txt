project(mris_calc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_calc mris_calc.c)
add_help(mris_calc mris_calc.help.xml)
target_link_libraries(mris_calc utils)

install(TARGETS mris_calc DESTINATION bin)
