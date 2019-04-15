project(mris_calc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_calc mris_calc.cpp)
add_help(mris_calc mris_calc.help.xml)
target_link_libraries(mris_calc utils)

install(TARGETS mris_calc DESTINATION bin)
