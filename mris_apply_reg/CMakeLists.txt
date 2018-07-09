project(mris_apply_reg)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_apply_reg mris_apply_reg.c)
add_help(mris_apply_reg mris_apply_reg.help.xml)
target_link_libraries(mris_apply_reg utils)

install(TARGETS mris_apply_reg DESTINATION bin)
