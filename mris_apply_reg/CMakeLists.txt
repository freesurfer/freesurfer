project(mris_apply_reg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_apply_reg mris_apply_reg.cpp)
add_help(mris_apply_reg mris_apply_reg.help.xml)
target_link_libraries(mris_apply_reg utils)

install(TARGETS mris_apply_reg DESTINATION bin)
