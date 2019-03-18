project(mris_make_template)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_make_template mris_make_template.cpp)
target_link_libraries(mris_make_template utils)

install(TARGETS mris_make_template DESTINATION bin)
