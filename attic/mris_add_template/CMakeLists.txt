project(mris_add_template)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_add_template mris_add_template.cpp)
target_link_libraries(mris_add_template utils)

install(TARGETS mris_add_template DESTINATION bin)
