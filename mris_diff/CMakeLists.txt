project(mris_diff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_diff mris_diff.c)
target_link_libraries(mris_diff utils)

install(TARGETS mris_diff DESTINATION bin)
