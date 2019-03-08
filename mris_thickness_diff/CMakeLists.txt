project(mris_thickness_diff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_thickness_diff mris_thickness_diff.cpp)
target_link_libraries(mris_thickness_diff utils)

install(TARGETS mris_thickness_diff DESTINATION bin)
