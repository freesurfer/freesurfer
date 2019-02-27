project(mris_talairach)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_talairach mris_talairach.cpp)
target_link_libraries(mris_talairach utils)

install(TARGETS mris_talairach DESTINATION bin)
