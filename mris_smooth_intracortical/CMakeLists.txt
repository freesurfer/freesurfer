project(mris_smooth_intracortical)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_smooth_intracortical mris_smooth_intracortical.cpp)
target_link_libraries(mris_smooth_intracortical utils)

install(TARGETS mris_smooth_intracortical DESTINATION bin)
