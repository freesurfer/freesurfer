project(mris_average_curvature)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_average_curvature mris_average_curvature.cpp)
target_link_libraries(mris_average_curvature utils)

install(TARGETS mris_average_curvature DESTINATION bin)
