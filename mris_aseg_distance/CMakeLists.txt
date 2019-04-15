project(mris_aseg_distance)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_aseg_distance mris_aseg_distance.cpp)
target_link_libraries(mris_aseg_distance utils)

install(TARGETS mris_aseg_distance DESTINATION bin)
