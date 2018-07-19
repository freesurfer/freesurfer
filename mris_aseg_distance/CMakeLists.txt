project(mris_aseg_distance)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_aseg_distance mris_aseg_distance.c)
target_link_libraries(mris_aseg_distance utils)

install(TARGETS mris_aseg_distance DESTINATION bin)
