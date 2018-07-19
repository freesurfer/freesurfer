project(mris_congeal)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_congeal mris_congeal.c)
target_link_libraries(mris_congeal utils)
install(TARGETS mris_congeal DESTINATION bin)
