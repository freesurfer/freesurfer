project(mris_congeal)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_congeal mris_congeal.cpp)
target_link_libraries(mris_congeal utils)
install(TARGETS mris_congeal DESTINATION bin)
