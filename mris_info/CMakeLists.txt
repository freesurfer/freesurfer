project(mris_info)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_info mris_info.cpp)
target_link_libraries(mris_info utils)
install(TARGETS mris_info DESTINATION bin)
