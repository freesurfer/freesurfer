project(mris_twoclass)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_twoclass mris_twoclass.c)
target_link_libraries(mris_twoclass utils)

install(TARGETS mris_twoclass DESTINATION bin)
