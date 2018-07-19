project(mris_reverse)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_reverse mris_reverse.c)
target_link_libraries(mris_reverse utils)

install(TARGETS mris_reverse DESTINATION bin)
