project(mris_entropy)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_entropy mris_entropy.c)
target_link_libraries(mris_entropy utils)

install(TARGETS mris_entropy DESTINATION bin)
