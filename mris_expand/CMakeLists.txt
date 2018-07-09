project(mris_expand)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_expand mris_expand.c)
target_link_libraries(mris_expand utils)

install(TARGETS mris_expand DESTINATION bin)
