project(mris_w_to_curv)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_w_to_curv mris_w_to_curv.c)
target_link_libraries(mris_w_to_curv utils)

install(TARGETS mris_w_to_curv DESTINATION bin)
