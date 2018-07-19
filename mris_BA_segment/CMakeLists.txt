project(mris_BA_segment)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_BA_segment mris_BA_segment.c)
target_link_libraries(mris_BA_segment utils)

install(TARGETS mris_BA_segment DESTINATION bin)
