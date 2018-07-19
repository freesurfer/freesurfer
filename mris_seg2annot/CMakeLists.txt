project(mris_seg2annot)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_seg2annot mris_seg2annot.c)
target_link_libraries(mris_seg2annot utils)

install(TARGETS mris_seg2annot DESTINATION bin)
