project(mris_seg2annot)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_seg2annot mris_seg2annot.cpp)
target_link_libraries(mris_seg2annot utils)

install(TARGETS mris_seg2annot DESTINATION bin)
