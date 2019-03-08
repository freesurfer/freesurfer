project(mris_BA_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_BA_segment mris_BA_segment.cpp)
target_link_libraries(mris_BA_segment utils)

install(TARGETS mris_BA_segment DESTINATION bin)
