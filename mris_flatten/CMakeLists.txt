project(mris_flatten)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_flatten mris_flatten.cpp)
target_link_libraries(mris_flatten utils)
install(TARGETS mris_flatten DESTINATION bin)

