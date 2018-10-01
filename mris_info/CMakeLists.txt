project(mris_info)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_info mris_info.cpp)
target_link_libraries(mris_info utils)
install(TARGETS mris_info DESTINATION bin)
