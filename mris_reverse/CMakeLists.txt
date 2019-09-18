project(mris_reverse)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_reverse mris_reverse.cpp)
target_link_libraries(mris_reverse utils)

install(TARGETS mris_reverse DESTINATION bin)
