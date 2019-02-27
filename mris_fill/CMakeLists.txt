project(mris_fill)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_fill mris_fill.cpp)
target_link_libraries(mris_fill utils)

install(TARGETS mris_fill DESTINATION bin)
