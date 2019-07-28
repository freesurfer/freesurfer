project(mris_w_to_curv)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_w_to_curv mris_w_to_curv.cpp)
target_link_libraries(mris_w_to_curv utils)

install(TARGETS mris_w_to_curv DESTINATION bin)
