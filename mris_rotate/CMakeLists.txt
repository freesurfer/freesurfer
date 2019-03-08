project(mris_rotate)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_rotate mris_rotate.cpp)
target_link_libraries(mris_rotate utils)

install(TARGETS mris_rotate DESTINATION bin)
