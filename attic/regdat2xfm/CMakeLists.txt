project(regdat2xfm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(regdat2xfm regdat2xfm.cpp)
target_link_libraries(regdat2xfm utils)

install(TARGETS regdat2xfm DESTINATION bin)
