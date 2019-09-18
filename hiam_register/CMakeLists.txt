project(hiam_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(hiam_register hiam_register.cpp)
target_link_libraries(hiam_register utils)

install(TARGETS hiam_register DESTINATION bin)
