project(label2flat)

include_directories(${FS_INCLUDE_DIRS})

add_executable(label2flat label2flat.cpp)
target_link_libraries(label2flat utils)

install(TARGETS label2flat DESTINATION bin)
