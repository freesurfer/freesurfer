project(dngtester)

include_directories(${FS_INCLUDE_DIRS})

# always debug
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# dngtester
add_executable(dngtester dngtester.c)
target_link_libraries(dngtester utils)

# cpptester
include_directories(${ITK_INCLUDE_DIRS})
add_executable(cpptester cpptester.cpp)
target_link_libraries(cpptester utils ${ITK_LIBRARIES})
