project(dngtester)

include_directories(${FS_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include)

# dngtester
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
add_executable(dngtester dngtester.c)
target_link_libraries(dngtester utils fsgdf)

# cpptester
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
if(ITK_FOUND)
  include_directories(${ITK_INCLUDE_DIR})
  add_executable(cpptester cpptester.cpp)
  target_link_libraries(cpptester utils fsgdf ${ITK_LIBRARIES})
endif()
