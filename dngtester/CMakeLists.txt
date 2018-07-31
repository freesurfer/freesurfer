project(dngtester)

include_directories(${FS_INCLUDE_DIRS})

# dngtester
add_executable(dngtester dngtester.c)
target_link_libraries(dngtester utils fsgdf)

# cpptester
if(ITK_FOUND)
  include_directories(${ITK_INCLUDE_DIR})
  add_executable(cpptester cpptester.cpp)
  target_link_libraries(cpptester utils fsgdf ${ITK_LIBRARIES})
endif()
