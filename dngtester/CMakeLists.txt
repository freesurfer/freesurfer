project(dngtester)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom)

# always debug
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# dngtester
add_executable(dngtester dngtester.cpp)
target_link_libraries(dngtester utils)

# cpptester
include_directories(${ITK_INCLUDE_DIRS})
add_executable(cpptester cpptester.cpp)
target_link_libraries(cpptester PRIVATE utils ${ITK_LIBRARIES})

if(ARMADILLO_FOUND)
  target_include_directories(cpptester PRIVATE ${ARMADILLO_INCLUDE_DIRS})
  target_link_libraries(cpptester PRIVATE ${ARMADILLO_LIBRARIES})
endif()
