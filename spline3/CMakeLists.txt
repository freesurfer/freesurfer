project(spline3)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(Spline3_test Spline3_test.cpp)
install(TARGETS Spline3_test DESTINATION bin)

