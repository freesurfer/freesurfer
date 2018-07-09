project(regdat2xfm)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(regdat2xfm regdat2xfm.cpp)
target_link_libraries(regdat2xfm utils)

install(TARGETS regdat2xfm DESTINATION bin)
