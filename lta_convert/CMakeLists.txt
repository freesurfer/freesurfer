project(lta_convert)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(lta_convert lta_convert.cpp)
add_help(lta_convert lta_convert.help.xml)
target_link_libraries(lta_convert utils)

install(TARGETS lta_convert DESTINATION bin)
