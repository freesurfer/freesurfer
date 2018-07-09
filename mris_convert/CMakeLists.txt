project(mris_convert)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_convert mris_convert.c)
target_link_libraries(mris_convert utils)

install(TARGETS mris_convert DESTINATION bin)
