project(histo_synthesize)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(histo_synthesize histo_synthesize.c)
target_link_libraries(histo_synthesize utils)

install(TARGETS histo_synthesize DESTINATION bin)
