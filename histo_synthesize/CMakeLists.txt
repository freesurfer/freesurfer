project(histo_synthesize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(histo_synthesize histo_synthesize.cpp)
target_link_libraries(histo_synthesize utils)

install(TARGETS histo_synthesize DESTINATION bin)
