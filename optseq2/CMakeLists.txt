project(optseq2)

include_directories(${FS_INCLUDE_DIRS})

add_executable(optseq2 optseq2.c)
target_link_libraries(optseq2 utils)

install(TARGETS optseq2 DESTINATION bin)
