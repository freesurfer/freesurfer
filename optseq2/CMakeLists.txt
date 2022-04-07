project(optseq2)

include_directories(${FS_INCLUDE_DIRS})

add_executable(optseq2 optseq2.cpp)
target_link_libraries(optseq2 utils)

add_test_script(NAME optseq2_test SCRIPT test.sh)

install(TARGETS optseq2 DESTINATION bin)
