project(optseq2)

include_directories(${FS_INCLUDE_DIRS})

add_executable(optseq2 optseq2.c)
target_link_libraries(optseq2 utils)

add_test_script(NAME optseq2_test SCRIPT test_optseq2.csh DEPENDS optseq2)
configure_file(test_data.tar.gz . COPYONLY)

install(TARGETS optseq2 DESTINATION bin)
