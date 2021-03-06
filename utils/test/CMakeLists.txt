add_custom_target(utils_tests)

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(filename ${files})
  get_filename_component(stripped_filename ${filename} NAME)
  string(REPLACE ".cpp" "" name ${stripped_filename})

  # fmarching3dnband and heap_mesh are not tests, they are needed in mris_diff_utils_test
  if(name STREQUAL "mris_diff_utils_test")
    set(additional "fmarching3dnband.cpp" "heap_mesh.cpp")
  endif()

  # TODO: eitherr fix tixtest or remove it alltogether, it appears to be unnecessary
  if(name STREQUAL "fmarching3dnband"
     OR name STREQUAL "heap_mesh"
     OR name STREQUAL "tixtest"
     )
    continue()
  endif()
  add_executable(${name} ${stripped_filename} ${additional})
  add_dependencies(utils_tests ${name})
  target_include_directories(${name} PUBLIC ${gtest_SOURCE_DIR}/include ${FS_INCLUDE_DIRS} ${ANN_INCLUDE_DIR}
                                            "/usr/local/opt/cppunit/include"
                             )
  target_link_libraries(${name}
                        utils
                        gtest
                        xml2
                        m
                        ${ANN_LIBRARIES}
                        "/usr/local/opt/cppunit/lib/libcppunit.dylib"
                        )

  add_test(${name} ${name})
  set_property(TEST ${name} PROPERTY LABELS UNIT)

  if(CMAKE_BUILD_TYPE STREQUAL "coverage")
    set_property(TARGET ${name}
                 PROPERTY ADDITIONAL_CLEAN_FILES
                          "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}.dir/${name}.cpp.gcno;${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}.dir/${name}.cpp.gcda"
                 )

  endif()
endforeach()

add_test_script(NAME
                utils_test
                SCRIPT
                test.sh
                DEPENDS
                test_TriangleFile_readWrite
                topology_test
                testcolortab
                test_c_nr_wrapper
                extest
                inftest
                tiff_write_image
                sc_test
                sse_mathfun_test
                )

add_subdirectories(mriBuildVoronoiDiagramFloat MRIScomputeBorderValues mrishash mriSoapBubbleFloat)
