add_executable(testcolortab EXCLUDE_FROM_ALL testcolortab.c)
target_link_libraries(testcolortab utils)

add_executable(test_c_nr_wrapper EXCLUDE_FROM_ALL test_c_nr_wrapper.c)
target_link_libraries(test_c_nr_wrapper utils)

add_executable(extest EXCLUDE_FROM_ALL extest.cpp)
target_link_libraries(extest utils)

add_executable(inftest EXCLUDE_FROM_ALL inftest.cpp)
target_link_libraries(inftest utils)

add_executable(tiff_write_image EXCLUDE_FROM_ALL tiff_write_image.c)
target_link_libraries(tiff_write_image utils)

add_executable(sc_test EXCLUDE_FROM_ALL sc_test.c)
target_link_libraries(sc_test utils)

add_test_executable(sse_mathfun_test EXCLUDE_FROM_ALL sse_mathfun_test.c)
target_link_libraries(sse_mathfun_test m)

add_test_script(NAME utils_test SCRIPT test.py
  DEPENDS
  testcolortab
  test_c_nr_wrapper
  extest
  inftest
  tiff_write_image
  sc_test
)

add_subdirectories(
  mriBuildVoronoiDiagramFloat
  MRIScomputeBorderValues
  mrishash
  mriSoapBubbleFloat
)