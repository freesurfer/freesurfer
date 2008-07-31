project(utilscpp-test)
include_directories(${utilscpp-test_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(utilscpp-test_SRCS
ctest.c
)

SET (UTILSCPP_TEST_LIBS utils rgb dicom hipsstubs log image unix utilscpp ${LIBXML2_LIBRARIES} ${TIFF_LIBRARIES} ${JPEG_LIBRARIES} ${EXPAT_LIBRARIES} ${MINC_LIBRARIES} ${ZLIB_LIBRARIES} ${VNL_LIBRARIES} ${VNL_ALGO_LIBRARIES} ${V3P_LIBRARIES} ${VCL_LIBRARIES} ${VPL_LIBRARIES} ${TCL_LIBRARIES} ${TK_LIBRARIES} ${VTK_LIBRARIES})
add_executable(utilscpp-test ${utilscpp-test_SRCS})
target_link_libraries(utilscpp-test ${UTILSCPP_TEST_LIBS})
install(TARGETS utilscpp-test DESTINATION bin)	

