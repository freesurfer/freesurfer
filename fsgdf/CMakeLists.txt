project(fsgdf)
include_directories(${fsgdf_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(fsgdf_SRCS
fsgdf.c
fsgdf_wrap.c
FsgdfPlot.cxx
)

SET (test-fsgdf-libs
fsgdf
${FS_LIBS}
)

add_library(fsgdf STATIC ${fsgdf_SRCS})
add_executable(test-fsgdf test.c)
target_link_libraries(test-fsgdf ${test-fsgdf-libs})
install(TARGETS test-fsgdf DESTINATION bin)	

