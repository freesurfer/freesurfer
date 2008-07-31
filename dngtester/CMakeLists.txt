project(dngtester)
include_directories(${dngtester_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(dngtester_SRCS
dngtester.c
)

SET (dngtester_LIBS
fsgdf
${FS_LIBS}
)


add_executable(dngtester ${dngtester_SRCS})
target_link_libraries(dngtester ${dngtester_LIBS})
install(TARGETS dngtester DESTINATION bin)	

