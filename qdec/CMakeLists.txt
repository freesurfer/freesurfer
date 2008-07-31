project(qdec)
include_directories(${qdec_SOURCE_DIR}
${qdec_SOURCE_DIR}/../qdecproject
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}
${vtkutils_SOURCE_DIR}
${vtkfsio_SOURCE_DIR}
) 

FIND_LIBRARY(TIX_LIBRARY
NAMES blt24 BLT24 blt23 BLT23 BLT blt
PATHS ${TCLTK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(BLT_LIBRARY
NAMES tix83 tix82 tix81 tix80 tix
PATHS ${TCLTK_POSSIBLE_LIB_PATHS}
)


FILE (GLOB qdecfilelist *.cxx)

SET(qdec_SRCS
${qdecfilelist}
)

SET (QDECLIBS
qdecproject 
vtkfsio
vtkutils
fsgdf
${TIX_LIBRARY}
${BLT_LIBRARY} 
${FS_LIBS}
)

add_executable(qdec ${qdec_SRCS})
target_link_libraries(qdec ${QDECLIBS})
install(TARGETS qdec DESTINATION bin)	

