project(qdec_glmfit)
include_directories(${qdec_glmfit_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${qdec_glmfit_SOURCE_DIR}/../qdecproject
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(qdec_glmfit_SRCS
qdec_glmfit.cxx
)

SET (QDECLIBS
qdecproject 
${FS_LIBS}
)

add_executable(qdec_glmfit ${qdec_glmfit_SRCS})
target_link_libraries(qdec_glmfit ${QDECLIBS})
install(TARGETS qdec_glmfit DESTINATION bin)	

