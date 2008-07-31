project(svm-lib)
include_directories(${svm-lib_SOURCE_DIR}
${svm-lib_SOURCE_DIR}/include
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(svm-lib_SRCS
svm-data-param.cpp
svm-kernel.cpp
svm-lib-c.cpp
svm-lib.cpp
svm-model.cpp
svm-param.cpp
svm-pr-loqo.c
svm-similarity.cpp
svm-sketch.cpp
)


add_library(svm-lib STATIC ${svm-lib_SRCS})
