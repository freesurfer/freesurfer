project(dicom)
include_directories(${dicom_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(dicom_SRCS
condition.c
ctnthread.c
dcm.c
dcmcond.c
dcmdict.c
dcmsupport.c
lst.c
)


add_library(dicom STATIC ${dicom_SRCS})

