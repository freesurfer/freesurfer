project(mrisp_paint)
include_directories(${mrisp_paint_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mrisp_paint_SRCS
mrisp_paint.c
)


add_executable(mrisp_paint ${mrisp_paint_SRCS})
target_link_libraries(mrisp_paint ${FS_LIBS})
install(TARGETS mrisp_paint DESTINATION bin)	

