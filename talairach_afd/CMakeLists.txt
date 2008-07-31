project(talairach_afd)
include_directories(${talairach_afd_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(talairach_afd_SRCS
talairach_afd.c
)


add_executable(talairach_afd ${talairach_afd_SRCS})
target_link_libraries(talairach_afd ${FS_LIBS})
install(TARGETS talairach_afd DESTINATION bin)	

