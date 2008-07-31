project(tkmedit)
include_directories(${tkmedit_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(tkmedit_SRCS
tkmDisplayArea.c
tkmedit.c
tkmFunctionalVolume.c
tkmMeditWindow.c
)


add_executable(tkmedit ${tkmedit_SRCS})
target_link_libraries(tkmedit ${FS_LIBS})
install(TARGETS tkmedit DESTINATION bin)	

add_subdirectory(test)
