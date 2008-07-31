project(mris_spherical_wavelets)
include_directories(${mris_spherical_wavelets_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_spherical_wavelets_SRCS
mris_spherical_wavelets.c
mris_spherical_wavelets.cpp
)


add_executable(mris_spherical_wavelets ${mris_spherical_wavelets_SRCS})
target_link_libraries(mris_spherical_wavelets ${FS_LIBS})
install(TARGETS mris_spherical_wavelets DESTINATION bin)	


