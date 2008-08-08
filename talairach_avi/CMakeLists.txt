project(talairach_avi Fortran C)
include_directories(${talairach_avi_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET (RMS_SRCS cgauss3d.c fftsol.f imgpad.f polfit.f matopr.f eigen.f param6opr.f)
SET (TRX_SRCS rec.c Getifh.c endianio.c)
SET (ACT_SRCS conc.c)

SET (TALAIRACH_FFILES
	eigen.f 
	ffind.f 
	fftsol.f 
	fimgreg.f 
	ft4imgn.f 
	ft4imgo.f 
	ft4ixyz.f 
	imgpad.f 
	imgvalm.f 
	imgvalx.f 
	matopr.f 
	param6opr.f 
	polfit.f 
	set_rnan.f 
	spline3dvgh.f 
	t4_sub.f 
	to_711-2B.f
)

SET (analyzeto4dfp_SOURCES 
analyzeto4dfp.c 
${TRX_SRCS} 
cflip.c
)
SET (analyzeto4dfp_LDADD ${LIBS_ICC} ${LIB_G2C_A})

SET (4dfptoanalyze_SOURCES 
4dfptoanalyze.c 
${TRX_SRCS} 
cflip.c 
Inithdr.c)
SET (4dfptoanalyze_LDADD ${LIBS_ICC} ${LIB_G2C_A})

SET (flip_4dfp_SOURCES 
flip_4dfp.c ${TRX_SRCS} 
cflip.c)
SET (flip_4dfp_LDADD ${LIBS_ICC} ${LIB_G2C_A})

add_executable(analyzeto4dfp ${analyzeto4dfp_SOURCES})
target_link_libraries(analyzeto4dfp ${analyzeto4dfp_LDADD})
install(TARGETS analyzeto4dfp DESTINATION bin)	

add_executable(4dfptoanalyze ${4dfptoanalyze_SOURCES})
target_link_libraries(4dfptoanalyze ${4dfptoanalyze_LDADD})
install(TARGETS 4dfptoanalyze DESTINATION bin)	

add_executable(flip_4dfp ${flip_4dfp_SOURCES})
target_link_libraries(flip_4dfp ${flip_4dfp_LDADD})
install(TARGETS flip_4dfp DESTINATION bin)	


SET (gauss_4dfp_SOURCES gauss_4dfp.c ${TRX_SRCS} ${RMS_SRCS} ${ACT_SRCS} cgauss3dd.c)
SET (gauss_4dfp_LDADD ${LIBS_ICC} ${FCLIBS} ${LIB_G2C_A})

add_executable(gauss_4dfp ${gauss_4dfp_SOURCES})
target_link_libraries(gauss_4dfp ${gauss_4dfp_LDADD})
install(TARGETS gauss_4dfp DESTINATION bin)	

SET (zero_lt_4dfp_SOURCES zero_lt_4dfp.c ${TRX_SRCS})
SET (zero_lt_4dfp_LDADD ${LIBS_ICC} ${LIB_G2C_A})

add_executable(zero_lt_4dfp ${zero_lt_4dfp_SOURCES})
target_link_libraries(zero_lt_4dfp ${zero_lt_4dfp_LDADD})
install(TARGETS zero_lt_4dfp DESTINATION bin)	

SET (LIN_SRCS t4_sub.f t4scale.c)
SET (TALAIRACH_IMG_SRCS fimgreg.f imgvalm.f imgvalx.f ffind.f)

SET (imgreg_4dfp_SOURCES imgreg_4dfp.c ${TRX_SRCS} ${RMS_SRCS} 
	${LIN_SRCS} ${TALAIRACH_IMG_SRCS} cflip.c)
SET (imgreg_4dfp_LDADD ${LIBS_ICC} ${FCLIBS} ${LIB_G2C_A})
add_executable(imgreg_4dfp ${imgreg_4dfp_SOURCES})
target_link_libraries(imgreg_4dfp ${imgreg_4dfp_LDADD})
install(TARGETS imgreg_4dfp DESTINATION bin)	

SET (T4IMGS_SRCS ft4imgo.f to_711-2B.f ft4imgn.f set_rnan.f)
SET (T4A_SRCS ft4ixyz.f)

SET (t4imgs_4dfp_SOURCES t4imgs_4dfp.c ${TRX_SRCS} ${RMS_SRCS} ${LIN_SRCS} 
	${TALAIRACH_IMG_SRCS} ${T4IMGS_SRCS} ${T4A_SRCS} cflip.c spline3dvgh.f)
SET (t4imgs_4dfp_LDADD ${LIBS_ICC} ${FCLIBS} ${LIB_G2C_A})
add_executable(t4imgs_4dfp ${t4imgs_4dfp_SOURCES})
target_link_libraries(t4imgs_4dfp ${t4imgs_4dfp_LDADD})
install(TARGETS t4imgs_4dfp DESTINATION bin)	


SET (compute_vox2vox_SOURCES compute_vox2vox.c ${TRX_SRCS} ${RMS_SRCS} ${LIN_SRCS} 
	${TALAIRACH_IMG_SRCS} cflip.c set_rnan.f ft4imgo.f )
SET (compute_vox2vox_LDADD ${LIBS_ICC} ${FCLIBS} ${LIB_G2C_A})
add_executable(compute_vox2vox ${compute_vox2vox_SOURCES})
target_link_libraries(compute_vox2vox ${compute_vox2vox_LDADD})
install(TARGETS compute_vox2vox DESTINATION bin)	



SET (ifh2hdr_SOURCES ifh2hdr.c Inithdr.c ${TRX_SRCS})
SET (ifh2hdr_LDADD ${LIBS_ICC} ${LIB_G2C_A})
add_executable(ifh2hdr ${ifh2hdr_SOURCES})
target_link_libraries(ifh2hdr ${ifh2hdr_LDADD})
install(TARGETS ifh2hdr DESTINATION bin)	



