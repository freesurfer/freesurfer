project(talairach_avi)
include_directories(${talairach_avi_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(talairach_avi_SRCS
4dfptoanalyze.c
analyzeto4dfp.c
cflip.c
cgauss3d.c
cgauss3dd.c
compute_vox2vox.c
conc.c
eigen.f
endianio.c
ffind.f
fftsol.f
fimgreg.f
flip_4dfp.c
ft4imgn.f
ft4imgo.f
ft4ixyz.f
gauss_4dfp.c
Getifh.c
ifh2hdr.c
imgpad.f
imgreg_4dfp.c
imgvalm.f
imgvalx.f
Inithdr.c
matopr.f
param6opr.f
polfit.f
rec.c
set_rnan.f
spline3dvgh.f
t4_sub.f
t4imgs_4dfp.c
t4scale.c
to_711-2B.f
zero_lt_4dfp.c
)


add_executable(talairach_avi ${talairach_avi_SRCS})
target_link_libraries(talairach_avi ${FS_LIBS})
install(TARGETS talairach_avi DESTINATION bin)	
