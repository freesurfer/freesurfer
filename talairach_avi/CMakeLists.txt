project(talairach_avi)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(ICC_LIBRARIES irc svml imf)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ax=CORE-AVX512 -qopenmp -132 -nofor_main")
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -fno-second-underscore -fallow-invalid-boz")
  else()
    if(HOST_OS MATCHES "CentOS9stream")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -fno-second-underscore -fallow-invalid-boz")
    else()  
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -fno-second-underscore")
    endif()
  endif()
endif()

install_configured(
  talairach_avi
  mpr2mni305
  brec
  brec.awk
  t4img_4dfp
  tal_QC_AZS
  DESTINATION bin
)

add_test_script(NAME talairach_avi_test SCRIPT test.sh)

set(TRX_SOURCES rec.c Getifh.c endianio.c)

add_executable(analyzeto4dfp analyzeto4dfp.c cflip.c ${TRX_SOURCES})
target_link_libraries(analyzeto4dfp ${ICC_LIBRARIES} )
install(TARGETS analyzeto4dfp DESTINATION bin)

add_executable(4dfptoanalyze 4dfptoanalyze.c cflip.c Inithdr.c ${TRX_SOURCES})
target_link_libraries(4dfptoanalyze ${ICC_LIBRARIES} )
install(TARGETS 4dfptoanalyze DESTINATION bin)

add_executable(flip_4dfp flip_4dfp.c cflip.c ${TRX_SOURCES})
target_link_libraries(flip_4dfp ${ICC_LIBRARIES})
install(TARGETS flip_4dfp DESTINATION bin)

add_executable(zero_lt_4dfp zero_lt_4dfp.c ${TRX_SOURCES})
target_link_libraries(zero_lt_4dfp ${ICC_LIBRARIES})
install(TARGETS zero_lt_4dfp DESTINATION bin)

add_executable(ifh2hdr ifh2hdr.c Inithdr.c ${TRX_SOURCES})
target_link_libraries(ifh2hdr ${ICC_LIBRARIES})
install(TARGETS ifh2hdr DESTINATION bin)

set(RMS_SOURCES
  cgauss3d.c
  fftsol.f
  imgpad.f
  polfit.f
  matopr.f
  eigen.f
  param6opr.f
)

add_definitions(-DHAVE_GFORTRAN)

# this prevents linking to gfortran dynmically
if(APPLE)
  get_filename_component(FORTRAN_REALPATH "${CMAKE_Fortran_COMPILER}" REALPATH)
  get_filename_component(FORTRAN_BASE_DIR "${FORTRAN_REALPATH}/../.." REALPATH)
  if (NOT APPLE_ARM64)
     file(GLOB_RECURSE STATIC_QUADMATH FOLLOW_SYMLINKS "${FORTRAN_BASE_DIR}/lib/*/libquadmath.a")
  endif()
  if(STATIC_QUADMATH)
    set(CMAKE_Fortran_FLAGS "-w -static-libgfortran -lgfortran -lgcc -lgcc_eh -lSystem -nodefaultlibs ${STATIC_QUADMATH} ${CMAKE_Fortran_FLAGS}")
  else()
    message(WARNING "Not linking FORTRAN executables with static libquadmath")
  endif()
else()
  if(NO_STATIC_GFORTRAN)
    message(WARNING "Not linking FORTRAN executables with static libgfortran")
  else()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgfortran")
  endif() 
endif()

add_executable(gauss_4dfp gauss_4dfp.c cgauss3dd.c conc.c ${TRX_SOURCES} ${RMS_SOURCES})
target_link_libraries(gauss_4dfp ${ICC_LIBRARIES})
install(TARGETS gauss_4dfp DESTINATION bin)

set(LIN_SOURCES t4_sub.f t4scale.c)
set(IMG_SOURCES fimgreg.f imgvalm.f imgvalx.f ffind.f)

add_executable(imgreg_4dfp imgreg_4dfp.c ${TRX_SOURCES} ${RMS_SOURCES}
               ${LIN_SOURCES} ${IMG_SOURCES} cflip.c)
target_link_libraries(imgreg_4dfp ${ICC_LIBRARIES})
install(TARGETS imgreg_4dfp DESTINATION bin)

add_executable(t4imgs_4dfp t4imgs_4dfp.c ${TRX_SOURCES} ${RMS_SOURCES}
               ${LIN_SOURCES} ${IMG_SOURCES} cflip.c spline3dvgh.f
               ft4imgo.f to_711-2B.f ft4imgn.f set_rnan.f ft4ixyz.f)
target_link_libraries(t4imgs_4dfp ${ICC_LIBRARIES})
install(TARGETS t4imgs_4dfp DESTINATION bin)

add_executable(compute_vox2vox compute_vox2vox.c ${TRX_SOURCES} ${RMS_SOURCES}
               ${LIN_SOURCES} ${IMG_SOURCES} cflip.c set_rnan.f ft4imgo.f)
target_link_libraries(compute_vox2vox ${ICC_LIBRARIES})
install(TARGETS compute_vox2vox DESTINATION bin)

install(FILES talairach_avi.help.xml DESTINATION docs/xml)

install_symlinks(
  3T18yoSchwartzReactN32_as_orig.4dfp.hdr
  3T18yoSchwartzReactN32_as_orig.4dfp.ifh
  3T18yoSchwartzReactN32_as_orig.4dfp.img
  3T18yoSchwartzReactN32_as_orig.4dfp.img.rec
  3T18yoSchwartzReactN32_as_orig.4dfp.mat
  3T18yoSchwartzReactN32_as_orig.lst
  711-2B_as_mni_average_305_mask.4dfp.hdr
  711-2B_as_mni_average_305_mask.4dfp.ifh
  711-2B_as_mni_average_305_mask.4dfp.img
  711-2B_as_mni_average_305_mask.4dfp.img.rec
  711-2C_as_mni_average_305.4dfp.hdr
  711-2C_as_mni_average_305.4dfp.ifh
  711-2C_as_mni_average_305.4dfp.img
  711-2C_as_mni_average_305.4dfp.img.rec
  711-2C_as_mni_average_305.4dfp.mat
  TYPE files
  DESTINATION average
)

# not required for minimal recon-all
if(NOT MINIMAL)

  install_symlinks(
    711-2B_as_mni_average_305.4dfp.hdr
    711-2B_as_mni_average_305.4dfp.ifh
    711-2B_as_mni_average_305.4dfp.img
    711-2B_as_mni_average_305.4dfp.img.rec
    mni_average_305.4dfp.hdr
    mni_average_305.4dfp.ifh
    mni_average_305.4dfp.img
    mni_average_305.4dfp.img.rec
    RLB700_atlas_as_orig.4dfp.hdr
    RLB700_atlas_as_orig.4dfp.ifh
    RLB700_atlas_as_orig.4dfp.img
    RLB700_atlas_as_orig.4dfp.img.rec
    SVIP_Adult_Comp_N24_as_orig.4dfp.hdr
    SVIP_Adult_Comp_N24_as_orig.4dfp.ifh
    SVIP_Adult_Comp_N24_as_orig.4dfp.img
    SVIP_Adult_Comp_N24_as_orig.4dfp.img.rec
    SVIP_Adult_Control_N12_as_orig.4dfp.hdr
    SVIP_Adult_Control_N12_as_orig.4dfp.ifh
    SVIP_Adult_Control_N12_as_orig.4dfp.img
    SVIP_Adult_Control_N12_as_orig.4dfp.img.rec
    SVIP_Adult_Patient_N12_as_orig.4dfp.hdr
    SVIP_Adult_Patient_N12_as_orig.4dfp.ifh
    SVIP_Adult_Patient_N12_as_orig.4dfp.img
    SVIP_Adult_Patient_N12_as_orig.4dfp.img.rec
    SVIP_Child_Comp_N24_as_orig.4dfp.hdr
    SVIP_Child_Comp_N24_as_orig.4dfp.ifh
    SVIP_Child_Comp_N24_as_orig.4dfp.img
    SVIP_Child_Comp_N24_as_orig.4dfp.img.rec
    SVIP_Full_Comp_N48_as_orig.4dfp.hdr
    SVIP_Full_Comp_N48_as_orig.4dfp.ifh
    SVIP_Full_Comp_N48_as_orig.4dfp.img
    SVIP_Full_Comp_N48_as_orig.4dfp.img.rec
    SVIP_Young_Control_N12_as_orig.4dfp.hdr
    SVIP_Young_Control_N12_as_orig.4dfp.ifh
    SVIP_Young_Control_N12_as_orig.4dfp.img
    SVIP_Young_Control_N12_as_orig.4dfp.img.rec
    SVIP_Young_Patient_N12_as_orig.4dfp.hdr
    SVIP_Young_Patient_N12_as_orig.4dfp.ifh
    SVIP_Young_Patient_N12_as_orig.4dfp.img
    SVIP_Young_Patient_N12_as_orig.4dfp.img.rec
    TYPE files
    DESTINATION average
  )

endif()

