project(utils)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/include/dicom
  ${CMAKE_SOURCE_DIR}/include/NrrdIO
  ${XML2_INCLUDE_DIR}
  ${TIFF_INCLUDE_DIR}
  ${MINC_INCLUDE_DIR}
  ${ITK_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}/GL
)

set(SOURCES
  afni.c
  annotation.c
  art1.c
  artmap.c
  aseg_edit_svm.c
  autoencoder.c
  backprop.c
  bfileio.c
  box.c
  Bruker.c
  chklc.c
  chronometer.c
  cluster.c
  class_array.c
  cma.c
  cmat.c
  cmdargs.c
  colortab.c
  connectcomp.c
  corio.c
  cpputils.cpp
  ctrpoints.c
  cvector.c
  dcdflib.c
  dct.c
  density.c
  diag.c
  DICOMRead.c
  diffuse.c
  dmatrix.c
  dti.c
  dtk.fs.c
  error.c
  evschutils.c
  fcd.c
  fftutils.c
  field_code.c
  filter.c
  fio.c
  flash.c
  fmarchmesh.c
  fmriutils.c
  fsenv.c
  fsglm.c
  fsPrintHelp.c
  fsinit.c
  gcaboundary.c
  gca.c
  gcamorph.c
  gcarray.c
  gclass.c
  gcsa.c
  getdelim.c
  geos.c
  getline.c
  gifti_io.c
  gifti_local.c
  gifti_xml.c
  gtm.c
  gw_ic2562.c
  gw_utils.c
  handle.c
  heap.c
  hippo.c
  histo.c
  joint_histo.c
  ic163842.c
  ic2562.c
  ic642.c
  ic42.c
  ic12.c
  icosahedron.c
  image.c
  imageio.c
  imautils.c
  intensity_eig.c
  ipmpar.c
  kernel.c
  kinput.c
  label.c
  lpafile.c
  machine.c
  MARS_DT_Boundary.c
  matfile.c
  matrix.c
  mgh_filter.c
  mgh_matrix.c
  min_heap.c
  mincutils.c
  minmaxrc.c
  morph.c
  mosaic.c
  mri2.c
  mri.c
  mriBSpline.c
  mriclass.c
  mri_conform.c
  mricurv.c
  mriflood.c
  mriFunctionalDataAccess.c
  mriHeadPointList.c
  mrihisto.c
  mri_identify.c
  MRIio_old.c
  mriio.c
  mri_level_set.c
  mrifilter.c
  mrimorph.c
  mrinorm.c
  mripolv.c
  mriprob.c
  mrisbiorthogonalwavelets.c
  mrisegment.c
  mriset.c
  mrishash.c
  mrisp.c
  MRISrigidBodyAlignGlobal.c
  mriSurface.c
  mrisurf.c
  mrisutils.c
  mri_tess.c
  mri_topology.c
  mri_transform.c
  mriTransform.c
  mriVolume.c
  mrivoxel.c
  nifti1_io.c
  offset.c
  path.c
  pdf.c
  pgmstubs.c
  prime.c
  Progname.c
  queue.c
  randomfields.c
  rbf.c
  rbm.c
  realm.c
  rec.c
  region.c
  registerio.c
  resample.c
  rescale.c
  retinotopy.c
  rfa.c
  rfutils.c
  rforest.c
  romp_support.c
  selxavgio.c
  sig.c
  signa.c
  stats.c
  stc.c
  surfcluster.c
  svm.c
  tags.c
  talairachex.c
  thread.c
  timer.c
  transform.c
  tritri.c
  tukey.c
  utils.c
  version.c
  vertexRotator.c
  vlabels.c
  volcluster.c
  voxlist.c
  x3DList.c
  xDebug.c
  xGLutWindow.c
  xGrowableArray.c
  xList.c
  xSparseVolume.c
  xTypes.c
  xUndoList.c
  xUtilities.c
  xVoxel.c
  znzlib.c
  fs_cost_function.cpp
  fs_lbfgs.cpp
  fs_lbfgs_observer.cpp
  fs_lbfgs_subject.cpp
  fs_powell.cpp
  mri_fastmarching.cpp
  mriio_nrrd.c
  mriio_nrrd_itk.cpp
  numerics.cpp
  chronometerpp.cpp
  gcamorphtestutils.cpp
  gcautils.cpp
  gcalinearnode.cpp
  gcalinearprior.cpp
  cmat.c
  mris_compVolFrac.c
  gcamcomputeLabelsLinearCPU.cpp
)

# add cephes sources
set(SOURCES ${SOURCES}
  cephes/bdtr.c
  cephes/chdtr.c
  cephes/const.c
  cephes/fdtr.c
  cephes/gamma.c
  cephes/igam.c
  cephes/igami.c
  cephes/incbet.c
  cephes/incbi.c
  cephes/mtherr.c
  cephes/ndtr.c
  cephes/ndtri.c
  cephes/polevl.c
  cephes/stdtr.c
)

# nrrdio-required definitions
add_definitions(-DTEEM_DIO=0 -DTEEM_QNANHIBIT=1)

# check for big endian
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  add_definitions(-DTEEM_ENDIAN=4321)
else()
  add_definitions(-DTEEM_ENDIAN=1234)
endif()

# non 32-bit
# todo: some logic needs to be setup to determine this
add_definitions(-DTEEM_32BIT=0)

# add NrrdIO sources
set(SOURCES ${SOURCES}
  NrrdIO/754.c
  NrrdIO/formatPNG.c
  NrrdIO/accessors.c
  NrrdIO/formatPNM.c
  NrrdIO/array.c
  NrrdIO/formatText.c
  NrrdIO/arraysNrrd.c
  NrrdIO/formatVTK.c
  NrrdIO/axis.c
  NrrdIO/gzio.c
  NrrdIO/biffbiff.c
  NrrdIO/keyvalue.c
  NrrdIO/comment.c
  NrrdIO/methodsNrrd.c
  NrrdIO/defaultsNrrd.c
  NrrdIO/miscAir.c
  NrrdIO/dio.c
  NrrdIO/mop.c
  NrrdIO/encodingAscii.c
  NrrdIO/parseAir.c
  NrrdIO/encodingBzip2.c
  NrrdIO/parseNrrd.c
  NrrdIO/encoding.c
  NrrdIO/encodingGzip.c
  NrrdIO/encodingHex.c
  NrrdIO/read.c
  NrrdIO/encodingRaw.c
  NrrdIO/reorder.c
  NrrdIO/endianAir.c
  NrrdIO/endianNrrd.c
  NrrdIO/sane.c
  NrrdIO/enum.c
  NrrdIO/simple.c
  NrrdIO/enumsNrrd.c
  NrrdIO/string.c
  NrrdIO/format.c
  NrrdIO/subset.c
  NrrdIO/formatEPS.c
  NrrdIO/write.c
  NrrdIO/formatNRRD.c
)

add_library(utils STATIC ${SOURCES})

# the following libraries might get built locally via ExternalProject,
# so we should make they're always built first
add_external_dependencies(utils xml2 tiff minc itk)

target_link_libraries(utils
  rgb
  hipsstubs
  fsgdf
  dicom
  log
  unix
  crypt
  ${ZLIB_LIBRARIES}
  ${XML2_LIBRARY}
  ${JPEG_LIBRARY}
  ${TIFF_LIBRARY}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARY}
  ${MINC_LIBRARIES}
  ${NETCDF_LIBRARY}
  ${ITK_LIBRARIES}
)

if(NOT APPLE)
  target_link_libraries(utils rt)
endif()


# utils binaries

set(GIFTI_SOURCES gifti_io.c gifti_xml.c nifti1_io.c znzlib.c)

# gifti_tool
add_executable(gifti_tool gifti_tool.c ${GIFTI_SOURCES})
target_link_libraries(gifti_tool ${EXPAT_LIBRARY} ${ZLIB_LIBRARIES})
add_external_dependencies(gifti_tool expat)
install(TARGETS gifti_tool DESTINATION bin)

# gifti test
add_executable(gifti_test gifti_test.c ${GIFTI_SOURCES})
target_link_libraries(gifti_test ${EXPAT_LIBRARY} ${ZLIB_LIBRARIES})
add_external_dependencies(gifti_test expat)
install(TARGETS gifti_test DESTINATION bin)

# xmlToHtml
add_executable(xmlToHtml xmlToHtml.c)
target_link_libraries(xmlToHtml ${XML2_LIBRARY} ${ZLIB_LIBRARIES})
add_external_dependencies(xmlToHtml xml2)
install(TARGETS xmlToHtml DESTINATION bin)

# fsPrintHelp
add_executable(fsPrintHelp fsPrintHelp.c)
target_compile_definitions(fsPrintHelp PUBLIC BUILD_MAIN)
target_link_libraries(fsPrintHelp ${XML2_LIBRARY} ${ZLIB_LIBRARIES})
add_external_dependencies(fsPrintHelp xml2)
install(TARGETS fsPrintHelp DESTINATION bin)
install(PROGRAMS fsPrintHelp DESTINATION bin RENAME fsman)
install(PROGRAMS fsPrintHelp DESTINATION bin RENAME fshelp)

