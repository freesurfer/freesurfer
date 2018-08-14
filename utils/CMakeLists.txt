project(utils)

include_directories(
  ${FS_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include/dicom
  ${CMAKE_SOURCE_DIR}/include/NrrdIO
  ${XML2_INCLUDE_DIR}
  ${JPEG_INCLUDE_DIR}
  ${TIFF_INCLUDE_DIR}
  ${MINC_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}/GL
  ${X11_INCLUDE_DIR}
  ${ITK_INCLUDE_DIR}
  ${Expat_INCLUDE_DIR}
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
  chronometerpp.cpp
  class_array.c
  cluster.c
  cma.c
  cmat.c
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
  fs_cost_function.cpp
  fs_lbfgs.cpp
  fs_lbfgs_observer.cpp
  fs_lbfgs_subject.cpp
  fs_powell.cpp
  fsenv.c
  fsglm.c
  fsinit.c
  fsPrintHelp.c
  gca.c
  gcaboundary.c
  gcalinearnode.cpp
  gcalinearprior.cpp
  gcamcomputeLabelsLinearCPU.cpp
  gcamorph.c
  gcamorphtestutils.cpp
  gcarray.c
  gcautils.cpp
  gclass.c
  gcsa.c
  geos.c
  getdelim.c
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
  ic12.c
  ic163842.c
  ic2562.c
  ic42.c
  ic642.c
  icosahedron.c
  image.c
  imageio.c
  imautils.c
  intensity_eig.c
  ipmpar.c
  joint_histo.c
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
  mri.c
  mri2.c
  mri_conform.c
  mri_fastmarching.cpp
  mri_identify.c
  mri_level_set.c
  mri_tess.c
  mri_topology.c
  mri_transform.c
  mriBSpline.c
  mriclass.c
  mricurv.c
  mrifilter.c
  mriflood.c
  mriFunctionalDataAccess.c
  mriHeadPointList.c
  mrihisto.c
  mriio.c
  mriio_nrrd.c
  mriio_nrrd_itk.cpp
  MRIio_old.c
  mrimorph.c
  mrinorm.c
  mripolv.c
  mriprob.c
  mris_compVolFrac.c
  mrisbiorthogonalwavelets.c
  mrisegment.c
  mriset.c
  mrishash.c
  mrisp.c
  MRISrigidBodyAlignGlobal.c
  mrisurf.c
  mriSurface.c
  mrisutils.c
  mriTransform.c
  mriVolume.c
  mrivoxel.c
  nifti1_io.c
  numerics.cpp
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
  rforest.c
  rfutils.c
  romp_support.c
  selxavgio.c
  sig.c
  signa.c
  stats.c
  stc.c
  surfcluster.c
  surfgrad.c
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

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  add_definitions(-DTEEM_32BIT=0)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
  add_definitions(-DTEEM_32BIT=1)
endif()

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

target_link_libraries(utils
  rgb
  hipsstubs
  dicom
  log
  unix
  ${ZLIB_LIBRARIES}
  ${XML2_LIBRARIES}
  ${JPEG_LIBRARIES}
  ${TIFF_LIBRARIES}
  ${MINC_LIBRARIES}
  ${NetCDF_LIBRARIES}
  ${ITK_LIBRARIES}
)

if(NOT APPLE)
  target_link_libraries(utils crypt rt)
endif()


# utils binaries

set(GIFTI_SOURCES gifti_io.c gifti_xml.c nifti1_io.c znzlib.c)

# gifti_tool
add_executable(gifti_tool gifti_tool.c ${GIFTI_SOURCES})
target_link_libraries(gifti_tool ${Expat_LIBRARIES} ${ZLIB_LIBRARIES})
install(TARGETS gifti_tool DESTINATION bin)

# gifti test
add_executable(gifti_test gifti_test.c ${GIFTI_SOURCES})
target_link_libraries(gifti_test ${Expat_LIBRARIES} ${ZLIB_LIBRARIES})
install(TARGETS gifti_test DESTINATION bin)

# xmlToHtml
add_executable(xmlToHtml xmlToHtml.c)
target_link_libraries(xmlToHtml ${XML2_LIBRARIES} ${ZLIB_LIBRARIES})
install(TARGETS xmlToHtml DESTINATION bin)

# fsPrintHelp
add_executable(fsPrintHelp fsPrintHelp.c)
target_compile_definitions(fsPrintHelp PRIVATE BUILD_MAIN)
target_link_libraries(fsPrintHelp ${XML2_LIBRARIES} ${ZLIB_LIBRARIES})
install(TARGETS fsPrintHelp DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/fsPrintHelp DESTINATION bin RENAME fsman)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/fsPrintHelp DESTINATION bin RENAME fshelp)

add_subdirectory(test)
