project(anatomicuts)

if(VTK_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Code
    ${CMAKE_SOURCE_DIR}/freeview/track_io
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )

  include(${ITK_USE_FILE})  
  include(${VTK_USE_FILE})

  if(NOT APPLE)
    add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedefs  -Wno-maybe-uninitialized )
  else()
    # Mac OS 10.13.6 native clang compiler does not support -Wno-maybe-uninitialized
    add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedefs -Wno-inconsistent-missing-override )
  endif()

  set(TRACKIO ../freeview/track_io/TrackIO.cpp)

  set(ITK_LIBRARIES "ITKCommon;ITKIOImageBase;ITKIOBMP;ITKIOBioRad;ITKIOBruker;ITKIOGDCM;ITKIOGE;ITKIONRRD;ITKIOTIFF;ITKIOVTK;ITKIOStimulate;ITKIONIFTI;ITKIOMINC;ITKIOGIPL;ITKIOHDF5;ITKIOJPEG;ITKIOLSM;ITKIOMRC;ITKIOMeta;ITKIOPNG")
  if(HOST_OS MATCHES "CentOS9stream" OR (HOST_OS MATCHES "Ubuntu22" AND HOST_OS_GCC_VERSION MATCHES "11"))
    set(ITK_LIBRARIES "${ITK_LIBRARIES};ITKIOJPEG2000")
  endif()

#AnatomiCuts correspondences
  add_executable(dmri_AnatomiCuts AnatomiCuts.cxx ${TRACKIO})
  target_link_libraries(dmri_AnatomiCuts utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_AnatomiCuts DESTINATION bin)

#AnatomiCuts
  add_executable(dmri_match AnatomiCuts_correspondences.cxx ${TRACKIO})
  target_link_libraries(dmri_match utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_match DESTINATION bin)

#TRK tools for header
  add_executable(trk_tools trk_tools.cxx ${TRACKIO})
  target_link_libraries(trk_tools utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS trk_tools DESTINATION bin)

#streamline filtering
  add_executable(streamlineFilter StreamlineFilter.cxx ${TRACKIO})
  target_link_libraries(streamlineFilter utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS streamlineFilter DESTINATION bin)
  
#dmri_stats_ac
  add_executable(dmri_stats_ac dmri_stats_ac.cxx ${TRACKIO})
  target_link_libraries(dmri_stats_ac utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_stats_ac DESTINATION bin)

#orientation planes
  add_executable(testOrientationPlanesFromParcellation testOrientationFromParcellationFilter.cxx ) 
  target_link_libraries(testOrientationPlanesFromParcellation ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS testOrientationPlanesFromParcellation DESTINATION bin)

#Save histograms
  add_executable(dmri_saveHistograms SaveHistograms.cxx ${TRACKIO}) 
  target_link_libraries(dmri_saveHistograms ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_saveHistograms DESTINATION bin)

#dmri_groupByEndpoints
  add_executable(dmri_groupByEndpoints dmri_groupByEndpoints.cxx ${TRACKIO}) 
  target_link_libraries(dmri_groupByEndpoints utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_groupByEndpoints DESTINATION bin)

#dmri_extractSurfaceMeasurements
  # also needs ITKIOTransformHDF5;ITKIOTransformMatlab;ITKIOTransformInsightLegacy
  add_executable(dmri_extractSurfaceMeasurements dmri_extractSurfaceMeasurements.cxx ${TRACKIO}) 
  target_link_libraries(dmri_extractSurfaceMeasurements utils nifti "${ITK_LIBRARIES};ITKIOTransformHDF5;ITKIOTransformMatlab;ITKIOTransformInsightLegacy" ${VTK_LIBRARIES})
  install(TARGETS dmri_extractSurfaceMeasurements DESTINATION bin)

#dmri_coloredFA
  add_executable(dmri_coloredFA dmri_coloredFA.cxx ${TRACKIO})
  target_link_libraries(dmri_coloredFA utils nifti ${OPENGL_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_coloredFA DESTINATION bin)

#dmri_projectEndPoints
  add_executable(dmri_projectEndPoints dmri_projectEndPoints.cxx ${TRACKIO})
  target_link_libraries(dmri_projectEndPoints utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS dmri_projectEndPoints DESTINATION bin)

#dmri_neighboringRegions
  # only needs ITKCommon
  add_executable(dmri_neighboringRegions dmri_neighboringRegions.cxx ${TRACKIO})
  target_link_libraries(dmri_neighboringRegions utils nifti "ITKCommon" ${VTK_LIBRARIES})
  install(TARGETS dmri_neighboringRegions DESTINATION bin)

endif()

install(PROGRAMS dmri_ac.sh DESTINATION bin)

install_pyscript(
	diffusionUtils
	anatomiCutsUtils
	dmri_violinPlots
)
