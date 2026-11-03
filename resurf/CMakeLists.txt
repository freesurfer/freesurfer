project(resurf)

if(ITK_FOUND AND VTK_FOUND)
  include(${ITK_USE_FILE})
  if(HOST_OS MATCHES "Rocky9" OR HOST_OS MATCHES "Ubuntu24")
     include_directories(${VTK_INCLUDE_DIRS})
  else()
     include(${VTK_USE_FILE})
  endif()

  if(HOST_OS MATCHES "Ubuntu24")
     if(NOT UBUNTU_ARM64)
        set(UB24_ITK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(UB24_VTK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        # set(UB24_COMMON_LINK_OPTIONS "-Wl,--unresolved-symbols=ignore-in-shared-libs -Wl,-rpath,/usr/lib/${UBUNTU_GCC_LIBDIR}/lib:/usr/lib/jvm/default-java/lib:/usr/lib/jvm/default-java/lib/server:/usr/lib/${UBUNTU_GCC_LIBDIR}/hdf5/serial")
        set(UB24_COMMON_LINK_OPTIONS "-Wl,--unresolved-symbols=ignore-in-shared-libs")
        set(EXT "5.3.so")
     else()
        # don't have native itk libs on arm64
        set(UB24_ITK_BASE "/usr/pubsw/packages/itk/5.4_ubuntu24/lib")
        set(UB24_VTK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        # set(UB24_COMMON_LINK_OPTIONS "-Wl,--unresolved-symbols=ignore-in-shared-libs -Wl,-rpath,${UB24_ITK_BASE}:${UB24_VTK_BASE}/lib:/usr/lib/jvm/default-java/lib:/usr/lib/jvm/default-java/lib/server:/usr/lib/${UBUNTU_GCC_LIBDIR}/hdf5/serial")
        set(UB24_COMMON_LINK_OPTIONS "-Wl,--unresolved-symbols=ignore-in-shared-libs")
        set(EXT "5.4.so")
     endif()

     # reset and do not use default glob that pulls in every library
     set(ITK_LIBRARIES "")

     set(MRIS_COMMON_ITK_LIBS libITKIOBMP-${EXT} libITKIOBioRad-${EXT} libITKIOBruker-${EXT} libITKIOGDCM-${EXT} libITKIOGE-${EXT} libITKIOGIPL-${EXT} libITKIOHDF5-${EXT} libITKIOJPEG-${EXT} libITKIOJPEG2000-${EXT} libITKIOLSM-${EXT} libITKIOMINC-${EXT} libITKIOMRC-${EXT} libITKIOMeta-${EXT} libITKIONIFTI-${EXT} libITKIONRRD-${EXT} libITKIOPNG-${EXT} libITKIOStimulate-${EXT} libITKIOTIFF-${EXT} libITKIOVTK-${EXT} libITKNrrdIO-${EXT} libITKniftiio-${EXT} libITKznz-${EXT} libitkdouble-conversion-${EXT} libitkopenjpeg-${EXT} libitksys-${EXT} libitkv3p_netlib-${EXT} libitkvnl-${EXT} libitkvnl_algo-${EXT})

     set(MRI_VESSEL_SEGMENT_ITK_LIBS ${MRIS_COMMON_ITK_LIBS} libITKCommon-${EXT} libITKIOImageBase-${EXT} libITKLabelMap-${EXT})
     set(MRIS_CURVATURE2IMAGE_ITK_LIBS libitkv3p_netlib-${EXT} libitkvnl-${EXT} libitkvnl_algo-${EXT})
     set(MRIS_EXTRACT_VALUES_ITK_LIBS ${MRIS_COMMON_ITK_LIBS})
     set(MRIS_IMAGE2VTK_ITK_LIBS ${MRIS_COMMON_ITK_LIBS} libITKIOMeshBYU-${EXT} libITKIOMeshFreeSurfer-${EXT} libITKIOMeshGifti-${EXT} libITKIOMeshOBJ-${EXT} libITKIOMeshOFF-${EXT} libITKIOMeshVTK-${EXT} libITKgiftiio-${EXT})

     # no entries for MRIS_MESH_SUBDIVIDE_ITK_LIBS
     set(MRIS_MULTIMODAL_ITK_LIBS libITKCommon-${EXT} libITKNrrdIO-${EXT} libITKniftiio-${EXT} libITKznz-${EXT} libitkdouble-conversion-${EXT} libitkopenjpeg-${EXT} libitkv3p_netlib-${EXT} libitkvnl-${EXT} libitkvnl_algo-${EXT})

     set(MRIS_MULTIMODAL_SURFACE_PLACEMENT_ITK_LIBS ${MRIS_COMMON_ITK_LIBS})
     set(MRIS_PROFILECLUSTERING_ITK_LIBS ${MRIS_COMMON_ITK_LIBS} libITKCommon-${EXT} libITKIOImageBase-${EXT} libITKStatistics-${EXT})
     set(MRIS_SURF2VTK_ITK_LIBS libITKCommon-${EXT} libITKNrrdIO-${EXT} libITKniftiio-${EXT} libITKznz-${EXT} libitkdouble-conversion-${EXT} libitkopenjpeg-${EXT} libitkv3p_netlib-${EXT} libitkvnl-${EXT} libitkvnl_algo-${EXT})

     # reset and do not use default glob that pulls in every library
     set(VTK_LIBRARIES "")

     set(MRIS_COMMON_VTK_LIBS_1 libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1)
     set(MRIS_COMMON_VTK_LIBS_2 libvtkCommonCore-9.1.so.1)

     # no entries for MRI_VESSEL_SEGMENT_VTK_LIBS
     set(MRIS_CURVATURE2IMAGE_VTK_LIBS ${MRIS_COMMON_VTK_LIBS_1})
     set(MRIS_EXTRACT_VALUES_VTK_LIBS ${MRIS_COMMON_VTK_LIBS_2})
     set(MRIS_IMAGE2VTK_VTK_LIBS ${MRIS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkCommonMisc-9.1.so.1 libvtkFiltersCore-9.1.so.1 libvtkFiltersModeling-9.1.so.1 libvtkImagingCore-9.1.so.1 libvtkImagingGeneral-9.1.so.1 libvtkImagingMorphological-9.1.so.1)
     # no entries for MRIS_MESH_SUBDIVIDE_VTK_LIBS
     set(MRIS_MULTIMODAL_VTK_LIBS ${MRIS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersCore-9.1.so.1 libvtkFiltersModeling-9.1.so.1 libvtkFiltersPoints-9.1.so.1)
     set(MRIS_MULTIMODAL_SURFACE_PLACEMENT_VTK_LIBS ${MRIS_COMMON_VTK_LIBS_1} libvtkIOCore-9.1.so.1 libvtkIOLegacy-9.1.so.1)
     set(MRIS_PROFILECLUSTERING_VTK_LIBS ${MRIS_COMMON_VTK_LIBS_2}) 
     # no entries for MRIS_SURF2VTK_VTK_LIBS
  else()
     set(ITK_LIBRARIES "ITKCommon;ITKIOImageBase;ITKIOBMP;ITKIOBioRad;ITKIOBruker;ITKIOGDCM;ITKIOGE;ITKIONRRD;ITKIOTIFF;ITKIOVTK;ITKIOStimulate;ITKIOMINC;ITKIOGIPL;ITKIOHDF5;ITKIOJPEG;ITKIOLSM;ITKIOMRC;ITKIOMeta;ITKIOPNG;ITKIONIFTI")
     if( (HOST_OS MATCHES "Ubuntu24" AND NOT UB24_SYSLIBS) OR ((HOST_OS MATCHES "CentOS9" OR HOST_OS MATCHES "CentOS9stream" OR HOST_OS MATCHES "Rocky9" OR HOST_OS MATCHES "Ubuntu22") AND HOST_OS_GCC_VERSION MATCHES "11"))
       set(ITK_LIBRARIES "${ITK_LIBRARIES};ITKIOJPEG2000;ITKIOMeshVTK")
     elseif(HOST_OS MATCHES "Rocky8" AND HOST_OS_GCC_VERSION MATCHES "10")
       set(ITK_LIBRARIES "${ITK_LIBRARIES};ITKIOJPEG2000;ITKIOMeshVTK")
     endif()
  endif()

  set(TRACKIO ../freeview/track_io/TrackIO.cpp)

  add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-unused-variable -Wno-unused-local-typedefs -Wno-reorder)
  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Code
    ${CMAKE_SOURCE_DIR}/anatomicuts/Code
    ${CMAKE_SOURCE_DIR}/freeview/track_io
    ${CMAKE_SOURCE_DIR}/freeview
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    /space/erebus/2/users/vsiless/code/commong/alglib/src
  )

  add_executable(mris_multimodal mris_multimodal.cxx)
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_multimodal PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_multimodal PRIVATE ${UB24_VTK_BASE})
     target_link_options(mris_multimodal PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_MULTIMODAL_ITK_LIBS} ${MRIS_MULTIMODAL_VTK_LIBS})
        target_link_libraries(mris_multimodal ${LIB})
     endforeach()

     target_link_libraries(mris_multimodal utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(mris_multimodal utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS mris_multimodal DESTINATION bin)


  # also needs ITKStatistics
  add_executable(mris_profileClustering mris_profileClustering.cxx)
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_profileClustering PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_profileClustering PRIVATE ${UB24_VTK_BASE})
     target_link_options(mris_profileClustering PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_PROFILECLUSTERING_ITK_LIBS} ${MRIS_PROFILECLUSTERING_VTK_LIBS})
        target_link_libraries(mris_profileClustering ${LIB})
     endforeach()

     target_link_libraries(mris_profileClustering utils nifti)
  else()
     target_link_libraries(mris_profileClustering utils nifti ${OPENGL_LIBRARIES} "${ITK_LIBRARIES};ITKStatistics" ${VTK_LIBRARIES})
  endif() 
  install(TARGETS mris_profileClustering DESTINATION bin)

 
  add_executable(mris_surf2vtk mris_surf2vtk.cxx)
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_surf2vtk PRIVATE ${UB24_ITK_BASE})
     target_link_options(mris_surf2vtk PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_SURF2VTK_ITK_LIBS})
        target_link_libraries(mris_surf2vtk ${LIB})
     endforeach()

     target_link_libraries(mris_surf2vtk utils nifti)
  else()
     target_link_libraries(mris_surf2vtk utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS mris_surf2vtk DESTINATION bin)


  add_executable(mris_image2vtk mris_image2vtk.cxx) # Code/vtkNIFTIImageReader.cxx)
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_image2vtk PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_image2vtk PRIVATE ${UB24_VTK_BASE})
     target_link_options(mris_image2vtk PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_IMAGE2VTK_ITK_LIBS} ${MRIS_IMAGE2VTK_VTK_LIBS})
        target_link_libraries(mris_image2vtk ${LIB})
     endforeach()

     target_link_libraries(mris_image2vtk utils nifti)
  elseif( (HOST_OS MATCHES "Ubuntu24" AND NOT UB24_SYSLIBS ) OR ((HOST_OS MATCHES "CentOS9" OR HOST_OS MATCHES "CentOS9stream" OR HOST_OS MATCHES "Rocky9" OR HOST_OS MATCHES "Ubuntu22") AND HOST_OS_GCC_VERSION MATCHES "11"))
    # also needs ITKIOMeshOFF;ITKIOMeshOBJ;ITKIOMeshGifti;ITKIOMeshFreeSurfer;ITKIOMeshBYU for ITK5.4
    target_link_libraries(mris_image2vtk utils nifti "${ITK_LIBRARIES};ITKIOMeshOFF;ITKIOMeshOBJ;ITKIOMeshGifti;ITKIOMeshFreeSurfer;ITKIOMeshBYU" ${VTK_LIBRARIES} )
  elseif(HOST_OS MATCHES "Rocky8" AND HOST_OS_GCC_VERSION MATCHES "10")
    target_link_libraries(mris_image2vtk utils nifti "${ITK_LIBRARIES};ITKIOMeshOFF;ITKIOMeshOBJ;ITKIOMeshGifti;ITKIOMeshFreeSurfer;ITKIOMeshBYU" ${VTK_LIBRARIES} )
  else()
    target_link_libraries(mris_image2vtk utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES} )
  endif()
  install(TARGETS mris_image2vtk DESTINATION bin)


  #file(GLOB ALGLIB_SRC /space/erebus/2/users/vsiless/code/commong/alglib/src/*.cpp)
  #add_executable(sigmoidFit sigmoidFit.cxx ${ALGLIB_SRC})
  #target_link_libraries(sigmoidFit utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  #install(TARGETS sigmoidFit DESTINATION bin)

  ## also needs ITKLabelMap
  ## Does not build on all OS's ??
  if(HOST_OS MATCHES "Ubuntu24")
     add_executable(mri_vessel_segment mri_vessel_segment.cxx)
     target_link_directories(mri_vessel_segment PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mri_vessel_segment PRIVATE ${UB24_VTK_BASE})
     target_link_options(mri_vessel_segment PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRI_VESSEL_SEGMENT_ITK_LIBS} ${MRI_VESSEL_SEGMENT_VTK_LIBS})
        target_link_libraries(mri_vessel_segment ${LIB})
     endforeach()

     target_link_libraries(mri_vessel_segment utils nifti)
     install(TARGETS mri_vessel_segment DESTINATION bin)
  elseif(NOT ((HOST_OS MATCHES "CentOS9" OR HOST_OS MATCHES "CentOS9stream" OR HOST_OS MATCHES "Rocky9" OR (HOST_OS MATCHES "Ubuntu24" AND NOT UB24_SYSLIBS) OR HOST_OS MATCHES "Ubuntu22") AND HOST_OS_GCC_VERSION MATCHES "11"))
     # fails if some vtk packages are not installed/present
     add_executable(mri_vessel_segment mri_vessel_segment.cxx)
     target_link_libraries(mri_vessel_segment utils nifti "${ITK_LIBRARIES};ITKLabelMap" ${VTK_LIBRARIES})
     if(LINK_SHARED_AS_NEEDED)
        target_link_options(mri_vessel_segment PRIVATE "LINKER:--as-needed")
     endif()
     install(TARGETS mri_vessel_segment DESTINATION bin)
  endif()


  add_executable(mris_multimodal_surface_placement mris_multimodal_surface_placement.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_multimodal_surface_placement PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_multimodal_surface_placement PRIVATE ${UB24_VTK_BASE})
     target_link_options(mris_multimodal_surface_placement PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_MULTIMODAL_SURFACE_PLACEMENT_ITK_LIBS} ${MRIS_MULTIMODAL_SURFACE_PLACEMENT_VTK_LIBS})
        target_link_libraries(mris_multimodal_surface_placement ${LIB})
     endforeach()

     target_link_libraries(mris_multimodal_surface_placement utils nifti)
  else()
     target_link_libraries(mris_multimodal_surface_placement utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS mris_multimodal_surface_placement DESTINATION bin)


  add_executable(mris_curvature2image mris_curvature2image.cxx)
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_curvature2image PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_curvature2image PRIVATE ${UB24_VTK_BASE})
     target_link_options(mris_curvature2image PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_CURVATURE2IMAGE_ITK_LIBS} ${MRIS_CURVATURE2IMAGE_VTK_LIBS})
        target_link_libraries(mris_curvature2image ${LIB})
     endforeach()

     target_link_libraries(mris_curvature2image utils nifti)
  else()
     target_link_libraries(mris_curvature2image utils nifti ${VTK_LIBRARIES})
   endif()
  install(TARGETS mris_curvature2image DESTINATION bin)


  add_executable(mris_extract_values mris_extract_values.cxx)
  if(HOST_OS MATCHES "Ubuntu24")
     target_link_directories(mris_extract_values PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_extract_values PRIVATE ${UB24_VTK_BASE})
     target_link_options(mris_extract_values PRIVATE ${UB24_COMMON_LINK_OPTIONS})
     foreach(LIB ${MRIS_EXTRACT_VALUES_ITK_LIBS} ${MRIS_EXTRACT_VALUES_VTK_LIBS})
        target_link_libraries(mris_extract_values ${LIB})
     endforeach()

     target_link_libraries(mris_extract_values utils nifti)
  else()
     target_link_libraries(mris_extract_values utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS mris_extract_values DESTINATION bin)

  if(LINK_SHARED_AS_NEEDED)
     target_link_options(mris_multimodal PRIVATE "LINKER:--as-needed")
     target_link_options(mris_profileClustering PRIVATE "LINKER:--as-needed")
     target_link_options(mris_surf2vtk PRIVATE "LINKER:--as-needed")
     target_link_options(mris_image2vtk PRIVATE "LINKER:--as-needed")
     # target_link_options(mri_vessel_segment PRIVATE "LINKER:--as-needed")
     target_link_options(mris_multimodal_surface_placement PRIVATE "LINKER:--as-needed")
     target_link_options(mris_curvature2image PRIVATE "LINKER:--as-needed")
     target_link_options(mris_extract_values PRIVATE "LINKER:--as-needed")
  endif()

endif()
