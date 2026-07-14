project(mris_mesh_subdivide)

if(VTK_FOUND AND ITK_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

  add_compile_options(-Wno-deprecated)

  add_executable(mris_mesh_subdivide mris_mesh_subdivide.cxx)

  if(HOST_OS MATCHES "Ubuntu24")
     # set(UB24_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
     if(NOT UBUNTU_ARM64)
        set(UB24_ITK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(UB24_VTK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(EXT "5.3.so")
     else()
        # don't have native itk sysem libs on arm64
        set(UB24_ITK_BASE "/usr/pubsw/packages/itk/5.4_ubuntu24/lib")
        set(UB24_VTK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(EXT "5.4.so")
     endif()

     # reset and do not use default glob that pulls in every library
     set(VTK_LIBRARIES "")
     set(MRIS_MESH_SUBDIVIDE_VTK_LIBS libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1 libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersModeling-9.1.so.1)

     # reset and do not use default glob that pulls in every library
     set(ITK_LIBRARIES "")
     set(MRIS_MESH_SUBDIVIDE_ITK_LIBS libitkv3p_netlib-${EXT} libitkvnl-${EXT} libitkvnl_algo-${EXT})

     target_link_directories(mris_mesh_subdivide PRIVATE ${UB24_ITK_BASE})
     target_link_directories(mris_mesh_subdivide PRIVATE ${UB24_VTK_BASE})
     foreach(LIB ${MRIS_MESH_SUBDIVIDE_ITK_LIBS} ${MRIS_MESH_SUBDIVIDE_VTK_LIBS})
        target_link_libraries(mris_mesh_subdivide ${LIB})
     endforeach()

     target_link_libraries(mris_mesh_subdivide utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(mris_mesh_subdivide utils ${VTK_LIBRARIES})
  endif()

  if(LINK_SHARED_AS_NEEDED)
     target_link_options(mris_mesh_subdivide PRIVATE "LINKER:--as-needed")
  endif()

  install(TARGETS mris_mesh_subdivide DESTINATION bin)

endif()
