project(lineprof)

if(ITK_FOUND AND VTK_FOUND AND PETSC_FOUND)

  include_directories(
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    ${PETSC_INCLUDE_DIR}
  )

  # temporary macro to deal with BoundingBoxType
  # todo: this should be removed once we updgrade everything
  # if(NOT APPLE_ARM64)
  add_definitions(-DHAVE_ITK45)
  # endif()

  set(SOURCES
    Tracer.cpp
    Preprocessor.cpp
    PetscSolver.cpp
    LineProf.cpp
  )

  add_compile_options(-Wno-deprecated)

  add_library(lineprof STATIC ${SOURCES})

  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
  endif()

  # if(FREEVIEW_LINEPROF)
  if(APPLE_ARM64 OR (HOST_OS MATCHES "Ubuntu22" AND HOST_OS_GCC_VERSION MATCHES "11"))
    message(STATUS "For HOST_OS ${HOST_OS} skipping lineprof")
  else()
     add_test_executable(lineprof_test LineProfTest.cpp)
     target_link_libraries(lineprof_test
       lineprof
       ${ZLIB_LIBRARIES}
       ${VTK_LIBRARIES}
       ${PETSC_LIBRARIES}
       ${LAPACK_LIBRARIES}
       ${BLAS_LIBRARIES}
       ${GFORTRAN_LIBRARIES}
       ${ITK_LIBRARIES}
       ${QUADMATH_LIBRARIES}
     )
  endif()

endif()
