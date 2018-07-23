project(lineprof)

if(ITK_FOUND AND VTK_FOUND AND PETSC_FOUND)

  include_directories(
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    ${PETSC_INCLUDE_DIR}
  )

  # temporary macro to deal with BoundingBoxType
  # this should be removed once we updgrade everything
  add_definitions(-DHAVE_ITK45)

  set(SOURCES
    Tracer.cpp
    Preprocessor.cpp
    PetscSolver.cpp
    LineProf.cpp
  )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

  add_library(lineprof STATIC ${SOURCES})
  target_link_libraries(lineprof ${PETSC_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES})

endif()
