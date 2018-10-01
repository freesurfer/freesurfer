project(qdecproject)

include_directories(${FS_INCLUDE_DIRS})

set(SOURCES
  QdecContrast.cxx
  QdecDataTable.cxx
  QdecFactor.cxx
  QdecGlmDesign.cxx
  QdecGlmFit.cxx
  QdecGlmFitResults.cxx
  QdecProject.cxx
  QdecSubject.cxx
)

add_library(qdecproject STATIC ${SOURCES})
target_link_libraries(qdecproject utils)
