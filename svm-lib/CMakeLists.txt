project(svm)

set(SOURCES
  svm-data-param.cpp
  svm-kernel.cpp
  svm-lib-c.cpp
  svm-lib.cpp
  svm-model.cpp
  svm-param.cpp
  svm-similarity.cpp
  svm-sketch.cpp
)

if(SUPPRESS_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-value")
endif()

add_library(svm STATIC ${SOURCES})
