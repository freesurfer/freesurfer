project(svm)

add_library(svm STATIC
  svm-data-param.cpp
  svm-kernel.cpp
  svm-lib-c.cpp
  svm-lib.cpp
  svm-model.cpp
  svm-param.cpp
  svm-similarity.cpp
  svm-sketch.cpp
)
