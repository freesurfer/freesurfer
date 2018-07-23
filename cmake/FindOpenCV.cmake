# OpenCS Find Module

if(NOT OpenCV_DIR)
  set(OpenCV_DIR ${FS_PACKAGES_DIR}/opencv/2.2.0)
endif()

# find the include dir
find_path(OpenCV_INCLUDE_DIR HINTS ${OpenCV_DIR} NAMES opencv opencv2 PATH_SUFFIXES include)

# find the libraries
foreach(LIB opencv_core opencv_imgproc opencv_highgui opencv_ml)
  find_library(tmp HINTS ${OpenCV_DIR} NAMES ${LIB} PATH_SUFFIXES lib)
  set(OpenCV_LIBRARIES ${OpenCV_LIBRARIES} ${tmp})
  unset(tmp CACHE)  # this is necessary for find_library to work (plus it clears it from the cache)
endforeach()

find_package_handle_standard_args(OpenCV DEFAULT_MSG OpenCV_INCLUDE_DIR OpenCV_LIBRARIES)
