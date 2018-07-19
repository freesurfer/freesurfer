# TCLTKTIXBLT Find Module

if(NOT TCLTKTIXBLT_DIR)
  set(TCLTKTIXBLT_DIR ${FS_PACKAGES_DIR}/tcltktixblt/8.4.6)
endif()

find_path(TCLTKTIXBLT_INCLUDE_DIR PATHS ${TCLTKTIXBLT_DIR} NAMES tcl.h PATH_SUFFIXES include NO_DEFAULT_PATH)

# find the (many) libraries
foreach(LIB tix8.1.8.4 tk8.4 tcl8.4)
  find_library(tmp PATHS ${TCLTKTIXBLT_DIR} NAMES ${LIB} PATH_SUFFIXES lib NO_DEFAULT_PATH)
  set(TCLTKTIXBLT_LIBRARIES ${TCLTKTIXBLT_LIBRARIES} ${tmp})
  unset(tmp CACHE)  # this is necessary for find_library to work (plus it clears it from the cache)
endforeach()

find_package_handle_standard_args(TCLTKTIXBLT DEFAULT_MSG TCLTKTIXBLT_INCLUDE_DIR TCLTKTIXBLT_LIBRARIES)

# make sure the tcltktixblt path gets linked from the lib directory during install
symlink(${TCLTKTIXBLT_DIR} ${CMAKE_INSTALL_PREFIX}/lib/tcltktixblt)
