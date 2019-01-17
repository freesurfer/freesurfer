# Xt Find Module

find_path(Xt_INCLUDE_DIR HINTS /usr/X11 NAMES GL/gl.h PATH_SUFFIXES include)
find_library(Xt_LIBRARY HINTS /usr/X11 NAMES Xt PATH_SUFFIXES lib)
find_package_handle_standard_args(Xt DEFAULT_MSG Xt_INCLUDE_DIR Xt_LIBRARY)
set(Xt_LIBRARIES ${Xt_LIBRARY})
