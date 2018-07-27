# WXWidgets Find Module

# WXWidgets is only used by mris_decimate_gui, so it's not too necessary that we supply
# it as an external package for those that don't have it installed. It's also a very messy
# package to deal with, and it has dependencies all over the place. Because of that, this module
# is really only cutomized for a martinos build

# set up the default wxwidgets path on martinos machines 
set(MARTINOS_WXWidgets /usr/pubsw/packages/wxWidgets/wxGTK-2.8.9)
if(NOT WXWidgets_DIR AND EXISTS ${MARTINOS_WXWidgets})
  set(WXWidgets_DIR ${MARTINOS_WXWidgets})
endif()

# find the include dir
find_path(WXWidgets_INCLUDE_DIR HINTS ${WXWidgets_DIR} NAMES wx/wx.h PATH_SUFFIXES include/wx-2.8)

# find the wx gtk gl library
find_library(WXLIB HINTS ${WXWidgets_DIR} NAMES libwx_gtk2_gl-2.8.a PATH_SUFFIXES lib)

find_package_handle_standard_args(WXWidgets DEFAULT_MSG WXWidgets_INCLUDE_DIR WXLIB)

# the wx-config command is quite helpful in locating library dependencies, so we'll use the output
# of this when we link (but we'll use our own jpeg, tiff, and expat libraries)
set(SED_CMD "| sed 's/[^ ]*jpeg[^ ]*//g' | sed 's/[^ ]*tiff[^ ]*//g' | sed 's/[^ ]*expat[^ ]*//g'")
execute_process(COMMAND bash -c "${WXWidgets_DIR}/bin/wx-config --libs ${SED_CMD}" OUTPUT_VARIABLE WXCONFIG_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND bash -c "${WXWidgets_DIR}/bin/wx-config --libs gl ${SED_CMD}" OUTPUT_VARIABLE WXCONFIG_LIBSGL OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(WXCONFIG_LIBS)
separate_arguments(WXCONFIG_LIBSGL)
set(WXWidgets_LIBRARIES ${WXCONFIG_LIBS} ${WXCONFIG_LIBSGL} ${WXLIB})

# get wx flags
execute_process(COMMAND bash -c "${WXWidgets_DIR}/bin/wx-config --cxxflags" OUTPUT_VARIABLE WX_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(WXWidgets_FLAGS "${WX_CXX_FLAGS} -D__WXGTK20__")

# add some extra include dirs that wxwidgets requires
set(WXWidgets_INCLUDE_DIRS ${WXWidgets_INCLUDE_DIR}
  /usr/include/gtk-2.0
  /usr/lib/gtk-2.0/include
  /usr/lib64/gtk-2.0/include
  /usr/include/gdk-pixbuf-2.0
  /usr/include/atk-1.0
  /usr/include/cairo
  /usr/include/pango-1.0
)
