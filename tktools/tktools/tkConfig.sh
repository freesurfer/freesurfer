# tkConfig.sh --
# 
# This shell script (for sh) is generated automatically by Tk's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for Tk extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.
#
# RCS: @(#) $Id: tkConfig.sh.in,v 1.8 2002/07/28 03:23:08 mdejong Exp $

# Tk's version number.
TK_VERSION='8.4'
TK_MAJOR_VERSION='8'
TK_MINOR_VERSION='4'
TK_PATCH_LEVEL='.6'

# -D flags for use with the C compiler.
TK_DEFS=' -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -D_LARGEFILE64_SOURCE=1 -DTCL_WIDE_INT_IS_LONG=1 -DPEEK_XCLOSEIM=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TIME_H=1 -DTIME_WITH_SYS_TIME=1 -DHAVE_PW_GECOS=1 '

# Flag, 1: we built a shared lib, 0 we didn't
TK_SHARED_BUILD=1

# This indicates if Tk was build with debugging symbols
TK_DBGX=

# The name of the Tk library (may be either a .a file or a shared library):
TK_LIB_FILE='libtk8.4.so'

# Additional libraries to use when linking Tk.
TK_LIBS='-L/usr/X11R6/lib64 -lX11 -ldl  -lieee -lm'

# Top-level directory in which Tcl's platform-independent files are
# installed.
TK_PREFIX='/usr/pubsw/packages/tcltktixblt/current'

# Top-level directory in which Tcl's platform-specific files (e.g.
# executables) are installed.
TK_EXEC_PREFIX='/usr/pubsw/packages/tcltktixblt/current'

# -I switch(es) to use to make all of the X11 include files accessible:
TK_XINCLUDES='-I/usr/X11R6/include'

# Linker switch(es) to use to link with the X11 library archive.
TK_XLIBSW='-L/usr/X11R6/lib64 -lX11'

# -l flag to pass to the linker to pick up the Tcl library
TK_LIB_FLAG='-ltk8.4'

# String to pass to linker to pick up the Tk library from its
# build directory.
TK_BUILD_LIB_SPEC='-L/usr/pubsw/packages/tcltktixblt/current/src/tk8.4.6/unix -ltk8.4'

# String to pass to linker to pick up the Tk library from its
# installed directory.
TK_LIB_SPEC='-L/usr/pubsw/packages/tcltktixblt/current/lib -ltk8.4'

# Location of the top-level source directory from which Tk was built.
# This is the directory that contains a README file as well as
# subdirectories such as generic, unix, etc.  If Tk was compiled in a
# different place than the directory containing the source files, this
# points to the location of the sources, not the location where Tk was
# compiled.
TK_SRC_DIR='/usr/pubsw/packages/tcltktixblt/current/src/tk8.4.6'

# Needed if you want to make a 'fat' shared library library
# containing tk objects or link a different wish.
TK_CC_SEARCH_FLAGS='-Wl,-rpath,${LIB_RUNTIME_DIR}:/usr/X11R6/lib64'
TK_LD_SEARCH_FLAGS='-Wl,-rpath,${LIB_RUNTIME_DIR}:/usr/X11R6/lib64'

# The name of the Tk stub library (.a):
TK_STUB_LIB_FILE='libtkstub8.4.a'

# -l flag to pass to the linker to pick up the Tk stub library
TK_STUB_LIB_FLAG='-ltkstub8.4'

# String to pass to linker to pick up the Tk stub library from its
# build directory.
TK_BUILD_STUB_LIB_SPEC='-L/usr/pubsw/packages/tcltktixblt/current/src/tk8.4.6/unix -ltkstub8.4'

# String to pass to linker to pick up the Tk stub library from its
# installed directory.
TK_STUB_LIB_SPEC='-L/usr/pubsw/packages/tcltktixblt/current/lib -ltkstub8.4'

# Path to the Tk stub library in the build directory.
TK_BUILD_STUB_LIB_PATH='/usr/pubsw/packages/tcltktixblt/current/src/tk8.4.6/unix/libtkstub8.4.a'

# Path to the Tk stub library in the install directory.
TK_STUB_LIB_PATH='/usr/pubsw/packages/tcltktixblt/current/lib/libtkstub8.4.a'
