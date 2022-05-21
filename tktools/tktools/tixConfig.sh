# tixConfig.sh --
# $Id: tixConfig.sh.in,v 1.1 2000/11/03 02:28:58 idiscovery Exp $
#
# 
# This shell script (for sh) is generated automatically by Tix's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for Tcl extensions so that they don't have to figure this all
# out for themselves.
#
# The information in this file is specific to a single platform.
#
# SCCS: @(#) tclConfig.sh.in 1.20 97/07/01 11:40:19

# String to pass to linker to pick up the Tcl library from its
# build directory.
TIX_BUILD_LIB_SPEC='-L/usr/pubsw/packages/tcltktixblt/current/src/tix-8.1.4/unix/tk8.4 -ltix8.1.8.4'
