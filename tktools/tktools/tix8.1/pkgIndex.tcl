# Tcl package index file, version 1.0
#
# $Id: pkgIndex.tcl.in,v 1.1.2.1 2002/11/11 07:36:50 idiscovery Exp $
#

package ifneeded Tix 8.1.8.4 [list load "[file join [file dirname $dir] libtix8.1.8.4.so]" Tix]
package ifneeded Tixsam 8.1.8.4 [list load "[file join [file dirname $dir] libtixsam8.1.8.4.so]" Tixsam]
package ifneeded wm_default 1.0 [list source [file join $dir pref WmDefault.tcl]]
