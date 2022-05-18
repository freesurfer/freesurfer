# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: Verify.tcl,v 1.2.2.1 2001/11/03 07:26:10 idiscovery Exp $
#
# Verify.tcl --
#
#	Config option verification routines.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tixVerifyBoolean {val} {
    return [tixGetBoolean $val]
}

proc tixVerifyDirectory {val} {
    if {![file isdir $val]} {
	error "\"$val\" is not a directory"
    }
    return $val
}

