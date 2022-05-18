# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: Compat.tcl,v 1.1.1.1.2.1 2001/11/03 06:43:50 idiscovery Exp $
#
# Compat.tcl --
#
# 	This file wraps around many incompatibilities from Tix 3.6
#	to Tix 4.0.
#
#	(1) "box" to "Box" changes
#	(2) "DlgBtns" to "ButtonBox" changes
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


proc tixDlgBtns {args} {
    return [eval tixButtonBox $args]
}

proc tixStdDlgBtns {args} {
    return [eval tixStdButtonBox $args]
}

proc tixCombobox {args} {
    return [eval tixComboBox $args]
}

proc tixFileSelectbox {args} {
    return [eval tixFileSelectBox $args]
}

proc tixScrolledListbox {args} {
    return [eval tixScrolledListBox $args]
}

proc tixInit {args} {
    eval tix config $args
    puts stderr "tixInit no longer needed for this version of Tix"
}
