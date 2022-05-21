# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: FileCmpt.tcl,v 1.1.1.1.2.1 2001/11/03 06:43:50 idiscovery Exp $
#
# FileCmpt.tcl --
#
#	File access portibility routines.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


# Internal file names
# (1) Idempotent: [tixFileIntName $intName] == $intName
# (2) Does not contain "~", "..", "."
# (3) All DOS type C:foo will be translated to absoulte path such as
#     /\C:\windows\foo
# (4) Does not contail trailing "/" or "\\" characters
#

proc tixFileResolveName {nativeName {defParent ""}} {
    if {$defParent != ""} {
	return [tixNativeName [tixFileIntName $nativeName [tixFileIntName $defParent]]]
    } else {
        return [tixNativeName [tixFileIntName $nativeName]]
    }
}

proc tixNSubFolder {parent sub} {
    return [tixNativeName [tixSubFolder \
	[tixFileIntName $parent] [tixFileIntName $sub]]]
}
