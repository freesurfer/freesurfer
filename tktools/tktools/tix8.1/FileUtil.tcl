# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: FileUtil.tcl,v 1.1.1.1.2.1 2001/11/03 06:43:50 idiscovery Exp $
#
# FileUtil.tcl ---
#
#
#	Utility functions for filename handling.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tixResolveDir {dir} {
    set dir [tixFile tildesubst $dir]
    set dir [tixFile trimslash $dir]
    
    if {$dir == "/"} {
	return $dir
    }

    if {[string index $dir 0] != "/"} {
	# Isn't an absolute path
	#
	set appPWD [pwd]
	catch {
	    cd $dir
	    set dir [pwd]
	}
	cd $appPWD
	return $dir
    }

    set names [split $dir "/"]

    # Get rid of all "."
    set n /
    foreach name [lrange $names 1 end] {
	if {[string compare "." $name]} {
	    lappend n $name
	}
    }
    if {$n == "/"} {
	return /
    }

    # Get rid of all ".."
    #
    set list [tixCompressDotDot $n 0]

    if {$list == "/"} {
	return /
    }

    # General case
    #
    set dir ""
    foreach sub [lrange $list 1 end] {
	append dir /$sub
    }
    return $dir
}

proc tixCompressDotDot {list i} {
    set done 0

    while {1} {
	if {$i >= [llength $list]} {
	    return $list
	}

	if {[lindex $list $i] != ".."} {
	    incr i
	    continue
	}

	# We encounter a ".."
	#
	if {$i == 0} {
	    # Can't handle this
	    #
	    return ""
	}
	if {$i == 1} {
	    set l [lindex $list 0]
	    set list [concat $l [lrange $list 2 end]]
	    continue
	}

	set l [lrange $list 0 [expr $i-2]]
	set list [concat $l [lrange $list [expr $i+1] end]]
	incr i -1
    }
}
