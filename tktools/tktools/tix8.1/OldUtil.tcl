# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: OldUtil.tcl,v 1.3.2.1 2001/11/03 07:16:33 idiscovery Exp $
#
# OldUtil.tcl -
#
#	This is an undocumented file.
#	   Are these features used in Tix : NO.
#	   Should I use these features    : NO.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc setenv {name args} {
    global env

    if {[llength $args] == 1} {
        return [set env($name) [lindex $args 0]]
    } else {
        if {[info exists env($ename)] == 0} {
            bgerror "Error in setenv: "
                    "environment variable \"$name\" does not exist"
        } else {
            return $env($name)
        }
    }
}
#----------------------------------------------------------------------
#
#
#           U T I L I T Y   F U N C T I O N S  F O R   T I X 
#
#
#----------------------------------------------------------------------

# RESET THE STRING IN THE ENTRY
proc tixSetEntry {entry string} {
    set oldstate [lindex [$entry config -state] 4]
    $entry config -state normal
    $entry delete 0 end
    $entry insert 0 $string
    $entry config -state $oldstate
}

# GET THE FIRST SELECTED ITEM IN A LIST
proc tixListGetSingle {lst} {
    set indices [$lst curselection]
    if {$indices != "" } {
	return [$lst get [lindex $indices 0]]
    } else {
	return ""
    }
}

#----------------------------------------------------------------------
# RECORD A DIALOG'S POSITION AND RESTORE IT THE NEXT TIME IT IS OPENED
#----------------------------------------------------------------------
proc tixDialogRestore {w {flag -geometry}} {
    global tixDPos

    if {[info exists tixDPos($w)]} {
	if {![winfo ismapped $w]} {
	    wm geometry $w $tixDPos($w)
	    wm deiconify $w
	}
    } elseif {$flag == "-geometry"} {
	update
	set tixDPos($w) [winfo geometry $w]
    } else {
	update
	set tixDPos($w) +[winfo rootx $w]+[winfo rooty $w]
    }
}
#----------------------------------------------------------------------
# RECORD A DIALOG'S POSITION AND RESTORE IT THE NEXT TIME IT IS OPENED
#----------------------------------------------------------------------
proc tixDialogWithdraw {w {flag -geometry}} {
    global tixDPos

    if {[winfo ismapped $w]} {
	if {$flag == "-geometry"} {
	    set tixDPos($w) [winfo geometry $w]
	} else {
	    set tixDPos($w) +[winfo rootx $w]+[winfo rooty $w]
	}
	wm withdraw $w
    }
}
#----------------------------------------------------------------------
# RECORD A DIALOG'S POSITION AND RESTORE IT THE NEXT TIME IT IS OPENED
#----------------------------------------------------------------------
proc tixDialogDestroy {w {flag -geometry}} {
    global tixDPos

    if {[winfo ismapped $w]} {
	if {$flag == "-geometry"} {
	    set tixDPos($w) [winfo geometry $w]
	} else {
	    set tixDPos($w) +[winfo rootx $w]+[winfo rooty $w]
	}
    }
    destroy $w
}

# Obsolete
#
proc tixQueryAppResource {name class default} {

    set value [option get . $name $class]
    if {$value == ""} {
	return $default
    } else {
	return $value
    }    
}
proc tixCreateToplevel {w {type -mapped}} {
    upvar #0 $w data

    toplevel $w
    wm minsize $w 0 0
    if {$type == "-withdrawn"} {
	wm withdraw $w
    }

    bind $w <Destroy>    [bind Toplevel <Destroy>]
    bind $w <Map>        [bind Toplevel <Map>]
    bind $w <Unmap>      [bind Toplevel <Unmap>]
    bind $w <Visibility> [bind Toplevel <Visibility>]
    bind $w <Destroy>    "+_tixToplevelDestroy $w"
    bind $w <Map>        "+_tixToplevelMap $w"
    bind $w <Unmap>      "+_tixToplevelUnmap $w"
    bind $w <Visibility> "+_tixToplevelVisibility $w"
}

proc _tixToplevelDestroy {w} {
    upvar #0 $w data

    unset data
}

proc _tixToplevelUnmap {w} {
    upvar #0 $w data

    foreach dlg $data(dialogs) {
	set data($dlg,geom) [winfo geometry $dlg]
	wm withdraw $dlg
    }
}

proc _tixToplevelMap {w} {
    upvar #0 $w data

    foreach dlg $data(dialogs) {
	wm geometry $dlg $data($dlg,geom)
	wm deiconify $dlg
    }
}

proc _tixToplevelVisibility {w} {
    upvar #0 $w data

    foreach dlg $data(dialogs) {
	raise $dlg $w
    }
}

proc tixCreateDialogShell {w {type -mapped}} {
    toplevel $w
    set parent [winfo parent $w]
    upvar #0 $parent data

    wm minsize $w 0 0
    wm withdraw $w
    update
    mwm transfor $w [winfo parent $w]
    lappend data(dialogs) $w
    bind $w <Destroy> "_tixDialogDestroy $w"

    if {$type != "-withdrawn"} {
	wm deiconify $w
    }
}

proc _tixDialogDestroy {w} {
    set parent [winfo parent $w]
    upvar #0 $parent data

    catch {unset $data($w,geom)}
}


proc _tixInitMainWindow {w} {
    upvar #0 $w data

    set data(dialogs) ""

    bind $w <Destroy>    +[bind Toplevel <Destroy>]
    bind $w <Map>        +[bind Toplevel <Map>]
    bind $w <Unmap>      +[bind Toplevel <Unmap>]
    bind $w <Visibility> +[bind Toplevel <Visibility>]
    bind $w <Destroy>    "+_tixToplevelDestroy $w"
    bind $w <Map>        "+_tixToplevelMap $w"
    bind $w <Unmap>      "+_tixToplevelUnmap $w"
    bind $w <Visibility> "+_tixToplevelVisibility $w"
}

# The "mwm" command comes from tkmwm, a cousin package of Tix
# If this wish does not support mwm, the following line prevent code
# that uses "mwm" from breaking.
#
if {[info commands mwm] == ""} {
    proc mwm {args} {}
}

#----------------------------------------------------------------------
# Automatically initialization call
#----------------------------------------------------------------------

# This has been disabled

if 0 {
    _tixInitMainWindow .
}

