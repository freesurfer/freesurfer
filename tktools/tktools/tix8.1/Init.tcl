# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: Init.tcl,v 1.3.2.6 2002/12/02 04:02:51 idiscovery Exp $
#
# Init.tcl --
#
#	Initializes the Tix library and performes version checking to ensure
#	the Tcl, Tk and Tix script libraries loaded matches with the binary
#	of the respective packages.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

global tix_library
if {![info exists tix_library]} {
    set msg "Unable to find the Tix library directory TIX_LIBRARY"
    append msg "This probably means that Tix wasn't installed properly."
    error $msg
}
global auto_path
if {![tixStrEq $tix_library ""]} {
    lappend auto_path $tix_library [file join  $tix_library pref]
}
lappend auto_path [file dir [info name]]

proc __tixError {errorMsg} {
    error [concat $errorMsg \
       "Please check your TIX_LIBRARY environment variable and " \
       "your Tix installation."]
}

proc __tixInit {} {
    global tix tixPriv env tix_version tix_patchLevel tk_version tix_library
    global tcl_version

    if {[info exists tix(initialized)]} {
	return
    }
    if {[info command "@scope"] != ""} {
	set hasItcl 1
    } else {
	set hasItcl 0
    }

    # STEP 0: Version checking using the Tcl7.5 package mechanism. This is not
    #	      done if we are linked to Tcl 7.4.
    #
    if {[string compare [info command package] ""]} {
	if {![string comp [info command tixScriptVersion] ""] && 
		![auto_load tixScriptVersion]} {
	    __tixError [concat "Cannot determine version of Tix script " \
		"library. Requires version $tix_version."]
	}

	if {!$hasItcl} {
	    set pkgVersion  $tix_version.$tcl_version
	} else {
	    # The extra .1 indicates that the Tix binary is specially
	    # compiled for Itcl. This is necessary for the "package
	    # require" command to load in the correct shared library
	    # file.
	    set pkgVersion  $tix_version.$tcl_version.1
	}

	package provide Tix $pkgVersion
	if {[tixStrEq $tix_library ""]} {
	    package provide Tixsam $pkgVersion
	}
    }

    # STEP 1: Version checking
    #
    #
    package require -exact Tix $tix_version.$tcl_version

    # STEP 2: Initialize file compatibility modules
    #
    #
    if {[info exists tixPriv(isWindows)]} {
	tixInitFileCmpt:Win
    } elseif {[info exists env(WINDOWS_EMU_DEBUG)]} {
	tixInitFileCmpt:Win
	tixWinFileEmu
    } else {
	tixInitFileCmpt:Unix
    }

    # STEP 3: Initialize the Tix application context
    #
    #

    tixAppContext tix

    # STEP 4: Initialize the bindings for widgets that are implemented in C
    #
    #
    if {[string compare [info command tixHList] ""]} {
	tixHListBind
    }
    if {[string compare [info command tixTList] ""]} {
	tixTListBind
    }
    if {[string compare [info command tixGrid]  ""]} {
	tixGridBind
    }
    tixComboBoxBind
    tixControlBind
    tixFloatEntryBind
    tixLabelEntryBind
    tixScrolledGridBind
    tixScrolledListBoxBind

    global tcl_platform tcl_interactive
    if {$tcl_platform(platform) == "windows"} {
	if {[info exists tcl_interactive] && !$tcl_interactive && \
		[info exists env(TIX_CONSOLE)] && $env(TIX_CONSOLE) != "0"} {
	    # On Windows, initialize the console even if there is no script.
	    # The problem here is that frozen/wrapped exes never have a script.
	    # To invoke this,  simply set the environment variable TIX_CONSOLE
	    # to 1 if you want the console shown, and -1 if you want it hidden.
	    # after idle tixConsoleEvalAppend $tcl_interactive
	    set tcl_interactive 1
	    if {$env(TIX_CONSOLE) == "-1"} {after idle catch {console hide}}
	} else {
	    # To invoke this,  simply set the environment variable TIX_CONSOLE
	    # to 1 if you want the console shown, and -1 if you want it hidden.
	    # after idle tixConsoleEvalAppend $tcl_interactive
	    if {[info exists env(TIX_CONSOLE)] && $env(TIX_CONSOLE) == "-1"} {
		# You *must* use after idle
		after idle catch {console hide}
	    }
	}
    }

    # In the past, the interactive initialization file was inconsistent,
    # and on Windows, $env(HOME) is undefined or most users don't even
    # know where there HOME is (Profiles\User\Application Data\)!
    # So a site wide initialization file tixwishrc.tcl is now used,
    # which must be in the same directory as the executable. To restore
    # the past behaviour, simply add the following line to that file:
    #  if {[file isfile [set file ~/.tixwishrc]]} {source $file}

    set bindir [file dirname [info nameofexe]]
    if {[file isfile [set file [file join $bindir tixwishrc.tcl]]]} {
	global tcl_rcFileName
	set tcl_rcFileName $file
    }

    rename __tixError ""
    rename __tixInit ""
}

# tixWidgetClassEx --
#
#       This procedure is similar to tixWidgetClass, except it
#       performs a [subst] on the class declaration before evaluating
#       it. This gives us a chance to specify platform-specific widget
#       default without using a lot of ugly double quotes.
#
#       The use of subst'able entries in the class declaration should
#       be restrained to widget default values only to avoid producing
#       unreadable code.
#
# Arguments:
# name -	The name of the class to declare.
# classDecl -	Various declarations about the class. See documentation
#               of tixWidgetClass for details.

proc tixWidgetClassEx {name classDecl} {
    tixWidgetClass $name [uplevel [list subst $classDecl]]
}


