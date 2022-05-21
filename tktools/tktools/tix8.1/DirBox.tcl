# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: DirBox.tcl,v 1.2.2.1 2001/11/03 07:29:16 idiscovery Exp $
#
# DirBox.tcl --
#
#	Implements the tixDirSelectBox widget.
#
# 	   - overrides the -browsecmd and -command options of the
#	     HList subwidget
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

tixWidgetClass tixDirSelectBox {
    -classname TixDirSelectBox
    -superclass tixPrimitive
    -method {
    }
    -flag {
	-command -disablecallback -value
    }
    -configspec {
	{-command command Command ""}
	{-disablecallback disableCallback DisableCallback 0 tixVerifyBoolean}
	{-label label Label "Directory:"}
	{-value value Value ""}
    }
    -forcecall {
	-value -label
    }
    -default {
	{*combo*listbox.height 		5}
	{*combo.label.anchor		w}
	{*combo.labelSide		top}
	{*combo.history			true}
	{*combo.historyLimit		20}
    }
}

proc tixDirSelectBox:InitWidgetRec {w} {
    upvar #0 $w data
    tixChainMethod $w InitWidgetRec
}

proc tixDirSelectBox:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget
    set data(w:dircbx) [tixFileComboBox $w.dircbx]
    set data(w:dirlist)  [tixDirList $w.dirlist]

    pack $data(w:dircbx) -side top -fill x -padx 4 -pady 2
    pack $data(w:dirlist) -side top -fill both -expand yes -padx 4 -pady 2

    if {![string comp $data(-value) ""]} {
	set data(-value) [tixFSPWD]
    }
}

proc tixDirSelectBox:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    $data(w:dircbx) config -command "tixDirSelectBox:Cmd-DirCbx $w"
    $data(w:dirlist) config -command "tixDirSelectBox:Cmd-DirList $w"\
	-browsecmd "tixDirSelectBox:Browse-DirList $w"
}

#----------------------------------------------------------------------
# Incoming event: User
#----------------------------------------------------------------------

# User activates the FileComboBox
#
#
proc tixDirSelectBox:Cmd-DirCbx {w args} {
    upvar #0 $w data

    set fInfo [tixEvent value]
    set path [lindex $fInfo 0]

    if {![file exists $path]} {
	tk_dialog .tix_error "" "Directory \"$path\" does not exist." \
	    error 0 Ok
	$data(w:dircbx) config \
	    -text [tixFSDisplayName [tixFSNormDir $data(-value)]] \
	    -directory $data(-value)
	return

	#
	# The following code is not used because directories cannot be created
	# on Windows
	#

	# 1.1 Check for validity. The pathname cannot contain invalid chars
	#
	if {![tixFSIsValid $path]} {
	    tk_dialog .tix_error "Error" \
		"\"$path\" is not a valid directory name" \
		error 0 Ok
	    $data(w:dircbx) config \
		-text [tixFSDisplayName [tixFSNormDir $data(-value)]] \
		-directory $data(-value)
	    return
	}

	# 1.2 Prompt for creation
	#
	set choice [tk_dialog .tix_error "" \
	    "Directory \"$path\" does not exist. Do you want to create it?" \
	    question 1 Yes No]
	if {$choice == 1} {
	    $data(w:dircbx) config \
		-text [tixFSDisplayName [tixFSNormDir $data(-value)]] \
		-directory $data(-value)
	    return
	} else {
	    if {![tixFSCreateDirs $path]} {
		tk_dialog .tix_error "Error" \
		    "Cannot create directory \"$path\". Permission denied" \
		    error 0 Ok
		$data(w:dircbx) config \
		    -text [tixFSDisplayName [tixFSNormDir $data(-value)]] \
		    -directory $data(-value)
		return
	    }
	    tixDirSelectBox:SetValue $w $path 1 1
	}
    } elseif {![file isdirectory $path]} {
	# 2.1: Can't choose a non-directory file
	#
	tk_dialog .tix_error "Error" \
	    "\"$path\" is not a directory." \
	    error 0 Ok
	$data(w:dircbx) config \
	    -text [tixFSDisplayName [tixFSNormDir $data(-value)]] \
	    -directory $data(-value)
	return
    } else {
	# OK. It is an existing directory
	#
	tixDirSelectBox:SetValue $w $path 1 1
    }
}

# User activates the dir list
#
#
proc tixDirSelectBox:Cmd-DirList {w args} {
    upvar #0 $w data

    set dir $data(-value)
    catch {
	set dir [tixEvent flag V]
    }
    set dir [tixFSNormDir $dir]
    tixDirSelectBox:SetValue $w $dir 0 0
}

# User browses the dir list
#
#
proc tixDirSelectBox:Browse-DirList {w args} {
    upvar #0 $w data

    set dir $data(-value)
    catch {
	set dir [tixEvent flag V]
    }
    set dir [tixFSNormDir $dir]
    tixDirSelectBox:SetValue $w $dir 0 0
}

#----------------------------------------------------------------------
# Incoming event: Application
#----------------------------------------------------------------------
proc tixDirSelectBox:config-value {w value} {
    upvar #0 $w data
    set value [tixFSNormDir $value]

    tixDirSelectBox:SetValue $w $value 1 1
    return $value
}

proc tixDirSelectBox:config-label {w value} {
    upvar #0 $w data

    $data(w:dircbx) subwidget combo config -label $value
}

#----------------------------------------------------------------------
#
#			Internal functions
#
#----------------------------------------------------------------------

# Arguments:
#	callback:Bool	Should we invoke the the -command.
# 	setlist:Bool	Should we set the -value of the DirList subwidget.
#
proc tixDirSelectBox:SetValue {w dir callback setlist} {
    upvar #0 $w data

    set data(-value) $dir
    $data(w:dircbx) config -text [tixFSDisplayName $dir] \
	-directory [tixFSDisplayName $dir] 
    if {$setlist && [file isdirectory $dir]} {
	tixSetSilent $data(w:dirlist) $dir
    }

    if {$callback} {
	if {!$data(-disablecallback) && ![tixStrEq $data(-command) ""]} {
	    set bind(specs) {%V}
	    set bind(%V)    $data(-value)

	    tixEvalCmdBinding $w $data(-command) bind $data(-value)
	}
    }
}
