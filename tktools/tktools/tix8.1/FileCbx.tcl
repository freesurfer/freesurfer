# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: FileCbx.tcl,v 1.3.2.1 2001/11/03 06:43:50 idiscovery Exp $
#
# tixFileCombobox --
#
#	A combobox widget for entering file names, directory names, file
#	patterns, etc.
#
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.

# tixFileComboBox displays and accepts the DOS pathnames only. It doesn't
# recognize UNC file names or Tix VPATHS.
#
tixWidgetClass tixFileComboBox {
    -classname TixFileComboBox
    -superclass tixPrimitive
    -method {
	invoke
    }
    -flag {
	-command -defaultfile -directory -text
    }
    -forcecall {
	-directory
    }
    -configspec {
	{-defaultfile defaultFile DefaultFile ""}
	{-directory directory Directory ""}
	{-command command Command ""}
	{-text text Text ""}
    }
    -default {
    }
}

proc tixFileComboBox:InitWidgetRec {w} {
    upvar #0 $w data

    tixChainMethod $w InitWidgetRec

    if {![string comp $data(-directory) ""]} {
	set data(-directory) [tixFSPWD]
    }
}

proc tixFileComboBox:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget
    set data(w:combo) [tixComboBox $w.combo -editable true -dropdown true]
    pack $data(w:combo) -expand yes -fill both
}

proc tixFileComboBox:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings
    $data(w:combo) config -command [list tixFileComboBox:OnComboCmd $w]
}

proc tixFileComboBox:OnComboCmd {w args} {
    upvar #0 $w data

    set text [string trim [tixEvent value]]

    set fInfo [tixFSNorm [tixFSVPath $data(-directory)] \
	$text $data(-defaultfile) "" errorMsg]
    if {[info exists errorMsg]} {

    } else {
	tixSetSilent $data(w:combo) [lindex $fInfo 0]
	if {[string compare $data(-command) ""]} {
	    set bind(specs) {%V}
	    set bind(%V)    $fInfo
	    tixEvalCmdBinding $w $data(-command) bind $fInfo
	}
    }
}

proc tixFileComboBox:config-text {w val} {
    upvar #0 $w data

    tixSetSilent $data(w:combo) $val
}

proc tixFileComboBox:config-directory {w val} {
    upvar #0 $w data

    set data(-directory) [tixFSNormDir $val]
    return $data(-directory)
}

proc tixFileComboBox:invoke {w} {
    upvar #0 $w data

    $data(w:combo) invoke
}


