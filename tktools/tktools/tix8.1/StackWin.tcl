# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: StackWin.tcl,v 1.1.1.1.2.1 2001/11/03 07:23:17 idiscovery Exp $
#
# StackWin.tcl --
#
#	Similar to NoteBook but uses a Select widget to represent the pages.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

tixWidgetClass tixStackWindow {
    -classname TixStackWindow
    -superclass tixVStack
    -method {
    }
    -flag {
    }
    -configspec {
    }
}

proc tixStackWindow:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    set data(w:tabs) [tixSelect $w.tabs]

    # We can't use the packer because it will conflict with the
    # geometry management of the VStack widget.
    #
    tixManageGeometry $data(w:tabs) "tixVStack:ClientGeomProc $w"
}

proc tixStackWindow:add {w child args} {
    upvar #0 $w data

    set ret [eval tixChainMethod $w add $child $args]

    # Find out the -label option
    #
    tixForEach {flag value} $args {
	if {$flag == "-label"} {
	    set label $value
	}
    }

    $data(w:tabs) add $child -command "$w raise $child" \
	-text $label

    return $ret
}

proc tixStackWindow:raise {w child} {
    upvar #0 $w data

    $data(w:tabs) config -value $child

    tixChainMethod $w raise $child
}

proc tixStackWindow:Resize {w} {
    upvar #0 $w data

    # We have to take care of the size of the tabs so that 
    #
    set tW [winfo reqwidth  $data(w:tabs)]
    set tH [winfo reqheight $data(w:tabs)]

    tixMoveResizeWindow $data(w:tabs) $data(-ipadx) $data(-ipady) $tW $tH
    tixMapWindow $data(w:tabs)

    set data(pad-y1) [expr $tH + $data(-ipadx)]
    set data(minW)   [expr $tW + 2 * $data(-ipadx)]
    set data(minH)   [expr $tH + 2 * $data(-ipady)]

    # Now that we know data(pad-y1), we can chain the call
    #
    tixChainMethod $w Resize
}
