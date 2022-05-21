# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: NoteBook.tcl,v 1.3.2.3 2002/01/24 10:08:58 idiscovery Exp $
#
# NoteBook.tcl --
#
#	tixNoteBook: NoteBook type of window.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

tixWidgetClass tixNoteBook {
    -classname TixNoteBook
    -superclass tixVStack
    -method {
    }
    -flag {
    }
    -configspec {
	{-takefocus takeFocus TakeFocus 0 tixVerifyBoolean} 
    }
    -default {
	{.nbframe.tabPadX	8}
	{.nbframe.tabPadY	5}
	{.nbframe.borderWidth	2}
	{*nbframe.relief	raised}
    }
}

proc tixNoteBook:InitWidgetRec {w} {
    upvar #0 $w data

    tixChainMethod $w InitWidgetRec

    set data(pad-x1) 0
    set data(pad-x2) 0
    set data(pad-y1) 20
    set data(pad-y2) 0
}

proc tixNoteBook:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget
    
    set data(w:top) [tixNoteBookFrame $w.nbframe -slave 1 -takefocus 1]
    set data(w:nbframe) $data(w:top)

    bind $data(w:top) <ButtonPress-1> "tixNoteBook:MouseDown $w %x %y"
    bind $data(w:top) <ButtonRelease-1> "tixNoteBook:MouseUp $w %x %y"

    bind $data(w:top) <B1-Motion> "tixNoteBook:MouseDown $w %x %y"

    bind $data(w:top) <Left>  "tixNoteBook:FocusNext $w prev"
    bind $data(w:top) <Right> "tixNoteBook:FocusNext $w next"

    bind $data(w:top) <Return> "tixNoteBook:SetFocusByKey $w"
    bind $data(w:top) <space>  "tixNoteBook:SetFocusByKey $w"
}

#----------------------------------------------------------------------
# Public methods
#----------------------------------------------------------------------
proc tixNoteBook:add {w child args} {
    upvar #0 $w data

    set ret [eval tixChainMethod $w add $child $args]

    set new_args ""
    tixForEach {flag value} $args {
	if {$flag != "-createcmd" && $flag != "-raisecmd"} {
	    lappend new_args $flag
	    lappend new_args $value
	}
    }

    eval $data(w:top) add $child $new_args

    return $ret
}

proc tixNoteBook:raise {w child} {
    upvar #0 $w data

    tixChainMethod $w raise $child

    if {[$data(w:top) pagecget $child -state] == "normal"} {
	$data(w:top) activate $child
    }
}

proc tixNoteBook:delete {w child} {
    upvar #0 $w data

    tixChainMethod $w delete $child
    $data(w:top) delete $child
}

#----------------------------------------------------------------------
# Private methods
#----------------------------------------------------------------------
proc tixNoteBook:Resize {w} {
    upvar #0 $w data

    # We have to take care of the size of the tabs so that 
    #
    set rootReq [$data(w:top) geometryinfo]
    set tW [lindex $rootReq 0]
    set tH [lindex $rootReq 1]

    set data(pad-x1) 2 
    set data(pad-x2) 2
    set data(pad-y1) [expr $tH + $data(-ipadx) + 1]
    set data(pad-y2) 2
    set data(minW)   [expr $tW]
    set data(minH)   [expr $tH]

    # Now that we know data(pad-y1), we can chain the call
    #
    tixChainMethod $w Resize
}

proc tixNoteBook:MouseDown {w x y} {
    upvar #0 $w data

    focus $data(w:top)

    set name [$data(w:top) identify $x $y]
    $data(w:top) focus $name
    set data(w:down) $name
}

proc tixNoteBook:MouseUp {w x y} {
    upvar #0 $w data

    #it could happen (using the tk/menu) that a MouseUp
    #proceeds without a MouseDown event!!
    if {! [info exists data(w:down)] || ! [info exists data(w:top)]} {
	return
    }
	
    set name [$data(w:top) identify $x $y]

    if {$name != "" && $name == $data(w:down) && [$data(w:top) pagecget $name -state] == "normal" } {
        $data(w:top) activate $name
        tixCallMethod $w raise $name
    } else {
        $data(w:top) focus ""
    }
}


#----------------------------------------------------------------------
#
# Section for keyboard bindings
#
#----------------------------------------------------------------------

proc tixNoteBook:FocusNext {w dir} {
    upvar #0 $w data

    if {[$data(w:top) info focus] == ""} {
	set name [$data(w:top) info active]
	$data(w:top) focus $name

	if {$name != ""} {
	    return
	}
    } else {
	set name [$data(w:top) info focus$dir]
 	$data(w:top) focus $name
   }
}

proc tixNoteBook:SetFocusByKey {w} {
    upvar #0 $w data

    set name [$data(w:top) info focus]

    if {$name != "" && [$data(w:top) pagecget $name -state] == "normal"} {
	tixCallMethod $w raise $name
	$data(w:top) activate $name
    }
}

#----------------------------------------------------------------------
# Automatic bindings for alt keys
#----------------------------------------------------------------------
proc tixNoteBookFind {w char} {
    set char [string tolower $char]

    foreach child [winfo child $w] {
	if {![winfo ismapped $w]} {
	    continue
	}
	switch [winfo class $child] {
	    {Toplevel} {
		continue
	    }
	    TixNoteBook {
		set nbframe [$child subwidget nbframe]
		foreach page [$nbframe info pages] {
		    set char2 [string index [$nbframe pagecget $page -label] \
			[$nbframe pagecget $page -underline]]
		    if {([string compare $char [string tolower $char2]] == 0)||
			($char == "")} {
			if {[$nbframe pagecget $page -state] != "disabled"} {
			    return [list $child $page]
			}
		    }
		}
	    }
	}
	# Well, this notebook doesn't match with the key, but maybe
	# it contains a "subnotebook" that will match ..
	set match [tixNoteBookFind $child $char]
	if {$match != ""} {
	    return $match
	}
    }
    return ""
}

proc tixTraverseToNoteBook {w char} {
    if {$char == ""} {
	return 0
    }
    if {![winfo exists $w]} {
	return 0
    }
    set list [tixNoteBookFind [winfo toplevel $w] $char]
    if {$list != ""} {
	[lindex $list 0] raise [lindex $list 1]
	return 1
    }
    return 0
}

#----------------------------------------------------------------------
# Set default class bindings
#----------------------------------------------------------------------

bind all <Alt-KeyPress> "+tixTraverseToNoteBook %W %A"
bind all <Meta-KeyPress> "+tixTraverseToNoteBook %W %A"

