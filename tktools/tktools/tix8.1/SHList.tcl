# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: SHList.tcl,v 1.3.2.2 2001/12/09 02:54:02 idiscovery Exp $
#
# SHList.tcl --
#
#	This file implements Scrolled HList widgets
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

tixWidgetClass tixScrolledHList {
    -classname TixScrolledHList
    -superclass tixScrolledWidget
    -method {
    }
    -flag {
	 -highlightbackground -highlightcolor -highlightthickness
    }
    -configspec {
	{-highlightbackground -highlightBackground HighlightBackground #d9d9d9}
	{-highlightcolor -highlightColor HighlightColor black}
	{-highlightthickness -highlightThickness HighlightThickness 2}
    }
    -default {
	{.scrollbar			auto}
	{*f1.borderWidth		1}
	{*hlist.Background		#c3c3c3}
	{*hlist.highlightBackground	#d9d9d9}
	{*hlist.relief			sunken}
	{*hlist.takeFocus		1}
	{*Scrollbar.takeFocus		0}
    }
    -forcecall {
	-highlightbackground -highlightcolor -highlightthickness
    }
}

proc tixScrolledHList:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    set data(pw:f1) \
	[frame $w.f1 -takefocus 0]
    set data(w:hlist) \
	[tixHList $w.f1.hlist -bd 0 -takefocus 1 -highlightthickness 0]

    pack $data(w:hlist) -in $data(pw:f1) -expand yes -fill both -padx 0 -pady 0

    set data(w:hsb) \
	[scrollbar $w.hsb -orient horizontal -takefocus 0]
    set data(w:vsb) \
	[scrollbar $w.vsb -orient vertical -takefocus 0]
    
    set data(pw:client) $data(pw:f1)
}

proc tixScrolledHList:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    $data(w:hlist) config \
	-xscrollcommand "$data(w:hsb) set"\
	-yscrollcommand "$data(w:vsb) set"\
	-sizecmd [list tixScrolledWidget:Configure $w]

    $data(w:hsb) config -command "$data(w:hlist) xview"
    $data(w:vsb) config -command "$data(w:hlist) yview"

}

#----------------------------------------------------------------------
#
#		option configs
#----------------------------------------------------------------------
proc tixScrolledHList:config-takefocus {w value} {
    upvar #0 $w data
  
    $data(w:hlist) config -takefocus $value
}	

proc tixScrolledHList:config-highlightbackground {w value} {
    upvar #0 $w data

    $data(pw:f1) config -highlightbackground $value
}

proc tixScrolledHList:config-highlightcolor {w value} {
    upvar #0 $w data

    $data(pw:f1) config -highlightcolor $value
}

proc tixScrolledHList:config-highlightthickness {w value} {
    upvar #0 $w data

    $data(pw:f1) config -highlightthickness $value
}


#----------------------------------------------------------------------
#
#		Widget commands
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#
#		Private Methods
#----------------------------------------------------------------------
# virtual
#
proc tixScrolledHList:RepackHook {w} {
    upvar #0 $w data

if 0 {
    if {[tixGetBoolean [$data(w:hlist) cget -header]]} {
	set data(vsbPadY) [winfo height $data(w:hlist).tixsw:header]
    } else {
	set data(vsbPadY) 0
    }

    puts $data(vsbPadY)\ $data(w:hlist).tixsw:header
}
    tixChainMethod $w RepackHook
}
#----------------------------------------------------------------------
# virtual functions to query the client window's scroll requirement
#----------------------------------------------------------------------
proc tixScrolledHList:GeometryInfo {w mW mH} {
    upvar #0 $w data

    set extra [expr [$w.f1 cget -bd]+[$w.f1 cget -highlightthickness]]

    set mW [expr $mW - $extra*2]
    set mH [expr $mH - $extra*2]

    if {$mW < 1} {
	set mW 1
    }
    if {$mH < 1} {
	set mH 1
    }

    return [$data(w:hlist) geometryinfo $mW $mH]
}
