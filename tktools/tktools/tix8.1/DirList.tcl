# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: DirList.tcl,v 1.3.2.1 2001/11/03 06:43:50 idiscovery Exp $
#
# DirList.tcl --
#
#	Implements the tixDirList widget.
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

tixWidgetClass tixDirList {
    -classname TixDirList
    -superclass tixScrolledHList
    -method {
	chdir
    }
    -flag {
	 -browsecmd -command -dircmd -disablecallback 
	 -root -rootname -showhidden -value
    }
    -configspec {
	{-browsecmd browseCmd BrowseCmd ""}
	{-command command Command ""}
	{-dircmd dirCmd DirCmd ""}
	{-disablecallback disableCallback DisableCallback 0 tixVerifyBoolean}
	{-root root Root ""}
	{-rootname rootName RootName ""}
	{-showhidden showHidden ShowHidden 0 tixVerifyBoolean}
	{-value value Value ""}
    }
    -default {
	{.scrollbar			auto}
	{*borderWidth			1}
	{*hlist.background		#c3c3c3}
	{*hlist.indent			7}
	{*hlist.relief			sunken}
	{*hlist.height			10}
	{*hlist.width			20}
	{*hlist.padX			2}
	{*hlist.padY			0}
	{*hlist.wideSelection		0}
	{*hlist.drawBranch		0}
	{*hlist.highlightBackground	#d9d9d9}
	{*hlist.itemType		imagetext}
	{*hlist.takeFocus		1}
    }
    -forcecall {
	-value
    }
}

# Important data members:
#
# data(vpath)
#	The currently selected vpath. This internal variable is useful on
#	the Win95 platform, where an directory may correspond to more than
#	one node in the hierarchy. For example, C:\Windows\Desktop\Foo
#	can appead as "Desktop\Foo" and
#	"Desktop\My Computer\C:\Windows\Desktop\Foo". This variable tells us
#	which icon should we show given the same DOS pathname.
#

proc tixDirList:InitWidgetRec {w} {
    upvar #0 $w data

    tixChainMethod $w InitWidgetRec
}

proc tixDirList:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    $data(w:hlist) config \
	-separator [tixFSSep] \
	-selectmode "single"

    # We must creat an extra copy of these images to avoid flashes on
    # the screen when user changes directory
    #
#    set data(images) [image create compound -window $data(w:hlist)]
#    $data(images) add image -image [tix getimage act_fold]
#    $data(images) add image -image [tix getimage folder]
#    $data(images) add image -image [tix getimage openfold]
}

proc tixDirList:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    $data(w:hlist) config \
	-browsecmd [list tixDirList:Browse $w] \
	-command [list tixDirList:Command $w]

    if {[tixStrEq $data(-value) ""]} {
	set data(-value) [tixFSPWD]
    }
    if [catch {
	set data(vpath) [tixFSVPath [tixFSNormDir $data(-value)]]
    }] {
	set data(vpath) [tixFSVPath [tixFSNormDir [tixFSPWD]]]
    }
}

#----------------------------------------------------------------------
# Incoming-Events
#----------------------------------------------------------------------
proc tixDirList:Browse {w args} {
    upvar #0 $w data

    uplevel #0 set TRANSPARENT_GIF_COLOR [$data(w:hlist) cget -bg]
    set vpath [tixEvent flag V]
    set value [$data(w:hlist) info data $vpath]

    tixDirList:HighLight $w $vpath

    set data(vpath)  $vpath
    set data(-value) $value

    tixDirList:CallBrowseCmd $w $data(-value)
}

proc tixDirList:Command {w args} {
    upvar #0 $w data

    set vpath [tixEvent value]
    set value [$data(w:hlist) info data $vpath]
    set data(-value) $value

    tixDirList:LoadDir $w [tixFSNormDir $value] $vpath
    tixDirList:HighLight $w $vpath

    set data(vpath) $vpath
    tixDirList:CallCommand $w $data(-value)
}

#----------------------------------------------------------------------
# Outgoing-Events
#----------------------------------------------------------------------

proc tixDirList:CallBrowseCmd {w value} {
    upvar #0 $w data

    if {$data(-browsecmd) != ""} {
	set bind(specs) "%V"
	set bind(%V) $value
	tixEvalCmdBinding $w $data(-browsecmd) bind $value
    }
}

proc tixDirList:CallCommand {w value} {
    upvar #0 $w data

    if {$data(-command) != "" && !$data(-disablecallback)} {
	set bind(specs) "%V"
	set bind(%V) $value
	tixEvalCmdBinding $w $data(-command) bind $value
    }
}

#----------------------------------------------------------------------
# 		Directory loading
#----------------------------------------------------------------------
proc tixDirList:LoadDir {w {npath ""} {vpath ""}} {
    upvar #0 $w data

    tixBusy $w on $data(w:hlist)

    $data(w:hlist) delete all

    if {![string compare $npath ""]} {
	set npath [tixFSNormDir $data(-value)]
	set vpath [tixFSVPath $npath]
    }

    tixDirList:ListHierachy $w $npath $vpath
    tixDirList:ListSubDirs $w $npath $vpath

    tixWidgetDoWhenIdle tixBusy $w off $data(w:hlist)
}

proc tixDirList:ListHierachy {w dir vpath} {
    upvar #0 $w data
    uplevel #0 set TRANSPARENT_GIF_COLOR [$data(w:hlist) cget -bg]

    foreach p [tixFSSplit $vpath] {
	set vpath [lindex $p 0]
	set text  [lindex $p 1]
	set path  [lindex $p 2]

	$data(w:hlist) add $vpath -text $text -data $path \
	    -image [tix getimage openfold]
    }
}

proc tixDirList:ListSubDirs {w dir vpath} {
    upvar #0 $w data
    uplevel #0 set TRANSPARENT_GIF_COLOR [$data(w:hlist) cget -bg]

    $data(w:hlist) entryconfig $vpath \
	-image [tix getimage act_fold]

    foreach ent [tixFSListDir $vpath 1 0 0 $data(-showhidden)] {
	set vp   [lindex $ent 0]
	set name [lindex $ent 1]
	set path [lindex $ent 2]

	$data(w:hlist) add $vp -text $name -data $path \
	    -image [tix getimage folder]
    }
}

proc tixDirList:SetValue {w npath vpath {flag ""}} {
    upvar #0 $w data

    if {![string compare $flag reload] ||
	![$data(w:hlist) info exists $vpath]} {
    	tixDirList:LoadDir $w $npath $vpath
    }

    tixDirList:HighLight $w $vpath

    set data(-value) [tixFSDisplayName $npath]
    set data(vpath) $vpath
    tixDirList:CallCommand $w $data(-value)
}

proc tixDirList:HighLight {w vpath} {
    upvar #0 $w data

    if {![tixStrEq $data(vpath) $vpath]} {
	set old $data(vpath)

	if {[$data(w:hlist) info exists $old]} {
	    # Un-highlight the originally selected entry by changing its
	    # folder image

	    if {[$data(w:hlist) info children $old] == ""} {
		$data(w:hlist) entryconfig $old\
		    -image [tix getimage folder]
	    } else {
		$data(w:hlist) entryconfig $old\
		    -image [tix getimage openfold]
	    }
	}
    }

    # Highlight the newly selected entry
    #
    $data(w:hlist) entryconfig $vpath \
	-image [tix getimage act_fold]
    $data(w:hlist) anchor set $vpath
    $data(w:hlist) select clear
    $data(w:hlist) select set $vpath
    $data(w:hlist) see $vpath
}

#----------------------------------------------------------------------
# Config options
#----------------------------------------------------------------------
proc tixDirList:config-value {w value} {
    upvar #0 $w data

    tixDirList:chdir $w $value
    return $data(-value)
}

proc tixDirList:config-showhidden {w value} {
    upvar #0 $w data

    tixWidgetDoWhenIdle tixDirList:LoadDir $w
}

#----------------------------------------------------------------------
# Public methods
#----------------------------------------------------------------------
proc tixDirList:chdir {w value} {
    upvar #0 $w data

    set path [tixFSNormDir $value]
    tixDirList:SetValue $w $path [tixFSVPath $path]
}
