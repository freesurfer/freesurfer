# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: DirTree.tcl,v 1.2.2.2 2002/01/24 10:08:58 idiscovery Exp $
#
# DirTree.tcl --
#
#	Implements directory tree for Unix file systems
#
#       What the indicators mean:
#
#	(+): There are some subdirectories in this directory which are not
#	     currently visible.
#	(-): This directory has some subdirectories and they are all visible
#
#      none: The dir has no subdirectori(es).
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

tixWidgetClass tixDirTree {
    -classname TixDirTree
    -superclass tixVTree
    -method {
	activate chdir refresh
    }
    -flag {
	-browsecmd -command -directory -disablecallback -showhidden -value
    }
    -configspec {
	{-browsecmd browseCmd BrowseCmd ""}
	{-command command Command ""}
	{-disablecallback disableCallback DisableCallback 0 tixVerifyBoolean}
	{-showhidden showHidden ShowHidden 0 tixVerifyBoolean}
	{-value value Value ""}
    }
    -alias {
	{-directory -value}
    }
    -default {
	{.scrollbar			auto}
	{*Scrollbar.takeFocus           0}
	{*borderWidth                   1}
	{*hlist.indicator               1}
	{*hlist.background              #c3c3c3}
	{*hlist.drawBranch              1}
	{*hlist.height                  10}
	{*hlist.highlightBackground      #d9d9d9}
	{*hlist.indent                  20}
	{*hlist.itemType                imagetext}
	{*hlist.padX                    3}
	{*hlist.padY                    0}
	{*hlist.relief                  sunken}
	{*hlist.takeFocus               1}
	{*hlist.wideSelection           0}
	{*hlist.width                   20}
    }
}

proc tixDirTree:InitWidgetRec {w} {
    upvar #0 $w data

    tixChainMethod $w InitWidgetRec

    if {$data(-value) == ""} {
	global env
	if {[info exists env(PWD)]} {
	    set data(-value) $env(PWD)
	} else {
	    set data(-value) [pwd]
	}
    }

    tixDirTree:SetDir $w [tixFileIntName $data(-value)]
}

proc tixDirTree:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget
    tixDoWhenMapped $w "tixDirTree:StartUp $w"

    $data(w:hlist) config \
	-separator [tixDirSep] \
	-selectmode "single" -drawbranch 1

    # We must creat an extra copy of these images to avoid flashes on
    # the screen when user changes directory
    #
    set data(images) [image create compound -window $data(w:hlist)]
    $data(images) add image -image [tix getimage act_fold]
    $data(images) add image -image [tix getimage folder]
    $data(images) add image -image [tix getimage openfold]
}

proc tixDirTree:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

# %% do I still need this?
#   bind $data(w:hlist) <3> "tixDirTree:DeleteSib $w %x %y"
}

# This procedure is supposed to "trim" the directory tree view to
# just the current directory and its ancestors.
#
#proc tixDirTree:DeleteSib {w x y} {
#    upvar #0 $w data
#
#    set ent [$data(w:hlist) nearest $y]
#
#    if {$ent != ""} {
#	$data(w:hlist) anchor set $ent
#
#	for {set e $ent} {$e != "/"} {set e [$data(w:hlist) info parent $e]} {
#	    $data(w:hlist) delete siblings $e
#	}
#	tixDirTree:Browse $w $ent
#    }
#}

# %% This functions needs to be optimized
#
#
proc tixDirTree:HasSubDir {w dir} {
    upvar #0 $w data

    if {[tixListDir $dir 1 0 0 $data(-showhidden)] != ""} {
	return 1
    } else {
	return 0
    }
}


# Add one dir into the parent directory, sorted alphabetically
#
proc tixDirTree:AddToList {w dir parent name image} {
    upvar #0 $w data

    set added 0
    foreach sib [$data(w:hlist) info children $parent] {
	if {[string compare $dir $sib] < 0} {
	    $data(w:hlist) add $dir -before $sib -text $name -image $image
	    set added 1
	    break
	}
    }
    if {!$added} {
	$data(w:hlist) add $dir -text $name -image $image
    }

    if {[tixDirTree:HasSubDir $w $dir]} {
	tixVTree:SetMode $w $dir open
    }
}

# Add $dir and all ancestors of $dir into the HList widget
#
#
proc tixDirTree:AddAncestors {w dir} {
    upvar #0 $w data

    set path ""
    set parent ""
    foreach name [tixFileSplit $dir] {
	set path [tixSubFolder $path $name]
	if {![$data(w:hlist) info exists $path]} {
	    tixDirTree:AddToList $w $path $parent [tixFileDisplayName $path] \
		[tix getimage openfold]
	}
	set parent $path
    }
}

# Add all the sub directories of $dir into the HList widget
#
#
proc tixDirTree:ListDirs {w dir} {
    upvar #0 $w data
    uplevel #0 set TRANSPARENT_GIF_COLOR [$data(w:hlist) cget -bg]

    tixBusy $w on $data(w:hlist)

    foreach name [tixListDir $dir 1 0 0 $data(-showhidden)] {
	set subdir [tixSubFolder $dir $name]
	if {![$data(w:hlist) info exists $subdir]} {
	    tixDirTree:AddToList $w $subdir $dir [tixFileDisplayName $subdir] \
		[tix getimage folder]
	}
    }

    tixWidgetDoWhenIdle tixBusy $w off $data(w:hlist)
}

proc tixDirTree:LoadDir {w dir {mode toggle}} {
    if {![winfo exists $w]} {
	return
    }

    upvar #0 $w data
    uplevel #0 set TRANSPARENT_GIF_COLOR [$data(w:hlist) cget -bg]

    # Add the directory and set it to the active directory
    #
    if {![$data(w:hlist) info exists $dir]} {
	tixDirTree:AddAncestors $w $dir
    }
    $data(w:hlist) entryconfig $dir -image [tix getimage act_fold]

    if {$mode == "toggle"} {
	if {[$data(w:hlist) info children $dir] == ""} {
	    set mode expand
	} else {
	    set mode flatten
	}
    }

    if {$mode == "expand"} {
	tixDirTree:ListDirs $w $dir
	if {[$data(w:hlist) info children $dir] == ""} {
	    tixVTree:SetMode $w $dir none
	} else {
	    tixVTree:SetMode $w $dir close
	}
    } else {
	$data(w:hlist) delete offsprings $dir
	tixVTree:SetMode $w $dir open
    }
}

proc tixDirTree:ToggleDir {w value mode} {
    upvar #0 $w data

    tixDirTree:LoadDir $w $value $mode
    tixDirTree:CallCommand $w
}

proc tixDirTree:CallCommand {w} {
    upvar #0 $w data

    if {$data(-command) != "" && !$data(-disablecallback)} {
	set bind(specs) {%V}
	set bind(%V)    $data(-value)

	tixEvalCmdBinding $w $data(-command) bind $data(-value)
    }
}

proc tixDirTree:CallBrowseCmd {w ent} {
    upvar #0 $w data

    if {$data(-browsecmd) != "" && !$data(-disablecallback)} {
	set bind(specs) {%V}
	set bind(%V)    $data(-value)

	tixEvalCmdBinding $w $data(-browsecmd) bind [list $data(-value)]
    }
}

proc tixDirTree:StartUp {w} {
    if {![winfo exists $w]} {
	return
    }

    upvar #0 $w data

    tixDirTree:LoadDir $w $data(i-directory)
}

proc tixDirTree:ChangeDir {w value {forced 0}} {
    upvar #0 $w data

    if {!$forced && $data(i-directory) == $value} {
	return
    }
    uplevel #0 set TRANSPARENT_GIF_COLOR [$data(w:hlist) cget -bg]

    if {!$forced && [$data(w:hlist) info exists $value]} {
	# Set the old directory to "non active"
	#
	if {[$data(w:hlist) info exists $data(i-directory)]} {
	    $data(w:hlist) entryconfig $data(i-directory) \
		-image [tix getimage folder]
	}

	$data(w:hlist) entryconfig $value  \
		-image [tix getimage act_fold]

    } else {
	if {$forced} {
	    if {[$data(w:hlist) info children $value] == ""} {
		set mode flatten
	    } else {
		set mode expand
	    }
	} else {
	    set mode toggle
	}
	tixDirTree:LoadDir $w $value $mode
	tixDirTree:CallCommand $w
    }
    tixDirTree:SetDir $w $value
}


proc tixDirTree:SetDir {w intName} {
    upvar #0 $w data

    set data(i-directory) $intName
    set data(-value)  [tixNativeName $intName]
}

#----------------------------------------------------------------------
#
# Virtual Methods
#
#----------------------------------------------------------------------
proc tixDirTree:OpenCmd {w ent} {
    tixDirTree:ToggleDir $w $ent expand
    tixDirTree:ChangeDir $w $ent
    tixDirTree:CallBrowseCmd $w $ent
}

proc tixDirTree:CloseCmd {w ent} {
    tixDirTree:ToggleDir $w $ent flatten
    tixDirTree:ChangeDir $w $ent
    tixDirTree:CallBrowseCmd $w $ent
}

proc tixDirTree:Command {w B} {
    upvar #0 $w data
    upvar $B bind

    set ent [tixEvent flag V]
    tixChainMethod $w Command $B

    if {$data(-command) != ""} {
	tixEvalCmdBinding $w $data(-command) bind $ent
    }
}

# This is a virtual method
#
proc tixDirTree:BrowseCmd {w B} {
    upvar #0 $w data
    upvar $B bind
    
    set ent [tixEvent flag V]

#    if {[$data(w:hlist) indicator exist $ent] && 
#	[$data(w:hlist) info children $ent] == ""} {
#	
#	tixVTree:Activate $w $ent open
#   }

    if {[string index $ent 0] != "/"} {
        # This is a hack because %V may have been modified by
	# callbrowsecmd ....
        set ent [tixFileIntName $ent]
    } 
    tixDirTree:ChangeDir $w $ent
    tixDirTree:CallBrowseCmd $w $ent
}

#----------------------------------------------------------------------
#
# Public Methods
#
#----------------------------------------------------------------------
proc tixDirTree:chdir {w value} {
    tixDirTree:ChangeDir $w [tixFileIntName $value]
}

proc tixDirTree:refresh {w {dir ""}} {
    upvar #0 $w data

    if {$dir == ""} {
	set dir $data(-value)
    }

    tixDirTree:ChangeDir $w [tixFileIntName $dir] 1


    # Delete any stale directories that no longer exist
    #
    foreach sub [$data(w:hlist) info children [tixFileIntName $dir]] {
	if {![file exists [tixNativeName $sub]]} {
	    $data(w:hlist) delete entry $sub
	}
    }
}

proc tixDirTree:config-directory {w value} {
    tixDirTree:ChangeDir $w [tixFileIntName $value]
}
