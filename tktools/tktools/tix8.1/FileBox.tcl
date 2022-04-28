# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: FileBox.tcl,v 1.3.2.2 2002/01/24 10:08:58 idiscovery Exp $
#
# FileBox.tcl --
#
#	Implements the File Selection Box widget.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


# ToDo
#   (1)	If user has entered an invalid directory, give an error dialog
#

tixWidgetClass tixFileSelectBox {
    -superclass tixPrimitive
    -classname  TixFileSelectBox
    -method {
	filter invoke
    }
    -flag {
	-browsecmd -command -dir -directory -disablecallback
	-grab -pattern -selection -value
    }
    -configspec {
	{-browsecmd browseCmd BrowseCmd ""}
	{-command command Command ""}
	{-directory directory Directory ""}
	{-disablecallback disableCallback DisableCallback 0 tixVerifyBoolean}
	{-grab grab Grab global}
	{-pattern pattern Pattern *}
	{-value value Value ""}
    }
    -alias {
	{-selection -value}
	{-dir -directory}
    }
    -forcecall {
	-value
    }
    -default {
	{.relief			raised}
	{*filelist*Listbox.takeFocus	true}
	{.borderWidth 			1}
	{*Label.anchor			w}
	{*Label.borderWidth		0}
	{*TixComboBox*scrollbar		auto}
	{*TixComboBox*Label.anchor	w}
	{*TixScrolledListBox.scrollbar	auto}
	{*Listbox.exportSelection	false}
	{*directory*Label.text  	"Directories:"}
	{*directory*Label.underline	0}
	{*file*Label.text		"Files:"}
	{*file*Label.underline		2}
	{*filter.label			"Filter:"}
	{*filter*label.underline	3}
	{*filter.labelSide		top}
	{*selection.label		"Selection:"}
	{*selection*label.underline	0}
	{*selection.labelSide		top}
    }
}


proc tixFileSelectBox:InitWidgetRec {w} {
    upvar #0 $w data
    global env

    tixChainMethod $w InitWidgetRec

    if {$data(-directory) == ""} {
	set data(-directory) [pwd]
    }
    if {$data(-pattern) == ""} {
	set data(-pattern) [tixFilePattern allFiles]
    }

    tixFileSelectBox:SetPat $w [tixFileIntName $data(-pattern)]
    tixFileSelectBox:SetDir $w [tixFileIntName $data(-directory)]

    set data(flag)      0
    set data(fakeDir)   0
}

#----------------------------------------------------------------------
#		Construct widget
#----------------------------------------------------------------------
proc tixFileSelectBox:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    set frame1 [tixFileSelectBox:CreateFrame1 $w]
    set frame2 [tixFileSelectBox:CreateFrame2 $w]
    set frame3 [tixFileSelectBox:CreateFrame3 $w]

    pack $frame1 -in $w -side top -fill x
    pack $frame3 -in $w -side bottom -fill x
    pack $frame2 -in $w -side top -fill both -expand yes
}

proc tixFileSelectBox:CreateFrame1 {w} {
    upvar #0 $w data

    frame $w.f1 -border 10
    tixComboBox $w.f1.filter -editable 1\
	-command [list $w filter] -anchor e \
	-options {
	    slistbox.scrollbar auto
	    listbox.height 5
	    label.anchor w
	}
    set data(w:filter) $w.f1.filter

    pack $data(w:filter) -side top -expand yes -fill both
    return $w.f1
}

proc tixFileSelectBox:CreateFrame2 {w} {
    upvar #0 $w data

    tixPanedWindow $w.f2 -orientation horizontal
    #     THE LEFT FRAME
    #-----------------------
    set dir [$w.f2 add directory -size 120]
    $dir config -relief flat
    label $dir.lab
    set data(w:dirlist) [tixScrolledListBox $dir.dirlist\
		       -scrollbar auto\
		       -options {listbox.width 4 listbox.height 6}]

    pack $dir.lab -side top -fill x -padx 10
    pack $data(w:dirlist) -side bottom -expand yes -fill both -padx 10

    #     THE RIGHT FRAME
    #-----------------------
    set file [$w.f2 add file -size 160]
    $file config -relief flat
    label $file.lab
    set data(w:filelist) [tixScrolledListBox $file.filelist \
		       -scrollbar auto\
		       -options {listbox.width 4 listbox.height 6}]

    pack $file.lab -side top -fill x -padx 10
    pack $data(w:filelist) -side bottom -expand yes -fill both -padx 10

    return $w.f2
}

proc tixFileSelectBox:CreateFrame3 {w} {
    upvar #0 $w data

    frame $w.f3 -border 10
    tixComboBox $w.f3.selection -editable 1\
	-command [list tixFileSelectBox:SelInvoke $w] \
	-anchor e \
	-options {
	    slistbox.scrollbar auto
	    listbox.height 5
	    label.anchor w
	}

    set data(w:selection) $w.f3.selection

    pack $data(w:selection) -side top -fill both

    return $w.f3
}

proc tixFileSelectBox:SelInvoke {w args} {
    upvar #0 $w data

    set event [tixEvent type]

    if {$event != "<FocusOut>" && $event != "<Tab>"} {
	$w invoke
    }
}

proc tixFileSelectBox:SetValue {w value} {
    upvar #0 $w data

    set data(i-value) $value
    set data(-value)  [tixNativeName $value 0]
}

proc tixFileSelectBox:SetDir {w value} {
    upvar #0 $w data

    set data(i-directory) $value
    set data(-directory)  [tixNativeName $value]
}

proc tixFileSelectBox:SetPat {w value} {
    upvar #0 $w data

    set data(i-pattern) $value
    set data(-pattern)  [tixNativeName $value 0]
}


#----------------------------------------------------------------------
#                           BINDINGS
#----------------------------------------------------------------------

proc tixFileSelectBox:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    tixDoWhenMapped $w "tixFileSelectBox:FirstMapped $w"

    $data(w:dirlist) config \
	-browsecmd [list tixFileSelectBox:SelectDir $w] \
	-command   "tixFileSelectBox:InvokeDir $w"

    $data(w:filelist) config \
	-browsecmd [list tixFileSelectBox:SelectFile $w] \
	-command   "tixFileSelectBox:InvokeFile $w"
}

#----------------------------------------------------------------------
#                           CONFIG OPTIONS
#----------------------------------------------------------------------
proc tixFileSelectBox:config-directory {w value} {
    upvar #0 $w data

    if {$value == ""} {
	set value [pwd]
    }
    tixFileSelectBox:SetDir $w [tixFileIntName $value]
    tixFileSelectBox:SetFilter $w $data(i-directory) $data(i-pattern)
    $w filter

    return $data(-directory)
}

proc tixFileSelectBox:config-pattern {w value} {
    upvar #0 $w data

    if {$value == ""} {
	set value [tixFilePattern allFiles]
    }

    tixFileSelectBox:SetPat $w [tixFileIntName $value]
    tixFileSelectBox:SetFilter $w $data(i-directory) $data(i-pattern)

    # Returning a value means we have overridden the value and updated
    # the widget record ourselves.
    #
    return $data(-pattern)
}

proc tixFileSelectBox:config-value {w value} {
    upvar #0 $w data

    tixFileSelectBox:SetValue $w [tixFileIntName $value]
    tixSetSilent $data(w:selection) $value

    return $data(-value)
}

#----------------------------------------------------------------------
#                    PUBLIC METHODS
#----------------------------------------------------------------------
proc tixFileSelectBox:filter {w args} {
    upvar #0 $w data

    $data(w:filter) popdown
    tixFileSelectBox:InterpFilter $w
    tixFileSelectBox:LoadDir $w
}

proc tixFileSelectBox:invoke {w args} {
    upvar #0 $w data

    if {[$data(w:selection) cget -value] !=
	[$data(w:selection) cget -selection]} {
	    # this will in turn call "invoke" again ...
	    #
	    $data(w:selection) invoke
	    return
    }
    
    # record the filter
    #
    set filter [tixFileSelectBox:InterpFilter $w]
    $data(w:filter) addhistory $filter

    # record the selection
    #
    set userInput [string trim [$data(w:selection) cget -value]]
    tixFileSelectBox:SetValue $w \
	[tixFileIntName $userInput $data(i-directory)]
    $data(w:selection) addhistory $data(-value)

    $data(w:filter) align
    $data(w:selection)  align

    if {$data(-command) != "" && !$data(-disablecallback)} {
	set bind(specs) "%V"
	set bind(%V) $data(-value)
	tixEvalCmdBinding $w $data(-command) bind $data(-value)
    }
}

#----------------------------------------------------------------------
#                    INTERNAL METHODS
#----------------------------------------------------------------------
# InterpFilter:
#	Interprets the value of the w:filter widget. 
#
# Side effects:
#	Changes the fields data(-directory) and data(-pattenn) 
#
proc tixFileSelectBox:InterpFilter {w {filter ""}} {
    upvar #0 $w data

    if {$filter == ""} {
	set filter [$data(w:filter) cget -selection]
	if {$filter == ""} {
	    set filter [$data(w:filter) cget -value]
	}
    }

    set i_filter [tixFileIntName $filter]

    if {[file isdir $filter]} {
	tixFileSelectBox:SetDir $w $i_filter
	tixFileSelectBox:SetPat $w [tixFilePattern allFiles]
    } else {
	set nDir [file dir $filter]
	if {$nDir == "" || $nDir == "."} {
	    tixFileSelectBox:SetDir $w [tixFileIntName $data(i-directory)]
	} else {
	    tixFileSelectBox:SetDir $w [tixFileIntName $nDir]
	}
	tixFileSelectBox:SetPat $w [tixFileIntName [file tail $filter]]
    }

    tixFileSelectBox:SetFilter $w $data(i-directory) $data(i-pattern)

    return $data(filter)
}

proc tixFileSelectBox:SetFilter {w dir pattern} {
    upvar #0 $w data

    set data(filter) [tixSubFolder $dir $pattern]
    tixSetSilent $data(w:filter) [tixNativeName $data(filter)]
}

proc tixFileSelectBox:LoadDirIntoLists {w} {
    upvar #0 $w data

    $data(w:dirlist)  subwidget listbox delete 0 end
    $data(w:filelist) subwidget listbox delete 0 end

    set dir $data(i-directory)

    # (1) List the directories
    #
    set isDrive 0
    catch {
	set nDir [tixNativeName $dir]
	if {[llength [file split $nDir]] == 1} {
	    set isDrive 1
	}
    }
    foreach name [tixListDir $dir 1 0 1 1] {
	if {![string compare ".." $name]} {
	    if $isDrive {
		continue
	    }
	}
	$data(w:dirlist) subwidget listbox insert end $name
    }

    # (2) List the files
    #
    # %% UNIX'ISM:
    # If the pattern is "*" force glob to list the .* files.
    # However, since the user might not
    # be interested in them, shift the listbox so that the "normal" files
    # are seen first
    #
    # NOTE: if we pass $pat == "" but with $showHidden set to true,
    #       tixListDir will list "* .*" in Unix. See the comment on top of
    #	    the tixListDir code.
    #
    if {[string compare $data(i-pattern) *] == 0} {
	set pat ""
    } else {
	set pat $data(i-pattern)
    }

    set top 0
    foreach name [tixListDir $dir 0 1 0 0 $pat] {
	$data(w:filelist) subwidget listbox insert end $name
	if {[string match .* $name]} {
	    incr top
	}
    }

    $data(w:filelist) subwidget listbox yview $top
}

proc tixFileSelectBox:LoadDir {w} {
    upvar #0 $w data

    tixBusy $w on [$data(w:dirlist) subwidget listbox]

    tixFileSelectBox:LoadDirIntoLists $w
    tixFileSelectBox:MkDirMenu $w

    if {[$data(w:dirlist) subwidget listbox size] == 0} {
	# fail safe, just in case the user has inputed an errnoeuos
	# directory
	$data(w:dirlist) subwidget listbox insert 0 ".."
    }

    tixWidgetDoWhenIdle tixBusy $w off [$data(w:dirlist) subwidget listbox]
}

# %% unimplemented
#
proc tixFileSelectBox:MkDirMenu {w} {
    upvar #0 $w data
}

# User single clicks on the directory listbox
#
proc tixFileSelectBox:SelectDir {w} {
    upvar #0 $w data

    if {$data(fakeDir) > 0} {
	incr data(fakeDir) -1
	$data(w:dirlist) subwidget listbox select clear 0 end
	$data(w:dirlist) subwidget listbox activate -1
	return
    }

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    set subdir [tixListboxGetCurrent [$data(w:dirlist) subwidget listbox]]
    if {$subdir == ""} {
	set subdir "."
    }

    tixFileSelectBox:SetFilter $w \
	[tixFileIntName [tixSubFolder $data(i-directory) $subdir]] \
	$data(i-pattern)
    set data(flag) 0
}

proc tixFileSelectBox:InvokeDir {w} {
    upvar #0 $w data

    set theDir [$data(w:dirlist) subwidget listbox get active]

    tixFileSelectBox:SetDir $w [tixFileIntName \
	[tixSubFolder $data(i-directory) $theDir]]

    $data(w:dirlist) subwidget listbox select clear 0 end

    tixFileSelectBox:SetFilter $w $data(i-directory) $data(i-pattern)
    tixFileSelectBox:InterpFilter $w [tixNativeName $data(filter)]

    tixFileSelectBox:LoadDir $w

    if {![tixEvent match <Return>]} {
	incr data(fakeDir) 1
    }
}

proc tixFileSelectBox:SelectFile {w} {
    upvar #0 $w data

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    # Reset the "Filter:" box to the current directory:
    #	
    $data(w:dirlist) subwidget listbox select clear 0 end
    tixFileSelectBox:SetFilter $w $data(i-directory) $data(i-pattern)

    # Now select the file
    #
    set selected [tixListboxGetCurrent [$data(w:filelist) subwidget listbox]]
    if {$selected != ""} {
	# Make sure that the selection is not empty!
	#
	tixFileSelectBox:SetValue $w \
	    [tixFileIntName [tixSubFolder $data(i-directory) $selected]]
	tixSetSilent $data(w:selection) $data(-value)

	if {$data(-browsecmd) != ""} {
	    tixEvalCmdBinding $w $data(-browsecmd) "" $data(-value)
	}
    }
    set data(flag) 0
}

proc tixFileSelectBox:InvokeFile {w} {
    upvar #0 $w data

    set selected [tixListboxGetCurrent [$data(w:filelist) subwidget listbox]]
    if {$selected  != ""} {
	$w invoke
    }
}

# This is only called the first this fileBox is mapped -- load the directory
#
proc tixFileSelectBox:FirstMapped {w} {
    if {![winfo exists $w]} {
	return
    }

    upvar #0 $w data

    tixFileSelectBox:SetFilter $w $data(i-directory) $data(i-pattern)
    tixFileSelectBox:LoadDir $w
    $data(w:filter) align
}


#----------------------------------------------------------------------
#
#
#              C O N V E N I E N C E   R O U T I N E S 
#
#
#----------------------------------------------------------------------

# This is obsolete. Use the widget tixFileSelectDialog instead
#
#
proc tixMkFileDialog {w args} {
    set option(-okcmd)    ""
    set option(-helpcmd)  ""

    tixHandleOptions option {-okcmd -helpcmd} $args

    toplevel $w
    wm minsize $w 10 10

    tixStdDlgBtns $w.btns
    
    if {$option(-okcmd) != ""} {
	tixFileSelectBox $w.fsb -command "wm withdraw $w; $option(-okcmd)"
    } else {
	tixFileSelectBox $w.fsb -command "wm withdraw $w"
    }

    $w.btns button ok     config -command "$w.fsb invoke"
    $w.btns button apply  config -command "$w.fsb filter" -text Filter
    $w.btns button cancel config -command "wm withdraw $w"

    if {$option(-helpcmd) == ""} {
	$w.btns button help config -state disabled
    } else {
	$w.btns button help config -command $option(-helpcmd)
    }
    wm protocol $w WM_DELETE_WINDOW "wm withdraw $w"
    pack $w.btns  -side bottom -fill both
    pack $w.fsb   -fill both -expand yes

    return $w.fsb
}


