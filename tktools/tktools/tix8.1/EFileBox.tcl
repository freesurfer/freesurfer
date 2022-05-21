# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: EFileBox.tcl,v 1.3.2.1 2001/11/03 07:29:54 idiscovery Exp $
#
# EFileBox.tcl --
#
#	Implements the Extended File Selection Box widget.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


#
# ToDo
#   (1)	If user has entered an invalid directory, give an error dialog
#

tixWidgetClass tixExFileSelectBox {
    -classname TixExFileSelectBox
    -superclass tixPrimitive
    -method {
	filter invoke
    }
    -flag {
	-browsecmd -command -dialog -dir -dircmd -directory 
	-disablecallback -filetypes -pattern -selection -showhidden -value
    }
    -forcecall {
	-filetypes
    }
    -configspec {
	{-browsecmd browseCmd BrowseCmd ""}
	{-command command Command ""}
	{-dialog dialog Dialog ""}
	{-dircmd dirCmd DirCmd ""}
	{-directory directory Directory ""}
	{-disablecallback disableCallback DisableCallback 0 tixVerifyBoolean}
	{-filetypes fileTypes FileTypes ""}
	{-pattern pattern Pattern *}
	{-showhidden showHidden ShowHidden 0 tixVerifyBoolean}
	{-value value Value ""}
    }
    -alias {
	{-dir -directory}
	{-selection -value}
    }

    -default {
	{*dir.label 			{Directories:}}
	{*dir.editable 			true}
	{*dir.history 			true}
	{*dir*listbox.height 		5}
	{*file.label  			Files:}
	{*file.editable 		true}
	{*file.history 			false}
	{*file*listbox.height 		5}
	{*types.label 			{List Files of Type:}}
	{*types*listbox.height 		3}
	{*TixComboBox.labelSide 	top}
	{*TixComboBox*Label.anchor 	w}
	{*dir.label.underline 		0}
	{*file.label.underline		0}
	{*types.label.underline 	14}
	{*TixComboBox.anchor 		e}
	{*TixHList.height 		7}
	{*filelist*listbox.height 	7}
	{*hidden.wrapLength 		3c}
	{*hidden.justify 		left}
    }
}

proc tixExFileSelectBox:InitWidgetRec {w} {
    upvar #0 $w data
    global env

    tixChainMethod $w InitWidgetRec

    if {$data(-directory) == ""} {
	global env

	if {[info exists env(PWD)]} {
	    set data(-directory) $env(PWD)
	} else {
	    set data(-directory) [pwd]
	}
    }
    set data(oldDir)    ""
    set data(flag)      0
}


#----------------------------------------------------------------------
#		Construct widget
#----------------------------------------------------------------------
proc tixExFileSelectBox:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    # listbox frame
    set lf [frame $w.lf]

    # The pane that contains the two listboxes
    #
    set pane  [tixPanedWindow $lf.pane -orientation horizontal]
    set dpane [$pane add 1 -size 160]
    set fpane [$pane add 2 -size 160]

    $dpane config -relief flat
    $fpane config -relief flat

    # The File List Pane
    #
    set data(w:file)  [tixComboBox $fpane.file\
	-command "tixExFileSelectBox:Cmd-FileCombo $w"\
	-prunehistory true \
	-options { \
	    label.anchor w \
	}]
    set data(w:filelist) [tixScrolledListBox $fpane.filelist \
	-command "tixExFileSelectBox:Cmd-FileList $w 1"\
	-browsecmd "tixExFileSelectBox:Cmd-FileList $w 0"]
    pack $data(w:file)  -padx 8 -pady 4 -side top -fill x
    pack $data(w:filelist) -padx 8 -pady 4 -side top -fill both -expand yes

    # The Directory Pane
    #
    set data(w:dir)   [tixComboBox $dpane.dir \
	-command "tixExFileSelectBox:Cmd-DirCombo $w"\
	-prunehistory true \
	-options { \
	    label.anchor w \
	}]
    set data(w:dirlist) [tixDirList  $dpane.dirlist \
	-command "tixExFileSelectBox:Cmd-DirList $w"\
	-browsecmd "tixExFileSelectBox:Browse-DirList $w"]
    pack $data(w:dir)   -padx 8 -pady 4 -side top -fill x
    pack $data(w:dirlist) -padx 8 -pady 4 -side top -fill both -expand yes

    # The file types listbox
    #
    set data(w:types) [tixComboBox $lf.types\
	-command "tixExFileSelectBox:Cmd-TypeCombo $w" \
	-options { \
	    label.anchor w \
	}]

    pack $data(w:types)  -padx 12 -pady 4 -side bottom -fill x -anchor w
    pack $pane -side top -padx 4 -pady 4 -expand yes -fill both

    # Buttons to the right
    #
    set bf [frame $w.bf]
    set data(w:ok)     [button $bf.ok -text Ok -width 6 \
	-underline 0 -command [list tixExFileSelectBox:Ok $w]]
    set data(w:cancel) [button $bf.cancel -text Cancel -width 6 \
	-underline 0 -command [list tixExFileSelectBox:Cancel $w]]
    set data(w:hidden) [checkbutton $bf.hidden -text "Show Hidden Files"\
	-underline 0\
       	-variable [format %s(-showhidden) $w] -onvalue 1 -offvalue 0\
	-command [list tixExFileSelectBox:SetShowHidden $w]]

    pack $data(w:ok) $data(w:cancel) $data(w:hidden)\
	-side top -fill x -padx 6 -pady 3

    pack $bf -side right -fill y -pady 6
    pack $lf -side left -expand yes -fill both

    tixDoWhenMapped $w "tixExFileSelectBox:Map $w"

    if {$data(-filetypes) == ""} {
	$data(w:types) config -state disabled
    }
}


#----------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------
proc tixExFileSelectBox:config-showhidden {w value} {
    upvar #0 $w data

    set data(-showhidden) $value
    tixExFileSelectBox:SetShowHidden $w
}

# Update both DirList and {file list and dir combo}
#
#
proc tixExFileSelectBox:config-directory {w value} {
    upvar #0 $w data

    if {![tixIsAbsPath $value]} {
	return $data(-directory)
    }

    set data(-directory) [tixFSAbsPath $value]
    tixSetSilent $data(w:dirlist) $data(-directory) 
    tixSetSilent $data(w:dir) $data(-directory) 
    tixWidgetDoWhenIdle tixExFileSelectBox:LoadFiles $w reload

    return $data(-directory) 
}

proc tixExFileSelectBox:config-filetypes {w value} {
    upvar #0 $w data

    $data(w:types) subwidget listbox delete 0 end

    foreach name [array names data] {
	if {[string match type,* $name]} {
	    catch {unset data($name)}
	}
    }

    if {$value == ""} {
	$data(w:types) config -state disabled
    } else {
	$data(w:types) config -state normal
    
	foreach type $value {
	    $data(w:types) insert end [lindex $type 1]
	    set data(type,[lindex $type 1]) [lindex $type 0]
	}
	tixSetSilent $data(w:types) ""
    }
}

#----------------------------------------------------------------------
# MISC Methods
#----------------------------------------------------------------------
proc tixExFileSelectBox:SetShowHidden {w} {
    upvar #0 $w data

    $data(w:dirlist) config -showhidden $data(-showhidden)

    tixWidgetDoWhenIdle tixExFileSelectBox:LoadFiles $w reload
}

# User activates the dir combobox
#
#
proc tixExFileSelectBox:Cmd-DirCombo {w args} {
    upvar #0 $w data

    set dir [tixEvent flag V]
    if {![tixIsAbsPath $dir]} {
	return
    }
    set dir [tixFSAbsPath $dir]

    if {![file isdirectory $dir]} {
	return
    }

    $data(w:dirlist) config -value $dir
    set data(-directory) $dir
}

# User activates the dir list
#
#
proc tixExFileSelectBox:Cmd-DirList {w args} {
    upvar #0 $w data

    set dir $data(-directory)
    catch {
	set dir [tixEvent flag V]
    }
    set dir [tixFSAbsPath $dir]

    tixSetSilent $data(w:dir) $dir
    set data(-directory) $dir

    tixWidgetDoWhenIdle tixExFileSelectBox:LoadFiles $w noreload
}

# User activates the dir list
#
#
proc tixExFileSelectBox:Browse-DirList {w args} {
    upvar #0 $w data

    set dir [tixEvent flag V]
    tixExFileSelectBox:Cmd-DirList $w $dir
}

proc tixExFileSelectBox:IsPattern {w string} {
    foreach char [split $string ""] {
	if {$char == "*" || $char == "?" || $char == "\{"  || $char == "\[" } {
	    return 1
	}
    }
    return 0
}

proc tixExFileSelectBox:Cmd-FileCombo {w value} {
    upvar #0 $w data

    if {[tixEvent type] == "<Return>"} {
	tixExFileSelectBox:Ok $w
    }
}

proc tixExFileSelectBox:Ok {w} {
    upvar #0 $w data

    set value [string trim [$data(w:file) subwidget entry get]]
    if {$value == ""} {
	set value $data(-pattern)
    }
    tixSetSilent $data(w:file) $value

    if {[tixExFileSelectBox:IsPattern $w $value]} {
	set data(-pattern) $value
	tixWidgetDoWhenIdle tixExFileSelectBox:LoadFiles $w reload
    } else {
	if {[tixIsAbsPath $value]} {
	    set intName [tixFileIntName $value]
	} else {
	    set intName [tixSubFolder [tixFileIntName $data(-directory)] \
	    	[tixFileIntName $value]]
	}
	set data(-value) [tixNativeName $intName]
	tixExFileSelectBox:Invoke $w
    }
}

proc tixExFileSelectBox:Cancel {w} {
    upvar #0 $w data

    if {$data(-dialog) != ""} {
	eval $data(-dialog) popdown
    }
}

proc tixExFileSelectBox:Invoke {w} {
    upvar #0 $w data

    # Save some old history
    #
    $data(w:dir)  addhistory [$data(w:dir) cget -value]
    $data(w:file) addhistory $data(-pattern)
    $data(w:file) addhistory $data(-value)
    if {$data(-dialog) != ""} {
	eval $data(-dialog) popdown
    }
    if {$data(-command) != "" && !$data(-disablecallback)} {
	set bind(specs) "%V"
	set bind(%V) $data(-value)
	tixEvalCmdBinding $w $data(-command) bind $data(-value)
    }
}

proc tixExFileSelectBox:Cmd-FileList {w invoke args} {
    upvar #0 $w data

    set index [lindex [$data(w:filelist) subwidget listbox curselection] 0]
    if {$index == ""} {
	set index 0
    }

    set file [$data(w:filelist) subwidget listbox get $index]
    tixSetSilent $data(w:file) $file

    set data(-value) [tixNativeName [tixSubFolder \
	[tixFileIntName $data(-directory)] [tixFileIntName $file]]]

    if {$invoke == 1} {
	tixExFileSelectBox:Invoke $w
    } else {
	if {$data(-browsecmd) != ""} {
	    tixEvalCmdBinding $w $data(-browsecmd) "" $data(-value)
	}
    }
}

proc tixExFileSelectBox:Cmd-TypeCombo {w args} {
    upvar #0 $w data

    set value [tixEvent flag V]

    if {[info exists data(type,$value)]} {
	set data(-pattern) $data(type,$value)
	tixSetSilent $data(w:file) $data(-pattern)
	tixWidgetDoWhenIdle tixExFileSelectBox:LoadFiles $w reload
    }
}

proc tixExFileSelectBox:LoadFiles {w flag} {
    upvar #0 $w data

    if {$flag != "reload" && $data(-directory) == $data(oldDir)} {
	return
    }

    if {![winfo ismapped [winfo toplevel $w]]} {
	tixDoWhenMapped [winfo toplevel $w] \
	    "tixExFileSelectBox:LoadFiles $w $flag"
	return
    }

    set listbox [$data(w:filelist) subwidget listbox]
    $listbox delete 0 end

    set data(-value) ""

    tixBusy $w on [$data(w:dirlist) subwidget hlist]

    # wrap in a catch so you can't get stuck in a Busy state
    if {[catch {
	set intDir [tixFileIntName $data(-directory)]
	foreach name [tixListDir $intDir 0 1 0 \
		$data(-showhidden) $data(-pattern)] {
	    
	    $listbox insert end [tixFileDisplayName [tixSubFolder $intDir $name]]
	}

	if {$data(oldDir) != $data(-directory)} {
	    # Otherwise if the user has already selected a file and then presses
	    # "show hidden", the selection won't be wiped out.
	    tixSetSilent $data(w:file) $data(-pattern)
	}
    } err]} {
	tixDebug "tixExFileSelectBox:LoadFiles error for $w\n$err"
    }
    set data(oldDir) $data(-directory)

    tixWidgetDoWhenIdle tixBusy $w off [$data(w:dirlist) subwidget hlist]
}

#
# Called when thd listbox is first mapped
proc tixExFileSelectBox:Map {w} {
    if {![winfo exists $w]} {
	return
    }
    upvar #0 $w data

    set bind(specs) "%V"
    set bind(%V) $data(-value)
    tixEvalCmdBinding $w bind \
	"tixExFileSelectBox:Cmd-DirList $w" $data(-directory)
}

#----------------------------------------------------------------------
# Public commands
#
#----------------------------------------------------------------------
proc tixExFileSelectBox:invoke {w} {
    tixExFileSelectBox:Invoke $w
}

proc tixExFileSelectBox:filter {w} {
    tixExFileSelectBox:LoadFiles $w reload
}

