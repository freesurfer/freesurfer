# tkUtils.tcl (tku)

# tkuMakeMenu isMenuButton "Menu Name" {item...}
# item = { command   "Item Name" command                [group_name] }
# item = { radio     "Item Name" command variable value [group_name] }
# item = { check     "Item Name" command variable       [group_name] }
# item = { cascade   "Item Name" {item...}              [group_name] }
# item = { separator }

set kLabelFont -*-lucida-medium-r-normal-*-12-*-*-*-*-*-*-*
set kNormalFont -*-lucida-medium-r-normal-*-12-*-*-*-*-*-*-*

proc tkuLabelFont {} {
    global kLabelFont
    return $kLabelFont
}

proc tkuNormalFont {} {
    global kNormalFont
    return $kNormalFont
}

# tkuMakeActiveLabel
# -variable : variable for the value (required)
# -label : label preceding the value
# -width : width of the variable entry
proc tkuMakeActiveLabel { ifwTop args } {
    global kLabelFont

    # set default arguments for all fields
    set aArgs(-label) ""
    set aArgs(-width) 0
    set aArgs(-height) 0
    set aArgs(-font) $kLabelFont
    set aArgs(-labelwidth) 0
    set aArgs(-labelheight) 0

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-variable} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeActiveLabel: no $arg specified"
	    return
	}
    }

    frame $ifwTop

    if { $aArgs(-label) != "" } {
	label $ifwTop.lw \
		-font $aArgs(-font) \
		-text $aArgs(-label) \
		-width $aArgs(-labelwidth) \
		-height $aArgs(-labelheight) \
		-anchor e

	pack $ifwTop.lw \
		-side left \
		-expand yes \
		-fill x \
		-anchor e \
		-padx 5
    }

#    entry $ifwTop.ew \
#      -textvariable $aArgs(-variable) \
#      -width $aArgs(-width) \
#      -state disabled \
#      -relief flat

    label $ifwTop.ew \
	    -textvariable $aArgs(-variable) \
	    -font $aArgs(-font) \
	    -width $aArgs(-width) \
	    -height $aArgs(-height) \
	    -relief flat \
	    -anchor w

    pack $ifwTop.ew \
	    -side left \
	    -expand yes \
	    -fill x \
	    -anchor w
}

# tkuMakeNormalLabel
# -label : text for the label (reqd)
# -width : hard width of the label or 0 for flexible
# -wrap : the length at which to wrap text
# -font : font to use
# -anchor : side to anchor (n, e, se, etc)
proc tkuMakeNormalLabel { ifwTop args } {
    global kNormalFont

    # set default arguments for all fields
    set aArgs(-wrap) 0
    set aArgs(-width) 0
    set aArgs(-font) $kNormalFont
    set aArgs(-justify) left

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-label} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeNormalLabel: no $arg specified"
	    return
	}
    }

    frame $ifwTop
    
    label $ifwTop.lw \
	-font $aArgs(-font) \
	-text $aArgs(-label) \
	-width $aArgs(-width) \
	-justify left \
	-anchor $aArgs(-anchor) \
	-wraplength $aArgs(-wrap)

    if { $aArgs(-wrap) != 0 } {
	$ifwTop.lw configure -wraplength $aArgs(-wrap)
    }

    pack $ifwTop.lw \
      -side left \
      -anchor w
}

proc tkuMakeEntry { ifwTop args } {
    global kLabelFont
    
    # set default arguments for all fields
    set aArgs(-label) ""
    set aArgs(-width) 0
    set aArgs(-font) $kLabelFont
    set aArgs(-labelwidth) 0
    set aArgs(-setCmd) ""

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-variable} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeEntry: no $arg specified"
	    return
	}
    }

    frame $ifwTop

    if { $aArgs(-label) != "" } {
	
	label $ifwTop.lwLabel \
		-width $aArgs(-labelwidth) \
		-text $aArgs(-label) \
		-font $aArgs(-font)
	
	pack $ifwTop.lwLabel \
		-side left \
		-anchor w
    }

    entry $ifwTop.ewEntry \
      -textvariable $aArgs(-variable) \
      -width $aArgs(-width)
    
    pack $ifwTop.ewEntry \
      -side right \
      -anchor e \
      -expand yes \
      -fill x

    if { $aArgs(-setCmd) != "" } {
	bind $ifwTop.ewEntry <Return> "$aArgs(-setCmd) [set $aArgs(-variable)]"
    }

}


# tkuMakeMenu
# -menu : location of the new menu button widget (required)
# -label : the name of the menu button
# -items : list of items
proc tkuMakeMenu { args } {

    # Set menu items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-menu -label -items} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeMenu: no $arg specified"
	    return
	}
    }

    # for every letter in menu name..
    # if this letter is not in the underline list...
    # add the letter to the list.
    # underline this index.
    # if no letters left, underline nothing.
    
    menubutton $aArgs(-menu) \
	-text $aArgs(-label) \
	-menu $aArgs(-menu).mw
    #      -underline $nUnderline
    
    # start an underline list for this menu
    
    tkuAddMenuItemsToMenu $aArgs(-menu).mw $aArgs(-items)
}

proc tkuAddMenuItemsToMenu { isMenu ilMenuItems } {
    
    menu $isMenu
    set nItemNum 1
    
    foreach lItem $ilMenuItems {
	
	set sType [lindex $lItem 0]
	
	set bProcessed 0
	
	if { [string compare $sType "command"] == 0 } {
	    
	    # for every letter in item name..
	    # if this letter is not in the local underline list...
	    # add the letter to the list.
	    # underline this index.
	    # if no letters left, underline nothing.
	    
	    $isMenu add command \
		-label [lindex $lItem 1] \
		-command [lindex $lItem 2]
	    #              -underline $nUnderline
	    
	    set sGroupName [lindex $lItem 3]
	    if { [string compare $sGroupName ""] != 0 } {
		tkuAddItemToMenuGroup $sGroupName $isMenu $nItemNum
	    }
	    
	    set bProcessed 1
	}
	if { [string compare $sType "radio" ] == 0 } {
	    $isMenu add radio \
		-label [lindex $lItem 1] \
		-command [lindex $lItem 2] \
		-variable [lindex $lItem 3] \
		-value [lindex $lItem 4]
	    
	    set sGroupName [lindex $lItem 5]
	    if { [string compare $sGroupName ""] != 0 } {
		tkuAddItemToMenuGroup $sGroupName $isMenu $nItemNum
	    }
	    
	    set bProcessed 1
	}
	if { [string compare $sType "check" ] == 0 } {
	    $isMenu add check \
		-label [lindex $lItem 1] \
		-command [lindex $lItem 2] \
		-variable [lindex $lItem 3]
	    
	    set sGroupName [lindex $lItem 4]
	    if { [string compare $sGroupName ""] != 0 } {
		tkuAddItemToMenuGroup $sGroupName $isMenu $nItemNum
	    }
	    
	    set bProcessed 1
	}
	if { [string compare $sType "separator"] == 0 } {
	    $isMenu add separator
	    set bProcessed 1
	}
	if { [string compare $sType "cascade"] == 0 } {
	    $isMenu add cascade \
		-label [lindex $lItem 1] \
		-menu $isMenu.cmw$nItemNum
	    
	    set lCascadeItems [lindex $lItem 2]
	    tkuAddMenuItemsToMenu $isMenu.cmw$nItemNum $lCascadeItems
	    
	    set sGroupName [lindex $lItem 3]
	    if { [string compare $sGroupName ""] != 0 } {
		tkuAddItemToMenuGroup $sGroupName $isMenu $nItemNum
	    }
	    set bProcessed 1
	}
	
	if { $bProcessed == 0 } {
	    puts "Error!!!!! $sType not recognized"
	}
	
	if { $bProcessed == 1 } {
	    incr nItemNum
	}
    }
}

proc tkuAddItemToMenuGroup { isGroupName ifwMenuObject inMenuItemNum } {
    
    global glMenuGroups
    
    # add this menu / item pair to the list
    lappend glMenuGroups($isGroupName) "$ifwMenuObject $inMenuItemNum"
}

proc tkuSetMenuItemGroupStatus { isGroupName ibEnable } {
    
    global glMenuGroups
    
    # for each menu / item pair in the list
    foreach lMenuItemPair $glMenuGroups($isGroupName) {
	
	# first item is a menu button
	set mbwMenu [lindex $lMenuItemPair 0]
	
	# second item is a list of items
	set nMenuItem [lindex $lMenuItemPair 1]
	
	if { $ibEnable == 0 } {
	    if { [catch {$mbwMenu entryconfigure $nMenuItem -state disabled} sResult] } {
		set sType [$mbwMenu type $nMenuItem]
		puts "error, $isGroupName: $mbwMenu $nMenuItem $sType\n\t$sResult"
		
	    }
	} else {
	    if { [catch {$mbwMenu entryconfigure $nMenuItem -state normal} sResult] } {
		puts "error, $isGroupName: $mbwMenu $nMenuItem\n\t$sResult"
	    }
	}
    }
}



# replaces the percent symbols in a string with a substitution string:
# i.e. tkuDoSubPercent { %s1 "Hello %s1" "World" }
# returns "Hello World"
proc tkuDoSubPercent { isPercent isString isSubstitution } {
    if {![string match %* $isPercent]} {
	return $isString;
    }
    regsub -all {\\|&} $isSubstitution {\\\0} isSubstitution
    regsub -all $isPercent $isString $isSubstitution sResult
    return $sResult
}

proc tkuDoFileDlog { args } {

    global gDialog

    set lFields {1 2 3 4 5}
    set knWidth 400
    
    # set default arguments for all fields
    set aArgs(-title) "No title"
    set aArgs(-okCmd) "puts \$sFileName"
    set aArgs(-cancelCmd) ""
    foreach nField $lFields {
	set aArgs(-prompt$nField) ""
	set aArgs(-type$nField) "file"
	set aArgs(-note$nField) ""
	set aArgs(-default$nField) ""
	set sFileName$nField ""
    }
    
    # allow passed options to override defaults
    array set aArgs $args

    # dialog name. make it the name of the title, subbing dashes for spaces.
    regsub -all { } .wwDlogFile$aArgs(-title) {-} wwDialog

    # do percent substitutions in ok command for %s1 thru %s5
    foreach nField $lFields {
	set aArgs(-okCmd) \
	    [tkuDoSubPercent %s$nField $aArgs(-okCmd) \$sFileName$nField]
    }
    
    # if we can bring up the dialog
    if { [tkuCreateDialog $wwDialog "$aArgs(-title)" {-borderwidth 10}] } {
	
	# for each field...
	foreach nField $lFields {
	    
	    # create a variable for this prompt. even if we don't use this
	    # field, we'll need it later (ugh)
	    set fwPrompt$nField  $wwDialog.fwPrompt$nField
	    
	    # if they didn't enter a prompt, skip this field
	    if { [string match "$aArgs(-prompt$nField)" ""] == 1 } {
		continue;
	    }
	    
	    # switch on the type for this field and create the approriate
	    # selecter widget. bind it to sFileName[1..5]
	    switch $aArgs(-type$nField) {
		file { 
		    set sFileName$nField ""
		    tkuMakeFileSelector [set fwPrompt$nField] \
			-text "$aArgs(-prompt$nField)" \
			-variable sFileName$nField \
			-getdefaultdir $aArgs(-default$nField) 
		}
		dir { 
		    set sFileName$nField "";
		    tkm_MakeDirectorySelector [set fwPrompt$nField] \
			"$aArgs(-prompt$nField)" sFileName$nField \
			$aArgs(-default$nField) 
		}
		text { 
		    set sFileName$nField $aArgs(-default$nField)
		    tkm_MakeEntry [set fwPrompt$nField] \
			"$aArgs(-prompt$nField)" sFileName$nField 
		}
		default { continue; }
	    }
	    
	    # if they requested a note, make one, otherwise just an empty
	    # frame.
	    set fwNote $wwDialog.fwNote$nField
	    if { [string match "$aArgs(-note$nField)" ""] == 0 } {
		tkuMakeNormalLabel $fwNote \
		    -label "$aArgs(-note$nField)"
	    } else {
		frame $fwNote
	    }
	    
	    # pack this prompt and note
	    pack [set fwPrompt$nField] $fwNote \
		-side top       \
		-expand yes     \
		-fill x         \
		-padx 5         \
		-pady 5
	}
	
	# create the buttons. bind the ok command and the cancel command.
	# before calling the ok command, update the file selector variables.
	# TODO: nicer way to write this command string?
	tkuMakeCancelOKButtons $wwDialog.fwButtons $wwDialog \
	    -okCmd  "catch { tkuUpdateFileSelectorVariable $fwPrompt1 };\
                     catch { tkuUpdateFileSelectorVariable $fwPrompt2 };\
                     catch { tkuUpdateFileSelectorVariable $fwPrompt3 };\
                     catch { tkuUpdateFileSelectorVariable $fwPrompt4 };\
                     catch { tkuUpdateFileSelectorVariable $fwPrompt5 };\
                     $aArgs(-okCmd)"

	pack $wwDialog.fwButtons
	
	# after the next idle, the window will be mapped. set the min
	# width to our width and the min height to the mapped height.
	after idle [format {
	    update idletasks
	    wm minsize %s %d [winfo reqheight %s]
	} $wwDialog $knWidth $wwDialog] 
    }
}

# ============================================================ FILE SELECTORS

# tkuMakeFile
# -variable : variable name of the entry value (required)
# -text : text for the prompt
# -getdefaultdir : a function that returns the default directory
# -defaultdir : the default directory
proc tkuMakeFileSelector { ifwTop args } {

    # set default arguments for all fields
    set aArgs(-text) "Choose a file:"

    # Set menu items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-variable} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeFileSelector: no $arg specified"
	    return
	}
    }

    # create the frame
    frame $ifwTop -width 200
    
    upvar $aArgs(-variable) theVar 

    tixFileEntry $ifwTop.few \
	-label $aArgs(-text) \
	-labelside top \
	-variable $aArgs(-variable) \
	-options {
	    entry.expand yes
	    entry.fill x
	}
    
    # set the value of the field to the value of the variable
    $ifwTop.few config -value $theVar
    
    # set the default location 
    if {[info exists aArgs(-getdefaultdir)]} {
	$ifwTop.few filedialog subwidget fsbox \
	    configure -directory [eval $aArgs(-getdefaultdir)]
    }
    if {[info exists aArgs(-defaultdir)]} {
	$ifwTop.few filedialog subwidget fsbox \
	    configure -directory $aArgs(-defaultdir)
    }

    pack $ifwTop.few \
      -side left \
      -expand yes \
      -fill x
}

proc tkuUpdateFileSelectorVariable { ifwTop } {

    $ifwTop.few update
}

# tkuMakeDirectorySelector
# -variable : variable name of the entry value (required)
# -text : text for the prompt
# -getdefaultdir : a function that returns the default directory
# -defaultdir : the default directory
proc tkuMakeDirectorySelector { isFrame args } {
    
    # set default arguments for all fields
    set aArgs(-text) "Choose a directory:"

    # Set menu items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-variable} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeDirectorySelector: no $arg specified"
	    return
	}
    }

    frame $isFrame -width 200
    
    upvar $aArgs(-variable) theVar 

    tixFileEntry $isFrame.few \
	    -label $aArgs(-text) \
	    -labelside top \
	    -variable $aArgs(-variable) \
	    -dialogtype tixDirSelectDialog \
	    -options {
	entry.expand yes
	entry.fill x
    }
    
    # set the value of the field to the value of the variable
    $isFrame.few config -value $theVar

    # set the default location 
    if { [info exists aArgs(-getdefaultdir)] } {
	[$isFrame.few filedialog subwidget dirbox] \
	       subwidget dirlist configure -value [eval $aArgs(-getdefaultdir)]
    }
    if {[info exists aArgs(-defaultdir)]} {
	$ifwTop.few filedialog subwidget fsbox \
		configure -directory $aArgs(-defaultdir)
    }
    
    pack $isFrame.few \
	    -side left \
	    -expand yes \
	    -fill x
}

proc tkuUpdateDirectorySelectorVariable { isFrame } {

    $isFrame.few update;
}

# ============================================================ BUTTON GROUPS

# tkuMakeCancelOKButtons
# -okCmd : command to execute when OK button is pressed
# -cancelCmd : command to execute when Cancel button is pressed
proc tkuMakeCancelOKButtons { ifwTop iwwTop args } {

    global kLabelFont

    # Set menu items and make sure we have the ones we require,
    array set aArgs $args

    set okCmd ""
    if {[info exists aArgs(-okCmd)]} { set okCmd $aArgs(-okCmd) }
    set cancelCmd ""
    if {[info exists aArgs(-cancelCmd)]} { set cancelCmd $aArgs(-cancelCmd) }

    frame $ifwTop
    
    button $ifwTop.bwOK \
      -text "OK" \
      -command "$okCmd; tkuCloseDialog $iwwTop"

    button $ifwTop.bwCancel \
      -text "Cancel" \
      -command "$cancelCmd; tkuCloseDialog $iwwTop"

    bind $iwwTop <Return> \
      "$ifwTop.bwOK flash; $ifwTop.bwOK invoke"
    bind $iwwTop <Escape> \
      "$ifwTop.bwCancel flash; $ifwTop.bwCancel invoke"

    pack $ifwTop.bwOK $ifwTop.bwCancel \
      -side right \
      -padx 5 \
      -pady 5
}


# ============================================================ DIALOG BOXES

proc tkuCreateDialog { iwwTop isTitle iArgs } {

    global gDialog
    
    if [winfo exists $iwwTop] {
	switch -- [wm state $iwwTop] {
	    normal {
		raise $iwwTop 
	    }
	    withdrawn - iconic {
		wm deiconify $iwwTop
		catch {
		    wm geometry $iwwTop $gDialog(geo,$iwwTop)
		}
	    }
	}
	return 0
    } else {
	eval {toplevel $iwwTop} $iArgs
	wm title $iwwTop $isTitle
	return 1
    }
}

proc tkuCloseDialog { iwwTop } {
    
    global gDialog
    
    catch {
	set gDialog(geo,$iwwTop) [wm geometry $iwwTop]
	#  wm withdraw $iwwTop
	destroy $iwwTop
    }
}

# tkuCreateMsgDlog
# -msg : message for the dlog (required)
# -variable : alternatively, the variable containing text for the msg dlog
# -parent : parent window on which to center
# -width : width of the dlog
# -height : height of the dlog
proc tkuCreateMsgDlog { args } {
    
    # set default arguments for all fields
    set aArgs(-width) 200
    set aArgs(-height) 50
    set aArgs(-parent) ""
    set aArgs(-variable) ""
    
    array set aArgs $args
    foreach arg {-msg} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuCreateMsgDlog: no $arg specified"
	    return
	}
    }

    if { $aArgs(-variable) != "" } {
	upvar #0 $aArgs(-variable) msgVariable
	puts "msgVariable = $msgVariable"
    }

    set width $aArgs(-width)
    set height $aArgs(-height)
    toplevel .msgDlog -width $width -height $height -borderwidth 0
    
    set x 0
    set y 0
    if { $aArgs(-parent) != "" } {

	set parent $aArgs(-parent)
	set x [expr [winfo x $parent] + ([winfo width $parent] / 2) - \
		($width / 2)]
	set y [expr [winfo y $parent] + ([winfo height $parent] / 2) - \
		($height / 2)]
    }

    frame .msgDlog.fw
    if { $aArgs(-variable) != "" } {
	label .msgDlog.fw.lw -textvariable msgVariable
    } else {
	label .msgDlog.fw.lw -text $aArgs(-msg)
    }

    pack .msgDlog.fw.lw
    pack .msgDlog.fw -expand yes -fill both

    wm geometry .msgDlog ${width}x${height}+${x}+${y}
    update idletasks

    return .msgDlog
}

proc tkuDestroyMsgDlog { iTop } {

    destroy $iTop
}
