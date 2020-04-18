##
## tkUtils.tcl
##
##
## Copyright (C) 2002-2011, CorTechs Labs, Inc. (La Jolla, CA) and
## The General Hospital Corporation (Boston, MA).
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer/CorTechs Software License Agreement' contained
## in the file 'license.cortechs.txt' found in the FreeSurfer distribution,
## and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferCorTechsLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##

# tkUtils.tcl (tku)

# tkuMakeMenu isMenuButton "Menu Name" {item...}
# item = { command   "Item Name" command                [group_name] }
# item = { radio     "Item Name" command variable value [group_name] }
# item = { check     "Item Name" command variable       [group_name] }
# item = { cascade   "Item Name" {item...}              [group_name] }
# item = { separator }

set kLabelFont -*-lucida-bold-r-normal-*-14-*-*-*-*-*-*-*
set kNormalFont -*-lucida-medium-r-normal-*-12-*-*-*-*-*-*-*
set kSmallFont -*-lucida-medium-r-normal-*-9-*-*-*-*-*-*-*

proc tkuLabelFont {} {
    global kLabelFont
    return $kLabelFont
}

proc tkuNormalFont {} {
    global kNormalFont
    return $kNormalFont
}

proc tkuSmallFont {} {
    global kSmallFont
    return $kSmallFont
}

proc tkuBalloonWait {} {
    return 500
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

    if { [info exists aArgs(-bg)] } {
	$ifwTop configure -bg $aArgs(-bg)
    }

    if { $aArgs(-label) != "" } {
	label $ifwTop.lw \
		-font $aArgs(-font) \
		-text $aArgs(-label) \
		-width $aArgs(-labelwidth) \
		-height $aArgs(-labelheight) \
		-anchor e
	
	if { [info exists aArgs(-bg)] } {
	    $ifwTop.lw configure -bg $aArgs(-bg)
	}

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
    
    if { [info exists aArgs(-bg)] } {
	$ifwTop.ew configure -bg $aArgs(-bg)
    }

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
    set aArgs(-anchor) w

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

proc tkuSetNormalLabelEnabled { ifwTop ibEnabled } {

    set state enabled
    if { $ibEnabled } { 
	set state normal 
    } else {
	set state disabled 
    }

    $ifwTop.lw configure -state $state
}

# The -notify action here is really cool. If you want to remind the
# user that they have to press return to trigger the -command action,
# use the -notify option. Whenever the value in the field is different
# from the initial value or the value when return was last pressed,
# the text willl become red. When return is pressed, the text goes
# back to black. HOWEVER, note that if you change the data in the
# variable linked to the entry outside of typing it into the entry,
# you must use the tkuRefreshEntryNotify command to tell the entry to
# update its saved value.

proc tkuMakeEntry { ifwTop args } {
    global kLabelFont
    global gaEntryValue
    
    # set default arguments for all fields
    set aArgs(-label) ""
    set aArgs(-width) 0
    set aArgs(-font) $kLabelFont
    set aArgs(-labelwidth) 0
    set aArgs(-labelanchor) w
    set aArgs(-command) ""
    set aArgs(-notify) 0

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
	    -anchor $aArgs(-labelanchor) \
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

    if { $aArgs(-notify) } {

	set gaEntryValue($ifwTop.ewEntry,variable) $aArgs(-variable)
	set gaEntryValue($ifwTop.ewEntry,value) [uplevel set $aArgs(-variable)]

	bind $ifwTop.ewEntry <KeyRelease> "tkuColorEntry $ifwTop.ewEntry"

	bind $ifwTop.ewEntry <Return> "$aArgs(-command); set gaEntryValue($ifwTop.ewEntry,value) \[set $aArgs(-variable)\]; tkuColorEntry $ifwTop.ewEntry"

    } else {

	bind $ifwTop.ewEntry <Return> "$aArgs(-command)"
    }

}

proc tkuSetEntryEnabled { ifwTop ibEnabled } {

    set state enabled
    if { $ibEnabled } { 
	set state normal 
    } else { 
	set state disabled 
    }

    catch { 
	$ifwTop.lwLabel configure -state $state
    }
    $ifwTop.ewEntry configure -state $state
}

proc tkuColorEntry { iew } {
    global gaEntryValue

    if { "$gaEntryValue($iew,value)" != "[$iew get]" } {
	$iew config -fg red
    } else {
	$iew config -fg black
    }
    
}

proc tkuRefreshEntryNotify { ifwTop } {
    global gaEntryValue

    set gaEntryValue($ifwTop.ewEntry,value) [uplevel set $gaEntryValue($ifwTop.ewEntry,variable)]

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
	
	switch $sType {

	    "command" - "radio" - "check" {
		
		set sNameAndAccel [lindex $lItem 1]
		set sName ""
		set sAccel ""
		if { [string first : $sNameAndAccel] != -1 } {
		    set sName [string range $sNameAndAccel \
				   0 [expr [string first : $sNameAndAccel] - 1]]
		    set sAccel [string range $sNameAndAccel \
				[expr [string first : $sNameAndAccel] + 1] end]
		} else {
		    set sName $sNameAndAccel
		}

		set sGroupName ""

		switch $sType { 
		    "command" {
			$isMenu add command \
			    -label $sName \
			    -command [lindex $lItem 2] \
			    -font [tkuNormalFont] \
			    -accelerator $sAccel

			set sGroupName [lindex $lItem 3]
		    }
		    "radio" {
			$isMenu add radio \
			    -label $sName \
			    -command [lindex $lItem 2] \
			    -variable [lindex $lItem 3] \
			    -value [lindex $lItem 4] \
			    -font [tkuNormalFont] \
			    -accelerator $sAccel
			
			set sGroupName [lindex $lItem 5]
		    }
		    "check" {
			$isMenu add check \
			    -label $sName \
			    -command [lindex $lItem 2] \
			    -variable [lindex $lItem 3] \
			    -font [tkuNormalFont] \
			    -accelerator $sAccel
			
			set sGroupName [lindex $lItem 4]
		    }
		}

		if { [string compare $sGroupName ""] != 0 } {
		    tkuAddItemToMenuGroup $sGroupName $isMenu $nItemNum
		}

		set bProcessed 1
	    }
	    "separator" {
		$isMenu add separator
		set bProcessed 1
	    }
	    "cascade" {
		$isMenu add cascade \
		    -label [lindex $lItem 1] \
		    -menu $isMenu.cmw$nItemNum \
		    -font [tkuNormalFont]      
		
		set lCascadeItems [lindex $lItem 2]
		tkuAddMenuItemsToMenu $isMenu.cmw$nItemNum $lCascadeItems
		
		set sGroupName [lindex $lItem 3]
		if { [string compare $sGroupName ""] != 0 } {
		    tkuAddItemToMenuGroup $sGroupName $isMenu $nItemNum
		}
		set bProcessed 1
	    }
	    default {}
	}

	if { $bProcessed == 0 } {
	    puts "Error!!!!! $sType not recognized"
	}
	
	if { $bProcessed == 1 } {
	    incr nItemNum
	}
    }
}

proc tkuSetMenuItemName { ifwMenu inIndex isName } {

    $ifwMenu.mw entryconfigure $inIndex -label $isName
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

# tkuMakeOptionMenu
# -entries : a list of text items to put into the menu (required)
# -command : command to call when an item is selected. appends the item index
# -menutext : the initial text for the menu button
# -label : an optional label to the left of the menu
# -labelwidth : optional width of the label
# -labelanchor : optional anchor of the label
# -font : optional font of the label
proc tkuMakeOptionMenu { ifwTop args } {

    # Set default arguments.
    set aArgs(-menutext) "Choose:"
    set aArgs(-label) ""
    set aArgs(-labelwidth) 0
    set aArgs(-labelanchor) e
    set aArgs(-font) [tkuNormalFont]

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-entries -command} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeLongOptionMenu: no $arg specified"
	    return
	}
    }
    
    frame $ifwTop

    if { $aArgs(-label) != "" } {
	
	label $ifwTop.lwLabel \
	    -anchor $aArgs(-labelanchor) \
	    -width $aArgs(-labelwidth) \
	    -text $aArgs(-label) \
	    -font $aArgs(-font)
	
	pack $ifwTop.lwLabel \
		-side left \
		-anchor e
    }

    # Create a raised frame for the button look.
    frame $ifwTop.fwmw \
	-relief raised \
	-border 2

    # Create the menu button that houses the menu.
    menubutton $ifwTop.fwmw.mbw \
	-text $aArgs(-menutext) \
	-menu $ifwTop.fwmw.mbw.mw \
	-indicatoron 1

    # Create the menu.
    set mw [menu $ifwTop.fwmw.mbw.mw]

    tkuPopulateOptionMenu $mw $mw.mbw $aArgs(-entries) $aArgs(-command)

    if { 0 } {
    # If we have more than 30 entries...
    set lEntries $aArgs(-entries)
    if { [llength $lEntries] > 30 } {

	# For each entry...
	set nEntry 0
	set nSubMenu 0
	while { $nEntry < [llength $lEntries] } {

	    # Get an entry 29 items down (or < 29 if we don't have
	    # that many items.
	    set nTopEntry [expr $nEntry + 29]
	    if { $nTopEntry >= [llength $lEntries] } {
		set nTopEntry [expr [llength $lEntries] - 1]
	    }

	    # Create a submenu. Add the submenu to the main menu,
	    # giving it a label consisting of the entry and the entry
	    # 29 items down.
	    menu $mw.mw$nSubMenu
	    $mw add cascade -menu $mw.mw$nSubMenu \
		-label "[lindex $lEntries $nEntry] -> [lindex $lEntries $nTopEntry]"

	    # Look at the entry 30 items from now.
	    incr nEntry 30
	    incr nSubMenu
	}
    }

    # For each entry...
    set nEntry 0
    foreach entry $lEntries {

	# If we have more than 30, we're adding it to one of our
	# submenus. Otherwise, we're adding to the main menu.
	if { [llength $lEntries] > 30 } {
	    set curMenu $mw.mw[expr $nEntry / 30]
	} else {
	    set curMenu $mw
	}

	# Add the item. The command sets the text of the menu button
	# to this entry's text, and calls the user command with the
	# entry index.
	$curMenu add command \
	    -command "$ifwTop.fwmw.mbw config -text [lindex $lEntries $nEntry]; $aArgs(-command) $nEntry" \
	    -label [lindex $lEntries $nEntry]

	incr nEntry
    }
}

    pack $ifwTop.fwmw.mbw \
	-side left \
	-anchor w
    
    pack $ifwTop.fwmw \
	-side left \
	-anchor w \
	-padx 5
    
}

proc tkuSetOptionMenuText { ifwTop {isText "Choose:"} } {

    # Catch the error if there is one (e.g. if this is not really an
    # option menu).
    set err [catch {
	$ifwTop.fwmw.mbw config -text $isText
    } sResult]
	
    if { $err != 0 } {
	puts "tkuSetOptionMenuText: Error setting text, maybe not really a menu?"
    }
}

proc tkuPopulateOptionMenu { iMenu iMenuButton ilEntries iCommand } {
    
    # Start fresh.
    $iMenu delete 0 end

    # If we have more than 30 entries...
    if { [llength $ilEntries] > 30 } {

	# For each entry...
	set nEntry 0
	set nSubMenu 0
	while { $nEntry < [llength $ilEntries] } {

	    # Get an entry 29 items down (or < 29 if we don't have
	    # that many items.
	    set nTopEntry [expr $nEntry + 29]
	    if { $nTopEntry >= [llength $ilEntries] } {
		set nTopEntry [expr [llength $ilEntries] - 1]
	    }

	    # Create a submenu. Add the submenu to the main menu,
	    # giving it a label consisting of the entry and the entry
	    # 29 items down.
	    # If it already exists, that's cool.
	    catch { menu $iMenu.mw$nSubMenu }
	    $iMenu.mw$nSubMenu delete 0 end
	    $iMenu add cascade -menu $iMenu.mw$nSubMenu \
		-label "[lindex $ilEntries $nEntry] -> [lindex $ilEntries $nTopEntry]"

	    # Look at the entry 30 items from now.
	    incr nEntry 30
	    incr nSubMenu
	}
    }

    # For each entry...
    set nEntry 0
    foreach entry $ilEntries {

	# If we have more than 30, we're adding it to one of our
	# submenus. Otherwise, we're adding to the main menu.
	if { [llength $ilEntries] > 30 } {
	    set curMenu $iMenu.mw[expr $nEntry / 30]
	} else {
	    set curMenu $iMenu
	}

	# Add the item. The command sets the text of the menu button
	# to this entry's text, and calls the user command with the
	# entry index.
	$curMenu add command \
	    -command "$iMenuButton config -text [lindex $ilEntries $nEntry]; $iCommand $nEntry" \
	    -label [lindex $ilEntries $nEntry]

	incr nEntry
    }
}

# tkuMakeCheckboxes
# -orientation : orientation of checkboxes (h, v)
# -font : overall font
# -labelwidth : overall label width
# -labelheight : overall label height
# -checkboxes : list of checkbox data (required)
#    -type : type of checkbox label (text, image)
#    -label : text label
#    -image : image label
#    -variable : variable to set (required)
#    -font : font for this checkbox label
#    -command : command to run when clicked
proc tkuMakeCheckboxes { ifwTop args } {

    # set default arguments for all fields
    set aArgs(-orientation) v
    set aArgs(-font) [tkuLabelFont]
    set aArgs(-labelwidth) 0
    set aArgs(-labelheight) 0

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-checkboxes} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeCheckboxes: no $arg specified"
	    return
	}
    }

    frame $ifwTop
    
    # for each checkbox...
    set nCheckbox 0
    foreach lCheckbox $aArgs(-checkboxes) {
	
	# Set arg items and make sure we have the ones we require,
	set aCheckbox(-type) text
	set aCheckbox(-label) ""
	set aCheckbox(-font) $aArgs(-font)
	set aCheckbox(-labelwidth) $aArgs(-labelwidth)
	set aCheckbox(-labelheight) $aArgs(-labelheight)
	set aCheckbox(-command) ""
	array set aCheckbox $lCheckbox
	foreach arg {-variable} {
	    if {![info exists aCheckbox($arg)]} {
    	     puts "tkuMakeCheckboxes: no $arg specified in checkbox $nCheckbox"
		return
	    }
	}

	# make names for the checkbox and the label.
	set cbw $ifwTop.cb$nCheckbox
	set lw  $ifwTop.lw$nCheckbox
	
	# text or image?
	switch $aCheckbox(-type) {
	    
	    text {
		
		# text. make a normal checkbox and label.
		checkbutton $cbw \
		    -variable $aCheckbox(-variable) \
		    -command $aCheckbox(-command)
		label $lw \
		    -font $aCheckbox(-font) \
		    -text $aCheckbox(-label)
		
		# if horizontal, pack all in the same row. if vertical.
		# pack the checkbox, than the label, in the
		# same row as the number of this checkbox.
		switch $aArgs(-orientation) {
		    h - x { 
			grid $cbw -column [expr 2 * $nCheckbox] -row 0
			grid $lw -column [expr 1 + [expr 2 *$nCheckbox]] \
			    -row 0 -sticky w
		    }
		    v - y { 
			grid $cbw -column 0 -row $nCheckbox
			grid $lw -column 1 -row $nCheckbox -sticky w
		    }
		}
	    }
	    
	    image { 
		# image. create a checkbox with an image label. no text label.
		checkbutton $cbw \
		    -image $aCheckbox(-image) \
		    -variable $aCheckbox(-variable) \
		    -command $aCheckbox(-command) \
		    -indicatoron false \
		    -selectcolor gray
		
		# if horizontal, pack in increasing columns. if vertical,
		# pack in increasing rows.
		switch $isDirection {
		    h - x { 
			grid $cbw -column $nCheckbox -row 0
		    }
		    v - y { 
			grid $cbw -column 0 -row $nCheckbox
		    }
		}
	    }
	    
	    default { continue }
	}
	
	incr nCheckbox
    }
    
    grid columnconfigure $ifwTop 0 -weight 0
    grid columnconfigure $ifwTop 1 -weight 1
}

# tkuMakeSliders
# -font : overall font
# -sliders : list of slider data (required)
#    -orientation : orientation of this slider (h, v)
#    -label : label to left of slider
#    -postlabel : label to right of slider
#    -font : font for this slider
#    -variable  : variable to set (required)
#    -min : min value
#    -max : max value
#    -resolution : numerical resolution of slider
#    -length : length of slider
#    -command : command to call when changed
#    -entry : whether to include a text entry
#    -entrywidth : width of entry
#    -limitentry : if 1, allows user to type numbers outside of the slider
#                  range into the entry field. NOTE that if this option is
#                  used, the -variable MUST be global to function properly.
proc tkuMakeSliders { ifwTop args } {

    # set default arguments for all fields
    set aArgs(-font) [tkuLabelFont]

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-sliders} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeSliders: no $arg specified"
	    return
	}
    }

    frame $ifwTop
    
    # for each slider...
    set nSlider 0
    foreach lSlider $aArgs(-sliders) {
	
	# Set arg items and make sure we have the ones we require,
	set aSlider(-orientation) h
	set aSlider(-label) ""
	set aSlider(-postlabel) ""
	set aSlider(-font) $aArgs(-font)
	set aSlider(-min) 0
	set aSlider(-max) 10
	set aSlider(-resolution) 1.0
	set aSlider(-length) 100
	set aSlider(-entry) 0
	set aSlider(-limitentry) 1
	set aSlider(-entrywidth) 0
	set aSlider(-command) ""
	array set aSlider $lSlider
	foreach arg {-variable} {
	    if {![info exists aSlider($arg)]} {
    	     puts "tkuMakeSlider: no $arg specified in slider $nSlider"
		return
	    }
	}

	# Make labels.
	tkuMakeNormalLabel $ifwTop.lw$nSlider -label $aSlider(-label)
	tkuMakeNormalLabel $ifwTop.lwPost$nSlider -label $aSlider(-postlabel)

	# Make the slider.
	set s $ifwTop.sw$nSlider
	scale $s \
	    -orient $aSlider(-orientation) \
	    -from $aSlider(-min) \
	    -to $aSlider(-max) \
	    -length $aSlider(-length) \
	    -resolution $aSlider(-resolution) \
	    -showvalue false

	# Save the command here because we may need to modify it.
	set cmd $aSlider(-command)

	# If we're limiting the entry, just assign the given variable
	# name to the slider.
	if { $aSlider(-limitentry) } {
	    $s config -variable ::$aSlider(-variable)
	} else {

	    # If not, we do something tricky. We create an additional
	    # variable, the variable name with -slider attached, and
	    # use that as the slider variable. We modify the slider
	    # command so that when that variable changes, we update
	    # the real variable. NOTE that we're using global
	    # variables here because that's the only thing we can
	    # trace with the trace command (see below).
	    $s config -variable ::$aSlider(-variable)-slider
	    set oldCmd $cmd
	    set cmd "set ::$aSlider(-variable) \${::$aSlider(-variable)-slider}; $oldCmd"

	    # Make an inital assignment to the variable since our
	    # slider no longer does this for us.
	    set ::$aSlider(-variable) $aSlider(-min)

	    # Now we need to update the slider's fake variable
	    # whenever our real one is updated. So we put a trace on
	    # the real variable so that when it's edited, we call a
	    # callback function.
	    trace variable ::$aSlider(-variable) w \
		"tkuSliderCallback $ifwTop.sw$nSlider ::$aSlider(-variable)-slider"
	}

	bind $ifwTop.sw$nSlider <ButtonRelease> $cmd
	bind $ifwTop.sw$nSlider <B1-Motion> $cmd

	if { $aSlider(-entry) } {
	    entry $ifwTop.ew$nSlider \
		-textvariable ::$aSlider(-variable) \
		-width $aSlider(-entrywidth) \
		-selectbackground green \
		-insertbackground black
	    bind $ifwTop.ew$nSlider <Return> $aSlider(-command)
	}

	# if horizontal, pack all in the same row. if vertical.
	# pack the checkbox, than the label, in the
	# same row as the number of this checkbox.
	grid $ifwTop.lw$nSlider     -column 0 -row $nSlider -sticky w
	grid $ifwTop.sw$nSlider     -column 1 -row $nSlider -sticky ew
	if { $aSlider(-entry) } {
	    grid $ifwTop.ew$nSlider -column 2 -row $nSlider -sticky w
	    grid $ifwTop.lwPost$nSlider -column 3 -row $nSlider -sticky w
	} else {
	    grid $ifwTop.lwPost$nSlider -column 2 -row $nSlider -sticky w
	}		
	
	incr nSlider
    }

    grid columnconfigure $ifwTop 0 -weight 0
    grid columnconfigure $ifwTop 1 -weight 1
    grid columnconfigure $ifwTop 2 -weight 0
    grid columnconfigure $ifwTop 3 -weight 0
}


proc tkuSliderCallback { iswSlider iSliderVar iEntryVar iArrayIndex iOp } {

    # We get passed the variable name and, if it's an array, the array
    # index. I don't know why. We want the value of this variable, so
    # we check to see if the name we got is an array or not, and find
    # the value appropriately.
    set entryValue 0
    if { [array exists ::$iEntryVar] } {
	upvar ::${iEntryVar}($iArrayIndex) entryVar
	set entryValue $entryVar
    } else {
	upvar ::${iEntryVar} entryVar
	set entryValue $entryVar
    }

    # Grab a reference to the slider variable.
    upvar $iSliderVar sliderVar

    # Get our min and max from the slider.
    set min [$iswSlider cget -from]
    set max [$iswSlider cget -to]
    
    # Position the slider according to the new value. If it's outside
    # of the slider range, just set the slider var to the max or min.
    if { $entryValue >= $min &&
	 $entryValue <= $max } {
	set sliderVar $entryValue
    }
    if { $entryValue < $min } {
	set sliderVar $min
    }
    if { $entryValue > $max } {
	set sliderVar $max
    }
}

# Based on new min and max, will change -from and -to of all sliders
# in this group and recalculate a decent resolution.
proc tkuUpdateSlidersRange { ifwTop iMin iMax {iIncrement -1} } {

    # See what our value range is and set a good resolution.
    set newResolution $iIncrement
    if { $newResolution == -1 } {
	set diff [expr $iMax - $iMin]
	if { $diff > 1000 } {
	    set newResolution 10
	} elseif { $diff > 100 } {
	    set newResolution 1
	} elseif { $diff > 10 } {
	    set newResolution .1
	} elseif { $diff > 1 } {
	    set newResolution .001
	} elseif { $diff > 0.1 } {
	    set newResolution .0001
	} elseif { $diff > 0.01 } {
	    set newResolution .00001
	} elseif { $diff > 0.001 } {
	    set newResolution .000001
	} elseif { $diff > 0.0001 } {
	    set newResolution .0000001
	} elseif { $diff > 0.00001 } {
	    set newResolution .00000001
	} elseif { $diff > 0.000001 } {
	    set newResolution .000000001
	} else {
	    set newResolution .0000000001
	}
    }

    set nSlider 0
    set err 0
    while { $err == 0 } {

	set err [catch { 
	    $ifwTop.sw$nSlider config -from $iMin -to $iMax \
		-resolution $newResolution}]

	incr nSlider
    }
}

proc tkuSetSlidersEnabled { ifwTop ibEnabled } {

    set state enabled
    if { $ibEnabled } { 
	set state normal 
    } else {
	set state disabled 
    }

    set nSlider 0
    set err 0
    while { $err == 0 } {

	catch { 
	    tkuSetNormalLabelEnabled $ifwTop.lw$nSlider $ibEnabled
	    tkuSetNormalLabelEnabled $ifwTop.lwPost$nSlider $ibEnabled
	    $ifwTop.ew$nSlider config -state $state
	}

	set err [catch { 
	    $ifwTop.sw$nSlider config -state $state } sResult]
	incr nSlider
    }
}

# tkuMakeToolbar
# -allowzero : whether or not the toolbar can have no buttons on
# -radio : whether or not the toolbar can only have one button on
# -variable : variable to set (required)
# -command : command to call when value changes
# -font : font for all buttons
# -buttons : list of button arrays
#    -type : type of button (text,image)
#    -name : name and value of this button
#    -label : text label for text button
#    -image : image for image button
#    -font : font for this button
#    -balloon : balloon tip for this button
proc tkuMakeToolbar { ifwTop args } {
    
    # set default arguments for all fields
    set aArgs(-allowzero) false
    set aArgs(-radio) true
    set aArgs(-command) ""
    set aArgs(-font) [tkuNormalFont]

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-variable -buttons} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeToolbar: no $arg specified"
	    return
	}
    }

    frame $ifwTop

    tixBalloon $ifwTop.balloon -initwait [tkuBalloonWait]
    
    tixSelect $ifwTop.tbw \
	-allowzero $aArgs(-allowzero) \
	-radio $aArgs(-radio) \
	-command $aArgs(-command) \
	-disablecallback true
    
    tkuEnableLater $ifwTop.tbw
    
    foreach lButton $aArgs(-buttons) {
	
	# set default arguments for all fields
	set aButton(-type) text
	set aButton(-label) ""
	set aButton(-font) $aArgs(-font)
	array set aButton $lButton
	foreach arg {-name} {
	    if {![info exists aButton($arg)]} {
		puts "tkuMakeToolbar: no $arg specified for toolbar button"
		return
	    }
	}

	switch $aButton(-type) {

	    text {
		$ifwTop.tbw  \
		    add $aButton(-name) \
		    -text $aButton(-label) \
		    -font $aButton(-font)
	    } 
	    image {
		$ifwTop.tbw \
		    add $aButton(-name) \
		    -image $aButton(-image)
	    }
	    default {
		puts "didn't handle type $aButton(-type)"
	    }
	} 

	if { [info exists aButton(-balloon)] } {
	    set button [$ifwTop.tbw subwidget $aButton(-name)]
	    $ifwTop.balloon bind $button \
		-balloonmsg $aButton(-balloon)
	}
    }
    
    $ifwTop.tbw config -variable $aArgs(-variable)
    
    pack $ifwTop.tbw
}

proc tkuEnableLater { ifwWidget } {
    global glDisabledWidgets
    lappend glDisabledWidgets $ifwWidget
}

proc tkuFinish { } {
    global glDisabledWidgets
    foreach widget $glDisabledWidgets {
	$widget config -disablecallback false
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

# tkuDoFileDlog
# -title : dlog title
# -okCmd : cmd to execute on OK
# -cancel : cmd to execute on cancel
# -promptN : prompt for item N (required to make this item active)
# -noteN : note for item N
# -typeN : type for item N (file,dir,text,checkbutton,note)
# -defaultfuncN : func to return default dir for item N
# -defaultdirN : default dir for item N
# -shortcutdirs :  shortcut list for all file and dir itmes
# -shortcutdirsN : shortcut list for item N (will override above)
# -defaultvalueN : initial value for item N

# USAGE:
# tkuDoFileDlog -title "Save File" \
#  -prompt1 "Save File As:" \
#  -note1 "The file name to write" \
#  -defaultfunc1 [list GetDefaultLocation SaveFile] \
#  -shortcutdirs {/usr /usr/bin ~/documents} \
#  -okCmd {SaveFile %s1; \
#          SetDefaultLocation SaveFile %s1}
proc tkuDoFileDlog { args } {

    global gDialog
    global sFileName1 sFileName2 sFileName3 sFileName4 sFileName5

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
		continue
	    }
	    
	    set sDefaultDirOption ""
	    set sDefaultDirValue ""
	    if { [info exists aArgs(-defaultfunc$nField)] } {
		set sDefaultDirOption -defaultfunc
		set sDefaultDirValue $aArgs(-defaultfunc$nField)
	    }
	    if { [info exists aArgs(-defaultdir$nField)] } {
		set sDefaultDirOption -defaultdir
		set sDefaultDirValue $aArgs(-defaultdir$nField)
	    }
	    set sShortcutDirsOption ""
	    set sShortcutDirsValue ""
	    if { [info exists aArgs(-shortcutdirs$nField)] } {
		set sShortcutDirsOption -shortcutdirs
		set sShortcutDirsValue $aArgs(-shortcutdirs$nField)
	    }
	    if { [info exists aArgs(-shortcutdirs)] } {
		set sShortcutDirsOption -shortcutdirs
		set sShortcutDirsValue $aArgs(-shortcutdirs)
	    }

	    if { [info exists aArgs(-defaultvalue$nField)] } {
		set sFileName$nField $aArgs(-defaultvalue$nField)
	    }
	    
	    # switch on the type for this field and create the approriate
	    # selecter widget. bind it to sFileName[1..5]
	    switch $aArgs(-type$nField) {
		file { 
		    tkuMakeFileSelector [set fwPrompt$nField] \
			-text "$aArgs(-prompt$nField)" \
			-variable sFileName$nField \
			$sDefaultDirOption $sDefaultDirValue \
			$sShortcutDirsOption $sShortcutDirsValue
		}
		dir { 
		    tkuMakeDirectorySelector [set fwPrompt$nField] \
			-text "$aArgs(-prompt$nField)" \
			-variable sFileName$nField \
			$sDefaultDirOption $sDefaultDirValue \
			$sShortcutDirsOption $sShortcutDirsValue
		}
		text { 
		    tkuMakeEntry [set fwPrompt$nField] \
			-label "$aArgs(-prompt$nField)" \
			-variable sFileName$nField 
		}
		checkbox { 
		    tkuMakeCheckboxes [set fwPrompt$nField] \
			-orientation h \
			-checkboxes [list \
                           [list -type text -label "$aArgs(-prompt$nField)" \
				-variable sFileName$nField] ]
		}
		menu { 
		    # Make sure we have a -menu argument.
		    if { ![info exists aArgs(-menu$nField)] } {
			puts "Missing -menu for field $nField"
			return
		    }

		    # This one is complicated. We have a frame with a
		    # normal label in it, then a menubutton. The
		    # menubutton has a menu which we fill from the
		    # -menu argument we should've gotten. That uses a
		    # callback which will set the text of the menu
		    # button as well as the variable name for this
		    # field.
		    set fwMenu [set fwPrompt$nField]
		    set lw     $fwMenu.lw
		    set fmbw   $fwMenu.fmbw
		    set mbw    $fmbw.mbw
		    set mw     $mbw.mw

		    frame $fwMenu

		    tkuMakeNormalLabel $lw \
			-label "$aArgs(-prompt$nField)"

		    frame $fmbw -relief raised -border 2

		    menubutton $mbw \
			-text "Choose" \
			-menu $mw \
			-indicatoron 1

		    pack $mbw

		    menu $mw
		    foreach lValueLabel $aArgs(-menu$nField) {
			set value [lindex $lValueLabel 0]
			set sLabel [lindex $lValueLabel 1]
			$mw add command \
			    -label $sLabel \
	    -command "FileDlogMenuCallback $nField $mbw \"$sLabel\" $value"
		    }

		    pack $lw $fmbw \
			-anchor w \
			-side left

		    # Unless we got one, default value is -1;
		    if { [info exists aArgs(-defaultitem$nField)] &&
			 $aArgs(-defaultitem$nField) != -1 } {
			$mw invoke $aArgs(-defaultitem$nField)
		    } else {
			set sFileName$nField -1
		    }
		}
		note { 
		    tkuMakeNormalLabel [set fwPrompt$nField] \
			-label "$aArgs(-prompt$nField)"
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

proc FileDlogMenuCallback { inField imbw isLabel iValue } {
    global sFileName$inField
    $imbw config -text "$isLabel"
    set sFileName$inField $iValue
}

		    

# ============================================================ FILE SELECTORS

# tkuMakeFileSelector
# -variable : variable to set with filename (required)
# -text : text prompt
# -command : command to execute when file is chosen
# -defaultfunc : function to call to get default dir
# -defaultdir : default dir
# -shortcutdirs : dirs to go in shortcut menu
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

    frame $ifwTop -width 200
    
    # the entry
    tixLabelEntry $ifwTop.ew \
	-label $aArgs(-text) \
	-labelside acrosstop \
	-options "entry.textVariable $aArgs(-variable) \
                  entry.expand yes \
                  entry.fill x \
                  label.font [tkuLabelFont]"
    
    [$ifwTop.ew subwidget entry] icursor end

    set sDefaultDirOption ""
    if { [info exists aArgs(-defaultfunc)] } {
	set sDefaultDirOption "-defaultfunc $aArgs(-defaultfunc)"
    }
    if { [info exists aArgs(-defaultdir)] } {
	set sDefaultDirOption "-defaultdir $aArgs(-defaultdir)"
    }
    set sShortcutDirsOption ""
    if { [info exists aArgs(-shortcutdirs)] } {
	set sShortcutDirsOption "-shortcutdirs $aArgs(-shortcutdirs)"
    }
    set sCommandOption ""
    if { [info exists aArgs(-command)] } {
	set sCommandOption "-command [list $aArgs(-command)]"
    }

    # the browse button
    button $ifwTop.bw \
	-text "Browse..." \
	-command "tkuBrowseFile -variable $aArgs(-variable) $sDefaultDirOption $sShortcutDirsOption $sCommandOption"
    
    # pack it in a grid
    grid $ifwTop.ew -column 0 -row 0 -sticky ew
    grid $ifwTop.bw -column 1 -row 0
    grid columnconfigure $ifwTop 0 -weight 1
    grid columnconfigure $ifwTop 1 -weight 0
}

# tkuBrowseFile
# -variable : variable name to set (required)
# -command : command to execute when file is chosen
# -defaultfunc : function to call to get default dir
# -defaultdir : default dir
# -shortcutdirs : dirs to go in shortcut menu
proc tkuBrowseFile { args } {
    global gPostHandleSelectFileCommand

    # Set menu items and make sure we have the ones we require,
    set aArgs(-command) ""    
    array set aArgs $args
    foreach arg {-variable} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuBrowseFile: no $arg specified"
	    return
	}
    }

    # create the dialog box if it doesn't already exist
    set wwDirDlog [tix filedialog tixFileSelectDialog]
    
    # set the default location. if it's actually returning a file
    # name, get just the directory portion instead.
    set fnDefaultDir ""
    if { [info exists aArgs(-defaultfunc)] } {
	set fnDefaultDir [$aArgs(-defaultfunc)]
    }
    if { [info exists aArgs(-defaultdir)] } {
	set fnDefaultDir $aArgs(-defaultdir)
    }
    if { $fnDefaultDir != "" } {
	if { [file isfile $fnDefaultDir] } {
	    set fnDefaultDir [file dirname $fnDefaultDir]
	}
	[$wwDirDlog subwidget fsbox] configure -directory $fnDefaultDir
    }

    # add shortcuts
    if { [info exists aArgs(-shortcutdirs)] } {
	foreach sDirectory $aArgs(-shortcutdirs) {
	    [[$wwDirDlog subwidget fsbox] subwidget filter] \
		appendhistory $sDirectory
	}
    }
    
    # when they click ok, call the tkuHandleSelectDirectory function,
    # passing in the variable from the parent dialog.
    $wwDirDlog config -command "tkuHandleSelectFile $aArgs(-variable)"
    
    $wwDirDlog popup

    set gPostHandleSelectFileCommand $aArgs(-command)
}

proc tkuHandleSelectFile { iVariable ifnFile } {
    global gPostHandleSelectFileCommand
    # set the variable.
    upvar $iVariable theVar 
    set theVar $ifnFile
    if { $gPostHandleSelectFileCommand != "" } {
	uplevel #0 {eval $gPostHandleSelectFileCommand}
    }
}

proc tkuUpdateFileSelectorVariable { ifwSelector } {
    $ifwSelector.ew update;
}

# tkuMakeDirectorySelector:
# -variable : variable to set with dirname (required)
# -text : text prompt
# -defaultfunc : function to call to get default dir
# -defaultdir : default dir
# -shortcutdirs : dirs to go in shortcut menu
proc tkuMakeDirectorySelector { ifwTop args } {
    
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

    frame $ifwTop -width 200
    
    # the entry
    tixLabelEntry $ifwTop.ew \
	-label $aArgs(-text) \
	-labelside acrosstop \
	-options "entry.textVariable $aArgs(-variable) \
                  entry.expand yes \
                  entry.fill x \
                  label.font [tkuLabelFont]"
    
    [$ifwTop.ew subwidget entry] icursor end

    set sDefaultDirOption ""
    if { [info exists aArgs(-defaultfunc)] } {
	set sDefaultDirOption "-defaultfunc $aArgs(-defaultfunc)"
    }
    if { [info exists aArgs(-defaultdir)] } {
	set sDefaultDirOption "-defaultdir $aArgs(-defaultdir)"
    }
    set sShortcutDirsOption ""
    if { [info exists aArgs(-shortcutdirs)] } {
	set sShortcutDirsOption "-shortcutdirs $aArgs(-shortcutdirs)"
    }

    # the browse button
    button $ifwTop.bw \
	-text "Browse..." \
	-command "tkuBrowseDirectory -variable $aArgs(-variable) $sDefaultDirOption $sShortcutDirsOption"
    
    # pack it in a grid
    grid $ifwTop.ew -column 0 -row 0 -sticky ew
    grid $ifwTop.bw -column 1 -row 0
    grid columnconfigure $ifwTop 0 -weight 1
    grid columnconfigure $ifwTop 1 -weight 0
}

# tkuBrowseDirectory
# -variable : variable name to set (required)
# -defaultfunc : function to call to get default dir
# -defaultdir : default dir
# -shortcutdirs : dirs to go in shortcut menu
proc tkuBrowseDirectory { args } {
    
    # Set menu items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-variable} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuBrowseDirectory: no $arg specified"
	    return
	}
    }

    # create the dialog box if it doesn't already exist
    set wwDirDlog .wwDirDlog
    if ![winfo exists $wwDirDlog] {
	tixDirSelectDialog $wwDirDlog
    }
    
    # set the default location. if it's actually returning a file
    # name, get just the directory portion instead.
    set fnDefaultDir ""
    if { [info exists aArgs(-defaultfunc)] } {
	set fnDefaultDir [$aArgs(-defaultfunc)]
    }
    if { [info exists aArgs(-defaultdir)] } {
	set fnDefaultDir $aArgs(-defaultdir)
    }
    if { $fnDefaultDir != "" } {
	if { [file isfile $fnDefaultDir] } {
	    set fnDefault [file dirname $fnDefaultDir]
	}
	[$wwDirDlog subwidget dirbox] configure -value $fnDefaultDir
    }
    
    # add shortcuts
    if { [info exists aArgs(-shortcutdirs)] } {
	foreach sDirectory $aArgs(-shortcutdirs) {
	    [[[$wwDirDlog subwidget dirbox] subwidget dircbx] subwidget combo]\
		appendhistory $sDirectory
	}
    }
    
    # when they click ok, call the tkuHandleSelectDirectory function,
    # passing in the variable from the parent dialog.
    $wwDirDlog config -command "tkuHandleSelectDirectory $aArgs(-variable)"
    
    $wwDirDlog popup
}

proc tkuHandleSelectDirectory { iVariable ifnDir } {
    # set the variable.
    upvar $iVariable theVar 
    set theVar $ifnDir
}

proc tkuUpdateDirectorySelectorVariable { ifwSelector } {
    $ifwSelector.ew update;
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

# tkuMakeCloseButton
# -closeCmd : command to execute when Cancel button is pressed
proc tkuMakeCloseButton { ifwTop iwwTop args } {

    global kLabelFont

    # Set menu items and make sure we have the ones we require,
    array set aArgs $args

    set closeCmd ""
    if {[info exists aArgs(-closeCmd)]} { set closeCmd $aArgs(-closeCmd) }

    frame $ifwTop
    
    button $ifwTop.bwClose \
      -text "Close" \
      -command "$closeCmd; tkuCloseDialog $iwwTop"

    bind $iwwTop <Escape> \
      "$ifwTop.bwClose flash; $ifwTop.bwClose invoke"

    pack $ifwTop.bwClose \
      -side right \
      -padx 5 \
      -pady 5
}

# tkuMakeApplyCloseButtons
# -applyCmd : command to execute when Apply button si pressed
# -applyLabel : option label for Apply button
# -closeCmd : command to execute when Close button is pressed
proc tkuMakeApplyCloseButtons { ifwTop iwwTop args } {

    global kLabelFont

    # Set menu items and make sure we have the ones we require,
    array set aArgs $args

    set closeCmd ""
    if {[info exists aArgs(-closeCmd)]} { set closeCmd $aArgs(-closeCmd) }
    set applyCmd ""
    if {[info exists aArgs(-applyCmd)]} { set applyCmd $aArgs(-applyCmd) }
    set sApplyLabel "Apply"
    if {[info exists aArgs(-applyLabel)]} { set sApplyLabel $aArgs(-applyLabel) }

    frame $ifwTop
    
    button $ifwTop.bwApply \
      -text $sApplyLabel \
      -command "$applyCmd"
    button $ifwTop.bwClose \
      -text "Close" \
      -command "$closeCmd; tkuCloseDialog $iwwTop"

    bind $iwwTop <Escape> \
      "$ifwTop.bwClose flash; $ifwTop.bwClose invoke"

    pack $ifwTop.bwClose $ifwTop.bwApply \
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

proc tkuErrorDlog { isMsg } {

    tk_messageBox -type ok \
      -icon error \
      -message $isMsg \
      -title "Error"
}

proc tkuFormattedErrorDlog { isTitle isMsg isDesc } {

    global gDialog

    set wwDialog .wwFormattedErrorDlog

    # try to create the dlog...
    if { [tkuCreateDialog $wwDialog "Error" {-borderwidth 10}] } {
	
	set fwText       $wwDialog.fwText
	set fwButtons    $wwDialog.fwButtons
	
	text $fwText -width 40 \
	    -height 10 \
	    -spacing3 10 \
	    -relief flat \
	    -wrap word
	$fwText insert end "Error: $isTitle \n" {tTitle}
	$fwText insert end "$isMsg \n" {tMsg}
	$fwText insert end "$isDesc \n" {tDesc}
	$fwText tag configure tTitle -font [tkuLabelFont]
	$fwText tag configure tMsg -font [tkuNormalFont]
	$fwText tag configure tDesc -font [tkuNormalFont]
	$fwText configure -state disabled
	
	# button.
	tkuMakeCloseButton $fwButtons $wwDialog
	
	pack $fwText $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc tkuAlertDlog { isMsg } {

    tk_messageBox -type ok \
      -icon info \
      -message $isMsg \
      -title "Note"
}

# tkuMakeColorPickers
# -font : font for all picker labels
# -pickers
#   -font : font for picker label
#   -label : label for picker
#   -redVariable : variable for red component from 0 - 255 (required)
#   -blueVariable : variable for blue component from 0 - 255 (required)
#   -greenVariable : variable for green component from 0 - 255 (required)
#   -command : command to execute when color is selected
set gColorPickerInfo(id) 0
proc tkuMakeColorPickers { ifwTop args } {
    global gColorPickerInfo

    # set default arguments for all fields
    set aArgs(-font) [tkuLabelFont]

    # Set arg items and make sure we have the ones we require,
    array set aArgs $args
    foreach arg {-pickers} {
	if {![info exists aArgs($arg)]} {
	    puts "tkuMakeColorPickers: no $arg specified"
	    return
	}
    }

    frame $ifwTop
    
    # for each picker...
    set nPicker 0
    foreach lPicker $aArgs(-pickers) {
	
	# Set arg items and make sure we have the ones we require,
	set aPicker(-font) $aArgs(-font)
	set aPicker(-label) ""
	set aPicker(-command) ""
	array set aPicker $lPicker
	foreach arg {-redVariable -blueVariable -greenVariable} {
	    if {![info exists aPicker($arg)]} {
    	     puts "tkuMakeColorPickers:: no $arg specified in picker $nPicker"
		return
	    }
	}
	
	upvar #0 $aPicker(-redVariable) red
	upvar #0 $aPicker(-greenVariable) green
	upvar #0 $aPicker(-blueVariable) blue

	set id $gColorPickerInfo(id)
	incr gColorPickerInfo(id)

	lappend gColorPickerInfo($ifwTop,idList) $id

	tkuMakeNormalLabel $ifwTop.fwLabel-$nPicker -label $aPicker(-label) \
	    -font $aPicker(-font)
	
	canvas $ifwTop.fwColor-$nPicker -height 20 -width 20
	
	button $ifwTop.bwChange-$nPicker \
	    -font $aPicker(-font) \
	    -text "Select color..." \
	    -command "tkuColorPicker_CreateWindow $id tkuColorPickerCallback"
	
	$ifwTop.fwColor-$nPicker create rectangle 0 0 20 20 \
	    -fill [format "#%.2x%.2x%.2x" $red $green $blue] \
	    -tag color
	
	set gColorPickerInfo($id,canvas) $ifwTop.fwColor-$nPicker
	set gColorPickerInfo($id,command) $aPicker(-command)
	set gColorPickerInfo($id,redVar) $aPicker(-redVariable)
	set gColorPickerInfo($id,greenVar) $aPicker(-greenVariable)
	set gColorPickerInfo($id,blueVar) $aPicker(-blueVariable)

	grid $ifwTop.fwLabel-$nPicker  -column 0 -row $nPicker -sticky w
	grid $ifwTop.fwColor-$nPicker  -column 1 -row $nPicker -sticky we
	grid $ifwTop.bwChange-$nPicker -column 2 -row $nPicker -sticky e

	incr nPicker
    }

    grid columnconfigure $ifwTop 0 -weight 0
    grid columnconfigure $ifwTop 1 -weight 1
    grid columnconfigure $ifwTop 2 -weight 0

}

proc tkuUpdateColorPickerValues { ifwTop } {
    global gColorPickerInfo

    if { ![info exists gColorPickerInfo($ifwTop,idList)] } { 
	puts "tkuUpdateColorPickerValues: no idList for $ifwTop"
	return
    }
    
    foreach id $gColorPickerInfo($ifwTop,idList) {
	
	upvar #0 $gColorPickerInfo($id,redVar) red
	upvar #0 $gColorPickerInfo($id,greenVar) green
	upvar #0 $gColorPickerInfo($id,blueVar) blue
	
	$gColorPickerInfo($id,canvas) delete color
	
	$gColorPickerInfo($id,canvas) create rectangle 0 0 20 20 \
	    -fill [format "#%.2x%.2x%.2x" $red $green $blue]
    }
}

proc tkuColorPickerCallback { iID iRed iGreen iBlue } {
    global gColorPickerInfo

    upvar #0 $gColorPickerInfo($iID,redVar) red
    upvar #0 $gColorPickerInfo($iID,greenVar) green
    upvar #0 $gColorPickerInfo($iID,blueVar) blue

    set red $iRed
    set green $iGreen
    set blue $iBlue

    set gColorPickerInfo(callback) $gColorPickerInfo($iID,command)
    uplevel #0 eval {$gColorPickerInfo(callback)}

    $gColorPickerInfo($iID,canvas) delete color

    $gColorPickerInfo($iID,canvas) create rectangle 0 0 20 20 \
	-fill [format "#%.2x%.2x%.2x" $iRed $iGreen $iBlue]
}

proc tkuColorPicker_CreatePicker { iwTop iID iCallbackFunction {izSquare 16}} {
    global gColorPickerCB

    # Picker config.
    set kzSquareX $izSquare
    set kzSquareY $izSquare
    set klSquaresX 25
    set klSquaresY 9
    
    # 216 colors
    set kcReds 4
    set kcGreens 4
    set kcBlues 8

    set gColorPickerCB $iCallbackFunction

    frame $iwTop

    set cwPicker       $iwTop.cwPicker

    # Create the canvas with the right width and height
    canvas $cwPicker \
	-width [expr $kzSquareX * $klSquaresX] \
	-height [expr $kzSquareY * $klSquaresY]

    # Find the color increments for each next square.
    set redInc   [expr 256 / $kcReds]
    set greenInc [expr 256 / $kcGreens]
    set blueInc  [expr 256 / $kcBlues]

    # Start off with color 0,0,0.
    set r 0
    set g 0
    set b 0

    # For each square...
    for { set nY 0 } { $nY < $klSquaresY } { incr nY } {
	for { set nX 0 } { $nX < $klSquaresX } { incr nX } {
	    
	    # Create a square in the current color, converting our
	    # numerical values to a hex color value. Also give it the
	    # common tag 'color' and a tag made out of its numerical
	    # color components.
	    $cwPicker create rectangle \
		[expr $nX * $kzSquareX] [expr $nY * $kzSquareY] \
		[expr ($nX+1) * $kzSquareX] [expr ($nY+1) * $kzSquareY] \
		-fill [format "#%.2x%.2x%.2x" $r $g $b] \
		-tags "color $r-$g-$b"

	    # Increment our numerical color values in the order of r,
	    # g, and then b. With the increments we have, we'll
	    # actually get to 256, but we really want a top value of
	    # 255, so check for the 256 and bump it down to
	    # 255. (Dividing with 255 when finding the increments
	    # doesn't give us the full 0-255 range.
	    incr r $redInc
	    if { $r == 256 } { set r 255 }
	    if { $r > 255 } {
		set r 0
		incr g $greenInc
		if { $g == 256 } { set g 255 }
		if { $g > 255 } {
		    set g 0
		    incr b $blueInc
		    if { $b == 256 } { set b 255 }
		    if { $b > 255 } {
			set b 0
		    }
		}
	    }
	}
    }

    # When we click on a color square, call the function.
    $cwPicker bind color <Button-1> "tkuColorPicker_HandleClick $iID %W"

    pack $cwPicker -fill both
}

proc tkuColorPicker_HandleClick { iID icwPicker } {
    global gColorPickerCB

    # Get the numerical tag from the current element, the one the
    # mouse is in. This is our r-g-b tag. Extract the numerical elements.
    set color [lindex [$icwPicker gettags current] 1]
    scan $color "%d-%d-%d" r g b

    # Detroy the color picker dlog.
    destroy .wwColorPicker-$iID

    # Call the callback function.
    $gColorPickerCB $iID $r $g $b
}

proc tkuColorPicker_CreateWindow { iID iCallbuckFunction } {

    # If it doesn't exist, make a new window, put a color picker
    # inside it with our given callback function, and pack it.
    if { ![winfo exists .wwColorPicker] } {

	toplevel .wwColorPicker-$iID
	wm title .wwColorPicker-$iID "Choose a color..."
	tkuMakeNormalLabel .wwColorPicker-$iID.lwInstructions \
	    -label "Chose a color..."
	tkuColorPicker_CreatePicker .wwColorPicker-$iID.cpwColor \
	    $iID $iCallbuckFunction
	pack .wwColorPicker-$iID.lwInstructions .wwColorPicker-$iID.cpwColor \
	    -side top

    } else {

	raise .wwColorPicker-$iID
    }
}




