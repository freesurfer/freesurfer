##
## tkm_wrappers.tcl
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

# tkm_MakeBigLabel fwFrame "Label Text"
# tkm_MakeSmallLabel fwFrame "Label Text"
# tkm_MakeNormalLabel fwFrame "Label Text"
# tkm_MakeActiveLabel fwFrame "Left Text" variable [width]

# tkm_MakeCheckboxes fwFrame direction=x|y {checkbox...}
# checkbox = { text "Label"     variable {function} ["Tool tip"] }
# checkbox = { image image_name variable {function} ["Tool tip"] }

# tkm_MakeRadioButtons fwFrame direction=x|y "Frame Title" variable {radio button...}
# radio button = { text  label      value cmd tooltip }
# radio button = { image image_name value cmd tooltip }
# tkm_MakeRadioButton fwFrame "Label Text" variable value [{function}]

# tkm_MakeToolbar fwFrame 0|1=radio variable {function} {button...}
# button = { text   name "Label"    ["Tool tip"] }
# button = { bitmap name image_name ["Tool tip"] }

# tkm_MakeEntry fwFrame "Left Text" variable [width] [{function}]

# tkm_MakeButtons fwFrame {button...} [direction=x|y]
# button = { text "Label" {function} "Tool tip" }
# button = { image image_name {function} "Tool tip" }

# tkm_MakeMenu isMenuButton "Menu Name" {item...}
# item = { command   "Item Name" command                [group_name] }
# item = { radio     "Item Name" command variable value [group_name] }
# item = { check     "Item Name" command variable       [group_name] }
# item = { cascade   "Item Name" {item...}              [group_name] }
# item = { separator }

# tkm_SetEnableGroupStatus group_name 0|1

# tkm_MakeEntryWithIncDecButtons fwFrame "Label" variable {function} step range
# tkm_MakeSlider fwFrame {"Left text" "Right text"} variable
#                min max length {function} 1|0=include_entry [resolution] 
#                [horizontal|vertical]

# tkm_MakeSliders fwFrame {slider...}
# slider = { {"Left text" "Right text"} variable min max length 
#            {function} 1|0=include_entry [resolution] [horizontal|vertical] }


# tkm_MakeSliders fwFrame {slider...}
# slider = { {"Left text" "Right text"} variable min max length {function]
#            include_entry=1|0 [resolution] }

# tkm_DoFileDlog { [option arg option arg...] }
# -title string:        the title of the dlog box
# -prompt[1..5] string: the prompt for the nth file
# -type[1..5] file|dir: the type of selector (file by default)
# -note[1..5] string:   a note to place under the prompt
# -okCmd string:        command to execute when okay button is hit. 
#                       %s[1..5] is replaced with the file name entered.
# -cancelCmd string:    command to execute when cancel button is hit


# tkm_MakeCancelApplyOKButtons fwFrame wwDlog {ok_fuction} [{cancel_function}]
# tkm_MakeCloseButton        fwFrame wwDlog [{close_function}]
# tkm_MakeApplyCloseButtons  fwFrame wwDlog {apply_function} [{close_function}]
#                            [Apply_button_text]
# tkm_MakeCancelOKButtons    fwFrame wwDlog {ok_fuction} [{cancel_function}]
# tkm_MakeApplyCloseButtons  fwFrame wwDlog {apply_function} [{close_function}]
#                            [Apply_button_text]
# tkm_MakeDialogButtons      fwFrame wwDlog {button}
# button = {type cmd [label]}
# types: apply, OK, help, close

# tkm_MakeFileSelector fwFrame "Prompt:" variable default_lcoation_func
# tkm_UpdateFileSelectorVariable fwFrame

# tkm_MakeDirectorySelector fwFrame "Prompt:" variable default_lcoation_func
# tkm_UpdateDirectorySelectorVariable fwFrame

set kNormalFont -*-lucida-medium-r-normal-*-12-*-*-*-*-*-*-*
set kSmallFont -*-lucida-medium-i-normal-*-10-*-*-*-*-*-*-*
set kLabelFont -*-lucida-bold-r-normal-*-14-*-*-*-*-*-*-*
set kHighlightBgColor white

catch {
    if { [string compare $env(TKM_FONT_SIZE) small] == 0 } {
	set kNormalFont -*-lucida-medium-r-normal-*-11-*-*-*-*-*-*-*
	set kSmallFont -*-lucida-medium-i-normal-*-10-*-*-*-*-*-*-*
	set kLabelFont -*-lucida-bold-r-normal-*-11-*-*-*-*-*-*-*
    }
}

catch {
    if { [string compare $env(TKM_FONT_SIZE) large] == 0 } {
	set kNormalFont -*-lucida-medium-r-normal-*-14-*-*-*-*-*-*-*
	set kSmallFont -*-lucida-medium-i-normal-*-10-*-*-*-*-*-*-*
	set kLabelFont -*-lucida-bold-r-normal-*-14-*-*-*-*-*-*-*
    }
}

catch {
    if { [string compare $env(TKM_FONT_FAMILY) helvetica] == 0 } {
	
	set kNormalFont -*-helvetica-medium-r-normal-*-12-*-*-*-*-*-*-*
	set kSmallFont -*-helvetica-medium-i-normal-*-10-*-*-*-*-*-*-*
	set kLabelFont -*-helvetica-bold-r-normal-*-14-*-*-*-*-*-*-*
	set kHighlightBgColor white
	
	catch {
	    if { [string compare $env(TKM_FONT_SIZE) small] == 0 } {
		set kNormalFont -*-helvetica-medium-r-normal-*-11-*-*-*-*-*-*-*
		set kSmallFont -*-helvetica-medium-i-normal-*-10-*-*-*-*-*-*-*
		set kLabelFont -*-helvetica-medium-r-normal-*-11-*-*-*-*-*-*-*
	    }
	}
	
	catch {
	    if { [string compare $env(TKM_FONT_SIZE) large] == 0 } {
		set kNormalFont -*-helvetica-medium-r-normal-*-14-*-*-*-*-*-*-*
		set kSmallFont -*-helvetica-medium-i-normal-*-10-*-*-*-*-*-*-*
		set kLabelFont -*-helvetica-medium-r-normal-*-14-*-*-*-*-*-*-*
	    }
	}
    }
}

set knBalloonWait 500

proc tkm_GetNormalFont {} { global kNormalFont; return $kNormalFont }
proc tkm_GetSmallFont {} { global kSmallFont; return $kSmallFont }
proc tkm_GetLabelFont {} { global kLabelFont; return $kLabelFont }

proc tkm_MakeBigLabel { isFrame isText {inWrapLength 0} } {
    
    global kLabelFont
    
    frame $isFrame
    
    tixLabelWidget $isFrame.label \
	-label $isText
    
    $isFrame.label subwidget label configure -font $kLabelFont -justify left
    if { $inWrapLength != 0 } {
	$isFrame.label subwidget label configure -wraplength $inWrapLength
    }
    
    pack $isFrame.label \
	-side left \
	-anchor w
}

proc tkm_MakeSmallLabel { isFrame isText {inWrapLength 0} } {
    
    global kSmallFont
    
    frame $isFrame
    
    tixLabelWidget $isFrame.label \
	-label $isText
    
    $isFrame.label subwidget label configure -font $kSmallFont -justify left
    
    if { $inWrapLength != 0 } {
	$isFrame.label subwidget label configure -wraplength $inWrapLength
    }
    
    pack $isFrame.label \
	-side left \
	-anchor w
}

proc tkm_MakeNormalLabel { isFrame isText {inWrapLength 0} } {
    
    global kNormalFont
    
    frame $isFrame
    
    tixLabelWidget $isFrame.label \
	-label $isText
    
    $isFrame.label subwidget label configure -font $kNormalFont -justify left
    if { $inWrapLength != 0 } {
	$isFrame.label subwidget label configure -wraplength $inWrapLength
    }
    
    pack $isFrame.label \
	-side left \
	-anchor w
    
}

proc tkm_MakeActiveLabel { isFrame isText isVariable {inWidth -1} } {
    
    global kNormalFont kHighlightBgColor
    global tk_version

    frame $isFrame
    
    if { $isText != "" } {
	tkm_MakeNormalLabel $isFrame.lw "$isText"
	
	pack $isFrame.lw \
	    -side left \
	    -anchor w
    }

    set entryState disabled
    if { $tk_version >= 8.4 } {
	set entryState readonly
    }

    entry $isFrame.ew \
	-textvariable $isVariable \
	-width $inWidth \
	-font $kNormalFont \
	-state $entryState \
	-relief flat \
	-foreground black
    
    pack $isFrame.ew \
	-side left \
	-anchor w \
	-expand yes \
	-fill x
}

proc tkm_MakeCheckboxes { isFrame isDirection ilCheckboxes } {
    
    global kLabelFont kNormalFont kHighlightBgColor
    global knBalloonWait
    
    frame $isFrame
    
    # { text label variable command tooltip }
    # { image image_name variable command tooltip }
    
    tixBalloon $isFrame.baw -initwait $knBalloonWait
    
    # for each checkbox...
    set nCheckbox 0
    foreach lCheckbox $ilCheckboxes {
	
	# grab the type.
	set sType [lindex $lCheckbox 0]
	
	# make names for the checkbox and the label.
	set cbw $isFrame.cb$nCheckbox
	set lw  $isFrame.lw$nCheckbox
	
	# text or image?
	switch $sType {
	    
	    text {
		
		# text. make a normal checkbox and label.
		checkbutton $cbw \
		    -variable [lindex $lCheckbox 2] \
		    -command [lindex $lCheckbox 3]
		label $lw \
		    -font $kNormalFont \
		    -text [lindex $lCheckbox 1]
		
		# if horizontal, pack all in the same row. if vertical.
		# pack the checkbox, than the label, in the
		# same row as the number of this checkbox.
		switch $isDirection {
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
		    -image [lindex $lCheckbox 1] \
		    -variable [lindex $lCheckbox 2] \
		    -command [lindex $lCheckbox 3] \
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
	
	# do we have a tool tip? if so, create one for the checkbox. try
	# to create one for the label, but if we did an image we don't have
	# one.
	if { [lindex $lCheckbox 4] != "" } {
	    $isFrame.baw bind $cbw \
		-balloonmsg [lindex $lCheckbox 4]
	    
	    catch { $isFrame.baw bind $lw \
			-balloonmsg [lindex $lCheckbox 4] }
	}
	
	incr nCheckbox
    }
    
    grid columnconfigure $isFrame 0 -weight 0
    grid columnconfigure $isFrame 1 -weight 1
}

proc tkm_MakeCheckbox { isFrame isType isTextOrImage iVariable iSetFunction {isTooltip ""} } {
    
    global kLabelFont kNormalFont kHighlightBgColor
    global knBalloonWait
    
    frame $isFrame
    
    if { $isType == "text" } {
	checkbutton $isFrame.cb \
	    -text $isTextOrImage \
	    -variable $iVariable \
	    -font $kNormalFont 
    } else {
	checkbutton $isFrame.cb \
	    -image $isTextOrImage \
	    -variable $iVariable \
	    -font $kNormalFont \
	    -indicatoron false \
	    -selectcolor gray
    }
    
    if { $isTooltip != "" } {
	tixBalloon $isFrame.baw -initwait $knBalloonWait
	$isFrame.baw bind $isFrame.cb -balloonmsg $isTooltip
    }
    
    bind $isFrame.cb <ButtonRelease-1> $iSetFunction
    
    pack $isFrame.cb \
	-anchor w
}

proc tkm_MakeRadioButtons { isFrame isDirection isTitle iVariable ilRadioButtons } {
    
    global kLabelFont kNormalFont kHighlightBgColor
    global knBalloonWait
    
    if { $isTitle != "" } {
	
	tixLabelFrame $isFrame \
	    -label $isTitle \
	    -labelside acrosstop \
	    -options { label.padX 5 }
	
	set fwMain [$isFrame subwidget frame]
    } else {
	set fwMain [frame $isFrame]
    }
    
    # { text  label      value cmd tooltip }
    # { image image_name value cmd tooltip }
    
    tixBalloon $fwMain.baw -initwait $knBalloonWait
    
    # for each radio button...
    set nRadioButton 0
    foreach lRadioButton $ilRadioButtons {
	
	# grab the type.
	set sType [lindex $lRadioButton 0]
	
	# make names for the radio and the label.
	set rbw $fwMain.rb$nRadioButton
	set lw  $fwMain.lw$nRadioButton
	
	# text or image?
	switch $sType {
	    
	    text {
		
		# text. make a normal radio and label.
		radiobutton $rbw              \
		    -command [lindex $lRadioButton 3]  \
		    -font $kNormalFont    \
		    -variable $iVariable  \
		    -relief flat          \
		    -value [lindex $lRadioButton 2]
		label $lw \
		    -font $kNormalFont \
		    -text [lindex $lRadioButton 1]
		
		# if horizontal, pack all in the same row. if vertical.
		# pack the checkbox, than the label, in the
		# same row as the number of this checkbox.
		switch $isDirection {
		    h - x { 
			grid $rbw -column [expr 2 * $nRadioButton] -row 0 \
			    -sticky e
			grid $lw -column [expr 1 + (2 * $nRadioButton)] \
			    -row 0 -sticky w
			grid columnconfigure $fwMain \
			    [expr 2 * $nRadioButton] -weight 0
			grid columnconfigure $fwMain \
			    [expr 1+ (2 * $nRadioButton)] -weight 1
		    }
		    v - y { 
			grid $rbw -column 0 -row $nRadioButton
			grid $lw -column 1 -row $nRadioButton -sticky w
		    }
		}
	    }
	    
	    image { 
		# image. create a checkbox with an image label. no text label.
		radiobutton $rbw              \
		    -image [lindex $lRadioButton 1] \
		    -variable $iVariable \
		    -command [lindex $lRadioButton 3] \
		    -indicatoron false \
		    -selectcolor gray \
		    -value [lindex $lRadioButton 2]
		
		# if horizontal, pack in increasing columns. if vertical,
		# pack in increasing rows.
		switch $isDirection {
		    h - x { 
			grid $rbw -column $nRadioButton -row 0
		    }
		    v - y { 
			grid $rbw -column 0 -row $nRadioButton
		    }
		}
	    }
	    
	    default { continue }
	}
	
	# do we have a tool tip? if so, create one for the checkbox. try
	# to create one for the label, but if we did an image we don't have
	# one.
	if { [lindex $lRadioButton 4] != "" } {
	    $fwMain.baw bind $rbw \
		-balloonmsg [lindex $lRadioButton 4]
	    
	    catch { $fwMain.baw bind $lw \
			-balloonmsg [lindex $lRadioButton 4] }
	}
	
	incr nRadioButton
    }
    
    switch $isDirection {
	v - y { 
	    grid columnconfigure $fwMain 0 -weight 0
	    grid columnconfigure $fwMain 1 -weight 1
	}
    }
}

proc tkm_MakeRadioButton { isFrame isText iVariable iValue {iCmd ""} } {
    
    global kLabelFont kNormalFont kHighlightBgColor 
    
    frame $isFrame
    
    radiobutton $isFrame.rbw        \
	-text $isText         \
	-command $iCmd        \
	-font $kNormalFont    \
	-variable $iVariable  \
	-relief flat          \
	-value $iValue
    
    pack $isFrame.rbw    \
	-side left \
	-anchor w
}

proc tkm_MakeEntry { isFrame isText iVariable {inWidth -1} {iSetFunction ""} } {
    
    global kLabelFont kNormalFont
    
    frame $isFrame
    
    if { $isText != "" } {
	
	label $isFrame.lwLabel \
	    -text $isText \
	    -font $kNormalFont
	
	pack $isFrame.lwLabel \
	    -side left \
	    -anchor w
    }
    
    entry $isFrame.ewEntry \
	-font $kNormalFont \
	-textvariable $iVariable \
	-width $inWidth \
	-selectbackground green \
	-insertbackground black
    
    pack $isFrame.ewEntry \
	-side right \
	-anchor e \
	-expand yes \
	-fill x
    
    if { $iSetFunction != "" } {
	bind $isFrame.ewEntry <Return> $iSetFunction
    }
    
}

proc tkm_MakeButtons { isFrame ilButtons {isDirection x} } {
    
    global kLabelFont kNormalFont
    global knBalloonWait
    
    frame $isFrame
    
    # { text label command tooltip }
    # { image image command tooltip }
    
    tixBalloon $isFrame.baw -initwait $knBalloonWait
    
    set nButton 0
    foreach lButton $ilButtons {
	
	set sType [lindex $lButton 0]
	
	if { [string compare $sType "text"] == 0 } {
	    
	    button $isFrame.bw$nButton \
		-font $kLabelFont \
		-text [lindex $lButton 1] \
		-command [lindex $lButton 2] 
	    
	} 
	if { [string compare $sType "image"] == 0 } {
	    
	    button $isFrame.bw$nButton \
		-image [lindex $lButton 1] \
		-command [lindex $lButton 2]
	    
	}
	
	switch $isDirection {
	    h - x {
		pack $isFrame.bw$nButton \
		    -side left \
		    -expand yes
	    }
	    v - y {
		pack $isFrame.bw$nButton \
		    -side top \
		    -expand yes \
		    -fill x \
		    -pady 2
	    }
	}
	
	if { [lindex $lButton 3] != "" } {
	    $isFrame.baw bind $isFrame.bw$nButton \
		-balloonmsg [lindex $lButton 3]
	}
	
	incr nButton
    }
}

proc tkm_MakeMenu { isMenuButton isMenuName ilMenuItems } {
    
    global kLabelFont kNormalFont glUnderlineList
    
    menubutton $isMenuButton \
	-text $isMenuName \
	-font $kNormalFont \
	-menu $isMenuButton.mw
    
    tkm_AddMenuItemsToMenu $isMenuButton.mw $ilMenuItems
}

# { command   "Item Name" command_to_execute                group_name } 
# { radio     "Item Name" command_to_execute variable value group_name }
# { check     "Item Name" command_to_execute variable       group_name }
# { cascade   "Item Name" {items list}                      group_name }
# { separator }

proc tkm_AddItemsToCascade { isCascadeMenuItem ilMenuItems } {
    
    tkm_AddMenuItemsToMenu $isCascadeMenuItem $ilMenuItems
}

proc tkm_AddMenuItemsToMenu { isMenu ilMenuItems } {
    
    global kLabelFont kNormalFont glUnderlineList
    
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
			    -font $kNormalFont \
			    -accelerator $sAccel

			set sGroupName [lindex $lItem 3]
		    }
		    "radio" {
			$isMenu add radio \
			    -label $sName \
			    -command [lindex $lItem 2] \
			    -variable [lindex $lItem 3] \
			    -value [lindex $lItem 4] \
			    -font $kNormalFont \
			    -accelerator $sAccel
			
			set sGroupName [lindex $lItem 5]
		    }
		    "check" {
			$isMenu add check \
			    -label $sName \
			    -command [lindex $lItem 2] \
			    -variable [lindex $lItem 3] \
			    -font $kNormalFont \
			    -accelerator $sAccel
			
			set sGroupName [lindex $lItem 4]
		    }
		}

		if { [string compare $sGroupName ""] != 0 } {
		    tkm_AddMenuItemToEnableGroup $sGroupName $isMenu $nItemNum
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
		    -font $kNormalFont      
		
		set lCascadeItems [lindex $lItem 2]
		tkm_AddMenuItemsToMenu $isMenu.cmw$nItemNum $lCascadeItems
		
		set sGroupName [lindex $lItem 3]
		if { [string compare $sGroupName ""] != 0 } {
		    tkm_AddMenuItemToEnableGroup $sGroupName $isMenu $nItemNum
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

proc tkm_MakeEntryWithIncDecButtons { isFrame isText iVariable iSetFunc ifStep {iRange {}} } {
    
    global kLabelFont kNormalFont kSmallFont
    global $iVariable
    global glDisabledWidgets
    
    frame $isFrame
    
    # Function to use if there is no set function provided. The new
    # value will be appended to this function when it changes, so in
    # effect this turns out to be: set IGNORE_ME $value
    if { [string length $iSetFunc] == 0 } {
	set iSetFunc "set IGNORE_ME"
    }
    
    # Make our entry command to set the value and button commands to
    # inc and dec the value. If they gave us a range, put the command
    # inside an if statement that checks the range.
    if { [llength $iRange] == 2 } {
	set min [lindex $iRange 0]
	set max [lindex $iRange 1]

	tkm_MakeEntry $isFrame.ew \
	    $isText $iVariable 4 "if { \[set $iVariable\] > $min && \[set $iVariable\] < $max } { $iSetFunc \[set $iVariable\] } else { tk_messageBox -icon error -message \"Value must be between $min and $max.\" }"
    
        # DNG added wrap-around
	button $isFrame.bwDec -text "-" \
	    -command "if { \[set $iVariable\] > $min } { incr $iVariable -$ifStep; $iSetFunc \[set $iVariable\] } else { set $iVariable $max;}" \
	    -padx 1 -pady 0
	
        # DNG added wrap-around
	button $isFrame.bwInc -text "+" \
	    -command "if { \[set $iVariable\] < $max } { incr $iVariable $ifStep; $iSetFunc \[set $iVariable\] } else { set $iVariable $min;}" \
	    -padx 0 -pady 0

    } else {

	tkm_MakeEntry $isFrame.ew \
	    $isText $iVariable 4 "$iSetFunc \[set $iVariable\]"
    
	button $isFrame.bwDec -text "-" \
	    -command "incr $iVariable -$ifStep; $iSetFunc \[set $iVariable\]" \
	    -padx 1 -pady 0
	
	button $isFrame.bwInc -text "+" \
	    -command "incr $iVariable $ifStep; $iSetFunc \[set $iVariable\]" \
	    -padx 0 -pady 0
    }

    pack $isFrame.ew $isFrame.bwDec $isFrame.bwInc \
	-anchor w -side left
    
}

# tkm_MakeSlider fwFrame {"prefix" "suffix"} var 0 100 50 {} 1 

proc tkm_MakeSlider { isFrame ilsText iVariable inMin inMax inLength iSetFunc ibIncludeEntry {ifResolution 1.0} } {
    
    frame $isFrame
    
    if { [lindex $ilsText 0] != "" } {
	
	tkm_MakeNormalLabel $isFrame.fwLabel [lindex $ilsText 0]
	pack $isFrame.fwLabel \
	    -side left    \
	    -anchor w     \
	    -expand no
    }
    
    scale $isFrame.sw -orient horizontal \
	-variable $iVariable         \
	-from $inMin                 \
	-to $inMax                   \
	-length $inLength            \
	-resolution $ifResolution    \
	-showvalue false
    bind $isFrame.sw <ButtonRelease> $iSetFunc
    bind $isFrame.sw <B1-Motion>     $iSetFunc
    
    pack $isFrame.sw    \
	-side left  \
	-anchor w   \
	-expand yes \
	-fill x
    
    if { [lindex $ilsText 1] != "" } {
	
	tkm_MakeNormalLabel $isFrame.fwLabel2 [lindex $ilsText 1]
	pack $isFrame.fwLabel2 \
	    -side left     \
	    -anchor w      \
	    -expand no
    }
    
    if { $ibIncludeEntry == 1 } {
	
	entry $isFrame.ew                \
	    -textvariable $iVariable \
	    -width 6                 \
	    -selectbackground green  \
	    -insertbackground black
	bind $isFrame.ew <Return> $iSetFunc
	
	pack $isFrame.ew   \
	    -side left \
	    -anchor w  \
	    -expand no
    }
}

proc tkm_MakeSliders { isFrame ilSliders } {
    
    frame $isFrame
    
    
    set nSlider 0
    foreach lSlider $ilSliders {
	
	set lText       [lindex $lSlider 0]
	set sLeftText   [lindex $lText 0]
	set sRightText  [lindex $lText 1]
	set sVariable   [lindex $lSlider 1]
	set fMin        [lindex $lSlider 2]
	set fMax        [lindex $lSlider 3]
	set nLength     [lindex $lSlider 4]
	set sFunction   [lindex $lSlider 5]
	set bEntry      [lindex $lSlider 6]
	set fResolution [lindex $lSlider 7]
	set sOrientation [lindex $lSlider 8]
	
	if { $fResolution == "" } {
	    set fResolution 1.0
	}
	if { $sOrientation == "" } {
	    set sOrientation horizontal
	}
	
	tkm_MakeNormalLabel $isFrame.fwLeft$nSlider  $sLeftText
	tkm_MakeNormalLabel $isFrame.fwRight$nSlider $sRightText
	
	scale $isFrame.sw$nSlider -orient $sOrientation \
	    -variable $sVariable         \
	    -from $fMin                  \
	    -to $fMax                    \
	    -length $nLength             \
	    -resolution $fResolution     \
	    -showvalue false
	bind $isFrame.sw$nSlider <ButtonRelease> $sFunction
	bind $isFrame.sw$nSlider <B1-Motion>     $sFunction
	
	grid $isFrame.fwLeft$nSlider   -column 0 -row $nSlider -sticky w
	grid $isFrame.sw$nSlider       -column 1 -row $nSlider -sticky we
	grid $isFrame.fwRight$nSlider  -column 2 -row $nSlider -sticky w
	
	if { $bEntry == 1 } {
	    
	    entry $isFrame.ew$nSlider        \
		-textvariable $sVariable \
		-width 6                 \
		-selectbackground green  \
		-insertbackground black
	    bind $isFrame.ew$nSlider <Return> $sFunction
	    
	    grid $isFrame.ew$nSlider -column 3 -row $nSlider
	}
	
	incr nSlider
    }
    
    grid columnconfigure $isFrame 0 -weight 0
    grid columnconfigure $isFrame 1 -weight 1
    grid columnconfigure $isFrame 2 -weight 0
    grid columnconfigure $isFrame 3 -weight 0
}


set gColorPickerInfo(id) 0
proc tkm_MakeColorPickers { isFrame ilPickers } {
    global kLabelFont
    global gColorPickerInfo

    frame $isFrame

    set nPicker 0
    foreach lPicker $ilPickers {

	set isLabel [lindex $lPicker 0]
	set iRedVar [lindex $lPicker 1]
	set iGreenVar [lindex $lPicker 2]
	set iBlueVar [lindex $lPicker 3]
	set iCmd [lindex $lPicker 4]

	upvar #0 $iRedVar red
	upvar #0 $iGreenVar green
	upvar #0 $iBlueVar blue
	
	set id $gColorPickerInfo(id)
	incr gColorPickerInfo(id)
	
	tkm_MakeNormalLabel $isFrame.fwLabel-$nPicker $isLabel
	
	canvas $isFrame.fwColor-$nPicker -height 20 -width 20
	
	button $isFrame.bwChange-$nPicker \
	    -font $kLabelFont \
	    -text "Change color..." \
	    -command "ColorPicker_CreateWindow $id tkm_ColorPickerCallback"
	
	$isFrame.fwColor-$nPicker create rectangle 0 0 20 20 \
	    -fill [format "#%.2x%.2x%.2x" $red $green $blue] \
	    -tag color
	
	set gColorPickerInfo($id,canvas) $isFrame.fwColor-$nPicker
	set gColorPickerInfo($id,command) $iCmd
	set gColorPickerInfo($id,redVar) $iRedVar
	set gColorPickerInfo($id,greenVar) $iGreenVar
	set gColorPickerInfo($id,blueVar) $iBlueVar

	grid $isFrame.fwLabel-$nPicker  -column 0 -row $nPicker -sticky w
	grid $isFrame.fwColor-$nPicker  -column 1 -row $nPicker -sticky we
	grid $isFrame.bwChange-$nPicker -column 2 -row $nPicker -sticky e

	incr nPicker
    }

    grid columnconfigure $isFrame 0 -weight 0
    grid columnconfigure $isFrame 1 -weight 1
    grid columnconfigure $isFrame 2 -weight 0
}

proc tkm_ColorPickerCallback { iID iRed iGreen iBlue } {
    global gColorPickerInfo

    upvar #0 $gColorPickerInfo($iID,redVar) red
    upvar #0 $gColorPickerInfo($iID,greenVar) green
    upvar #0 $gColorPickerInfo($iID,blueVar) blue

    set red $iRed
    set green $iGreen
    set blue $iBlue

    if { $gColorPickerInfo($iID,command) != "" } {
	$gColorPickerInfo($iID,command)
    }

    $gColorPickerInfo($iID,canvas) delete color

    $gColorPickerInfo($iID,canvas) create rectangle 0 0 20 20 \
	-fill [format "#%.2x%.2x%.2x" $iRed $iGreen $iBlue]
    
}

# replaces the percent symbols in a string with a substitution string:
# i.e. tkm_DoSubPercent { %s1 "Hello %s1" "World" }
# returns "Hello World"
proc tkm_DoSubPercent { isPercent isString isSubstitution } {
    if {![string match %* $isPercent]} {
	return $isString;
    }
    regsub -all {\\|&} $isSubstitution {\\\0} isSubstitution
    regsub -all $isPercent $isString $isSubstitution sResult
    return $sResult
}

# 
# puts up a standard file specification dlog with up to 5 entries
# options:
# -title string:        the title of the dlog box
# -prompt[1..5] string: the prompt for the nth file
# -type[1..5] file|dir: the type of selector (file by default)
# -note[1..5] string:   a note to place under the prompt
# -entry[1..5] string:  command that should eval to default field contents
# -default[1..5] string:command that should evaluate to a default dir
# -presets[1..5] list:  a list of preset directories to be available in the pdm
# -okCmd string:        command to execute when okay button is hit. 
# -cancelCmd string:    command to execute when cancel button is hit
#
proc tkm_DoFileDlog { ilArgs } {
    global sFileName1 sFileName2 sFileName3 sFileName4 sFileName5
    global gDialog
    
    set lFields {1 2 3 4 5}
    set knWidth 400
    
    # set default arguments for all fields
    set tArgs(-title) "No title"
    set tArgs(-okCmd) "puts \$sFileName"
    set tArgs(-cancelCmd) ""
    foreach nField $lFields {
	set tArgs(-prompt$nField) ""
	set tArgs(-type$nField) "file"
	set tArgs(-note$nField) ""
	set tArgs(-entry$nField) ""
	set tArgs(-default$nField) ""
	set tArgs(-presets$nField) ""
	set sFileName$nField ""
    }
    
    # allow passed options to override defaults
    array set tArgs $ilArgs
    
    # dialog name. make it the name of the title, subbing dashes for spaces.
    regsub -all { } .wwDlogFile$tArgs(-title) {-} wwDialog
    
    # do percent substitutions in ok command for %s1 thru %s5
    foreach nField $lFields {
	set tArgs(-okCmd) \
	    [tkm_DoSubPercent %s$nField $tArgs(-okCmd) \$sFileName$nField]
    }
    
    # if we can bring up the dialog
    if { [Dialog_Create $wwDialog "$tArgs(-title)" {-borderwidth 10}] } {
	
	# for each field...
	foreach nField $lFields {
	    
	    # create a variable for this prompt. even if we don't use this
	    # field, we'll need it later (ugh)
	    set fwPrompt$nField  $wwDialog.fwPrompt$nField
	    
	    # if they didn't enter a prompt, skip this field
	    if { [string match "$tArgs(-prompt$nField)" ""] == 1 } {
		continue;
	    }
	    
	    # switch on the type for this field and create the approriate
	    # selecter widget. bind it to sFileName[1..5]
	    switch $tArgs(-type$nField) {
		file { set sFileName$nField [eval $tArgs(-entry$nField)]; \
			   tkm_MakeFileSelector [set fwPrompt$nField] \
			   "$tArgs(-prompt$nField)" sFileName$nField \
			   $tArgs(-default$nField) $tArgs(-presets$nField)}
		dir { set sFileName$nField [eval $tArgs(-entry$nField)]; \
			  tkm_MakeDirectorySelector [set fwPrompt$nField] \
			  "$tArgs(-prompt$nField)" sFileName$nField \
			  $tArgs(-default$nField)  $tArgs(-presets$nField)}
		text { set sFileName$nField $tArgs(-entry$nField); \
			   tkm_MakeEntry [set fwPrompt$nField] \
			   "$tArgs(-prompt$nField)" sFileName$nField }
		default { continue; }
	    }
	    
	    # if they requested a note, make one, otherwise just an empty
	    # frame.
	    set fwNote $wwDialog.fwNote$nField
	    if { [string match "$tArgs(-note$nField)" ""] == 0 } {
		tkm_MakeSmallLabel $fwNote "$tArgs(-note$nField)" $knWidth
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
	set fwButtons $wwDialog.fwButtons
        tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    "catch { tkm_UpdateFileSelectorVariable $fwPrompt1 };\
    catch { tkm_UpdateFileSelectorVariable $fwPrompt2 };\
    catch { tkm_UpdateFileSelectorVariable $fwPrompt3 };\
    catch { tkm_UpdateFileSelectorVariable $fwPrompt4 };\
    catch { tkm_UpdateFileSelectorVariable $fwPrompt5 };\
    $tArgs(-okCmd)" \
	    "$tArgs(-cancelCmd)"
	
	# pack the buttons
	pack  $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
	
	# after the next idle, the window will be mapped. set the min
	# width to our width and the min height to the mapped height.
	after idle [format {
	    update idletasks
	    wm minsize %s %d [winfo reqheight %s]
	    wm geometry %s =%dx[winfo reqheight %s]
	} $wwDialog $knWidth $wwDialog $wwDialog $knWidth $wwDialog]
    }
}



proc tkm_MakeCancelApplyOKButtons { isFrame isTop iOKCmd {iCancelCmd ""} } {
    tkm_MakeDialogButtons $isFrame $isTop [list \
					       [list Apply "$iOKCmd"] \
					       [list Close "$iCancelCmd" "Cancel"] \
					       [list OK "$iOKCmd"] \
					      ]
}

proc tkm_MakeCloseButton { isFrame isTop {iCloseCmd ""} } {
    tkm_MakeDialogButtons $isFrame $isTop [list \
					       [list Close "$iCloseCmd"] \
					      ]
}

proc tkm_MakeApplyCloseButtons { isFrame isTop iApplyCmd {iCloseCmd ""} {isApplyText "Apply"}} {
    tkm_MakeDialogButtons $isFrame $isTop [list \
					       [list Apply "$iApplyCmd"] \
					       [list Close "$iCloseCmd"] \
					      ]
}

proc tkm_MakeCancelOKButtons { isFrame isTop iOKCmd {iCancelCmd ""} } {
    tkm_MakeDialogButtons $isFrame $isTop [list \
					       [list OK "$iOKCmd"] \
					       [list Close "$iCancelCmd" "Cancel"] \
					      ]
}

# tkm_MakeDialogButtons      fwFrame wwDlog {button}
# button = {type cmd [label]}
# types: apply, OK, help, close
proc tkm_MakeDialogButtons { isFrame isTop ilButtons } {
    
    set lButtonsToMake {}
    
    # Use these to keep track of whether or not we should bind certain
    # buttons. If we get the button in the list, we'll set its bind
    # index here to the button index we got.
    set nButton 0
    set nBindClose -1
    set nBindApply -1
    set nBindOK -1
    
    # Go through the list of buttons we got. First look for help
    # buttons, then close, then apply, then OK, so we can build a
    # button list to pass to tkm_MakeButtons later one. I use a switch
    # statement here because the syntax is easier.
    foreach lButton $ilButtons {
	set sType [lindex $lButton 0]
	switch $sType {
	    ok - OK {
		set sCommand "[lindex $lButton 1]; Dialog_Close $isTop"
		set sLabel "OK"
		if {[llength $lButton] == 3} {
		    set sLabel [lindex $lButton 2]
		}
		lappend lButtonsToMake [list text $sLabel $sCommand]
		set nBindOK $nButton
		incr nButton
	    }
	}
    }
    foreach lButton $ilButtons {
	set sType [lindex $lButton 0]
	switch $sType {
	    apply - Apply {
		set sCommand [lindex $lButton 1]
		set sLabel "Apply"
		if {[llength $lButton] == 3} {
		    set sLabel [lindex $lButton 2]
		}
		lappend lButtonsToMake [list text $sLabel $sCommand]
		set nBindApply $nButton
		incr nButton
	    }
	}
    }
    foreach lButton $ilButtons {
	set sType [lindex $lButton 0]
	switch $sType {
	    close - Close {
		set sCommand "[lindex $lButton 1]; Dialog_Close $isTop"
		set sLabel "Close"
		if {[llength $lButton] == 3} {
		    set sLabel [lindex $lButton 2]
		}
		lappend lButtonsToMake [list text $sLabel $sCommand]
		set nBindClose $nButton
		incr nButton
	    }
	}
    }
    foreach lButton $ilButtons {
	set sType [lindex $lButton 0]
	switch $sType {
	    help - Help {
		set sCommand [lindex $lButton 1]
		set sLabel "Help"
		if {[llength $lButton] == 3} {
		    set sLabel [lindex $lButton 2]
		}
		lappend lButtonsToMake [list text $sLabel $sCommand]
		incr nButton
	    }
	}
    }
    
    # Now pass the list we made to tkm_MakeButtons.
    tkm_MakeButtons $isFrame $lButtonsToMake
    
    # Now do the bindings. bind Return key to OK, space bar to Apply,
    # and Escape to Close. Note that we refer to the buttons by name,
    # and this is dependent on how tkm_MakeButtons names the buttons.
    if { $nBindOK > -1 } {
	bind $isTop <Return> \
	    "$isFrame.bw$nBindOK flash; $isFrame.bw$nBindOK invoke"
    }
    if { $nBindApply > -1 } {
	bind $isTop <space> \
	    "$isFrame.bw$nBindApply flash; $isFrame.bw$nBindApply invoke"
    }
    if { $nBindClose > -1 } {
	bind $isTop <Escape> \
	    "$isFrame.bw$nBindClose flash; $isFrame.bw$nBindClose invoke"
    }
}

# ================================================== FILE AND DIR SELECTORS

proc tkm_MakeFileSelector { isFrame isText iVariable {iDefaultFunc ""} {ilDirectories ""}} {
    
    frame $isFrame -width 200
    
    # the entry
    tixLabelEntry $isFrame.ew \
	-label $isText \
	-labelside acrosstop \
	-options "entry.textVariable $iVariable \
      entry.expand yes \
      entry.fill x"
    
    [$isFrame.ew subwidget entry] icursor end

    # the browse button
    button $isFrame.bw \
	-text "Browse..." \
	-command "tkm_BrowseFile $iVariable {$iDefaultFunc} {$ilDirectories}"
    
    # pack it in a grid
    grid $isFrame.ew -column 0 -row 0 -sticky ew
    grid $isFrame.bw -column 1 -row 0
    grid columnconfigure $isFrame 0 -weight 1
    grid columnconfigure $isFrame 1 -weight 0
}

proc tkm_BrowseFile { iVariable {iDefaultFunc ""} {ilDirectories ""} } {
    
    # create the dialog box if it doesn't already exist
    set wwDirDlog [tix filedialog tixFileSelectDialog]
    
    # set the default location. if it's actually returning a file
    # name, get just the directory portion instead.
    if { $iDefaultFunc != "" } {
	set fnDefault [eval $iDefaultFunc]
	if { [file isfile $fnDefault] } {
	    set fnDefault [file dirname $fnDefault]
	}
	[$wwDirDlog subwidget fsbox] configure -directory $fnDefault
    }
    
    foreach sDirectory $ilDirectories {
	[[$wwDirDlog subwidget fsbox] subwidget filter] \
	    appendhistory $sDirectory
    }
    
    # when they click ok, call the tkm_HandleSelectDirectory function,
    # passing in the variable from the parent dialog.
    $wwDirDlog config -command "tkm_HandleSelectFile $iVariable"
    
    $wwDirDlog popup
}

proc tkm_HandleSelectFile { iVariable iFile } {
    # set the variable.
    upvar $iVariable theVar 
    set theVar $iFile
}

proc tkm_UpdateFileSelectorVariable { isFrame } {
    $isFrame.ew update;
}

proc tkm_MakeDirectorySelector { isFrame isText iVariable {iDefaultFunc ""} {ilDirectories ""}} {
    
    frame $isFrame -width 200
    
    # the entry
    tixLabelEntry $isFrame.ew \
	-label $isText \
	-labelside acrosstop \
	-options "entry.textVariable $iVariable \
      entry.expand yes \
      entry.fill x"
    
    [$isFrame.ew subwidget entry] icursor end

    # the browse button
    button $isFrame.bw \
	-text "Browse..." \
	-command "tkm_BrowseDirectory $iVariable {$iDefaultFunc} {$ilDirectories}"
    
    # pack it in a grid
    grid $isFrame.ew -column 0 -row 0 -sticky ew
    grid $isFrame.bw -column 1 -row 0
    grid columnconfigure $isFrame 0 -weight 1
    grid columnconfigure $isFrame 1 -weight 0
}

proc tkm_BrowseDirectory { iVariable {iDefaultFunc ""} {ilDirectories ""} } {
    
    # create the dialog box if it doesn't already exist
    set wwDirDlog .wwDirDlog
    if ![winfo exists $wwDirDlog] {
	tixDirSelectDialog $wwDirDlog
    }
    
    # set the default location. if it's actually returning a file
    # name, get just the directory portion instead.
    if { $iDefaultFunc != "" } {
	set fnDefault [eval $iDefaultFunc]
	if { [file isfile $fnDefault] } {
	    set fnDefault [file dirname $fnDefault]
	}
	[$wwDirDlog subwidget dirbox] configure -value $fnDefault
    }
    
    foreach sDirectory $ilDirectories {
	[[[$wwDirDlog subwidget dirbox] subwidget dircbx] subwidget combo]\
	    appendhistory $sDirectory
    }
    
    # when they click ok, call the tkm_HandleSelectDirectory function,
    # passing in the variable from the parent dialog.
    $wwDirDlog config -command "tkm_HandleSelectDirectory $iVariable"
    
    $wwDirDlog popup
}

proc tkm_HandleSelectDirectory { iVariable iDir } {
    # set the variable.
    upvar $iVariable theVar 
    set theVar $iDir
}

proc tkm_UpdateDirectorySelectorVariable { isFrame } {
    $isFrame.ew update;
}

# =========================================================================

proc tkm_IsGroupEnabled { isGroupName } {
    global gbEnabledGroups

    if { ![info exists gbEnabledGroups($isGroupName)] } { 
	return false
    } else {
	return $gbEnabledGroups($isGroupName)
    }
}

proc tkm_AddMenuItemToEnableGroup { isGroupName ifwMenuObject inMenuItemNum } {
    global gbEnabledGroups
    global glEnableGroups

    if { ![info exists gbEnabledGroups($isGroupName)] } { 
	set gbEnabledGroups($isGroupName) false
    }
    
    # add this menu / item pair to the list
    lappend glEnableGroups($isGroupName) "menu $ifwMenuObject $inMenuItemNum"
}

proc tkm_AddCheckboxToEnableGroup { isGroupName ifwCheckbox } {
    global gbEnabledGroups
    global glEnableGroups
    
    if { ![info exists gbEnabledGroups($isGroupName)] } { 
	set gbEnabledGroups($isGroupName) false
    }
    
    # Add this checkbox
    lappend glEnableGroups($isGroupName) "checkbox $ifwCheckbox"
}

proc tkm_AddRadioButtonToEnableGroup { isGroupName ifwRadioButton } {
    global gbEnabledGroups
    global glEnableGroups
    
    if { ![info exists gbEnabledGroups($isGroupName)] } { 
	set gbEnabledGroups($isGroupName) false
    }
    
    # Add this radio button
    lappend glEnableGroups($isGroupName) "radiobutton $ifwRadioButton"
    $ifwRadioButton configure -state disabled
}

proc tkm_SetEnableGroupStatus { isGroupName ibEnable } {
    global gbEnabledGroups
    global glEnableGroups
    
    set gbEnabledGroups($isGroupName) $ibEnable
    
    # for each entry in the list
    foreach lEntry $glEnableGroups($isGroupName) {
	
	# first item is the type
	set sType [lindex $lEntry 0]
	switch $sType {

	    menu {
		# Second item is the menu name, third is the menu item.
		set mbwMenu [lindex $lEntry 1]
		set nMenuItem [lindex $lEntry 2]

		set sState normal
		if { $ibEnable == 0 } {
		    set sState disabled
		}

		if { [catch {
		    $mbwMenu entryconfigure $nMenuItem -state $sState
		} sResult] } {
		    puts "Error en/disabling $isGroupName: "
		    puts "\t$mbwMenu $nMenuItem to $sState"
		    puts "\t$sResult"
		}
	    }

	    checkbox - radiobutton {
		# Second item is widget name.
		set wWidget [lindex $lEntry 1]

		set sState normal
		if { $ibEnable == 0 } {
		    set sState disabled
		}

		if { [catch {
		    $wWidget configure -state $sState
		} sResult] } {
		    puts "Error en/disabling $isGroupName: "
		    puts "\t$wWidget to $sState"
		    puts "\t$sResult"
		}
	    }
	}
    }
}

proc tkm_MakeToolbar { isFrame ibRadio isVariable iCommand ilButtons } {
    
    global kNormalFont
    global knBalloonWait
    
    frame $isFrame
    
    if { $ibRadio == 1 } {
	set bAllowZero false
    } else {
	set bAllowZero true
    }
    
    tixSelect $isFrame.tbw \
	-allowzero $bAllowZero \
	-radio $ibRadio \
	-command $iCommand \
	-disablecallback true
    
    tkm_EnableLater $isFrame.tbw
    
    # { text   name text   tooltip }
    # { bitmap name bitmap tooltip }
    # note name is also the name of the subwidget so can't start with uppercase
    
    tixBalloon $isFrame.baw -initwait $knBalloonWait
    
    foreach lButton $ilButtons {
	
	set sType [lindex $lButton 0]
	
	if { [string compare $sType "text"] == 0 } {
	    set sName [lindex $lButton 1]
	    set sText [lindex $lButton 2]
	    $isFrame.tbw add $sName -text $sText \
		-font $kNormalFont
	} 
	if { [string compare $sType "image"] == 0 } {
	    set sName [lindex $lButton 1]
	    set sBitmap [lindex $lButton 2]
	    $isFrame.tbw add $sName -image $sBitmap
	    if { [lindex $lButton 3] != "" } {
		$isFrame.baw bind [$isFrame.tbw subwidget $sName] -balloonmsg [lindex $lButton 3]
	    }
	} 
    }
    
    $isFrame.tbw config -variable $isVariable
    
    pack $isFrame.tbw
}


proc tkm_EnableLater { ifwIn } {
    
    global glDisabledWidgets
    
    lappend glDisabledWidgets $ifwIn
}

proc tkm_Finish { } {
    
    global glDisabledWidgets
    
    # enabled the disabled widgets
    foreach select $glDisabledWidgets {
	$select config -disablecallback false
    }
}

set gfwProgressTime ""

proc tkm_MakeProgressDlog { isName isMessage } {
    
    #    global gfwProgressTime
    
    #    set wwDialog .wwProgressDlog
    
    # try to create the dlog...
    #    if { [Dialog_Create $wwDialog $isName {-borderwidth 10}] } {
    
    #  set fwLabel           $wwDialog.fwLabel
    #  set gfwProgressTime   $wwDialog.fwTime
    
    #  tkm_MakeBigLabel $fwLabel $isMessage
    #  tkm_MakeNormalLabel $gfwProgressTime "0 / 100"
    
    #  pack $fwLabel $gfwProgressTime -side top
    #    }
}

proc tkm_UpdateProgressDlog { isMessage inCurrent } {
    
    #    global gfwProgressTime
    
    #    [$gfwProgressTime.label subwidget label] configure -text "$inCurrent / 100"
}

proc tkm_DestroyProgressDlog { } {
    
    #    Dialog_Close .wwProgressDlog
}


proc Dialog_Create { iwwTop isTitle iArgs } {
    
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

proc Dialog_Wait { iwwTop iVarName {iFocus {}} } {
    
    upvar $iVarName var
    
    bind $iwwTop <Destroy> {list set $iVarName cancel}
    
    if {[string length $iFocus] == 0} {
	set focus $iwwTop
    }
    set $saveFocus [focus -displayof $iwwTop]
    focus $iFocus
    catch {
	tkwait visibilty $iwwTop
    }
    catch {
	grab $iwwTop
    }
    
    tkwait variable $iVarName
    catch {
	grab release $iwwTop
    }
    focus $saveFocus
}

proc Dialog_Close { iwwTop } {
    
    global gDialog
    
    catch {
	set gDialog(geo,$iwwTop) [wm geometry $iwwTop]
	#  wm withdraw $iwwTop
	destroy $iwwTop
    }
}

proc ColorPicker_CreatePicker { iwTop iID iCallbackFunction {izSquare 16} } {
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
    $cwPicker bind color <Button-1> "ColorPicker_HandleClick $iID %W"

    pack $cwPicker -fill both
}

proc ColorPicker_HandleClick { iID icwPicker } {
    global gColorPickerCB

    # Get the numerical tag from the current element, the one the
    # mouse is in. This is our r-g-b tag. Extract the numerical elements.
    set color [lindex [$icwPicker gettags current] 1]
    scan $color "%d-%d-%d" r g b

    # Detroy the color picker dlog.
    destroy .wwColorPicker

    # Call the callback function.
    $gColorPickerCB $iID $r $g $b
}

proc ColorPicker_CreateWindow { iID iCallbuckFunction } {

    # Make a new window, put a color picker inside it with our given
    # callback function, and pack it.
    toplevel .wwColorPicker 
    wm title .wwColorPicker "Choose a color..."
    tkm_MakeNormalLabel .wwColorPicker.lwInstructions "Chose a color..."
    ColorPicker_CreatePicker .wwColorPicker.cpwColor $iID $iCallbuckFunction
    pack .wwColorPicker.lwInstructions .wwColorPicker.cpwColor \
	-side top
}



