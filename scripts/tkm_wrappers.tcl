#! /usr/bin/tixwish

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

# tkm_MakeButtons fwFrame {button...}
# button = { text "Label" {function} "Tool tip" }
# button = { image image_name {function} "Tool tip" }

# tkm_MakeMenu isMenuButton "Menu Name" {item...}
# item = { command   "Item Name" command                [group_name] }
# item = { radio     "Item Name" command variable value [group_name] }
# item = { check     "Item Name" command variable       [group_name] }
# item = { cascade   "Item Name" {item...}              [group_name] }
# item = { separator }

# tkm_SetMenuItemGroupStatus group_name 0|1

# tkm_MakeEntryWithIncDecButtons fwFrame "Label" variable {function} step
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
# tkm_MakeCancelOKButtons    fwFrame wwDlog {ok_fuction} [{cancel_function}]

# tkm_MakeFileSelector fwFrame "Prompt:" variable default_lcoation_func
# tkm_UpdateFileSelectorVariable fwFrame

# tkm_MakeDirectorySelector fwFrame "Prompt:" variable default_lcoation_func
# tkm_UpdateDirectorySelectorVariable fwFrame

set kNormalFont -*-lucida-medium-r-normal-*-13-*-*-*-*-*-*-*
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

set knBalloonWait 500

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

    frame $isFrame

    if { $isText != "" } {
  tkm_MakeNormalLabel $isFrame.lw $isText

  pack $isFrame.lw \
    -side left \
    -anchor w
    }

    entry $isFrame.ew \
      -textvariable $isVariable \
      -width $inWidth \
      -font $kNormalFont \
      -state disabled \
      -relief flat

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
      grid $rbw -column [expr 2 * $nRadioButton] -row 0
      grid $lw -column [expr 1 + [expr 2 * $nRadioButton]] \
        -row 0 -sticky w
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
      -selectcolor gray

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
      $fwMain.baw bind $cbw \
        -balloonmsg [lindex $lRadioButton 4]

      catch { $fwMain.baw bind $lw \
        -balloonmsg [lindex $lRadioButton 4] }
  }

  incr nRadioButton
    }

    grid columnconfigure $fwMain 0 -weight 0
    grid columnconfigure $fwMain 1 -weight 1
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
  bind $isFrame.ewEntry <Return> "$iSetFunction [set $iVariable]"
    }

}

proc tkm_MakeButtons { isFrame ilButtons } {

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
      
      pack $isFrame.bw$nButton \
        -side left \
        -expand yes

  } 
  if { [string compare $sType "image"] == 0 } {

      button $isFrame.bw$nButton \
        -image [lindex $lButton 1] \
        -command [lindex $lButton 2]
      
      pack $isFrame.bw$nButton \
        -side left \
        -expand yes

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

    # for every letter in menu name..
    # if this letter is not in the underline list...
    # add the letter to the list.
    # underline this index.
    # if no letters left, underline nothing.

    menubutton $isMenuButton \
      -text $isMenuName \
      -font $kNormalFont \
      -menu $isMenuButton.mw
#      -underline $nUnderline

    # start an underline list for this menu

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

  if { [string compare $sType "command"] == 0 } {

      # for every letter in item name..
      # if this letter is not in the local underline list...
      # add the letter to the list.
      # underline this index.
      # if no letters left, underline nothing.

      $isMenu add command \
        -label [lindex $lItem 1] \
        -command [lindex $lItem 2] \
        -font $kNormalFont
#              -underline $nUnderline

      set sGroupName [lindex $lItem 3]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenu $nItemNum
      }
     
      set bProcessed 1
  }
  if { [string compare $sType "radio" ] == 0 } {
      $isMenu add radio \
        -label [lindex $lItem 1] \
        -command [lindex $lItem 2] \
        -variable [lindex $lItem 3] \
        -value [lindex $lItem 4] \
        -font $kNormalFont      

      set sGroupName [lindex $lItem 5]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenu $nItemNum
      }
      
      set bProcessed 1
  }
  if { [string compare $sType "check" ] == 0 } {
      $isMenu add check \
        -label [lindex $lItem 1] \
        -command [lindex $lItem 2] \
        -variable [lindex $lItem 3] \
        -font $kNormalFont      

      set sGroupName [lindex $lItem 4]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenu $nItemNum
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
        -menu $isMenu.cmw$nItemNum \
        -font $kNormalFont      

      set lCascadeItems [lindex $lItem 2]
      tkm_AddMenuItemsToMenu $isMenu.cmw$nItemNum $lCascadeItems

      set sGroupName [lindex $lItem 3]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenu $nItemNum
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

proc tkm_MakeEntryWithIncDecButtons { isFrame isText iVariable iSetFunc ifStep } {

    global kLabelFont kNormalFont
    global $iVariable
    global glDisabledWidgets

    frame $isFrame

    tixControl $isFrame.control \
      -command $iSetFunc \
      -label $isText \
      -variable $iVariable \
      -step $ifStep \
      -disablecallback true

    tkm_EnableLater $isFrame.control

    $isFrame.control subwidget label configure -font $kNormalFont

    pack $isFrame.control \
      -anchor w

}

proc nothing {} {

    label $isFrame.lwText \
      -text $isText \
      -font $kNormalFont

    entry $isFrame.ewEntry \
      -textvariable $iVariable \
      -width $inWidth
    bind $isFrame.ewEntry <Return> $iSetFunc

    button $isFrame.bwDec \
      -text "-" \
      -font $kNormalFont \
      -padx 2 \
      -pady 0 \
      -bd 0
    bind $isFrame.bwDec <ButtonRelease-1> $iDecFunc

    button $isFrame.bwInc \
      -text "+" \
      -font $kNormalFont \
      -padx 2 \
      -pady 0 \
      -bd 0
    bind $isFrame.bwInc <ButtonRelease-1> $iIncFunc

    pack $isFrame.lwText \
      -side left \
      -padx 2 \
      -anchor w

    pack $isFrame.bwInc $isFrame.bwDec $isFrame.ewEntry \
      -side right \
      -padx 2 \
      -anchor e
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

proc tkm_MakeColorPicker { isFrame isLabel iRedVar iGreenVar iBlueVar iCmd inWidth } {

    tixLabelFrame $isFrame \
      -label $isLabel \
      -labelside acrosstop \
      -options { label.padX 5 }
    
    set fwColorSub  [$isFrame subwidget frame]
    set fwColors    $fwColorSub.fwColors
    
    tkm_MakeSliders $fwColors [list\
      [list {"Red"} $iRedVar 0 1 $inWidth $iCmd 1 0.1 ] \
      [list {"Green"} $iGreenVar 0 1 $inWidth $iCmd 1 0.1 ] \
      [list {"Blue"} $iBlueVar 0 1 $inWidth $iCmd 1 0.1 ]]
    
    pack $fwColors \
      -side top \
      -anchor w \
      -expand yes \
      -fill x
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
# -okCmd string:        command to execute when okay button is hit. 
#                       %s[1..5] is replaced with the file name entered.
# -cancelCmd string:    command to execute when cancel button is hit
#
proc tkm_DoFileDlog { ilArgs } {

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
  set tArgs(-default$nField) ""
  set sFileName$nField ""
    }

    # allow passed options to override defaults
    array set tArgs $ilArgs

    # dialog name
    set wwDialog .wwDlogFile

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
    file { set sFileName$nField ""; \
      tkm_MakeFileSelector [set fwPrompt$nField] \
      "$tArgs(-prompt$nField)" sFileName$nField \
      $tArgs(-default$nField) }
    dir { set sFileName$nField ""; \
      tkm_MakeDirectorySelector [set fwPrompt$nField] \
      "$tArgs(-prompt$nField)" sFileName$nField \
      $tArgs(-default$nField) }
    text { set sFileName$nField $tArgs(-default$nField); \
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
  } $wwDialog $knWidth $wwDialog] 
    }
}



proc tkm_MakeCancelApplyOKButtons { isFrame isTop iOKCmd {iCancelCmd ""} } {

    global kLabelFont

    frame $isFrame
    
    button $isFrame.bwApply \
      -text "Apply" \
      -command "$iOKCmd" \
      -font $kLabelFont

    button $isFrame.bwOK \
      -text "OK" \
      -command "$iOKCmd; Dialog_Close $isTop" \
      -font $kLabelFont

    button $isFrame.bwCancel \
      -text "Cancel" \
      -command "$iCancelCmd; Dialog_Close $isTop" \
      -font $kLabelFont

    bind $isTop <Return> \
      "$isTop.fwButtons.bwOK flash; $isTop.fwButtons.bwOK invoke"
    bind $isTop <space> \
      "$isTop.fwButtons.bwApply flash; $isTop.fwButtons.bwApply invoke"
    bind $isTop <Escape> \
      "$isTop.fwButtons.bwCancel flash; $isTop.fwButtons.bwCancel invoke"
    pack $isFrame.bwOK $isFrame.bwApply $isFrame.bwCancel \
      -side right \
      -padx 5 \
      -pady 5

}

proc tkm_MakeCloseButton { isFrame isTop {iCloseCmd ""} } {

    global kLabelFont

    frame $isFrame

    button $isFrame.bwClose \
      -text "Close" \
      -command "$iCloseCmd; Dialog_Close $isTop" \
      -font $kLabelFont

    bind $isTop <Escape> \
      "$isTop.fwButtons.bwClose flash; $isTop.fwButtons.bwClose invoke"
    pack $isFrame.bwClose \
      -side right \
      -padx 5 \
      -pady 5

}

proc tkm_MakeApplyCloseButtons { isFrame isTop iApplyCmd {iCloseCmd ""} } {

    global kLabelFont

    frame $isFrame
    
    button $isFrame.bwApply \
      -text "Apply" \
      -command "$iApplyCmd" \
      -font $kLabelFont

    button $isFrame.bwClose \
      -text "Close" \
      -command "$iCloseCmd; Dialog_Close $isTop" \
      -font $kLabelFont

    bind $isTop <space> \
      "$isTop.fwButtons.bwApply flash; $isTop.fwButtons.bwApply invoke"
    bind $isTop <Escape> \
      "$isTop.fwButtons.bwClose flash; $isTop.fwButtons.bwClose invoke"
    pack $isFrame.bwApply $isFrame.bwClose \
      -side right \
      -padx 5 \
      -pady 5

}

proc tkm_MakeCancelOKButtons { isFrame isTop iOKCmd {iCancelCmd ""} } {

    global kLabelFont

    frame $isFrame
    
    button $isFrame.bwOK \
      -text "OK" \
      -command "$iOKCmd; Dialog_Close $isTop" \
      -font $kLabelFont

    button $isFrame.bwCancel \
      -text "Cancel" \
      -command "$iCancelCmd; Dialog_Close $isTop" \
      -font $kLabelFont

    bind $isTop <Return> \
      "$isTop.fwButtons.bwOK flash; $isTop.fwButtons.bwOK invoke"
    bind $isTop <Escape> \
      "$isTop.fwButtons.bwCancel flash; $isTop.fwButtons.bwCancel invoke"
    pack $isFrame.bwOK $isFrame.bwCancel \
      -side right \
      -padx 5 \
      -pady 5

}

proc tkm_MakeFileSelector { isFrame isText iVariable {iDefaultFunc ""}} {

    frame $isFrame -width 200
    
    upvar $iVariable theVar 

    tixFileEntry $isFrame.few \
      -label $isText \
      -labelside top \
      -variable $iVariable \
      -options {
         entry.expand yes
         entry.fill x
             }
    
    # set the value of the field to the value of the variable
    $isFrame.few config -value $theVar

    # set the default location 
    if { $iDefaultFunc != "" } {
  $isFrame.few filedialog subwidget fsbox configure -directory \
    [eval $iDefaultFunc]
    }

    pack $isFrame.few \
      -side left \
      -expand yes \
      -fill x
}

proc tkm_UpdateFileSelectorVariable { isFrame } {

    $isFrame.few update
}

proc tkm_MakeDirectorySelector { isFrame isText iVariable {iDefaultFunc ""}} {
    
    frame $isFrame -width 200
    
    upvar $iVariable theVar 

    tixFileEntry $isFrame.few \
      -label $isText \
      -labelside top \
      -variable $iVariable \
      -dialogtype tixDirSelectDialog \
      -options {
          entry.expand yes
          entry.fill x
            }
    
    # set the value of the field to the value of the variable
    $isFrame.few config -value $theVar

    # set the default location 
    if { $iDefaultFunc != "" } {
  [$isFrame.few filedialog subwidget dirbox] \
    subwidget dirlist configure -value [eval $iDefaultFunc]
    }
    
    pack $isFrame.few \
      -side left \
      -expand yes \
      -fill x
}

proc tkm_UpdateDirectorySelectorVariable { isFrame } {

    $isFrame.few update;
}

proc tkm_AddItemToMenuGroup { isGroupName ifwMenuObject inMenuItemNum } {
    
    global glMenuGroups

    # add this menu / item pair to the list
    lappend glMenuGroups($isGroupName) "$ifwMenuObject $inMenuItemNum"
}

proc tkm_SetMenuItemGroupStatus { isGroupName ibEnable } {

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
  puts "exists!"
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
