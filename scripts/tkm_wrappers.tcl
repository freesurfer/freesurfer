source $env(MRI_DIR)/lib/tcl/wrappers.tcl

set kNormalFont -*-lucida-medium-r-normal-*-13-*-*-*-*-*-*-*
set kSmallFont -*-lucida-medium-i-normal-*-10-*-*-*-*-*-*-*
set kLabelFont -*-lucida-bold-r-normal-*-14-*-*-*-*-*-*-*
set kHighlightBgColor white

proc tkm_MakeBigLabel { isTop isText } {

    global kLabelFont

    frame $isTop

    label $isTop.label \
      -text $isText \
      -font $kLabelFont

    pack $isTop.label \
      -side left \
      -anchor w
}

proc tkm_MakeSmallLabel { isTop isText } {

    global kSmallFont

    frame $isTop

    label $isTop.label \
      -text $isText \
      -font $kSmallFont

    pack $isTop.label \
      -side left \
      -anchor w
}

proc tkm_MakeNormalLabel { isTop isText } {

    global kNormalFont

    frame $isTop

    label $isTop.label \
      -text $isText \
      -font $kNormalFont

    pack $isTop.label \
      -side left \
      -anchor w

}

proc tkm_MakeActiveLabel { isTop isText isVariable {inWidth -1} } {

    global kNormalFont kHighlightBgColor

   frame $isTop

    label $isTop.lw \
      -text $isText \
      -font $kNormalFont

    entry $isTop.ew \
      -textvariable $isVariable \
      -width $inWidth \
      -font $kNormalFont \
      -state disabled \
      -highlightbackground $kHighlightBgColor \
      -relief flat

    pack $isTop.lw $isTop.ew \
      -side left \
      -anchor w
}

proc tkm_MakeCheckbox { isTop isText iVariable iSetFunction } {

    global kLabelFont kNormalFont kHighlightBgColor

    frame $isTop

    checkbutton $isTop.cb \
      -text $isText \
      -variable $iVariable \
      -font $kNormalFont \
      -highlightbackground $kHighlightBgColor

    bind $isTop.cb <ButtonRelease-1> $iSetFunction

    pack $isTop.cb \
      -anchor w
}

proc tkm_MakeRadioButton { isTop isText iVariable iValue {iCmd ""} } {

    global kLabelFont kNormalFont kHighlightBgColor 

    frame $isTop
    
    radiobutton $isTop.rbw        \
      -text $isText         \
      -command $iCmd        \
      -font $kNormalFont    \
      -variable $iVariable  \
      -relief flat          \
      -value $iValue        \
      -highlightbackground $kHighlightBgColor

    pack $isTop.rbw    \
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

    frame $isFrame

    set nButton 0
    foreach lButton $ilButtons {
  
  button $isFrame.bw$nButton \
    -font $kNormalFont \
    -text [lindex $lButton 0] \
    -command [lindex $lButton 1]
    
  pack $isFrame.bw$nButton \
    -side left

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

    menu $isMenuButton.mw
    set nItemNum 1

    foreach lItem $ilMenuItems {

  set sType [lindex $lItem 0]

  if { [string compare $sType "command"] == 0 } {

      # for every letter in item name..
      # if this letter is not in the local underline list...
      # add the letter to the list.
      # underline this index.
      # if no letters left, underline nothing.

      $isMenuButton.mw add command \
        -label [lindex $lItem 1] \
        -command [lindex $lItem 2] \
        -font $kNormalFont
#              -underline $nUnderline

      set sGroupName [lindex $lItem 3]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenuButton.mw $nItemNum
      }
      
  }
  if { [string compare $sType "radio" ] == 0 } {
      $isMenuButton.mw add radio \
        -label [lindex $lItem 1] \
        -command [lindex $lItem 2] \
        -variable [lindex $lItem 3] \
        -value [lindex $lItem 4] \
        -font $kNormalFont      

      set sGroupName [lindex $lItem 5]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenuButton.mw $nItemNum
      }
      
  }
  if { [string compare $sType "check" ] == 0 } {
      $isMenuButton.mw add check \
        -label [lindex $lItem 1] \
        -command [lindex $lItem 2] \
        -variable [lindex $lItem 3] \
        -font $kNormalFont      

      set sGroupName [lindex $lItem 4]
      if { [string compare $sGroupName ""] != 0 } {
    tkm_AddItemToMenuGroup $sGroupName $isMenuButton.mw $nItemNum
      }
      
  }
  if { [string compare $sType "separator"] == 0 } {
      $isMenuButton.mw add separator
  }

  incr nItemNum
    }
}

proc tkm_MakeEntryWithIncDecButtons { isFrame isText iVariable iSetFunc iDecFunc iIncFunc {inWidth 3} } {

    global kLabelFont kNormalFont
    global $iVariable

    frame $isFrame

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
    -width 4                 \
    -selectbackground green  \
    -insertbackground black
  bind $isFrame.ew <Return> $iSetFunc

  pack $isFrame.ew   \
    -side left \
    -anchor w  \
    -expand no
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

proc tkm_MakeFileSelector { isTop isText iVariable } {

    frame $isTop -width 200

    tkm_MakeNormalLabel $isTop.lw $isText

    tkm_MakeSmallLabel $isTop.lw2 "i.e. /home/myName/myDirectory/myFile.dat"
    tkm_MakeEntry $isTop.ew "" $iVariable
    
    pack $isTop.lw $isTop.lw2 $isTop.ew \
      -side top \
      -expand yes \
      -fill x \
      -anchor w
}

proc tkm_MakeDirectorySelector { isTop isText iVariable } {

    frame $isTop -width 200

    tkm_MakeNormalLabel $isTop.lw $isText

    tkm_MakeSmallLabel $isTop.lw2 "i.e. /home/myName/myDirectory"
    tkm_MakeEntry $isTop.ew "" $iVariable
    
    pack $isTop.lw $isTop.lw2 $isTop.ew \
      -side top \
      -expand yes \
      -fill x \
      -anchor w
}

proc tkm_AddItemToMenuGroup { isGroupName ifwMenuObject inMenuItemNum } {
    
    global glMenuGroups

    set glMenuGroups($isGroupName,0) $ifwMenuObject
    lappend glMenuGroups($isGroupName,1) $inMenuItemNum
}

proc tkm_SetMenuItemGroupStatus { isGroupName ibEnable } {

    global glMenuGroups

    # first item is a menu button
    set mbwMenu $glMenuGroups($isGroupName,0)
    
    # second item is a list of items
    foreach nMenuItem $glMenuGroups($isGroupName,1) {

  if { $ibEnable == 0 } {
      $mbwMenu entryconfigure $nMenuItem -state disabled
  } else {
      $mbwMenu entryconfigure $nMenuItem -state normal
  }
    }
}