source $env(MRI_DIR)/lib/tcl/wrappers.tcl

set kLabelFont   $pfont
set kNormalFont  $mfont
set kHighlightBgColor white

proc tkm_MakeBigLabel { isTop isText } {

    global kLabelFont

    frame $isTop

    label $isTop.label \
      -text $isText \
      -font $kLabelFont

    pack $isTop.label
}

proc tkm_MakeNormalLabel { isTop isText } {

    global kNormalFont

    frame $isTop

    label $isTop.label \
      -text $isText \
      -font $kNormalFont

    pack $isTop.label
}

proc tkm_MakeActiveLabel { isTop isText isVariable {inWidth -1} } {

    global kNormalFont kHighlightBgColor

    frame $isTop

    label $isTop.lw \
      -text $isText \
      -font $kNormalFont

    entry $isTop.ew                       \
      -textvariable $isVariable     \
      -width $inWidth               \
      -font $kNormalFont            \
      -state disabled               \
      -highlightbackground $kHighlightBgColor    \
      -relief flat 

    pack $isTop.lw $isTop.ew \
      -side left \
      -anchor w
}

proc tkm_MakeCheckbox { isTop isText iVariable iSetFunction } {

    global kLabelFont kNormalFont kHighlightBgColor
    global $iVariable

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

proc tkm_MakeEntry { isFrame isText iVariable {inWidth -1} {iSetFunction ""} } {

    global kLabelFont kNormalFont
    
    frame $isFrame

    label $isFrame.lwLabel \
      -text $isText \
      -font $kNormalFont

    entry $isFrame.ewEntry \
      -font $kNormalFont \
      -textvariable $iVariable \
      -width $inWidth \
      -selectbackground green \
      -insertbackground black

    if { $iSetFunction != "" } {
  bind $isFrame.ewEntry <Return> "$iSetFunction [set $iVariable]"
    }

    pack $isFrame.lwLabel \
      -side left \
      -anchor w

    pack $isFrame.ewEntry \
      -side right \
      -anchor e
}

proc tkm_MakeEntryWithIncDecButtons { isFrame isText iVariable iSetFunc iDecFunc iIncFunc } {

    global kLabelFont kNormalFont
    global $iVariable

    frame $isFrame

    label $isFrame.lwText \
      -text $isText \
      -font $kNormalFont

    entry $isFrame.ewEntry \
      -textvariable $iVariable \
      -width 3
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