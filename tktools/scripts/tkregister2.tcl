##
## tkregister2.tcl
##
##
## Copyright (C) 1996 Martin Sereno and Anders Dale
## Copyright (C) 2004-2011, CorTechs Labs, Inc. (La Jolla, CA) and
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
############################################################################

set program tkregister

### fonts from wrappers.tcl

### slider sizing
set sclenx 140
set scleny 120

### register vars set at startup
# $home $session $subject $registerdat $analysedat
# $tkrtitle

### literal file list,subdirs for setfile
set rgbdir rgb
set rgbfiles { rgb }

### others
set zrot 0.0
set xtrans 0.0
set ytrans 0.0
set xscale 100.0
set yscale 100.0
set zscale 100.0

### slice (slider!) limits
set cormax [expr $imnr1-1]
set cormin 0
set sagmax [expr $xdim/$zf-1]
set sagmin 0
set hormax [expr $ydim/$zf-1]
set hormin 0
set newimc 127
set newic 127
set newjc 127
set dontzoom FALSE

### slice planes
set cor 0
set hor 1
set sag 2

### overlay_modes
set target 1
set moveable 2

### events
set userok 0
set blinkflag FALSE
set initdelay $blinkdelay

### source standard widget wrapper code
source $env(FREESURFER_HOME)/lib/tcl/wrappers.tcl

############################################################################
proc setfile { varName value } {
  upvar $varName localvar
  upvar ${varName}abbrev localabbrev
  global home session subject              ;# setfile rgb ~/tmp/myrgb
  global rgbdir
  global rgbfiles 
  ## set subdir using name of var and premade lists
  if { [lsearch -exact $rgbfiles $varName] >= 0 } {
    set subdir $rgbdir
  } else {
    puts "bad arg: don't know how to set file type: $varName"
    puts "  setfile {$rgbfiles} value"
    prompt
    return
  }
  ## do expansions and make fullname
  if { [string range $value 0 0] == "/"} {
    set fullname $value
  } elseif { [string range $value 0 0] == "~"} {
    if { [string length $value] > 1 && [string range $value 1 1] != "/" } {
      puts "tkregister: can't reset subject home dir without restarting"
      prompt
      return
    }
    if { [string length $value] < 3 } {
      puts "tkregister: no filename specified for setfile"
      prompt
      return
    }
    set tildegone [string range $value 2 end]
    set subdir [file dirname $tildegone]      ;# overwrite, may be multilevel
    set filename [file tail $tildegone]
    if { $subdir == "." } {
      set fullname $home/$subject/$filename
    } else {
      set fullname $home/$subject/$subdir/$filename
    }
  } elseif { [string range $value 0 0] == "*"} {
    set stargone [string range $value 2 end]
    set fullname $session/$stargone
  } else {  ;# relative (guess session vs. subjects)
    if { $subdir == "." } {
      set fullname $session/$value
    } else {
      set fullname $session/$subdir/$value
    }
  }
  set localvar $fullname

  ### attempt to re-abbrev (first ~, then *, else set absolute)
  set homename $home/$subject
  set homelen [string length $homename]
  set endhome [incr homelen -1]
  set begintail [incr homelen 1]
  if { "$homename" == "[string range $fullname 0 $endhome]" } {
    set localabbrev ~[string range $fullname $begintail end]
    return
  }
  set sessionlen [string length $session]
  set endsession [incr sessionlen -1]
  set begintail [incr sessionlen 1]
  if { "$session" == "[string range $fullname 0 $endsession]" } {
    set localabbrev *[string range $fullname $begintail end]
    return
  }
  set localabbrev $fullname
}

proc resettransform { } {
  global zrot xtrans ytrans xscale yscale zscale
  set zrot 0.0
  set xtrans 0.0; set ytrans 0.0
  set xscale 100.0; set yscale 100.0; set zscale 100
}

proc fixfocus { varName index op } {
  global plane cor hor sag
  if {$plane==$cor} { focus .view.cor.bot.sc }
  if {$plane==$hor} { focus .view.hor.bot.sc }
  if {$plane==$sag} { focus .view.sag.bot.sc }
}

proc zoomcoords { varName index op } {  ;# trace nice, update real if changed
  global zf newimc newic newjc imc ic jc
  global dontzoom
  if {$dontzoom} { return }
  set imc [expr $newimc*$zf]
  set ic [expr $newic*$zf]
  set jc [expr $newjc*$zf]
}

proc unzoomcoords { plane } {  ;# update nice (stop loop)
  global zf newimc newic newjc imc ic jc
  global cor hor sag dontzoom
  set dontzoom TRUE
  set newimc [expr $imc/$zf]
  set newic [expr $ic/$zf]
  set newjc [expr $jc/$zf]
  set dontzoom FALSE
}

proc changeslice { dir plane } {
  global zf newimc newic newjc imc ic jc
  global cor hor sag
  if {$dir == "up"} { upslice }
  if {$dir == "down"} { downslice }
  if {$plane==$cor} { set newimc [expr $imc/$zf] }
  if {$plane==$hor} { set newic [expr $ic/$zf] }
  if {$plane==$sag} { set newjc [expr $jc/$zf] }
}

proc rotepi { angle } {
  global plane cor hor sag
  if {$plane==$cor} { rotate_brain_y [expr $angle*10.0] }
  if {$plane==$hor} { rotate_brain_z [expr -$angle*10.0] }
  if {$plane==$sag} { rotate_brain_x [expr $angle*10.0] }
}

proc transepi { dist axis } {
  global plane cor hor sag
  if {$axis=="x"} {
    if {$plane==$cor} { translate_brain_x [expr -$dist] }
    if {$plane==$hor} { translate_brain_x [expr -$dist] }
    if {$plane==$sag} { translate_brain_y $dist }
  }
  if {$axis=="y"} {
    if {$plane==$cor} { translate_brain_z $dist }
    if {$plane==$hor} { translate_brain_y $dist }
    if {$plane==$sag} { translate_brain_z $dist }
  }
}

proc scaleepi { factor axis } {
  global plane cor hor sag
  if {$axis=="x"} {
    if {$plane==$cor} { scale_brain_x [expr 100.0/$factor] }
    if {$plane==$hor} { scale_brain_x [expr 100.0/$factor] }
    if {$plane==$sag} { scale_brain_y [expr 100.0/$factor] }
  }
  if {$axis=="y"} {
    if {$plane==$cor} { scale_brain_z [expr 100.0/$factor] }
    if {$plane==$hor} { scale_brain_y [expr 100.0/$factor] }
    if {$plane==$sag} { scale_brain_z [expr 100.0/$factor] }
  }
}

proc comparebutton { } {
  global visible_mode last_visible_mode 
  global visible_plane last_visible_plane
  global target moveable overlay_mode
  if {$visible_plane != $last_visible_plane} { ;# change other plane to match
    if {$visible_mode == $target} { .mid.buff.mid.aMOVEABLE.bu invoke; return }
    if {$visible_mode == $moveable}  { .mid.buff.mid.aTARGET.bu invoke; return }
  }
  if {$visible_mode == $last_visible_mode} {
    if {$visible_mode == $target} { .mid.buff.mid.aMOVEABLE.bu invoke; return }
    if {$visible_mode == $moveable}  { .mid.buff.mid.aTARGET.bu invoke; return }
  }
  record_swapbuffers
}

proc showvisible { } {
  global visible_mode target moveable
  if {$visible_mode == $target} {
    .mid.buff.mid.aTARGET.bu config -relief sunken
    .mid.buff.mid.aMOVEABLE.bu config -relief raised
  }
  if {$visible_mode == $moveable} {
    .mid.buff.mid.aTARGET.bu config -relief raised
    .mid.buff.mid.aMOVEABLE.bu config -relief sunken
  }
}

proc findsendto { } {
  global fulltksurfer
  set fulltksurfer ""
  catch { set fulltksurfer [lrange [exec ps -af | grep /tksurfer] 7 7] }
}

proc testclose { } {
  global editedmatrix env registerdat
  if {$editedmatrix} {
    set resp [okclose $registerdat]
    if {$resp > 1} {write_reg}
    if {$resp > 0} {exit}
  } else {
    exit
  }
}

proc macro { } {  ;# restore old default
  global tkregisterinterface
  pack .head -before .view
  pack .view.main -before .view.cor -fill x
  pack .view.vals.left -before .view.vals.right -side left
  pack .mid -before .xform -side left
  pack .mid.buff.top -before .mid.buff.mid
  pack forget ".mid.buff.bot.aSAVE REG" .xform.tran.bot.la.ext
  pack forget .view.vals.bus
  .mid.buff.bot.aCOMPARE.bu config -pady 4 -padx 11
  pack .mid.buff.bot.aALIGN
  pack .mid.rgb -before .mid.buff
  pack .mid.scal -before .mid.rgb
  pack .xform.bot -after .xform.rot
  pack .xform.rot.bot -after .xform.rot.top
  set tkregisterinterface macro
}

proc mini { } {
  global tkregisterinterface ffontbb
  pack .view.vals.left -before .view.vals.right -side left
  pack .mid -before .xform -side left
  pack .mid.buff.top -before .mid.buff.mid
  pack .mid.scal -before .mid.buff -fill x
  pack .xform.tran.bot.la.ext
  pack forget .head .view.main .mid.rgb .xform.bot .mid.buff.top .xform.rot.bot
  pack forget .mid.buff.bot.aALIGN
  pack forget .view.vals.bus
  .mid.buff.bot.aCOMPARE.bu config -pady 3 -padx 4
  pack ".mid.buff.bot.aSAVE REG" -after .mid.buff.bot.aCOMPARE
  set tkregisterinterface mini
}

proc micro { } {
  global tkregisterinterface ffontbb
  pack forget .head .view.main .mid .xform.bot .mid.buff.top .xform.rot.bot
  pack forget .mid.scal .view.vals.left .xform.tran.bot.la.ext
  .mid.buff.bot.aCOMPARE.bu config -pady 3 -padx 4
  pack .view.vals.bus
  pack .xform.rot.bot -after .xform.rot.top
  set tkregisterinterface micro
}

proc mkalternates { } {
  global ffontbb ffontb sfont
  ### for mini
  set f .mid.buff.bot
  buttons $f "SAVE REG" {testreplace $registerdat write_reg} row 3 3
  "$f.aSAVE REG.bu" config -font $ffontbb
  pack forget ".mid.buff.bot.aSAVE REG"
  set f [frame .xform.tran.bot.la.ext]
  pack $f -side top
  frame $f.sp -height 55
  pack $f.sp -side top
  edlabval $f "zrot" 0 n 5 5 
  $f.zrot.e config -textvariable zrot -font $sfont
  $f.zrot.la config -font $sfont
  ### for micro
  set f [frame .view.vals.bus]
  pack $f -side left
  buttons $f "SAVE REG" {testreplace $registerdat write_reg} col 2 6
  "$f.aSAVE REG.bu" config -font $ffontb
  buttons $f "COMPARE" { } col 3 12
  $f.aCOMPARE.bu config -font $ffontbb
  bind $f.aCOMPARE.bu <ButtonPress-1> { comparebutton; set blinkflag TRUE }
  bind $f.aCOMPARE.bu <ButtonRelease-1> \
    { set blinkflag FALSE; set blinkdelay $initdelay; showvisible }
  bind $f.aCOMPARE.bu <ButtonRelease-3> { helpwin compare }
  pack forget .view.vals.bus
}

proc putaboveglwin { glwinx glwiny } {
  set geo [wm geometry .]
  set beg [expr [string first x $geo] + 1]
  set end [expr [string first + $geo] - 1]
  set intysize [string range $geo $beg $end]
  set intx $glwinx
  set inty [expr $glwiny - $intysize]
  set inty [expr $inty - 40]
  wm geometry . +${intx}+${inty}
}

############################################################################
wm title . "tkregister2: $tkrtitle"
wm geometry . +608+133    ;# +122+762 504x350+122+762
wm protocol . WM_DELETE_WINDOW testclose
wm resizable . 0 0

frame .head
pack .head -side top
  frame .head.pop
  pack .head.pop -side left
  frame .head.title
  pack .head.title -side left
  frame .head.env
  pack .head.env -side left

frame .view -borderwidth 1
pack .view -side left -anchor n
  # main
  frame .view.main -borderwidth 2 -relief groove
  pack .view.main -side top -fill x
    frame .view.main.reg
    pack .view.main.reg -side top
    frame .view.main.pnt
    pack .view.main.pnt -side top
  # three identical slice panels
  foreach v { cor sag hor } {
    frame .view.$v -borderwidth 2 -relief groove
    pack .view.$v -side top -fill x
      frame .view.$v.top
      pack .view.$v.top -side top
      frame .view.$v.bot
      pack .view.$v.bot -side top
  }
  # misc entries
  frame .view.vals
  pack .view.vals -side top
    frame .view.vals.left
    pack .view.vals.left -side left
    frame .view.vals.right
    pack .view.vals.right -side left

frame .mid -borderwidth 1
pack .mid -side left -anchor n
  # scale panel
  frame .mid.scal -borderwidth 2 -relief groove
  pack .mid.scal -side top -fill x
    frame .mid.scal.top
    pack .mid.scal.top -side top
    frame .mid.scal.bot
    pack .mid.scal.bot -side top
      frame .mid.scal.bot.la
      pack .mid.scal.bot.la -side right -anchor center
  # save rgb
  frame .mid.rgb -borderwidth 2 -relief groove
  pack .mid.rgb -side top -fill x
  # buffers
  frame .mid.buff -borderwidth 2 -relief groove
  pack .mid.buff -side top -fill x
    frame .mid.buff.top
    pack .mid.buff.top -side top
    frame .mid.buff.mid -borderwidth 1
    pack .mid.buff.mid -side top
    frame .mid.buff.bot
    pack .mid.buff.bot -side top

frame .xform -borderwidth 1
pack .xform -side left -anchor n
  # translate panel
  frame .xform.tran -borderwidth 2 -relief groove
  pack .xform.tran -side top
    frame .xform.tran.top
    frame .xform.tran.bot
    pack .xform.tran.top -side top 
    pack .xform.tran.bot -side top
      frame .xform.tran.bot.la
      pack .xform.tran.bot.la -side right -anchor center
  # rotate panel
  frame .xform.rot -borderwidth 2 -relief groove
  pack .xform.rot -side top
      frame .xform.rot.top
      pack .xform.rot.top -side top
      frame .xform.rot.bot -borderwidth 2
      pack .xform.rot.bot -side top -anchor center
  # last few
  frame .xform.bot -borderwidth 2
  pack .xform.bot -side top
    frame .xform.bot.a
    pack .xform.bot.a -side top
    frame .xform.bot.b
    pack .xform.bot.b -side top

############################################################################
### title 
set f .head.pop
buttons $f "POP GL" { pop_gl_window } row 0 5
set f .head.title
edlabval $f "scan" $session n 6 44
$f.scan.e config -font $ffontb -state disabled
$f.scan.e xview end
set f .head.env
buttons $f READENV { source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl; redraw } col 0 5

## main save panel
# save,read reg
set f .view.main.reg
buttons $f "SAVE REG" { testreplace $registerdat write_reg } row 4 5
"$f.aSAVE REG.bu" config -font $ffontbb
buttons $f "READ REG" { read_reg; \
     set overlay_mode $target; set updateflag TRUE; \
     set overlay_mode $moveable; set updateflag TRUE } row 2 5
# goto, save point
set f .view.main.pnt
#buttons $f "SEND PNT" { write_point } row 2 5
buttons $f "SEND PNT" { write_point; findsendto; \
         catch { send $fulltksurfer select_orig_vertex_coordinates } } row 2 5
buttons $f "GOTO PNT" { \
                goto_point; set updateflag TRUE; unzoomcoords $plane } row 2 5

### cor: button (suppress mode change after possible swap)
set f .view.cor.top
buttons $f CORONAL { set plane $cor; set overlay_mode $visible_mode; \
                     set updateflag TRUE; unzoomcoords $cor } row
edlabval $f "slice" 0 n 6 3
$f.slice.e config -textvariable newimc -font $sfont
bind $f.slice.e <Return> {set plane $cor; set overlay_mode $visible_mode; \
                          set updateflag TRUE; update idletasks}
### cor: slice
set f .view.cor.bot
scale $f.sc -from $cormin -to $cormax -length $sclenx -variable newimc \
   -orient horizontal -tickinterval 127 -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side top
bind $f.sc <ButtonRelease-1> {set plane $cor; set overlay_mode $visible_mode; \
                              set updateflag TRUE;update idletasks}
### sag: button
set f .view.sag.top
buttons $f SAGITTAL { set plane $sag; set overlay_mode $visible_mode; \
                      set updateflag TRUE; unzoomcoords $sag } row
edlabval $f "slice" 0 n 6 3
$f.slice.e config -textvariable newjc -font $sfont
bind $f.slice.e <Return> {set plane $sag; set overlay_mode $visible_mode; \
                          set updateflag TRUE; update idletasks}
### sag: slice
set f .view.sag.bot
scale $f.sc -from $sagmin -to $sagmax -length $sclenx -variable newjc \
   -orient horizontal -tickinterval 127 -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side top
bind $f.sc <ButtonRelease-1> {set plane $sag; set overlay_mode $visible_mode; \
                              set updateflag TRUE;update idletasks}
### hor: button
set f .view.hor.top
buttons $f HORIZONTAL { set plane $hor; set overlay_mode $visible_mode; \
                        set updateflag TRUE; unzoomcoords $hor } row
edlabval $f "slice" 0 n 6 3
$f.slice.e config -textvariable newic -font $sfont
bind $f.slice.e <Return> {set plane $hor; set overlay_mode $visible_mode; \
                          set updateflag TRUE; update idletasks}
### hor: slice
set f .view.hor.bot
scale $f.sc -from $hormin -to $hormax -length $sclenx -variable newic \
   -orient horizontal -tickinterval 127 -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side top
bind $f.sc <ButtonRelease-1> {set plane $hor; set overlay_mode $visible_mode; \
                              set updateflag TRUE;update idletasks}

### explicit plane focus
trace variable plane w fixfocus

### fix slice num on var update
trace variable newimc w zoomcoords
trace variable newic w zoomcoords
trace variable newjc w zoomcoords

### misc entries
# left
set f .view.vals.left
edlabval $f "contrast" 0 n 9 4
$f.contrast.e config -textvariable fsquash
bind $f.contrast.e <Return> { set_scale }
edlabval $f "midpoint" 0 n 9 4
$f.midpoint.e config -textvariable fthresh
bind $f.midpoint.e <Return> { set_scale }
# right
set f .view.vals.right
edlabval $f "fmov" 0 n 5 4
$f.fmov.e config -textvariable fscale_2
bind $f.fmov.e <Return> { set_scale; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
#buttons $f "BLUR" { blur 2.0 } col
checks $f "" "masktarg" maskflag col
$f.amasktarg.ck config -command { \
      if {$visible_mode != $overlay_mode} { set overlay_mode $visible_mode }; \
      set updateflag TRUE }

### scale epi ###
# title and horiz scale
set f .mid.scal.top
label $f.la -text "SCALE BRAIN (percent)" -font $ffontb -pady 0
pack $f.la -side top
scale $f.x -from 90 -to 110 -length $sclenx -variable xscale \
   -orient horizontal -tickinterval 20 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.x -side top
bind $f.x <ButtonRelease-1> { scaleepi $xscale x; resettransform; \
                              set overlay_mode $moveable; set updateflag TRUE }
# vertical scale
set f .mid.scal.bot
scale $f.y -from 110 -to 90 -length $scleny -variable yscale \
   -orient vertical -tickinterval 20 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.y -side left
bind $f.y <ButtonRelease-1> { scaleepi $yscale y; resettransform; \
                              set overlay_mode $moveable; set updateflag TRUE }
# entries
set f .mid.scal.bot.la
edlabval $f "x" 0 n 2 6
edlabval $f "y" 0 n 2 6
$f.x.e config -textvariable xscale -font $sfont
$f.y.e config -textvariable yscale -font $sfont
$f.x.la config -font $sfont
$f.y.la config -font $sfont
bind $f.x.e <Return> { scaleepi $xscale x; resettransform; \
                       set overlay_mode $moveable; set updateflag TRUE }
bind $f.y.e <Return> { scaleepi $yscale y; resettransform; \
                       set overlay_mode $moveable; set updateflag TRUE }

### save rgb button,field
set f .mid.rgb
setfile rgb $rgb  ;# make abbrev
#edlabval $f  "nm"   $rgb   n 3 19
#$f.nm.e config -textvariable rgbabbrev
entry $f.e -font $ffont -width 19 -textvariable rgbabbrev
$f.e config -selectbackground green -insertbackground black
pack $f.e -side top -fill x
buttons $f "SAVE RGB" { setfile rgb [.mid.rgb.e get]; \
                        testreplace $rgb save_rgb } col
### load buffers panel
# push
set f .mid.buff.top
label $f.title -text "Push on Front Buffer" -font $ffont -pady 0
pack $f.title -side top
set f .mid.buff.mid
buttons $f "TARGET" { set overlay_mode $target; set updateflag TRUE } row 2 7
buttons $f "MOVEABLE" { set overlay_mode $moveable; set updateflag TRUE} row 2 7
# why this works?? (ButtonRelease fails; showvisible fails in Press/Release/Cmd)
bind $f.aTARGET.bu <ButtonPress-1> { \
                         .mid.buff.mid.aTARGET.bu config -relief sunken; \
                         .mid.buff.mid.aMOVEABLE.bu config -relief raised }
bind $f.aMOVEABLE.bu <ButtonPress-1> { \
                         .mid.buff.mid.aTARGET.bu config -relief raised; \
                         .mid.buff.mid.aMOVEABLE.bu config -relief sunken }
# compare
set f .mid.buff.bot
buttons $f "COMPARE" { } row 4    ;# N.B.: button command on release
$f.aCOMPARE.bu config -font $ffontbb
bind $f.aCOMPARE.bu <ButtonPress-1> { comparebutton; set blinkflag TRUE } 
bind $f.aCOMPARE.bu <ButtonRelease-1> { set blinkflag FALSE; \
                                        set blinkdelay $initdelay; showvisible }
buttons $f "ALIGN" \
   { align_points; set overlay_mode $moveable; set updateflag TRUE } row 4 5
$f.aALIGN.bu config -font $ffontbb

### translate epi ###
# title and horiz scale
set f .xform.tran.top
label $f.la -text "TRANSLATE BRAIN (mm)" -font $ffontb -pady 0
pack $f.la -side top
scale $f.x -from 25 -to -25 -length $sclenx -variable xtrans \
   -orient horizontal -tickinterval 25 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.x -side top
bind $f.x <ButtonRelease-1> { transepi $xtrans x; resettransform; \
                              set overlay_mode $moveable; set updateflag TRUE }
# vertical scale
set f .xform.tran.bot
scale $f.y -from -25 -to 25 -length $scleny -variable ytrans \
   -orient vertical -tickinterval 25 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.y -side left
bind $f.y <ButtonRelease-1> { transepi $ytrans y; resettransform; \
                              set overlay_mode $moveable; set updateflag TRUE }
# entries
set f .xform.tran.bot.la
edlabval $f "x" 0 n 2 6
edlabval $f "y" 0 n 2 6
$f.x.e config -textvariable xtrans -font $sfont
$f.y.e config -textvariable ytrans -font $sfont
$f.x.la config -font $sfont
$f.y.la config -font $sfont
bind $f.x.e <Return> { transepi $xtrans x; resettransform; \
                       set overlay_mode $moveable; set updateflag TRUE }
bind $f.y.e <Return> { transepi $ytrans y; resettransform; \
                       set overlay_mode $moveable; set updateflag TRUE }

### rotate epi ###
set f .xform.rot.top
label $f.title -text "ROTATE BRAIN (deg)" -font $ffontb -pady 0
pack $f.title -side top
scale $f.z -from -30 -to 30 -length $sclenx -variable zrot \
    -orient horizontal -tickinterval 30 -showvalue false -font $sfont \
    -width 11 -resolution 0.2
pack $f.z -side top
bind $f.z <ButtonRelease-1> { rotepi $zrot; resettransform; \
                              set overlay_mode $moveable; set updateflag TRUE }
# entry
set f .xform.rot.bot
edlabval $f "zrot" 0 n 5 5
$f.zrot.e config -textvariable zrot -font $sfont
$f.zrot.la config -font $sfont
bind $f.zrot.e <Return> { rotepi $zrot; resettransform; \
                          set overlay_mode $moveable; set updateflag TRUE }
### sagittal mirror (TODO: fix blur hack for non-sgi)
set f .xform.bot.a
#buttons $f "SAGITTAL MIRROR" { mirror_brain; \
                    set overlay_mode $moveable; set updateflag TRUE } row
buttons $f "SAGMIRR" { mirror_brain; \
                    set overlay_mode $moveable; set updateflag TRUE } row
buttons $f "BLUR" { if { [exec uname] == "IRIX"} { blur 2.0 } } row

### blinktime
set f .xform.bot.b
edlabval $f "blinktime" 0 n 12 3
$f.blinktime.e config -textvariable blinktime

### for mini,micro
mkalternates

### update slice num's, etc
unzoomcoords $cor
unzoomcoords $sag
unzoomcoords $hor
set plane $cor
showvisible

############################################################################
### shortcut key bindings--to find keysyms: bind . <KeyPress> { puts %K }
# commands
bind . <Alt-A> { align_points; set overlay_mode $moveable; set updateflag TRUE }
bind . <Alt-Key-M> { {.xform.bot.aSAGITTAL MIRROR.bu} invoke }
bind . <Alt-r> { {.view.main.reg.aREAD REG.bu} invoke }
bind . <Alt-t> { {.view.main.reg.aSAVE REG.bu} invoke }
bind . <Alt-T> { {.view.main.reg.aSAVE REG.bu} invoke }
# planes
bind . <Alt-x> { set plane $sag; set updateflag TRUE }
bind . <Alt-y> { set plane $hor; set updateflag TRUE }
bind . <Alt-z> { set plane $cor; set updateflag TRUE }
# slices
bind . <Alt-Right> { changeslice up $plane; set updateflag TRUE }
bind . <Alt-Left> { changeslice down $plane; set updateflag TRUE }
# buffers
bind . <Alt-0> { comparebutton }
bind . <Alt-Key-1> { set overlay_mode $target; set updateflag TRUE } ;# 1-3 Key!
bind . <Alt-Key-2> { set overlay_mode $moveable; set updateflag TRUE }
# contrast
bind . <Alt-asterisk> { set fsquash [expr $fsquash * 1.1]; set_scale }
bind . <Alt-slash> { set fsquash [expr $fsquash / 1.1]; set_scale }
bind . <Alt-plus> { set fthresh [expr $fthresh + 0.05]; set_scale }
bind . <Alt-minus> { set fthresh [expr $fthresh - 0.05]; set_scale }
bind . <Alt-Up> { set fscale_2 [expr $fscale_2 * 1.5]; set updateflag TRUE }
bind . <Alt-Down> { set fscale_2 [expr $fscale_2 / 1.5]; set updateflag TRUE }
# translate
bind . <Alt-semicolon> {transepi -0.5 x; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-l>         {transepi 0.5 x;  set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-p>         {transepi -0.5 y; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-period>    {transepi 0.5 y;  set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-colon>     {transepi -2.5 x; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-L>         {transepi 2.5 x;  set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-P>         {transepi -2.5 y; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-greater>   {transepi 2.5 y;  set overlay_mode $moveable; \
                                                           set updateflag TRUE}
# rotate
bind . <Alt-braceleft>    {rotepi -0.5; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-braceright>   {rotepi 0.5;  set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-bracketleft>  {rotepi -5.0; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-bracketright> {rotepi 5.0;  set overlay_mode $moveable; \
                                                           set updateflag TRUE}
# scale
bind . <Alt-KP_Right> {scaleepi 102 x;   set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-KP_Left>  {scaleepi 98.04 x; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-KP_Up>    {scaleepi 102 y;   set overlay_mode $moveable; \
                                                           set updateflag TRUE}
bind . <Alt-KP_Down>  {scaleepi 98.04 y; set overlay_mode $moveable; \
                                                           set updateflag TRUE}
# interface size
bind . <Control-F1> { micro }
bind . <Control-F2> { mini }
bind . <Control-F3> { macro }

############################################################################
### right-click help
bind .mid.buff.bot.aALIGN.bu <ButtonRelease-3> { helpwin align }
bind .mid.buff.bot.aCOMPARE.bu <ButtonRelease-3> { helpwin compare }
bind .view.vals.right.amasktarg.ck <ButtonRelease-3> { helpwin masktarg }
bind .head.title.scan.e <ButtonRelease-3> { \
                   helpwin $env(HOME)/subjects/template/NOTES 80 30 other }

############################################################################
puts "tkregister.tcl: startup done"

fixcolors
if [info exists env(tkregisterinterface)] {
  if {$env(tkregisterinterface) == "macro"} { macro }
  if {$env(tkregisterinterface) == "mini"}  { mini }
} else {
  macro
  puts "tkregister.tcl: default macro interface (to change: macro,mini,micro)"
  puts "tkregister.tcl: or: setenv tkregisterinterface {macro,mini,micro}"
}

