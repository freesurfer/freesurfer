#! /usr/local/bin/tclsh7.4
############################################################################
#  Copyright (c) 1999 CorTechs - All rights reserved
############################################################################
set tk_strictMotif 1
set program tkmedit
set thisapp [tk appname]   ;# maybe multiple copies
### fonts from wrappers.tcl

set sclenx 150
set scleny 140

### medit vars set to abs paths at startup
# $home $session $subject $hemi $ext $dip $dec

### literal subdirs for setfile
set bemdir bem
set rgbdir rgb
set surfdir surf
set fsdir fs
set mridir mri

### file lists for setfile
set bemfiles { dip dec hpts htrans }
set rgbfiles { rgb }
set surffiles { insurf curv }
set fsfiles { fs fm }
set mrifiles { abs_imstem }

### default transform
set xrot 0
set yrot 0
set zrot 0
set xtrans 0
set ytrans 0

### default norm direction
set normdir 1

### default brush
set circleflag 1
set inplaneflag 1
set selectedpixval 1

### save panel
set userok 0

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
set parallelupdate TRUE

### slice planes
set cor 0
set hor 1
set sag 2
set changeplanelock 0

# kt
set gZoomCallNum 1
set g3DViewFlag 0

### source standard widget wrapper code
source $env(MRI_DIR)/lib/tcl/wrappers.tcl

############################################################################
proc setsession { name } {
  global session home subject
  if { [string range $name 0 0] == "~" } {
    if { [string length $name] == 1 } {
      set session $home/$subject
      return
    }
    if { [string length $name] < 3 || [string range $name 1 1] != "/"} {
      puts "bad tilde name for setsession"
      prompt
      return
    }
    set tildegone [string range $name 2 end]
    set subdir [file dirname $tildegone]
    set sessiontail [file tail $tildegone]
    if { $subdir == "." } {
      set session $home/$subject/$sessiontail
    } else {
      set session $home/$subject/$subdir/$sessiontail
    }
  } elseif { [string range $name 0 0] == "/" } {
    set session $name
  } else {
    puts "invalid name for setsession: must start with ~/ or /"
    prompt
    return
  }
}

proc setfile { varName value } {  ;# makes "varName"abbrev if it doesn't exist
  upvar $varName localvar               ;# setfile dec */bem/brain3d5-5.dec
  upvar ${varName}abbrev localabbrev    ;# setfile rgb ~/tmp/myrgb
  global home subject session imtype
  global bemdir rgbdir surfdir fsdir mridir
  global bemfiles rgbfiles surffiles fsfiles mrifiles
  ### set subdir using name of var and premade lists
  if { [lsearch -exact $bemfiles $varName] >= 0 } {
    set subdir $bemdir
  } elseif { [lsearch -exact $rgbfiles $varName] >= 0 } {
    set subdir $rgbdir
  } elseif { [lsearch -exact $surffiles $varName] >= 0 } {
    set subdir $surfdir
  } elseif { [lsearch -exact $fsfiles $varName] >= 0 } {
    set subdir $fsdir
  } elseif { [lsearch -exact $mrifiles $varName] >= 0 } {
    set subdir $mridir
  } else {
    puts "bad arg: don't know how to set file type: $varName"
    puts "  setfile {$bemfiles} value"
    puts "  setfile {$rgbfiles} value"
    puts "  setfile {$surffiles} value"
    puts "  setfile {$fsfiles} value"
    prompt
    return
  }
  ### do expansions and make fullname
  if { [string range $value 0 0] == "/"} {
    set fullname $value
  } elseif { [string range $value 0 0] == "~"} {
    if { [string length $value] > 1 && [string range $value 1 1] != "/" } {
      puts "tkmedit: can't reset subject home dir without restarting"
      prompt
      return
    }
    if { [string length $value] < 3 } {
      puts "tkmedit: no filename specified for setfile"
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
    if {$subdir == $bemdir} {
      set fullname $home/$subject/$subdir/$value
    } elseif { $subdir == $surfdir } {
      set fullname $home/$subject/$subdir/$value
    } elseif { $subdir == $fsdir } {
      set fullname $session/$subdir/$value
    } elseif { $subdir == $mridir } {
      if { $imtype == "local" } {
        set fullname $session/$value
      } else {
        set fullname $session/$subdir/$value
      }
    } elseif { $subdir == "." } {
      set fullname $session/$value
    } else {
      set fullname $session/$subdir/$value
    }
  }
  set localvar $fullname
  #puts $fullname

  ### attempt to re-abbrev (first ~, then *, else set absolute)
  set homename $home/$subject
  set homelen [string length $homename]
  set endhome [incr homelen -1]
  set begintail [incr homelen 1]
  if { $homename == [string range $fullname 0 $endhome] } {
    set localabbrev ~[string range $fullname $begintail end]
    return
  }
  set sessionlen [string length $session]
  set endsession [incr sessionlen -1]
  set begintail [incr sessionlen 1]
  if { $session == [string range $fullname 0 $endsession] } {
    set localabbrev *[string range $fullname $begintail end]
    return
  }
  set localabbrev $fullname
}

proc resettransform { } {
  global xrot yrot zrot xtrans ytrans
  set xrot 0; set yrot 0; set zrot 0
  set xtrans 0; set ytrans 0
}

proc fixdipdecname { varName index op } {
  global dip dec dip_spacing
  global dipabbrev decabbrev
  switch $op {
    w {    ;# on varName write
    ### read current abbrev in entry
    setfile dip $dipabbrev
    setfile dec $decabbrev
    set diptail [format "brain3d%d.dip" $dip_spacing]
    set dectail [format "brain3d%d.dec" $dip_spacing]
    ### work on absolute because file dirname does its own ~ expand!!
    set dip [file dirname $dip]/$diptail
    set dec [file dirname $dec]/$dectail
    ### update entry abbrev
    setfile dip $dip
    setfile dec $dec
    }
  }
}  

proc fixfocus { varName index op } {
  global plane cor hor sag changeplanelock
  if {$changeplanelock} { return }
  if {$plane==$cor} { focus .mri.main.left.view.pan.cor.bot.sc }
  if {$plane==$hor} { focus .mri.main.left.view.pan.hor.bot.sc }
  if {$plane==$sag} { focus .mri.main.left.view.pan.sag.bot.sc }
}

proc zoomcoords { varName index op } {  ;# trace nice, update real if changed
  global zf newimc newic newjc imc ic jc plane
  global dontzoom
  if {$dontzoom} { return }
  #set imc [expr $newimc*$zf]
  #set ic [expr $newic*$zf]
  #set jc [expr $newjc*$zf]

    # kt - to do coord conversions properly, call the VoxelToScreen cmd from
    # the c code. pass in our new* coords, which represent the coords on the
    # sliders, and use them to set our screen coords, used in the c code.
    # additionally, we need to flip ic here because the slider is in the
    # opposite orientation.
  # update - i guess we don't anymore. _shrug_
    #puts "zoom before: jc, ic, imc = $jc, $ic, $imc newjc, newic, newimc = $newjc, $newic, $newimc plane = $plane"
    # set gZoomCallNum [expr $gZoomCallNum+1]
    # set theTempNewIC [expr 255 - $newic]
    # puts "             tempic = $theTempNewIC"
    set jc [lindex [VoxelToScreen $newjc $newic $newimc $plane] 0]
    set ic [lindex [VoxelToScreen $newjc $newic $newimc $plane] 1]
    set imc [lindex [VoxelToScreen $newjc $newic $newimc $plane] 2]
    # set gZoomCallNum [expr $gZoomCallNum-1]
    #puts "zoom after:  jc, ic, imc = $jc, $ic, $imc newjc, newic, newimc = $newjc, $newic, $newimc plane = $plane"

    # end_kt
}

proc unzoomcoords { } {  ;# update nice (stop loop)
  global zf newimc newic newjc imc ic jc plane
  global dontzoom
  set dontzoom TRUE
    #set newimc [expr $imc/$zf]
    #set newic [expr $ic/$zf]
    #set newjc [expr $jc/$zf]
 
   # kt - do the same thing, other way around. set our slider coords
    # to the converted screen coords from our c code.
    #puts "unzoom before: jc, ic, imc = $jc, $ic, $imc newjc, newic, newimc = $newjc, $newic, $newimc plane = $plane"
    # set gZoomCallNum [expr $gZoomCallNum+1]
    
    set newjc [lindex [ScreenToVoxel $plane $jc $ic $imc] 0]
    set newic [lindex [ScreenToVoxel $plane $jc $ic $imc] 1]
    set newimc [lindex [ScreenToVoxel $plane $jc $ic $imc] 2]
    # set gZoomCallNum [expr $gZoomCallNum-1]
    #puts "unzoom after:  jc, ic, imc = $jc, $ic, $imc newjc, newic, newimc = $newjc, $newic, $newimc plane = $plane"

  set dontzoom FALSE

    # end_kt
}

# kt - this calls c code to save and restore the cursor before and after 
# setting the new plane. replaced all code that set the plane variable directly
# to use this function.
proc SetPlane { inNewPlane } {

    global plane jc ic imc

#    if { $plane == $inNewPlane } {
#  puts "SetPlane: new plane = plane, exiting.";
#  return;
#    }

    SaveCursorLocation;
    set plane $inNewPlane;
    RestoreCursorLocation;
    redraw;
    update idletasks;
    sendupdate;
}    

# kt - this does a similar thing except with switching between all3 view.
proc Set3DViewFlag { varName index op } {

    global all3flag g3DViewFlag

    SaveCursorLocation;
    set all3flag $g3DViewFlag;
    RestoreCursorLocation;
    redraw;
    update idletasks;
    sendupdate;
}


proc talupdate { varName index op } {
  global xtalairach ytalairach ztalairach
  coords_to_talairach
  set xtalairach $xtalairach   ;# touch tcl: trace doesn't report C updates
  set ytalairach $ytalairach
  set ztalairach $ztalairach
}

proc contupdate { varName index op } {  ;# trace zf, change binding
  global zf
  if {$zf < 3} {
    bind .mri.main.left.view.pan.cor.bot.sc <B1-Motion> { redraw }
    bind .mri.main.left.view.pan.sag.bot.sc <B1-Motion> { redraw }
    bind .mri.main.left.view.pan.hor.bot.sc <B1-Motion> { redraw }
  } else {
    bind .mri.main.left.view.pan.cor.bot.sc <B1-Motion> { }
    bind .mri.main.left.view.pan.sag.bot.sc <B1-Motion> { }
    bind .mri.main.left.view.pan.hor.bot.sc <B1-Motion> { }
  }
}

proc pixvaltitle { varName index op } {  ;# 0.2 msec
  global selectedpixval tkmeditinterface
  if {$tkmeditinterface == "micro"} {
    wm title . "pixelvalue = $selectedpixval"
  } else {
    wm title . "medit tools (pixelvalue = $selectedpixval)"
  }
}

proc sendupdate { } {
  global newimc newic newjc plane thisapp parallelupdate program
  set applist "{$program} {$program #2}"
  foreach app $applist {
    if { $app != $thisapp } {
      if {$parallelupdate} {
        catch { send -async $app "set plane $plane; set newimc $newimc; \
                                  set newic $newic; set newjc $newjc; redraw" }
      } else {
        catch { send -async $app "goto_point; set newimc $newimc; \
                                  set newic $newic; set newjc $newjc; redraw" }
      }
    }
  }
}

proc sendgoto { } {
  global plane thisapp parallelupdate program
  if {$parallelupdate} {
    set applist "{$program} {$program #2}"
    foreach app $applist {
      if { $app != $thisapp } {
        catch { send -async $app "goto_point; redraw; unzoomcoords" }
      }
    }
  }
}

proc gotoslice { slice } {
  global plane cor hor sag
  global newimc newic newjc
  if {$plane==$cor} { set newimc $slice }
  if {$plane==$hor} { set newic $slice }
  if {$plane==$sag} { set newjc $slice }
  redraw
}

proc changeslice { dir } {
  global zf newimc newic newjc imc ic jc
  global plane cor hor sag
  if {$dir == "up"} { upslice; redraw }
  if {$dir == "down"} { downslice; redraw }

    # kt - the c code handles the conversions now
    # if {$plane==$cor} { set newimc [expr $imc/$zf] }
    # if {$plane==$hor} { set newic [expr $ic/$zf] }
    # if {$plane==$sag} { set newjc [expr $jc/$zf] }
}

proc rotheadpts { angle } {
  global plane cor hor sag
  if {$plane==$cor} { rotate_brain_y [expr $angle*10.0] }
  if {$plane==$hor} { rotate_brain_z [expr -$angle*10.0] }
  if {$plane==$sag} { rotate_brain_x [expr $angle*10.0] }
}

proc transheadpts { dist axis } {
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

proc findsendto { } {
  global fulltksurfer
  set fulltksurfer ""
  catch { set fulltksurfer [lrange [exec ps -af | grep /tksurfer] 7 7] }
}

proc testclose { } {
  global editedimage abs_imstem
  if {$editedimage} {
    set resp [okclose ${abs_imstem}[format "%03d..." $editedimage]]
    if {$resp > 1} {write_images}
    if {$resp > 0} {exit}
  } else {
    exit
  }
}

proc histogram { } {
  global linearflag bwflag rgb
  #if { [exec uname] != "IRIX"} { return }
  if { [exec uname] != "IRIX" && [exec uname] != "IRIX64"} { return }
  set tmprgb /tmp/hist.rgb
  set linearflag TRUE
  set bwflag TRUE
  set rgbbak $rgb
  set rgb $tmprgb
  redraw
  save_rgb
  exec hist $tmprgb &   ;# else waits
  set linearflag FALSE
  set bwflag FALSE
  redraw
  exec rm $tmprgb  ;# unsafe
  set rgb $rgbbak
}

proc macro { } {
  global tkmeditinterface selectedpixval
  .mri.main.left.head.save.aSAVEIMG.bu config -text "SAVEIMG"
  .mri.main.left.head.pnt.aall.ck config -text "all"
  .mri.main.left.head.save.a3D.ck config -text "3D"
  pack forget .mri.main.left.head
  pack .mri.main.left.head -before .mri.main.left.view  ;# back to old nest
  set winlist { .meg .mri.main.right .mri.head }
  pack .mri.head             -before .mri.main -fill x
  pack .mri.main.right       -after .mri.main.left
  pack .mri.main.right.snorm -after .mri.main.right.fi
  pack .meg                  -after .mri
  set tkmeditinterface macro
  set selectedpixval $selectedpixval
}

proc mini+ { } {
  global tkmeditinterface selectedpixval
  .mri.main.left.head.save.aSAVEIMG.bu config -text "SAVEIMG"
  .mri.main.left.head.pnt.aall.ck config -text "all"
  .mri.main.left.head.save.a3D.ck config -text "3D"
  pack forget .mri.main.left.head
  pack .mri.main.left.head -before .mri.main.left.view
  pack .mri.head             -before .mri.main -fill x
  pack .mri.main.right       -after .mri.main.left
  pack .mri.main.right.snorm -after .mri.main.right.fi
  set winlist { .meg }
  foreach win $winlist { pack forget $win }
  set tkmeditinterface mini+
  set selectedpixval $selectedpixval
}

proc micro { } {
  global tkmeditinterface selectedpixval
  .mri.main.left.head.save.aSAVEIMG.bu config -text "SAVEIMG"
  .mri.main.left.head.save.a3D.ck config -text "3D"
  .mri.main.left.head.pnt.aall.ck config -text "all"
  pack forget .mri.main.left.head
  pack .mri.main.left.head -before .mri.main.left.view
  set winlist { .meg .mri.main.right .mri.head }
  foreach win $winlist { pack forget $win }
  set tkmeditinterface micro
  set selectedpixval $selectedpixval
}

proc mini { } {    ;# rearrange nesting
  global tkmeditinterface selectedpixval
  pack .mri.head       -before .mri.main
  pack .mri.main.right -after .mri.main.left
  pack forget .mri.main.left.head
  pack .mri.main.left.head
  pack .mri.main.right -in .mri.main.left.head    ;# redo nesting
  .mri.main.left.head.save.aSAVEIMG.bu config -text "SAVE IMGS"
  .mri.main.left.head.save.a3D.ck config -text "3Dbrush"
  .mri.main.left.head.pnt.aall.ck config -text "all3views"
  set winlist { .mri.head .meg .mri.main.right.snorm }
  foreach win $winlist { pack forget $win }
  set tkmeditinterface mini
  set selectedpixval $selectedpixval
}

############################################################################
wm title . "medit tools" ;# "tkmedit ($subject -- editable image: $imtype)"
wm geometry . +122+762    ;# +534+36
wm protocol . WM_DELETE_WINDOW testclose
wm resizable . 0 0

frame .mri
pack .mri -side left
  frame .mri.head
  pack .mri.head -side top -fill x
    frame .mri.head.pop
    pack .mri.head.pop -side left
    frame .mri.head.home
    pack .mri.head.home -side left
    frame .mri.head.env
    pack .mri.head.env -side left
  frame .mri.main
  pack .mri.main -side top
    # left column
    frame .mri.main.left
    pack .mri.main.left -side left
      frame .mri.main.left.head -borderwidth 2 -relief groove
      pack .mri.main.left.head -side top -fill x
        frame .mri.main.left.head.save
        pack .mri.main.left.head.save -side top
        frame .mri.main.left.head.pnt
        pack .mri.main.left.head.pnt -side top
      frame .mri.main.left.view
      pack .mri.main.left.view -side left
        frame .mri.main.left.view.pan
        pack .mri.main.left.view.pan -side top -fill y
        # three identical slice panels
        foreach v { cor sag hor } {
          frame .mri.main.left.view.pan.$v -borderwidth 2 -relief groove
          pack .mri.main.left.view.pan.$v -side top -fill x
            frame .mri.main.left.view.pan.$v.top
            pack .mri.main.left.view.pan.$v.top -side top
            frame .mri.main.left.view.pan.$v.bot
            pack .mri.main.left.view.pan.$v.bot -side top
        }
        # misc checks
        frame .mri.main.left.view.butt
        pack .mri.main.left.view.butt -side top
          frame .mri.main.left.view.butt.left -relief groove
          pack .mri.main.left.view.butt.left -side left
          frame .mri.main.left.view.butt.mid -width 1 -background darkgray 
          pack .mri.main.left.view.butt.mid -side left -fill y
          frame .mri.main.left.view.butt.right
          pack .mri.main.left.view.butt.right -side left
    # right column
    frame .mri.main.right ;# -borderwidth 2
    pack .mri.main.right -side left
      frame .mri.main.right.fi -borderwidth 2 -relief groove
      pack .mri.main.right.fi -side top
      frame .mri.main.right.snorm -borderwidth 2 -relief groove
      pack .mri.main.right.snorm -side top -fill x
        frame .mri.main.right.snorm.but
        pack .mri.main.right.snorm.but -side top
        frame .mri.main.right.snorm.bot
        pack .mri.main.right.snorm.bot -side top
          frame .mri.main.right.snorm.bot.lim
          pack .mri.main.right.snorm.bot.lim -side left
          frame .mri.main.right.snorm.bot.ffrac
          pack .mri.main.right.snorm.bot.ffrac -side left
          frame .mri.main.right.snorm.bot.dir
          pack .mri.main.right.snorm.bot.dir -side left
      frame .mri.main.right.wm ;# -borderwidth 2 -relief groove
      pack .mri.main.right.wm -side top -fill x
        frame .mri.main.right.wm.tru
        pack .mri.main.right.wm.tru -side top
          frame .mri.main.right.wm.tru.le
          pack .mri.main.right.wm.tru.le -side left
          frame .mri.main.right.wm.tru.mi
          pack .mri.main.right.wm.tru.mi -side left
          frame .mri.main.right.wm.tru.ri
          pack .mri.main.right.wm.tru.ri -side left
        frame .mri.main.right.wm.fil
        pack .mri.main.right.wm.fil -side top
          frame .mri.main.right.wm.fil.le
          pack .mri.main.right.wm.fil.le -side left
          frame .mri.main.right.wm.fil.ri
          pack .mri.main.right.wm.fil.ri -side left
        frame .mri.main.right.wm.thr
        pack .mri.main.right.wm.thr -side top
      frame .mri.main.right.im2 -borderwidth 2 -relief groove
      pack .mri.main.right.im2 -side top -fill both
      frame .mri.main.right.cmp
      pack .mri.main.right.cmp -side top -fill both

frame .meg
pack .meg -side left
  frame .meg.ti -borderwidth 1
  pack .meg.ti -side top
  frame .meg.dips -borderwidth 2 -relief groove
  pack .meg.dips -side top
    frame .meg.dips.spac
    pack .meg.dips.spac -side left
    frame .meg.dips.en
    pack .meg.dips.en -side left
  frame .meg.xform
  pack .meg.xform -side left
    frame .meg.xform.files
    pack .meg.xform.files -side top
    frame .meg.xform.pan
    pack .meg.xform.pan -side top
      # rotate plus fieldsign
      frame .meg.xform.pan.le
      pack .meg.xform.pan.le -side left -anchor n
        frame .meg.xform.pan.le.rg
        pack .meg.xform.pan.le.rg -side top -anchor n -fill x
          frame .meg.xform.pan.le.rg.rot -borderwidth 2 -relief groove
          pack .meg.xform.pan.le.rg.rot -side top -fill both
            frame .meg.xform.pan.le.rg.rot.top
            pack .meg.xform.pan.le.rg.rot.top -side top
            frame .meg.xform.pan.le.rg.rot.bot
            pack .meg.xform.pan.le.rg.rot.bot -side top -anchor center
        frame .meg.xform.pan.le.surf -borderwidth 2 -relief groove
        pack .meg.xform.pan.le.surf -side top -fill x
          frame .meg.xform.pan.le.surf.top
          pack .meg.xform.pan.le.surf.top -side top
            frame .meg.xform.pan.le.surf.top.le
            pack .meg.xform.pan.le.surf.top.le -side left
            frame .meg.xform.pan.le.surf.top.ri
            pack .meg.xform.pan.le.surf.top.ri -side left
          frame .meg.xform.pan.le.surf.bot
          pack .meg.xform.pan.le.surf.bot -side top
        frame .meg.xform.pan.le.ex
        pack .meg.xform.pan.le.ex -side top
          frame .meg.xform.pan.le.ex.le
          pack .meg.xform.pan.le.ex.le -side left
          frame .meg.xform.pan.le.ex.ri
          pack .meg.xform.pan.le.ex.ri -side left
      # translate panel
      frame .meg.xform.pan.tran -borderwidth 2 -relief groove
      pack .meg.xform.pan.tran -side left -anchor n -fill y
        frame .meg.xform.pan.tran.top
        pack .meg.xform.pan.tran.top -side top
        frame .meg.xform.pan.tran.bot
        pack .meg.xform.pan.tran.bot -side top
          frame .meg.xform.pan.tran.bot.la
          pack .meg.xform.pan.tran.bot.la -side right -anchor center
            frame .meg.xform.pan.tran.bot.la.ent
            pack .meg.xform.pan.tran.bot.la.ent -side top
            frame .meg.xform.pan.tran.bot.la.wm -borderwidth 2 -relief groove
            pack .meg.xform.pan.tran.bot.la.wm -side top
            frame .meg.xform.pan.tran.bot.la.en
            pack .meg.xform.pan.tran.bot.la.en -side top

############################################################################
### title 
set f .mri.head.pop
buttons $f "POP GL" { pop_gl_window } row 0 5
set f .mri.head.home
edlabval $f "home(~)" $home/$subject n 8 30
$f.home(~).e config -font $ffontb -state disabled
$f.home(~).e xview end
# readenv
set f .mri.head.env
buttons $f READENV { source $env(MRI_DIR)/lib/tcl/readenv.tcl; redraw } col 0 5
### main buttons (bigger font)
set f .mri.main.left.head.save
buttons $f "SAVEIMG" \
 { testreplace ${abs_imstem}[format "%03d" $editedimage] write_images } row 1 4
$f.aSAVEIMG.bu config -font $ffontbb

# brush

edlabval $f "rad" 0 n 4 2 row
$f.rad.e config -textvariable prad -font $sfont

checks $f "" "3D" inplaneflag row
$f.a3D.ck config -onvalue 0 -offvalue 1   ;# flip polarity

# point
set f .mri.main.left.head.pnt
# async so doesn't wait for answer from running offline tcl script
buttons $f "SEND PNT" { write_point; findsendto; \
  catch { send -async $fulltksurfer select_orig_vertex_coordinates } } row 2 2
buttons $f "GOTO PNT" {goto_point; redraw; unzoomcoords; sendgoto} row 2 2
checks $f "" "all" g3DViewFlag row
$f.aall.ck config -command { redraw } ;# kt- unzoom before all3

### cor: button (x=jc,y=imc,z=ic)
# on the top of the area...
set f .mri.main.left.view.pan.cor.top
# put a button that calls SetPlane
buttons $f "CORONAL" { SetPlane $cor;} row
# and an editable value
edlabval $f "yTal" 0 n 9 5 row
$f.yTal.e config -textvariable ytalairach -font $sfont
# when you press return in there, call this string of functions
bind $f.yTal.e <Return> \
  { SetPlane $cor; talairach_to_coords; redraw; unzoomcoords; sendupdate}
$f.yTal.la config -text "yTal P/A:"
# on the bottom...
set f .mri.main.left.view.pan.cor.bot
# make a slider linked to newimc
scale $f.sc -from $cormin -to $cormax -length $sclenx -variable newimc \
   -orient horizontal -tickinterval 127 -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side left
# when the button goes up, call SetPlane
bind $f.sc <ButtonRelease> { SetPlane $cor; }
# put another editable value
edlabval $f "none" 0 n 0 3 row
pack $f.none -side top   ;# repack to align with scale
# link it to newimc
$f.none.e config -textvariable newimc -font $sfont
# when you press return in it, call SetPlane
bind $f.none.e <Return> {SetPlane $cor;}

### sag: button
set f .mri.main.left.view.pan.sag.top
buttons $f SAGITTAL {SetPlane $sag;} row 2 8
edlabval $f "xTal" 0 n 9 5 row
$f.xTal.e config -textvariable xtalairach -font $sfont
bind $f.xTal.e <Return> \
  {SetPlane $sag; talairach_to_coords; redraw; unzoomcoords; sendupdate}
$f.xTal.la config -text "xTal R/L:"
### sag: slice
set f .mri.main.left.view.pan.sag.bot
scale $f.sc -from $sagmin -to $sagmax -length $sclenx -variable newjc \
   -orient horizontal -tickinterval 127 -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side left
bind $f.sc <ButtonRelease> { SetPlane $sag; redraw; }
#bind $f.sc <B1-Motion> { redraw }
edlabval $f "none" 0 n 0 3 row
pack $f.none -side top
$f.none.e config -textvariable newjc -font $sfont
bind $f.none.e <Return> {SetPlane $sag;}

### hor: button
set f .mri.main.left.view.pan.hor.top
buttons $f HORIZONTAL {SetPlane $hor;} row 2 6
edlabval $f "zTal" 0 n 9 5 row
$f.zTal.e config -textvariable ztalairach -font $sfont
bind $f.zTal.e <Return> \
  {SetPlane $hor; talairach_to_coords; redraw; unzoomcoords; sendupdate}
$f.zTal.la config -text "zTal I/S:"
### hor: slice
set f .mri.main.left.view.pan.hor.bot
scale $f.sc -from $hormin -to $hormax -length $sclenx -variable newic \
   -orient horizontal -tickinterval 127 -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side left
bind $f.sc <ButtonRelease> { SetPlane $hor; redraw; }
#bind $f.sc <B1-Motion> { redraw }
edlabval $f "none" 0 n 0 3 row
pack $f.none -side top
$f.none.e config -textvariable newic -font $sfont
bind $f.none.e <Return> {SetPlane $hor; }

### explicit plane focus
trace variable plane w fixfocus

### fix slice num on var update
trace variable newimc w zoomcoords
trace variable newic w zoomcoords
trace variable newjc w zoomcoords

### fix tal on nice update
trace variable newimc w talupdate
trace variable newic w talupdate
trace variable newjc w talupdate

### allow continuous update with zf=1
trace variable zf w contupdate

### current pixval
trace variable selectedpixval w pixvaltitle

### set all3 view
trace variable g3DViewFlag w Set3DViewFlag 

## trace the brush vars so we can call nice functions to do the updating
trace variable circleflag w SendBrushInfo
trace variable inplaneflag w SendBrushInfo
trace variable prad w SendBrushInfo

proc SendBrushInfo {  varName index op } {

    global circleflag inplaneflag prad

    SetBrushShape $circleflag
    SetBrush3DStatus [expr 1 - $inplaneflag]
    SetBrushRadius $prad
}

### misc buttons
set f .mri.main.left.view.butt.left
checks $f "" "overlay" do_overlay col
$f.aoverlay.ck config -command redraw
checks $f "" "surface" surfflag col 
$f.asurface.ck config -command redraw
set f .mri.main.left.view.butt.right
edlabval $f "contrast" 0 n 9 4
$f.contrast.e config -textvariable fsquash
bind $f.contrast.e <Return> { set_scale }
edlabval $f "midpoint" 0 n 9 4
$f.midpoint.e config -textvariable fthresh
bind $f.midpoint.e <Return> { set_scale }

set f .mri.main.right.fi
### surface
edlabval $f "surf" default r 6 15
setfile insurf $insurf  ;# make abbrev from C default
$f.surf.e config -textvariable insurfabbrev
$f.surf.br config -command { setfile insurf [.mri.main.right.fi.surf.e get]; \
                             read_binary_surf; redraw } -padx 7
bind .mri.main.right.fi.surf.e <Return> { .mri.main.right.fi.surf.br invoke}

### session (doesn't work anyway since 2nd imagedir not yet setfile'd)
#edlabval $f "sess" $session s 5 17
#$f.sess.e config -textvariable session
#$f.sess.bs config -command {setsession [.mri.main.right.fi.sess.e get]}

### output dataset
edlabval $f "outim" default w 6 15
setfile abs_imstem $abs_imstem   ;# make abbrev from C default
$f.outim.e config -textvariable abs_imstemabbrev
$f.outim.bw config -command { \
        setfile abs_imstem [.mri.main.right.fi.outim.e get]; \
        testreplace ${abs_imstem}[format "%03d" $editedimage] write_images } \
        -padx 7
bind .mri.main.right.fi.outim.e <Return> { .mri.main.right.fi.outim.bw invoke}

### save rgb button,field
edlabval $f "rgb" default  w 4 17
setfile rgb $rgb   ;# make abbrev
$f.rgb.e config -textvariable rgbabbrev
$f.rgb.bw config -command { setfile rgb [.mri.main.right.fi.rgb.e get]; \
                             testreplace $rgb save_rgb } -padx 7
### norm field
set f .mri.main.right.snorm.but
buttons $f TEST1 { norm_slice $normdir; \
                .mri.main.right.cmp.aCOMPARE.bu config -relief sunken } row 1 5
label $f.la1 -text " PieceWiseLinNorm " -font $ffontb
pack $f.la1 -side left
buttons $f ALL { norm_allslices $normdir } row 1 5
set f .mri.main.right.snorm.bot.lim
edlabval $f "lim3" 0 n 5 3
edlabval $f "lim2" 0 n 5 3
edlabval $f "lim1" 0 n 5 3
edlabval $f "lim0" 0 n 5 3
$f.lim3.e config -textvariable lim3
$f.lim2.e config -textvariable lim2
$f.lim1.e config -textvariable lim1
$f.lim0.e config -textvariable lim0
bind $f.lim3.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
bind $f.lim2.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
bind $f.lim1.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
bind $f.lim0.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
set f .mri.main.right.snorm.bot.ffrac
edlabval $f "ffrac3" 0 n 7 4
edlabval $f "ffrac2" 0 n 7 4
edlabval $f "ffrac1" 0 n 7 4
edlabval $f "ffrac0" 0 n 7 4
$f.ffrac3.e config -textvariable ffrac3
$f.ffrac2.e config -textvariable ffrac2
$f.ffrac1.e config -textvariable ffrac1
$f.ffrac0.e config -textvariable ffrac0
bind $f.ffrac3.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
bind $f.ffrac2.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
bind $f.ffrac1.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
bind $f.ffrac0.e <Return> { .mri.main.right.snorm.but.aTEST1.bu invoke }
set f .mri.main.right.snorm.bot.dir
radios $f "" "P/A" normdir 0 3 col
radios $f "" "I/S" normdir 1 3 col
radios $f "" "R/L" normdir 2 3 col
label $f.la1 -text "(radiol)" -font $sfont
pack $f.la1 -side top

### white matter test truncate
set f .mri.main.right.wm.tru.le
buttons $f "WMTRUNC" \
                {if {$truncflag} { set truncflag 0; \
                  .mri.main.right.wm.tru.le.aWMTRUNC.bu config -relief raised }\
                 else {set truncflag 1; \
                  .mri.main.right.wm.tru.le.aWMTRUNC.bu config -relief sunken};\
                 redraw } row 1 3
set f .mri.main.right.wm.tru.mi
edlabval $f "low" $white_lolim n 5 3
$f.low.e config -textvariable white_lolim
bind $f.low.e <Return> { set truncflag 0; \
                                 .mri.main.right.wm.tru.le.aWMTRUNC.bu invoke}
set f .mri.main.right.wm.tru.ri
edlabval $f "high" $white_hilim n 5 3
$f.high.e config -textvariable white_hilim
bind $f.high.e <Return> { set truncflag 0; \
                                  .mri.main.right.wm.tru.le.aWMTRUNC.bu invoke}
# wmfilter test
set f .mri.main.right.wm.fil
buttons $f "PLANEFILTER" { wmfilter_corslice } row 1 4
edlabval $f "grayhigh" $gray_hilim n 10 3
$f.grayhigh.e config -textvariable gray_hilim

# overlay thresh
set f .mri.main.right.wm.thr
edlabval $f fthr 0 n 5 3 row
edlabval $f fslope 0 n 7 3 row
edlabval $f fmid 0 n 5 3 row
$f.fthr.e config -textvariable f2thresh
$f.fslope.e config -textvariable fslope
$f.fmid.e config -textvariable fmid

### read 2nd image
set f .mri.main.right.im2
edlabval $f "secondimagedir" T1 r 13 8
$f.secondimagedir.br config -command { \
        read_second_images [.mri.main.right.im2.secondimagedir.e get]; \
        .mri.main.right.cmp.aCOMPARE.bu config -relief sunken }
bind $f.secondimagedir.e <Return> { \
                         .mri.main.right.im2.secondimagedir.br invoke }
$f.secondimagedir.la config -text "2nd imagedir:"

### image operation buttons
set f .mri.main.right.cmp
buttons $f "SAGMIRR" { mirror; redraw } row 2 5
buttons $f "HISTOGRAM" { histogram } row 1 4
buttons $f "COMPARE" {if {$drawsecondflag} { set drawsecondflag 0; \
                       .mri.main.right.cmp.aCOMPARE.bu config -relief raised } \
                      else {set drawsecondflag 1; \
                       .mri.main.right.cmp.aCOMPARE.bu config -relief sunken };\
                      redraw } row 2 4
$f.aCOMPARE.bu config -font $ffontbb

### title
set f .meg.ti
label $f.la -text "DIPOLE SAMPLE, MEG/EEG REGISTRATION, SURFPAINT" \
                                                      -font $ffontb -pady 0
pack $f.la -side top

### spacing and dipoles
set f .meg.dips.spac
label $f.la -text "3D DIPOLES" -font $ffontb -pady 2
pack $f.la -side top
edlabval $f "spacing" 0 n 8 2
$f.spacing.e config -textvariable dip_spacing
bind $f.spacing.e <Return> { fixdipdecname dip_spacing 0 w }
# dipoles
set f .meg.dips.en
edlabval $f  "dip"  default   w 4 19
edlabval $f  "dec"  default   w 4 19
setfile dip $dip    ;# make abbrevs
setfile dec $dec
$f.dip.e config -textvariable dipabbrev
$f.dec.e config -textvariable decabbrev
trace variable dip_spacing w fixdipdecname
# entries (don't abbrev $f inside { }: need curr $dipname but orig $f)
$f.dip.bw config -command { setfile dip [.meg.dips.en.dip.e get]; \
                                     testreplace $dip write_dipoles }
$f.dec.bw config -command { setfile dec [.meg.dips.en.dec.e get]; \
                                     testreplace $dec write_decimation }

### headpts read/write files
set f .meg.xform.files
edlabval $f "hpts"    default   r  7 22
edlabval $f "htrans"  default   rw 7 22
setfile hpts $hpts      ;# make abbrevs
setfile htrans $htrans
$f.hpts.e config -textvariable hptsabbrev
$f.htrans.e config -textvariable htransabbrev
# entries (can't abbrev $f inside { }: we need curr $dipname but orig $f)
$f.hpts.br config -command { setfile hpts [.meg.xform.files.hpts.e get]; \
                                        read_hpts; redraw }
$f.htrans.br config -command { setfile htrans [.meg.xform.files.htrans.e get]; \
                                        read_htrans; redraw }
$f.htrans.bw config -command { setfile htrans [.meg.xform.files.htrans.e get]; \
                                  testreplace $htrans write_htrans }
### rotate headpts
# title and horiz scale
set f .meg.xform.pan.le.rg.rot.top
label $f.title -text "ROTATE HDPTS (deg)" -font $ffontb -pady 0
pack $f.title -side top
scale $f.z -from 10 -to -10 -length $sclenx -variable zrot \
    -orient horizontal -tickinterval 10 -showvalue false -font $sfont \
    -width 11 -resolution 0.2
pack $f.z -side top
bind $f.z <ButtonRelease> { rotheadpts $zrot; resettransform; redraw }
# entry
set f .meg.xform.pan.le.rg.rot.bot
edlabval $f "zrot" 0 n 5 4
$f.zrot.e config -textvariable zrot -font $sfont -highlightthickness 1
$f.zrot.la config -font $sfont
bind $f.zrot.e <Return> { rotheadpts $zrot; resettransform; redraw }

### curv, fieldsign, linewidth
# label, linewidth
set f .meg.xform.pan.le.surf.top.le
label $f.title -text "SURFPAINT" -font $ffontb -pady 0
pack $f.title -side left
set f .meg.xform.pan.le.surf.top.ri
edlabval $f "linewidth" 0 n 10 3
$f.linewidth.e config -textvariable surflinewidth -highlightthickness 0
bind $f.linewidth.e <Return> { redraw }
# entries
set f .meg.xform.pan.le.surf.bot
edlabval $f "curv" 0 r 5 9
edlabval $f "fs" 0 r 3 11
edlabval $f "fm" 0 r 3 11
setfile curv $curv  ;# makes abbrevs
setfile fs $fs
setfile fm $fm
$f.curv.e config -textvariable curvabbrev
$f.fs.e config -textvariable fsabbrev
$f.fm.e config -textvariable fmabbrev
$f.curv.br config -command {setfile fs [.meg.xform.pan.le.surf.bot.curv.e get];\
                                          read_binary_curv; redraw }
$f.fs.br config -command { setfile fs [.meg.xform.pan.le.surf.bot.fs.e get];\
                                          read_fieldsign; redraw }
$f.fm.br config -command { setfile fm [.meg.xform.pan.le.surf.bot.fm.e get];\
                                          read_fsmask; redraw }
### extras
set f .meg.xform.pan.le.ex.le
buttons $f "SMOOTH 3D" { smooth_3d 1 } col 0 5
checks $f "" "fl" flossflag row
checks $f "" "sp" spackleflag row
$f.afl.ck config -highlightthickness 0
$f.asp.ck config -highlightthickness 0
set f .meg.xform.pan.le.ex.ri
checks $f "" "bwflag" bwflag col
checks $f "" "linearflg" linearflag col
$f.alinearflg.ck config -command redraw

### translate headpts
# title and horiz scale
set f .meg.xform.pan.tran.top
label $f.la -text "TRANSLATE HDPTS (mm)" -font $ffontb -pady 0
pack $f.la -side top
scale $f.x -from -25 -to 25 -length $sclenx -variable xtrans \
   -orient horizontal -tickinterval 25 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.x -side top
bind $f.x <ButtonRelease> { transheadpts $xtrans x; resettransform; redraw}
# vertical scale
set f .meg.xform.pan.tran.bot
scale $f.y -from 25 -to -25 -length $scleny -variable ytrans \
   -orient vertical -tickinterval 25 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.y -side left
bind $f.y <ButtonRelease> { transheadpts $ytrans y; resettransform; redraw}
# entries
set f .meg.xform.pan.tran.bot.la.ent
edlabval $f "x" 0 n 2 6
edlabval $f "y" 0 n 2 6
$f.x.e config -textvariable xtrans -font $sfont
$f.y.e config -textvariable ytrans -font $sfont
$f.x.la config -font $sfont
$f.y.la config -font $sfont
bind $f.x.e <Return> { transheadpts $xtrans x; resettransform; redraw}
bind $f.y.e <Return> { transheadpts $ytrans y; resettransform; redraw}
label $f.la -text " " -font $sfont  ;# space
pack  $f.la -side top

### update slice num's
unzoomcoords
set plane $cor
set zf $zf

############################################################################
### shortcut key bindings
# commands
bind . <Alt-f> { ".mri.main.left.head.pnt.aSEND PNT.bu" invoke }
bind . <Alt-r> { ".mri.main.left.head.pnt.aGOTO PNT.bu" invoke }
bind . <Alt-R> { setfile hpts [.meg.xform.files.hpts.e get]; read_hpts; \
                 setfile htrans [.meg.xform.files.htrans.e get]; \
                 read_htrans; redraw }
bind . <Alt-w> { .mri.main.left.head.save.aSAVEIMG.bu invoke }
bind . <Alt-W>  { setfile htrans [.meg.xform.files.htrans.e get]; \
                  testreplace $htrans write_htrans }
bind . <Alt-Key-M> { .main.left.view.butt.left.aoverlay.ck invoke }
bind . <Alt-m> { .mri.main.left.view.butt.left.aoverlay.ck invoke }
bind . <Alt-d> { .mri.main.left.view.butt.left.asurface.ck invoke }
# planes
bind . <Alt-y> { SetPlane $cor; }
bind . <Alt-x> { SetPlane $sag; }
bind . <Alt-z> { SetPlane $hor; }
# slices
bind . <Alt-Right> { changeslice up; sendupdate }
bind . <Alt-Left>  { changeslice down; sendupdate }
bind . <Alt-Up>    { changeslice up; sendupdate }
bind . <Alt-Down>  { changeslice down; sendupdate }
# contrast
bind . <Alt-asterisk> { set fsquash [expr $fsquash * 1.1]; set_scale }
bind . <Alt-slash> { set fsquash [expr $fsquash / 1.1]; set_scale }
bind . <Alt-plus> { set fthresh [expr $fthresh + 0.05]; set_scale }
bind . <Alt-minus> { set fthresh [expr $fthresh - 0.05]; set_scale }
# brush size
bind . <Alt-b> { set prad 0 }
bind . <Alt-B> { set prad [expr ($prad+1)%%5 ] }
bind . <Alt-v> { set pradtmp $prad; set prad $pradlast; set pradlast $pradtmp}
# translate small
bind . <Alt-semicolon> { transheadpts 0.5 x;  set overlay_mode 2; redraw }
bind . <Alt-l>         { transheadpts -0.5 x; set overlay_mode 2; redraw }
bind . <Alt-p>         { transheadpts 0.5 y;  set overlay_mode 2; redraw }
bind . <Alt-period>    { transheadpts -0.5 y; set overlay_mode 2; redraw }
# translate big
bind . <Alt-colon>     { transheadpts 2.5 x;  set overlay_mode 2; redraw }
bind . <Alt-L>         { transheadpts -2.5 x; set overlay_mode 2; redraw }
bind . <Alt-P>         { transheadpts 2.5 y;  set overlay_mode 2; redraw }
bind . <Alt-greater>   { transheadpts -2.5 y; set overlay_mode 2; redraw }
# rotate
bind . <Alt-braceleft>    { rotheadpts 5.0;  set overlay_mode 2; redraw }
bind . <Alt-braceright>   { rotheadpts -5.0; set overlay_mode 2; redraw }
bind . <Alt-bracketleft>  { rotheadpts 0.5;  set overlay_mode 2; redraw }
bind . <Alt-bracketright> { rotheadpts -0.5; set overlay_mode 2; redraw }
# second image
bind . <Alt-c> { .mri.main.right.cmp.aCOMPARE.bu invoke }
bind . <Alt-0> { .mri.main.right.cmp.aCOMPARE.bu invoke }
bind . <Alt-apostrophe>  { .mri.main.right.cmp.aCOMPARE.bu invoke }
# interface size
bind . <Control-F1> { micro }
bind . <Control-F2> { mini }
bind . <Control-F3> { mini+ }
bind . <Control-F4> { macro }

############################################################################
### right-click help
bind .mri.main.left.head.save.aSAVEIMG.bu <ButtonRelease-3> { helpwin saveimg }
bind .mri.main.left.head.save.rad.e <ButtonRelease-3> { helpwin brush }
bind .mri.main.left.head.save.a3D.ck <ButtonRelease-3> { helpwin inplane }
bind .mri.main.left.head.pnt.aall.ck <ButtonRelease-3> { helpwin all3 }
bind ".mri.main.left.head.pnt.aSEND PNT.bu" <ButtonRelease-3> {helpwin save_pnt}
bind ".mri.main.left.head.pnt.aGOTO PNT.bu" <ButtonRelease-3> {helpwin goto_pnt}
bind .mri.main.right.cmp.aSAGMIRR.bu <ButtonRelease-3> { helpwin sagmirror }
bind .mri.main.right.cmp.aHISTOGRAM.bu <ButtonRelease-3> { helpwin histogram }
bind ".meg.xform.pan.le.ex.le.aSMOOTH 3D.bu" <ButtonRelease-3> { helpwin smooth}
bind .mri.main.right.snorm.but.aTEST1.bu <ButtonRelease-3> { helpwin test1 }
bind .mri.main.right.snorm.but.aALL.bu <ButtonRelease-3> { helpwin all }
bind .mri.main.right.wm.tru.le.aWMTRUNC.bu <ButtonRelease-3> { helpwin wmtrunc }
bind .mri.main.right.wm.fil.aPLANEFILTER.bu <ButtonRelease-3> {helpwin wmfilter}
bind .mri.main.right.cmp.aCOMPARE.bu <ButtonRelease-3> { helpwin compare }
bind .mri.main.left.view.butt.right.contrast.e <ButtonRelease-3> \
                                                        { helpwin fsquash }
bind .mri.main.left.view.butt.right.midpoint.e <ButtonRelease-3> \
                                                        { helpwin fthresh }
bind .mri.main.right.im2.secondimagedir.br <ButtonRelease-3> { helpwin im2read }
bind .mri.main.right.im2.secondimagedir.e <ButtonRelease-3> { helpwin im2entry }
bind .mri.head.env.aREADENV.bu <ButtonRelease-3> { helpwin readenv }

############################################################################
puts "tkmedit.tcl: startup done"

# N.B.: w/command line script, these execute before gl window up
fixcolors
if [info exists env(tkmeditinterface)] {
  if {$env(tkmeditinterface) == "macro"} { macro }
  if {$env(tkmeditinterface) == "mini+"} { mini+ }
  if {$env(tkmeditinterface) == "mini"}  { mini }
  if {$env(tkmeditinterface) == "micro"} { micro }
} else {
  mini
  puts "tkmedit.tcl: default mini interface (to change: macro,mini+,mini,micro)"
  puts "tkmedit.tcl: or: setenv tkmeditinterface {macro,mini+,mini,micro}"
}

