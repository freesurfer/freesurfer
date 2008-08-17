if {[info exists env(LABEL)]} {
    labl_load $hemi.$env(LABEL)
}
if { $hemi == "lh" } {  rotate_brain_y 20
} else {  }
set gaLinkedVars(drawlabelflag) 1
set gaLinkedVars(labelstyle) 1
set truncphaseflag 1
set labeloutlinered 0
set labeloutlinegreen 0
set labeloutlineblue 1
set curv $env(CNAME)
read_binary_curv
#redraw
set drawlabelflag 1
labl_select -1
set rgb $env(HEMI).lateral.$env(SUFFIX).$env(SUBJECT).$env(ITER).rgb
save_rgb
puts "finished saving rgb $rgb"
#rotate_brain_y 180
#redraw
#set rgb $env(HEMI).medial.$env(SUFFIX).$env(SUBJECT).$env(ITER).rgb
#save_rgb
#puts "finished saving rgb $rgb"
exit
