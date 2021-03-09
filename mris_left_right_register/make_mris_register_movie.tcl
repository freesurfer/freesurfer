##
## make_mris_register_movie.tcl
##
##
## Copyright © 2021
## The General Hospital Corporation (Boston, MA). 
## All rights reserved.
##
## Distribution, usage and copying of this software is covered under the
## terms found in the License Agreement file named 'COPYING' found in the
## FreeSurfer source code root directory, and duplicated here:
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

#rotate_brain_y 190
set suffix sphere.left_right
set surf_dir $env(SUBJECTS_DIR)/$subject/surf
set_current_vertex_set 0
set fname $surf_dir/$hemi.sphere
read_surface_vertex_set 0 $fname
set target 1
set curv sno1_target_blur4.00
set curv $hemi.target${target}
set curv $hemi.target${target}
set curv $hemi.${suffix}.target${target}
read_binary_curv
redraw
set_current_vertex_set 0
set rgb $hemi.lateral.$subject.$suffix.target1.rgb
save_rgb
puts "finished saving rgb $rgb"
#rotate_brain_y 180
#redraw
#set rgb $hemi.medial.$subject.$suffix.target1.rgb
#save_rgb
#puts "finished saving rgb $rgb"
RestoreView
#set sulc inflated.H
#set cslope 17
#set sulc $hemi.0074blur2.00
#set sulc sno1_blur4.00
if {$target == 0} {
    set curv inflated.H
    read_binary_curv
}
if {$target == 1} {
    read_binary_sulc
}
set min_p 8
set max_p 10
set fname $surf_dir/$hemi.${suffix}
set fname $surf_dir/$hemi.${suffix}[format "%04d" 0]
read_surface_vertex_set 1 $fname
#mark_vertex 69625 1
for {set pno $min_p } { $pno <= $max_p} {incr pno 1} {
    RestoreView
		set_current_vertex_set 1
		set fname $surf_dir/$hemi.${suffix}[format "%04d" $pno]
		read_surface_vertex_set 1 $fname
		set rgb $hemi.lateral.$subject.${suffix}[format "%04d" $pno].rgb
		set_current_vertex_set 1
		redraw
		save_rgb
    puts "finished saving rgb $rgb"
#    rotate_brain_y 180
#    redraw
#		set rgb $hemi.medial.$subject.${suffix}[format "%04d" $pno].rgb
#		save_rgb
#    puts "finished saving rgb $rgb"
}
exit
