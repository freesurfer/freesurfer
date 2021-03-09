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

set suffix ec.reg
set suffix ec.tal.reg
set suffix reg
#rotate_brain_y 190
set suffix sphere.funcreg
set suffix sphere.noinf.ec.reg
set suffix sphere.ec.one.reg
set suffix sphere.ec.one.test.reg
set suffix sphere.ec.reg
set suffix sphere.new.reg
set surf_dir $env(SUBJECTS_DIR)/$subject/surf
set_current_vertex_set 0
set fname $surf_dir/$hemi.sphere
read_surface_vertex_set 0 $fname
set target 1
set curv sno1_target_blur4.00
set curv $hemi.target${target}
set curv $hemi.target
read_binary_curv
redraw
set_current_vertex_set 0
set rgb $hemi.$subject.$suffix.target1.rgb
save_rgb
#set sulc inflated.H
#set cslope 17
#set sulc $hemi.0074blur2.00
#set sulc sno1_blur4.00
read_binary_sulc
set min_p 90
set max_p 105
set fname $surf_dir/$hemi.${suffix}
set fname $surf_dir/$hemi.${suffix}[format "%04d" 0]
read_surface_vertex_set 1 $fname
mark_vertex 69625 1
for {set pno $min_p } { $pno <= $max_p} {incr pno 5} {
		set_current_vertex_set 1
		set fname $surf_dir/$hemi.${suffix}[format "%04d" $pno]
		read_surface_vertex_set 1 $fname
		set rgb $hemi.$subject.${suffix}[format "%04d" $pno].rgb
		redraw
		set_current_vertex_set 1
		redraw
		save_rgb
		puts "finished saving rgb $rgb"
}
exit
