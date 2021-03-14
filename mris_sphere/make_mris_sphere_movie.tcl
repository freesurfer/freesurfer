##
## make_mris_sphere_movie.tcl
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

set suffix hippo.smooth.sphere2
#rotate_brain_y 190
set surf_dir $env(SUBJECTS_DIR)/$subject/surf
set_current_vertex_set 0
redraw
set_current_vertex_set 0
read_binary_sulc
scale_brain 7
set min_p 0
set max_p 450
set fname $surf_dir/$hemi.${suffix}
set fname $surf_dir/$hemi.${suffix}[format "%04d" 0]
read_surface_vertex_set 1 $fname
set verticesflag 1
#mark_vertex 69625 1
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
