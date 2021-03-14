##
## make_mris_flatten_movie.tcl
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

set min_p 75
set max_p 90
set suffix occip.patch.new.flat
set offset .4
set surf_dir $env(SUBJECTS_DIR)/$subject/surf
read_binary_sulc
for {set pno $min_p } { $pno <= $max_p} {incr pno 5} {
		set patch $surf_dir/$hemi.${suffix}[format "%04d" $pno]
		read_binary_patch
    RestoreView
		set rgb $hemi.$subject.${suffix}[format "%04d" $pno].rgb
		redraw
		save_rgb
		puts "finished saving rgb $rgb"
}
#exit
