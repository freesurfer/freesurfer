set suffix ec.reg
set surf_dir $env(SUBJECTS_DIR)/$subject/surf
set_current_vertex_set 0
set fname $surf_dir/$hemi.sphere
read_surface_vertex_set 0 $fname
set curv $hemi.target
read_binary_curv
redraw
set_current_vertex_set 0
set rgb $hemi.$subject.sphere.$suffix.target1.rgb
save_rgb
read_binary_sulc
#set curv $hemi.blur4
#read_binary_curv
set min_p 0
set max_p 57
set fname $surf_dir/$hemi.sphere.${suffix}[format "%04d" 0]
read_surface_vertex_set 1 $fname
for {set pno $min_p } { $pno <= $max_p} {incr pno 1} {
		set_current_vertex_set 1
		set fname $surf_dir/$hemi.sphere.${suffix}[format "%04d" $pno]
		read_surface_vertex_set 1 $fname
		set rgb $hemi.$subject.sphere.${suffix}[format "%04d" $pno].rgb
		redraw
		set_current_vertex_set 1
		redraw
		save_rgb
		puts "finished saving rgb $rgb"
}
