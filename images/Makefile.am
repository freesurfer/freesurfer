#
# Makefile.am
# 

foodir=$(prefix)/lib/images

GIF=\
	icon_analysis_open.gif \
	icon_arrow_ccw.gif \
	icon_arrow_cw.gif \
	icon_arrow_down.gif \
	icon_arrow_expand_x.gif \
	icon_arrow_expand_y.gif \
	icon_arrow_left.gif \
	icon_arrow_right.gif \
	icon_arrow_rot_x_neg.gif \
	icon_arrow_rot_x_pos.gif \
	icon_arrow_rot_y_neg.gif \
	icon_arrow_rot_y_pos.gif \
	icon_arrow_rot_z_neg.gif \
	icon_arrow_rot_z_pos.gif \
	icon_arrow_shrink_x.gif \
	icon_arrow_shrink_y.gif \
	icon_arrow_up.gif \
	icon_aux_volume.gif \
	icon_brush_3d.gif \
	icon_brush_circle.gif \
	icon_brush_square.gif \
	icon_camera.gif \
	icon_color_label.gif \
	icon_colormap_grayscale.gif \
	icon_colormap_lut.gif \
	icon_color_picker.gif \
	icon_color_scalebar.gif \
	icon_control_point.gif \
	icon_cursor_goto.gif \
	icon_cursor_save.gif \
	icon_curv.gif \
	icon_curv_bw.gif \
	icon_cut_area.gif \
	icon_cut_clear.gif \
	icon_cut_closed_line.gif \
	icon_cut_line.gif \
	icon_cut_plane.gif \
	icon_disk.gif \
	icon_disk_multiview.gif \
	icon_disk_new_name.gif \
	icon_draw_line_closed.gif \
	icon_draw_line.gif \
	icon_edit_ctrlpts.gif \
	icon_edit_custom.gif \
	icon_edit_label.gif \
	icon_edit_parc.gif \
	icon_edit_volume.gif \
	icon_erase_label.gif \
	icon_erase_line.gif \
	icon_fill.gif \
	icon_fill_label.gif \
	icon_fill_roi.gif \
	icon_fill_volume.gif \
	icon_folder.gif \
	icon_home.gif \
	icon_label_filled.gif \
	icon_label_list.gif \
	icon_label_off.gif \
	icon_label_open.gif \
	icon_label_outline.gif \
	icon_label_save.gif \
	icon_label_table.gif \
	icon_label_to_marked.gif \
	icon_line_tool.gif \
	icon_linked_cursors.gif \
	icon_main_volume.gif \
	icon_marked_to_label.gif \
	icon_marker_crosshair.gif \
	icon_marker_diamond.gif \
	icon_matrix_erase.gif \
	icon_matrix_write.gif \
	icon_navigate.gif \
	icon_orientation_coronal.gif \
	icon_orientation_coronal_pdm.gif \
	icon_orientation_horizontal.gif \
	icon_orientation_horizontal_pdm.gif \
	icon_orientation_sagittal.gif \
	icon_orientation_sagittal_pdm.gif \
	icon_orientation_x.gif \
	icon_orientation_y.gif \
	icon_orientation_z.gif \
	icon_overlay.gif \
	icon_redraw.gif \
	icon_rotate_plane.gif \
	icon_seg_volume.gif \
	icon_selection_add.gif \
	icon_selection_subtract.gif \
	icon_snapshot_load.gif \
	icon_snapshot_save.gif \
	icon_stopwatch.gif \
	icon_surface_inflated.gif \
	icon_surface_main.gif \
	icon_surface_original.gif \
	icon_surface_pial.gif \
	icon_surface_white.gif \
	icon_table_open.gif \
	icon_view_31.gif \
	icon_view_mosaic.gif \
	icon_view_mosaic_pdm.gif \
	icon_view_multiple.gif \
	icon_view_multiple_pdm.gif \
	icon_view_single.gif \
	icon_view_single_pdm.gif \
	icon_zoom_in.gif \
	icon_zoom_out.gif

TIF=\
	icon_analysis_open.tif \
	icon_arrow_rot_x_neg.tif \
	icon_arrow_rot_x_pos.tif \
	icon_arrow_rot_y_neg.tif \
	icon_arrow_rot_y_pos.tif \
	icon_arrow_rot_z_neg.tif \
	icon_arrow_rot_z_pos.tif \
	icon_camera.tif \
	icon_colormap_grayscale.tif \
	icon_colormap_lut.tif \
	icon_curv_bw.tif \
	icon_disk_new_name.tif \
	icon_disk.tif \
	icon_disk_multiview.tif \
	icon_folder.tif \
	icon_home.tif \
	icon_label_open.tif \
	icon_label_save.tif \
	icon_marker_crosshair.tif \
	icon_matrix_erase.tif \
	icon_matrix_write.tif \
	icon_redraw.tif \
	icon_selection_add.tif \
	icon_selection_subtract.tif \
	icon_table_open.tif \
	icon_view_multiple.tif \
	icon_view_single.tif \
	icon_zoom_in.tif \
	icon_zoom_out.tif

foo_datafiles= $(GIF) $(TIF) \
	icon_zoom_out.png icon_label_list.xpm icon_fill_roi.png \
	icon_fill_volume.png icon_label_table.xpm

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
		$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
