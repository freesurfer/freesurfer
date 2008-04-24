## 
## Makefile.am 
##

dist_bin_SCRIPTS = mris_compute_lgi

foodir=$(prefix)/matlab
foo_DATA=\
	ComputeGeodesicProjection.m \
	compute_lgi.m \
	createMeshFacesOfVertex.m \
	dijk.m \
	find_corresponding_center_FSformat.m \
	freesurfer_fread3.m \
	freesurfer_read_surf.m \
	getFaceArea.m \
	getFacesArea.m \
	getMeshArea.m \
	getOrthogonalVector.m \
	getVerticesAndFacesInSphere.m \
	isInGeodesicROI.m \
	isVertexInRadius.m \
	MakeGeodesicOuterROI.m \
	make_outer_surface.m \
	make_roi_paths.m \
	mesh_adjacency.m \
	mesh_vertex_nearest.m \
	pred2path.m \
	PropagateGeodesic.m \
	read_normals.m \
	read_ROIlabel.m \
	redo_lgi.m \
	reorganize_verticeslist.m \
	SearchProjectionOnPial.m \
	transVertexToNormalAxisBase.m \
	write_lgi.m \
	write_path.m

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be \
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
