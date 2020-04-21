project(mris_compute_lgi)

install_configured(mris_compute_lgi DESTINATION bin)

install(FILES
  ComputeGeodesicProjection.m
  compute_lgi.m
  createMeshFacesOfVertex.m
  dijk.m
  find_corresponding_center_FSformat.m
  freesurfer_fread3.m
  freesurfer_read_surf.m
  getFaceArea.m
  getFacesArea.m
  getMeshArea.m
  getOrthogonalVector.m
  getVerticesAndFacesInSphere.m
  isInGeodesicROI.m
  isVertexInRadius.m
  MakeGeodesicOuterROI.m
  make_outer_surface.m
  make_roi_paths.m
  mesh_adjacency.m
  mesh_vertex_nearest.m
  pred2path.m
  PropagateGeodesic.m
  read_normals.m
  read_ROIlabel.m
  redo_lgi.m
  reorganize_verticeslist.m
  SearchProjectionOnPial.m
  transVertexToNormalAxisBase.m
  write_lgi.m
  write_path.m
  DESTINATION matlab
)
