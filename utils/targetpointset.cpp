/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 * targetpointset.cpp - this is a class that manages target point sets
 * for surfaces. A target point set allows the user to place points where 
 * the surface should be so that the surface placement optimization is
 * encouraged to place the surface there.
 */

#undef X
#include "mrisurf.h"
#include "surfcluster.h"

int SurfacePointSet::PruneByAngle(void)
{
  std::vector<int> vtxlist2; // all vertices
  std::vector<int> npervtxlist2; // number of target points that map to each vertex
  std::vector<std::vector<double>> txyzlist2; // target coords for all vertices
  std::vector<std::vector<double>> dxyz2; // delta from vertex to target
  std::vector<double> dist2; // angle between surface normal and vector to target
  std::vector<double> angle2; // angle between surface normal and vector to target
  if(m_debug) fprintf(m_debug_fp,"Prune Angle before n = %d\n",(int)vtxlist.size());
  for(int n=0; n < vtxlist.size(); n++){
    double a = angle[n];
    if((a < 90 && a > AngleDegThresh) || (a > 90 && (180-a) > AngleDegThresh)) continue;
    vtxlist2.push_back(vtxlist[n]);
    npervtxlist2.push_back(npervtxlist[n]);
    txyzlist2.push_back(txyzlist[n]);
    dxyz2.push_back(dxyz[n]);
    angle2.push_back(angle[n]);
    dist2.push_back(dist[n]);
  }
  vtxlist = vtxlist2;
  npervtxlist = npervtxlist2;
  txyzlist = txyzlist2;
  dxyz = dxyz2;
  angle = angle2;
  dist = dist2;
  if(m_debug) fprintf(m_debug_fp,"Prune Angle before n = %d\n",(int)vtxlist.size());
  fflush(m_debug_fp);
  return(0);
}

int SurfacePointSet::PruneByDistance(void)
{
  std::vector<int> vtxlist2; // all vertices
  std::vector<int> npervtxlist2; // number of target points that map to each vertex
  std::vector<std::vector<double>> txyzlist2; // target coords for all vertices
  std::vector<std::vector<double>> dxyz2; // delta from vertex to target
  std::vector<double> dist2; // angle between surface normal and vector to target
  std::vector<double> angle2; // angle between surface normal and vector to target
  if(m_debug) fprintf(m_debug_fp,"Prune Dist before n = %d\n",(int)vtxlist.size());
  for(int n=0; n < vtxlist.size(); n++){
    if(dist[n] < DistMmThresh) continue;
    vtxlist2.push_back(vtxlist[n]);
    npervtxlist2.push_back(npervtxlist[n]);
    txyzlist2.push_back(txyzlist[n]);
    dxyz2.push_back(dxyz[n]);
    angle2.push_back(angle[n]);
    dist2.push_back(dist[n]);
  }
  vtxlist = vtxlist2;
  npervtxlist = npervtxlist2;
  txyzlist = txyzlist2;
  dxyz = dxyz2;
  angle = angle2;
  dist = dist2;
  if(m_debug) fprintf(m_debug_fp,"Prune Dist before n = %d\n",(int)vtxlist.size());
  fflush(m_debug_fp);
  return(0);
}

/*!
  \fn DTK_TRACK_SET * SurfacePointSet::ConvertToTrack(int nsteps)
  \brief Creates a trackviz Track for each vertex to show where its target is. This can
  be loaded into freeview for debugging.
 */
DTK_TRACK_SET * SurfacePointSet::ConvertToTrack(int nsteps)
{
  DTK_TRACK_SET *dtkset;

  dtkset = (DTK_TRACK_SET *)calloc(sizeof(DTK_TRACK_SET), 1);
  dtkset->mri_template = mri;
  dtkset->hdr = (DTK_HDR *)calloc(sizeof(DTK_HDR), 1);

  strcpy(dtkset->hdr->id_string,"TRACK");
  strcpy(dtkset->hdr->voxel_order,"LAS");
  dtkset->mri_template = mri; // should make a copy
  dtkset->hdr->invert_x = 0;
  dtkset->hdr->voxel_size[0] = mri->xsize;
  dtkset->hdr->voxel_size[1] = mri->ysize;
  dtkset->hdr->voxel_size[2] = mri->zsize;
  dtkset->hdr->dim[0] = mri->width;
  dtkset->hdr->dim[1] = mri->height;
  dtkset->hdr->dim[2] = mri->depth;
  dtkset->hdr->n_scalars = 1;
  strcpy(dtkset->hdr->scalar_name[0],"intensity");
  dtkset->hdr->n_properties = 0;
  dtkset->hdr->n_count = vtxlist.size();
  dtkset->trk = (DTK_TRACK**) calloc(sizeof(DTK_TRACK*),dtkset->hdr->n_count);
  for(int k=0; k < dtkset->hdr->n_count; k++){
    dtkset->trk[k] = DTKallocTrack(nsteps, 1, 0);
    VERTEX *v = &(surf->vertices[vtxlist[k]]);
    double dx =  txyzlist[k][0] - v->x;
    double dy =  txyzlist[k][1] - v->y;
    double dz =  txyzlist[k][2] - v->z;
    for(int step=0; step < nsteps; step++){
      dtkset->trk[k]->c[step] = v->x + dx*step/(nsteps-1) + mri->c_r;
      dtkset->trk[k]->r[step] = v->y + dy*step/(nsteps-1) + mri->c_a;
      dtkset->trk[k]->s[step] = v->z + dz*step/(nsteps-1) + mri->c_s;
    }
  }

  return(dtkset);
}

/*!
  \fn int SurfacePointSet::WriteAsPatch(const char *fname, int ndil)
  \brief Creates a surface patch of the vertices targeted in the pointset
  for debugging. The patch can be dilated to get more context.
 */
int SurfacePointSet::WriteAsPatch(const char *fname, int ndil)
{
  std::vector<double> ripcopy;
  for(int n=0; n < surf->nvertices; n++){
    ripcopy.push_back(surf->vertices[n].ripflag);
    surf->vertices[n].ripflag = 1; // rip everything
  }
  for(int n=0; n < vtxlist.size(); n++) {
    surf->vertices[vtxlist[n]].ripflag = 0; // unrip
    if(ndil == 0) continue;
    // Dilate rip 
    SURFHOPLIST *shl = SetSurfHopList(vtxlist[n], surf, ndil);
    // For this vertex, go through each hop
    for(int h=0; h < ndil; h++){
      // For this hop, go through each neighbor
      for(int n=0; n < shl->nperhop[h]; n++){
	int nbrvno = shl->vtxlist[h][n];
	surf->vertices[nbrvno].ripflag = 0; // unrip neighbor
      }
    }
    SurfHopListFree(&shl);
  }
  int err = MRISwritePatch(surf, fname);
  // Reinstate ripflag
  for(int n=0; n < surf->nvertices; n++) surf->vertices[n].ripflag = ripcopy[n];
  return(err);
}

/*!
  \fn MRI *SurfacePointSet::MakeMask(void)
  \brief Creates a mask of the targeted vertices for debugging.
 */
MRI *SurfacePointSet::MakeMask(void)
{
  MRI *mri = MRIalloc(surf->nvertices,1,1,MRI_UCHAR);
  for(int n=0; n < vtxlist.size(); n++){
    int vno = vtxlist[n];
    MRIsetVoxVal(mri,vno,0,0,0,1);
  }
  return(mri);
}

/*!
  \fn double SurfacePointSet::CostAndGrad(double weight, int ComputeGradient)
  \brief Computes the cost (SSE) and gradient for use in the surface 
  placement optimizer. 
 */
double SurfacePointSet::CostAndGrad(double weight, int ComputeGradient)
{
  if(weight==0) return(0);
  this->MapPointSet(); // inefficient because computed twice (with and without grad)
  double sse = 0;
  for(int n=0; n < vtxlist.size(); n++){
    int vno = vtxlist[n];
    VERTEX *v = &(surf->vertices[vno]);
    sse += (dist[n]*dist[n]);
    if(ComputeGradient) {
      // cf mrisComputeTargetLocationTerm()
      double norm = dist[n];
      v->dx += weight * dxyz[n][0] / norm;
      v->dy += weight * dxyz[n][1] / norm;
      v->dz += weight * dxyz[n][2] / norm;
    }
    //printf("%3d %6d (%6.3f %6.3f %6.3f) (%6.3f %6.3f %6.3f) %6.4lf %6.4lf\n",
    //	   i,vno,x,y,z,v->x,v->y,v->z,mag,sse);
  }
  if(m_debug) fprintf(m_debug_fp,"  targetpointset C&G w=%g sse = %g\n",weight,sse);
  return(sse);
}

/*!
  \fn int SurfacePointSet::Print(FILE *fp)
  \brief Print info about the point set to a file
 */
int SurfacePointSet::Print(FILE *fp)
{
  for(int n=0; n < vtxlist.size(); n++){
    VERTEX *v = &(surf->vertices[vtxlist[n]]);
    fprintf(fp,"%3d %5d %3d [%8.4lf %8.4lf %8.4lf] [%8.4lf %8.4lf %8.4lf] %6.4f\n",
	    n,vtxlist[n],npervtxlist[n],
	    v->x,v->y,v->z,txyzlist[n][0],txyzlist[n][1],txyzlist[n][2],dist[n]);
  }
  return(0);
}

/*!
  \fn int SurfacePointSet::Print(FILE *fp)
  \brief Create a point set from the vertices
 */
fsPointSet SurfacePointSet::VerticesToPointSet(void)
{
  fsPointSet vertexps;
  vertexps.vox2ras = "tkreg";
  for(int n=0; n < vtxlist.size(); n++){
    fsPointSet::Point c;
    VERTEX *v = &(surf->vertices[vtxlist[n]]);
    c.x = v->x;
    c.y = v->y;
    c.z = v->z;
    c.index = n+1; // freeview point sets are 1-based
    c.value = dist[n]; // distance to the target point
    c.count = npervtxlist[n];
    vertexps.add(c);
  }
  return(vertexps);
}
/*!
  \fn int SurfacePointSet::Print(FILE *fp)
  \brief Create a point set from the target points
 */
fsPointSet SurfacePointSet::ConvertToPointSet(void)
{
  fsPointSet targetps;
  targetps.vox2ras = "tkreg";
  for(int n=0; n < vtxlist.size(); n++){
    fsPointSet::Point c;
    c.x = txyzlist[n][0];
    c.y = txyzlist[n][1];
    c.z = txyzlist[n][2];
    c.index = n+1; // freeview is 0-based
    c.value = vtxlist[n];
    c.count = npervtxlist[n];
    targetps.add(c);
  }
  return(targetps);
}

/*!
  \fn int SurfacePointSet::MapPointSet(void)
  \brief This is the function that takes the input target point set
  and creates a new point set, one point per vertex.
 */
int SurfacePointSet::MapPointSet(void)
{
  json PointSet = *pPointSet;
  double max_spacing;
  int max_vno;

  bool bTkReg = (PointSet["vox2ras"].get<std::string>() == std::string("tkreg"));
  json js_pts = PointSet["points"];

  // The computation of the hash may be inefficient here as the same
  // hash may have been computed as some other point in the surface placement process
  // This call to MHT has a memory leak. Could probably compute max spacing once
  MRIScomputeVertexSpacingStats(surf, NULL, NULL, &max_spacing, NULL, &max_vno, CURRENT_VERTICES);
  MRIS_HASH_TABLE *hash = MHTcreateVertexTable_Resolution(surf, CURRENT_VERTICES, max_spacing);

  // The conversion to TkReg could/should be done when pointset is read in,
  // but there are generally only a few points, so not too costly, and I
  // have not figured out how to change the vox2ras string
  if(!bTkReg && mri == NULL) mri = MRIallocFromVolGeom(&(surf->vg), MRI_UCHAR, 1,1);

  psvtxlist.clear();
  vtxlist.clear();
  npervtxlist.clear();
  txyzlist.clear();
  dist.clear();
  dxyz.clear();
  angle.clear();

  // First, map the primary point set
  for(int i = 0; i < PointSet["points"].size(); i++) {
    double x,y,z;
    float distance;
    x = js_pts[i]["coordinates"]["x"].get<float>();
    y = js_pts[i]["coordinates"]["y"].get<float>();
    z = js_pts[i]["coordinates"]["z"].get<float>();
    if(!bTkReg){ // convert to TkReg
      double sx, sy, sz;
      MRIRASToSurfaceRAS(mri, x, y, z, &sx, &sy, &sz);
      x = sx; y = sy; z = sz;
    }
    // Find the closest vertex to this point
    int vno = MHTfindClosestVertexNoXYZ(hash, surf, x, y, z, &distance);
    if (vno < 0){
      //printf("Failed to find closest vertex in hash, using brute force\n");
      vno = MRISfindClosestVertex(surf, x, y, z, &distance, CURRENT_VERTICES);
    }
    // Check whether this vertex is already in the list
    std::vector<int>::iterator itr = std::find(psvtxlist.begin(), psvtxlist.end(), vno);
    if(itr != psvtxlist.end()) {
      // This vertex is the closest vertex to two or more points. Match it to
      // the mean of those points.
      int m = itr - psvtxlist.begin();
      txyzlist[m][0] += x;
      txyzlist[m][1] += y;
      txyzlist[m][2] += z;
      npervtxlist[m]++;
      if(m_debug) fprintf(m_debug_fp,"PA %2d %5d %g %g %g\n",i,vno,x,y,z);
      continue;
    }
    // This vertex not in the list
    psvtxlist.push_back(vno);
    vtxlist.push_back(vno);
    npervtxlist.push_back(1);
    std::vector<double> xyz;
    xyz.push_back(x);
    xyz.push_back(y);
    xyz.push_back(z);
    txyzlist.push_back(xyz);
    if(m_debug) fprintf(m_debug_fp,"P  %2d %5d %g %g %g\n",i,vno,x,y,z);
    // Add this vertex to the hoplist vector if not there already
    int hit=0; 
    for(int n=0; n < m_shl.size(); n++){
      if(vno == m_shl[n]->cvtx){
	hit = 1;
	break;
      }
    }
    if(hit == 0){
      SURFHOPLIST *shl = SetSurfHopList(vno, surf, m_nhops);
      m_shl.push_back(shl);
    }
  } // end loop over points in point set
  MHTfree(&hash);

  // Compute the average position for those that have more than one point
  for(int n=0; n < psvtxlist.size(); n++){
    if(npervtxlist[n] > 1) for(int k=0; k<3; k++) txyzlist[n][k] /= npervtxlist[n];
  }

  // Go through the neighbors of the primary vertices and map them. This expands the
  // list of vertices that will be affected
  for(int k=0; k < psvtxlist.size(); k++){
    int vno0 = psvtxlist[k];
    SURFHOPLIST *shl = m_shl[k];
    // For this vertex, go through each hop
    for(int h=0; h < m_nhops; h++){
      // For this hop, go through each neighbor
      for(int n=0; n < shl->nperhop[h]; n++){
	int vno = shl->vtxlist[h][n]; // neigbor vertex no
	// Check whether this vertex is a primary, skip if so
	if( std::find(psvtxlist.begin(), psvtxlist.end(), vno) != psvtxlist.end() ) continue;
	// Determine whether it is already in the list
	std::vector<int>::iterator itr = std::find(vtxlist.begin(), vtxlist.end(), vno);
	int nlist = vtxlist.size();
	if(itr != vtxlist.end()) {
	  // Already in the list, so update it
	  int m = itr - vtxlist.begin();
	  for(int c=0; c<3; c++) txyzlist[m][c] += txyzlist[k][c]; // really? no += for vector?
	  npervtxlist[m]++;
	  if(m_debug) fprintf(m_debug_fp,"A %3d %2d %d %3d %6d %6d  [%5.2lf,%5.2lf,%5.2lf] %3d\n",
	  	 nlist,k,h,n,vno,vno0,txyzlist[k][0],txyzlist[k][1],txyzlist[k][2],m);
	  continue;
	}
	// Not in the list, so push it
	if(m_debug) fprintf(m_debug_fp,"S %3d %2d %d %3d %6d %6d [%5.2lf,%5.2lf,%5.2lf]\n",
	     nlist,k,h,n,vno,vno0,txyzlist[k][0],txyzlist[k][1],txyzlist[k][2]);
	vtxlist.push_back(vno);
	txyzlist.push_back(txyzlist[k]);
	npervtxlist.push_back(1);
      } // neighbor
    } // hop
  } // point/primary vertex

  // Compute average when npervtx is more than one
  for(int n=psvtxlist.size(); n < vtxlist.size(); n++){
    if(npervtxlist[n] > 1) for(int k=0; k<3; k++) txyzlist[n][k] /= npervtxlist[n];
  }

  // The vertices above may form a cluster with a hole. If so, fill the hole
  if(m_fill_holes){
    std::vector<double> valcopy;
    for(int vno=0; vno < surf->nvertices; vno++) {
      valcopy.push_back(surf->vertices[vno].val); // make a copy so can reinstate it below
      surf->vertices[vno].val = 1;
    }
    for(int n=0; n < vtxlist.size(); n++)        surf->vertices[vtxlist[n]].val = 0;
    SURFCLUSTERSUM *scs;
    int NClusters;
    scs = sclustMapSurfClusters(surf,0.5,-1,0,0,&NClusters,NULL,NULL);
    if(m_debug) fprintf(m_debug_fp,"NClusters = %d\n",NClusters);
    if(NClusters > 1) {
      for(int vno=0; vno < surf->nvertices; vno++) {
	VERTEX *v = &(surf->vertices[vno]);
	if(v->undefval <= 1) continue;
	// Find closest point in current target point set
	double dmin=10e10;
	int nmin=0;
	for(int n=0; n < vtxlist.size(); n++){
	  double d = powf(v->x-txyzlist[n][0],2) + powf(v->y-txyzlist[n][1],2) + powf(v->z-txyzlist[n][2],2);
	  if(d < dmin){
	    dmin = d;
	    nmin = n;
	  }
	}
	if(m_debug) fprintf(m_debug_fp,"H %6d %3d %3d  [%5.2lf,%5.2lf,%5.2lf]\n",
	     vno,v->undefval,nmin,txyzlist[nmin][0],txyzlist[nmin][1],txyzlist[nmin][2]);
	vtxlist.push_back(vno);
	txyzlist.push_back(txyzlist[nmin]);
	npervtxlist.push_back(1);
      }
    }
    free(scs);
    // Restore the val field
    for(int vno=0; vno < surf->nvertices; vno++)  surf->vertices[vno].val = valcopy[vno];
  }

  // compute distances and angles from vertex to target
  for(int n=0; n < vtxlist.size(); n++){
    VERTEX *v = &(surf->vertices[vtxlist[n]]);
    // dxyz Points from vertex to target
    double dx = txyzlist[n][0] - v->x;
    double dy = txyzlist[n][1] - v->y;
    double dz = txyzlist[n][2] - v->z;
    std::vector<double> dd = {dx,dy,dz};
    dxyz.push_back(dd);
    double d = sqrt(dx*dx+dy*dy+dz*dz);
    dist.push_back(d);
    // Compute the angle
    double c = (dx*v->nx + dy*v->ny+ dz*v->nz)/d;
    if(c>+1) c = +1;
    if(c<-1) c = -1;
    double a = acos(c)*180/M_PI;
    angle.push_back(a);
  }

  // Remove vertices whose target is an at extreme angle relative to
  // the surface normal. This may or may not be a good idea given that
  // the surface may be twisted in the area.
  if(m_prune_by_angle) this->PruneByAngle();
  if(m_prune_by_dist)  this->PruneByDistance();

  fflush(m_debug_fp);

  return(0);
}
SurfacePointSet::~SurfacePointSet(){
  if(m_debug) printf("SurfacePointSet::Destructor\n");
  if(mri) MRIfree(&mri);
  if(m_debug_fp != stdout && m_debug_fp != stderr) fclose(m_debug_fp);
  for(int n=0; n < m_shl.size(); n++) SurfHopListFree(&m_shl[n]);
  if(m_debug) printf("SurfacePointSet::Destructor done\n");
}
