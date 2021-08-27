/**
 * @brief routines for growing clusters on the surface
 *
 * routines for growing clusters on the surface
 * based on intensity thresholds and area threshold. Note: this
 * makes use of the undefval in the MRI_SURFACE structure.
 */
/*
 * Original Author: Doug Greve
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
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "matrix.h"
#include "mri.h"
#include "resample.h"
#include "timer.h"

// This must be included prior to volcluster.c (I think)
#define SURFCLUSTER_SRC
#include "surfcluster.h"
#undef SURFCLUSTER_SRC

#include "volcluster.h"

static int sclustCompare(const void *a, const void *b);

/* ------------------------------------------------------------
   sclustMapSurfClusters() - grows a clusters on the surface.  The
   cluster is a list of contiguous vertices that that meet the
   threshold criteria. The cluster does not exist as a list at this
   point. Rather, the clusters are mapped using using the undefval
   element of the MRI_SURF structure. If a vertex meets the cluster
   criteria, then undefval is set to the cluster number.
   ------------------------------------------------------------ */
SCS *sclustMapSurfClusters(MRI_SURFACE *Surf, float thmin, float thmax, int thsign, 
			   float minarea, int *nClusters, MATRIX *XFM, MRI *fwhmmap)
{
  SCS *scs, *scs_sorted;
  int vtx, vtx_inrange, vtx_clustno, CurrentClusterNo;
  float vtx_val, ClusterArea;
  int nVtxsInCluster;

  /* initialized all cluster numbers to 0 */
  for (vtx = 0; vtx < Surf->nvertices; vtx++) Surf->vertices[vtx].undefval = 0; /* overloads this elem of struct */

  /* Go through each vertex looking for one that meets the threshold
     criteria and has not been previously assigned to a cluster.
     When found, grow it out. */
  CurrentClusterNo = 1;
  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_val = Surf->vertices[vtx].val;
    vtx_clustno = Surf->vertices[vtx].undefval;
    vtx_inrange = clustValueInRange(vtx_val, thmin, thmax, thsign);

    if (vtx_clustno == 0 && vtx_inrange) {
      sclustGrowSurfCluster(CurrentClusterNo, vtx, Surf, thmin, thmax, thsign);
      if (minarea > 0) {
        /* If the cluster does not meet the area criteria, delete it */
        ClusterArea = sclustSurfaceArea(CurrentClusterNo, Surf, &nVtxsInCluster);
        if (ClusterArea < minarea) {
          sclustZeroSurfaceClusterNo(CurrentClusterNo, Surf);
          continue;
        }
      }
      CurrentClusterNo++;
    }
  }

  *nClusters = CurrentClusterNo - 1;
  if (*nClusters == 0) return (NULL);

  /* Get a summary of the clusters */
  scs = SurfClusterSummary(Surf, XFM, nClusters, fwhmmap);

  /* Sort the clusters by descending maxval */
  scs_sorted = SortSurfClusterSum(scs, *nClusters);

  if (Gdiag_no > 1) {
    printf("--- Surface Cluster Summary (unsorted) ---------------\n");
    DumpSurfClusterSum(stdout, scs, *nClusters);
    printf("---------- sorted ---------------\n");
    DumpSurfClusterSum(stdout, scs_sorted, *nClusters);
  }

  /* Remap the cluster numbers to match the sorted */
  sclustReMap(Surf, *nClusters, scs_sorted);

  free(scs);

  return (scs_sorted);
}
/* ------------------------------------------------------------
   sclustGrowSurfCluster() - grows a cluster on the surface from
   the SeedVtx. The cluster is a list of vertices that are
   contiguous with the seed vertex and that meet the threshold
   criteria. The cluster map itself is defined using the
   undefval of the MRI_SURF structure. If a vertex meets the
   cluster criteria, then undefval is set to the ClusterNo.
   The ClustNo cannot be 0.
   ------------------------------------------------------------ */
int sclustGrowSurfCluster(int ClusterNo, int SeedVtx, MRI_SURFACE *Surf, float thmin, float thmax, int thsign)
{
  int nbr, nbr_vtx, nbr_inrange, nbr_clustno;
  float nbr_val;

  if (ClusterNo == 0) {
    printf("ERROR: clustGrowSurfCluster(): ClusterNo is 0\n");
    return (1);
  }

  Surf->vertices[SeedVtx].undefval = ClusterNo;

  for (nbr = 0; nbr < Surf->vertices_topology[SeedVtx].vnum; nbr++) {
    nbr_vtx = Surf->vertices_topology[SeedVtx].v[nbr];
    nbr_clustno = Surf->vertices[nbr_vtx].undefval;
    if (nbr_clustno != 0) continue;
    nbr_val = Surf->vertices[nbr_vtx].val;
    if (fabs(nbr_val) < thmin) continue;
    nbr_inrange = clustValueInRange(nbr_val, thmin, thmax, thsign);
    if (!nbr_inrange) continue;
    sclustGrowSurfCluster(ClusterNo, nbr_vtx, Surf, thmin, thmax, thsign);
  }
  return (0);
}
/*----------------------------------------------------------------
  sclustSurfaceArea() - computes the surface area (in mm^2) of a
  cluster. Note:   MRIScomputeMetricProperties() must have been
  run on the surface.
  ----------------------------------------------------------------*/
float sclustSurfaceArea(int ClusterNo, MRI_SURFACE *Surf, int *nvtxs)
{
  int vtx, vtx_clusterno;
  float ClusterArea;

  *nvtxs = 0;
  ClusterArea = 0.0;
  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (vtx_clusterno != ClusterNo) continue;
    if (!Surf->group_avg_vtxarea_loaded)
      ClusterArea += Surf->vertices[vtx].area;
    else
      ClusterArea += Surf->vertices[vtx].group_avg_area;
    (*nvtxs)++;
  }

  if (Surf->group_avg_surface_area > 0 && !Surf->group_avg_vtxarea_loaded) {
    // In Dec 2008, a bug was found in this section of code.  The
    // above line read only: if(Surf->group_avg_surface_area > 0) This
    // caused the group vertex area ajustment to be applied twice
    // because the section of code immediately above has already
    // applied it if the group vertex area was already loaded (which
    // it always was). This caused the surface area to be too big (by
    // about 20-25% for fsaverage). This was fixed by adding:
    //   && !Surf->group_avg_vtxarea_loaded

    // This function is called by both mri_glmfit and mri_surfcluster.
    // To indicate this fix, a new global variable was created called
    // FixSurfClusterArea, which is set to 1. The mere presence of
    // this variable implies that this bug was fixed.  The variable
    // exists so that the CSD created by mri_glmfit can record the
    // fact that the cluster area is correct. When mri_surfcluster
    // reads in the CSD, the presense of this flag in the CSD file
    // will indicate that the area is correct. If it is not correct,
    // then mri_surfcluster will exit with error. This assures that an
    // old CSD file will not be used with the new mri_surfcluster.

    // This will not prevent new CSD files from being used with an old
    // version of mri_surfcluster. However, the CSD format was also
    // changed to be incompatible with the old CSD reader.

    // Always do this now (4/9/10)
    ClusterArea *= (Surf->group_avg_surface_area / Surf->total_area);
    // if (getenv("FIX_VERTEX_AREA") != NULL)
    // ClusterArea *= (Surf->group_avg_surface_area/Surf->total_area);
  }

  return (ClusterArea);
}
/*----------------------------------------------------------------
  float sclustWeight() - computes the cluster "weight", defined as
  the sum of the values in the cluster. If mri != NULL, the value
  is obtained from the mri structure. If mri==NULL, the the val
  field in Surf is used. If UseArea==1, then the value at a vertex
  is weighted by the area at the vertex.
  ----------------------------------------------------------------*/
float sclustWeight(int ClusterNo, MRI_SURFACE *Surf, MRI *mri, int UseArea)
{
  int vtx, vtx_clusterno;
  double ClusterWeight, h, vtxarea;

  ClusterWeight = 0.0;
  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (vtx_clusterno != ClusterNo) continue;
    if (mri == NULL)
      h = Surf->vertices[vtx].val;
    else
      h = MRIgetVoxVal(mri, vtx, 0, 0, 0);
    if (UseArea) {
      if (!Surf->group_avg_vtxarea_loaded)
        vtxarea = Surf->vertices[vtx].area;
      else
        vtxarea = Surf->vertices[vtx].group_avg_area;
      if (Surf->group_avg_surface_area > 0 && !Surf->group_avg_vtxarea_loaded)
        vtxarea *= (Surf->group_avg_surface_area / Surf->total_area);
      h *= vtxarea;
    }
    ClusterWeight += h;
  }
  return (ClusterWeight);
}

/*----------------------------------------------------------------
  sclustSurfaceMax() - returns the maximum intensity value of
  inside a given cluster and the vertex at which it occured.
----------------------------------------------------------------*/
float sclustSurfaceMax(int ClusterNo, MRI_SURFACE *Surf, int *vtxmax)
{
  int vtx, vtx_clusterno, first_hit;
  float vtx_val, vtx_val_max = 0;

  first_hit = 1;
  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (vtx_clusterno != ClusterNo) continue;

    vtx_val = Surf->vertices[vtx].val;
    if (first_hit) {
      vtx_val_max = vtx_val;
      *vtxmax = vtx;
      first_hit = 0;
      continue;
    }

    if (fabs(vtx_val) > fabs(vtx_val_max)) {
      vtx_val_max = vtx_val;
      *vtxmax = vtx;
    }
  }

  return (vtx_val_max);
}
/*----------------------------------------------------------------
  sclustSurfaceCentroid() - returns the centroid of a cluster.
----------------------------------------------------------------*/
int sclustSurfaceCentroid(const int ClusterNo, const MRI_SURFACE *Surf, double *xyz)
{
  int vtx, vtx_clusterno, nvtx;
  float xsum, ysum, zsum;
  nvtx = 0;
  xsum = 0;
  ysum = 0;
  zsum = 0;
  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (vtx_clusterno != ClusterNo) continue;
    xsum += Surf->vertices[vtx].x;
    ysum += Surf->vertices[vtx].y;
    zsum += Surf->vertices[vtx].z;
    nvtx++;
  }
  xyz[0] = xsum / nvtx;
  xyz[1] = ysum / nvtx;
  xyz[2] = zsum / nvtx;

  return (0);
}
/*----------------------------------------------------------------
  sclustZeroSurfaceClusterNo() - finds all the vertices with
  cluster number equal to ClusterNo and sets the cluster number
  to zero (cluster number is the undefval member of the surface
  structure). Nothing is done to the surface value. This function
  is good for pruning clusters that do not meet some other
  criteria (eg, area threshold).
  ----------------------------------------------------------------*/
float sclustZeroSurfaceClusterNo(int ClusterNo, MRI_SURFACE *Surf)
{
  int vtx, vtx_clusterno;

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (vtx_clusterno == ClusterNo) Surf->vertices[vtx].undefval = 0;
  }

  return (0);
}
/*----------------------------------------------------------------
  sclustZeroSurfaceNonClusters() - zeros the value of all the vertices
  that are not assocated with a cluster. The cluster number is the
  undefval member of the surface structure.
  ----------------------------------------------------------------*/
float sclustZeroSurfaceNonClusters(MRI_SURFACE *Surf)
{
  int vtx, vtx_clusterno;

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (vtx_clusterno == 0) Surf->vertices[vtx].val = 0.0;
  }

  return (0);
}
/*----------------------------------------------------------------
  sclustSetSurfaceClusterToClusterNo() - sets the value of a vertex to the
  cluster number. The cluster number is the undefval member of the
  surface structure.
  ----------------------------------------------------------------*/
float sclustSetSurfaceValToClusterNo(MRI_SURFACE *Surf)
{
  int vtx, vtx_clusterno;

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    Surf->vertices[vtx].val = vtx_clusterno;
  }

  return (0);
}
/*----------------------------------------------------------------
  sclustSetSurfaceClusterToCWP() - sets the value of a vertex to
  -log10(cluster-wise pvalue).
  ----------------------------------------------------------------*/
float sclustSetSurfaceValToCWP(MRI_SURFACE *Surf, SCS *scs)
{
  int vtx, vtx_clusterno;
  float val;

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;

    if (vtx_clusterno == 0)
      val = 0;
    else {
      val = scs[vtx_clusterno - 1].pval_clusterwise;
      if (val == 0.0)
        val = 50;
      else
        val = -log10(val);
      val = val * SIGN(scs[vtx_clusterno - 1].maxval);
    }
    Surf->vertices[vtx].val = val;
  }

  return (0);
}
/*----------------------------------------------------------------
  sclustCountClusters() - counts the number of clusters. Really
  just returns the largest cluster number, which will be the
  number of clusters if there are no holes.
  ----------------------------------------------------------------*/
float sclustCountClusters(MRI_SURFACE *Surf)
{
  int vtx, vtx_clusterno, maxclusterno;

  maxclusterno = 0;
  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    if (maxclusterno < vtx_clusterno) maxclusterno = vtx_clusterno;
  }
  return (maxclusterno);
}

/*----------------------------------------------------------------
  SurfClusterSummary() - (was "Fast") gives identical results as
    SurfClusterSummaryOld() but much, much faster.
  ----------------------------------------------------------------*/
SCS *SurfClusterSummary(MRI_SURFACE *Surf, MATRIX *T, int *nClusters, MRI *fwhmmap)
{
  int n, vtx, clusterno;
  SURFCLUSTERSUM *scs;
  MATRIX *xyz, *xyzxfm;
  float vtxarea, vtxval;
  int msecTime;
  double *weightvtx, *weightarea;  // to be consistent with orig code
  double fwhm;
  VERTEX *v;
  int ClusterUseAvgVertexArea=0;
  double avgvertexarea = 0, fwhmmean2 = 0;

  if(Surf->group_avg_vtxarea_loaded)
    avgvertexarea = Surf->group_avg_surface_area/Surf->nvertices;
  else
    avgvertexarea = Surf->total_area/Surf->nvertices;

  if(getenv("FS_CLUSTER_USE_AVG_VERTEX_AREA") != NULL){
    // When setenv FS_CLUSTER_USE_AVG_VERTEX_AREA 1, this computes the cluster
    // size as the vertex count * average area rather than the sum of the 
    // of each vertex area
    sscanf(getenv("FS_CLUSTER_USE_AVG_VERTEX_AREA"),"%d",&ClusterUseAvgVertexArea);
  }
  if(Gdiag_no > 0){
    printf("ClusterUseAvgVertexArea = %d, avgvertexarea = %g\n",ClusterUseAvgVertexArea,avgvertexarea);
    fflush(stdout);
  }

  Timer mytimer;

  *nClusters = sclustCountClusters(Surf);
  if (*nClusters == 0) return (NULL);

  xyz = MatrixAlloc(4, 1, MATRIX_REAL);
  xyz->rptr[4][1] = 1;
  xyzxfm = MatrixAlloc(4, 1, MATRIX_REAL);

  scs = (SCS *)calloc(*nClusters, sizeof(SCS));
  weightvtx = (double *)calloc(*nClusters, sizeof(double));
  weightarea = (double *)calloc(*nClusters, sizeof(double));

  if(fwhmmap){
    double fwhmsum=0, fwhmmean;
    int nhits=0;
    for (vtx = 0; vtx < Surf->nvertices; vtx++) {
      fwhm = MRIgetVoxVal(fwhmmap,vtx,0,0,0);
      if(fwhm > 0){
	fwhmsum += fwhm;
	nhits ++;
      }
    }
    fwhmmean = fwhmsum/nhits;
    if(Gdiag_no > 0) printf("fwhm mean = %g\n",fwhmmean);
    fwhmmean2 = (fwhmmean*fwhmmean);
  }

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    v = &(Surf->vertices[vtx]);
    clusterno = v->undefval;
    if (clusterno == 0) continue;

    n = clusterno - 1;
    scs[n].nmembers++;
    vtxval = v->val;

    // Initialize
    if (scs[n].nmembers == 1) {
      scs[n].maxval = vtxval;
      scs[n].vtxmaxval = vtx;
      weightvtx[n] = 0.0;
      weightarea[n] = 0.0;
      scs[n].cx = 0.0;
      scs[n].cy = 0.0;
      scs[n].cz = 0.0;
    }

    if(ClusterUseAvgVertexArea == 0){
      if (!Surf->group_avg_vtxarea_loaded)
	vtxarea = v->area;
      else
	vtxarea = v->group_avg_area;
    }
    else vtxarea = avgvertexarea; // effectively measure cluster size as vertex count

    // Convert to resels
    if(fwhmmap){
      fwhm = MRIgetVoxVal(fwhmmap,vtx,0,0,0);
      if(fwhm == 0) fwhm = 1.0; // invalid, not sure what to do 
      // Using fwhmmean2 here just provides a rescaling so that the
      // final cluster areas are reasonable. This might have a mild
      // effect on the distribution of cluster sizes when performing
      // non-stationary perm
      vtxarea *= (fwhmmean2)/(fwhm*fwhm);
    }

    scs[n].area += vtxarea;

    if (fabs(vtxval) > fabs(scs[n].maxval)) {
      scs[n].maxval = vtxval;
      scs[n].vtxmaxval = vtx;
    }
    weightvtx[n] += vtxval;
    weightarea[n] += (vtxval * vtxarea);
    scs[n].cx += v->x;
    scs[n].cy += v->y;
    scs[n].cz += v->z;
  }  // end loop over vertices

  for (n = 0; n < *nClusters; n++) {
    scs[n].clusterno = n + 1;
    scs[n].x = Surf->vertices[scs[n].vtxmaxval].x;
    scs[n].y = Surf->vertices[scs[n].vtxmaxval].y;
    scs[n].z = Surf->vertices[scs[n].vtxmaxval].z;
    scs[n].weightvtx = weightvtx[n];
    scs[n].weightarea = weightarea[n];
    scs[n].cx /= scs[n].nmembers;
    scs[n].cy /= scs[n].nmembers;
    scs[n].cz /= scs[n].nmembers;
    if (T != NULL) {
      xyz->rptr[1][1] = scs[n].x;
      xyz->rptr[2][1] = scs[n].y;
      xyz->rptr[3][1] = scs[n].z;
      MatrixMultiply(T, xyz, xyzxfm);
      scs[n].xxfm = xyzxfm->rptr[1][1];
      scs[n].yxfm = xyzxfm->rptr[2][1];
      scs[n].zxfm = xyzxfm->rptr[3][1];

      xyz->rptr[1][1] = scs[n].cx;
      xyz->rptr[2][1] = scs[n].cy;
      xyz->rptr[3][1] = scs[n].cz;
      MatrixMultiply(T, xyz, xyzxfm);
      scs[n].cxxfm = xyzxfm->rptr[1][1];
      scs[n].cyxfm = xyzxfm->rptr[2][1];
      scs[n].czxfm = xyzxfm->rptr[3][1];
    }
  }

  MatrixFree(&xyz);
  MatrixFree(&xyzxfm);
  free(weightvtx);
  free(weightarea);
  msecTime = mytimer.milliseconds();
  if (Gdiag_no > 0) printf("SurfClusterSumFast: n=%d, t = %g\n", *nClusters, msecTime / 1000.0);

  return (scs);
}

/*----------------------------------------------------------------*/
SCS *SurfClusterSummaryOld(MRI_SURFACE *Surf, MATRIX *T, int *nClusters)
{
  int n;
  SURFCLUSTERSUM *scs;
  MATRIX *xyz, *xyzxfm;
  double centroidxyz[3];
  int msecTime;
  const char *UFSS;

  // Must explicity "setenv USE_FAST_SURF_SMOOTHER 0" to turn off fast
  UFSS = getenv("USE_FAST_SURF_SMOOTHER");
  if (!UFSS) UFSS = "1";
  if (strcmp(UFSS, "0")) {
    scs = SurfClusterSummary(Surf, T, nClusters, NULL);
    return (scs);
  }
  if (Gdiag_no > 0) printf("SurfClusterSummary()\n");

  Timer mytimer;

  *nClusters = sclustCountClusters(Surf);
  if (*nClusters == 0) return (NULL);

  xyz = MatrixAlloc(4, 1, MATRIX_REAL);
  xyz->rptr[4][1] = 1;
  xyzxfm = MatrixAlloc(4, 1, MATRIX_REAL);

  scs = (SCS *)calloc(*nClusters, sizeof(SCS));

  for (n = 0; n < *nClusters; n++) {
    scs[n].clusterno = n + 1;
    scs[n].area = sclustSurfaceArea(n + 1, Surf, &scs[n].nmembers);
    scs[n].weightvtx = sclustWeight(n + 1, Surf, NULL, 0);
    scs[n].weightarea = sclustWeight(n + 1, Surf, NULL, 1);
    scs[n].maxval = sclustSurfaceMax(n + 1, Surf, &scs[n].vtxmaxval);
    scs[n].x = Surf->vertices[scs[n].vtxmaxval].x;
    scs[n].y = Surf->vertices[scs[n].vtxmaxval].y;
    scs[n].z = Surf->vertices[scs[n].vtxmaxval].z;
    sclustSurfaceCentroid(n + 1, Surf, &centroidxyz[0]);
    scs[n].cx = centroidxyz[0];
    scs[n].cy = centroidxyz[1];
    scs[n].cz = centroidxyz[2];
    if (T != NULL) {
      xyz->rptr[1][1] = scs[n].x;
      xyz->rptr[2][1] = scs[n].y;
      xyz->rptr[3][1] = scs[n].z;
      MatrixMultiply(T, xyz, xyzxfm);
      scs[n].xxfm = xyzxfm->rptr[1][1];
      scs[n].yxfm = xyzxfm->rptr[2][1];
      scs[n].zxfm = xyzxfm->rptr[3][1];

      xyz->rptr[1][1] = scs[n].cx;
      xyz->rptr[2][1] = scs[n].cy;
      xyz->rptr[3][1] = scs[n].cz;
      MatrixMultiply(T, xyz, xyzxfm);
      scs[n].cxxfm = xyzxfm->rptr[1][1];
      scs[n].cyxfm = xyzxfm->rptr[2][1];
      scs[n].czxfm = xyzxfm->rptr[3][1];
    }
  }

  MatrixFree(&xyz);
  MatrixFree(&xyzxfm);
  msecTime = mytimer.milliseconds();
  if (Gdiag_no > 0) printf("SurfClusterSum: n=%d, t = %g\n", *nClusters, msecTime / 1000.0);

  return (scs);
}
/*----------------------------------------------------------------*/
int DumpSurfClusterSum(FILE *fp, SCS *scs, int nClusters)
{
  int n;

  for (n = 0; n < nClusters; n++) {
    fprintf(fp,
            "%4d  %4d  %8.4f  %6d    %6.2f  %4d   %5.1f %5.1f %5.1f   "
            "%5.1f %5.1f %5.1f\n",
            n,
            scs[n].clusterno,
            scs[n].maxval,
            scs[n].vtxmaxval,
            scs[n].area,
            scs[n].nmembers,
            scs[n].x,
            scs[n].y,
            scs[n].z,
            scs[n].xxfm,
            scs[n].yxfm,
            scs[n].zxfm);
  }
  return (0);
}
/*----------------------------------------------------------------*/
SCS *SortSurfClusterSum(SCS *scs, int nClusters)
{
  SCS *scs_sorted;
  int n;

  scs_sorted = (SCS *)calloc(nClusters, sizeof(SCS));

  for (n = 0; n < nClusters; n++) memmove(&scs_sorted[n], &scs[n], sizeof(SCS));

  /* Note: scs_sorted.clusterno does not changed */
  qsort((void *)scs_sorted, nClusters, sizeof(SCS), sclustCompare);

  return (scs_sorted);
}

/*----------------------------------------------------------------
  sclustReMap() - remaps the cluster numbers (ie, undefval) in the
  surface structure based on the sorted surface cluster summary
  (SCS). It is assumed that the scs.clusterno in the sorted SCS
  is the cluster id that corresponds to the original cluster id.
  ----------------------------------------------------------------*/
int sclustReMap(MRI_SURFACE *Surf, int nClusters, SCS *scs_sorted)
{
  int vtx, c, cOld;
  int *Orig2Sorted;

  Orig2Sorted = (int *)calloc(nClusters, sizeof(int));

  for (c = 1; c <= nClusters; c++) {
    cOld = scs_sorted[c - 1].clusterno;
    Orig2Sorted[cOld - 1] = c;
    // printf("new = %3d old = %3d\n",c,scs_sorted[c-1].clusterno);
  }

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    c = Surf->vertices[vtx].undefval - 1;
    if (c < 0)
      Surf->vertices[vtx].undefval = 0;
    else
      Surf->vertices[vtx].undefval = Orig2Sorted[c];
  }

  // Change cluster numbers in table
  for (c = 1; c <= nClusters; c++) scs_sorted[c - 1].clusterno = c;

  free(Orig2Sorted);

  return (0);
}
/* Older, slower version of sclustReMap()*/
int sclustReMap0(MRI_SURFACE *Surf, int nClusters, SCS *scs_sorted)
{
  int vtx, vtx_clusterno, c;

  if (Gdiag_no > 1) {
    printf("sclustReMap:\n");
    for (c = 1; c <= nClusters; c++) printf("new = %3d old = %3d\n", c, scs_sorted[c - 1].clusterno);
  }

  for (vtx = 0; vtx < Surf->nvertices; vtx++) {
    vtx_clusterno = Surf->vertices[vtx].undefval;
    for (c = 1; c <= nClusters; c++) {
      if (vtx_clusterno == scs_sorted[c - 1].clusterno) {
        Surf->vertices[vtx].undefval = c;
        break;
      }
    }
  }

  return (0);
}

/*----------------------------------------------------------------*/
/*--------------- STATIC FUNCTIONS BELOW HERE --------------------*/
/*----------------------------------------------------------------*/

/*----------------------------------------------------------------
  sclustCompare() - compares two surface cluster summaries (for
  use with qsort().
  ----------------------------------------------------------------*/
static int sclustCompare(const void *a, const void *b)
{
  SCS sc1, sc2;

  sc1 = *((SURFCLUSTERSUM *)a);
  sc2 = *((SURFCLUSTERSUM *)b);

  if (sc1.pval_clusterwise < sc2.pval_clusterwise) return (-1);
  if (sc1.pval_clusterwise > sc2.pval_clusterwise) return (+1);

  if (sc1.area > sc2.area) return (-1);
  if (sc1.area < sc2.area) return (+1);

  if (fabs(sc1.maxval) > fabs(sc2.maxval)) return (-1);
  if (fabs(sc1.maxval) < fabs(sc2.maxval)) return (+1);

  return (0);
}

/*-------------------------------------------------------------------
  sclustMaxClusterArea() - returns the area of the cluster with the
  maximum area.
  -------------------------------------------------------------------*/
double sclustMaxClusterArea(SURFCLUSTERSUM *scs, int nClusters)
{
  int n;
  double maxarea;

  if (nClusters == 0) return (0);

  maxarea = scs[0].area;
  for (n = 0; n < nClusters; n++)
    if (maxarea < scs[n].area) maxarea = scs[n].area;
  return (maxarea);
}

/*-------------------------------------------------------------------
  sclustMaxClusterCount() - returns the area of the cluster with the
  maximum number of members (count)
  -------------------------------------------------------------------*/
int sclustMaxClusterCount(SURFCLUSTERSUM *scs, int nClusters)
{
  int n;
  int maxcount;

  if (nClusters == 0) return (0);

  maxcount = scs[0].nmembers;
  for (n = 0; n < nClusters; n++)
    if (maxcount < scs[n].nmembers) maxcount = scs[n].nmembers;
  return (maxcount);
}

/*-------------------------------------------------------------------
  sclustMaxClusterWeightVtx() - returns the weightvtx of the cluster
  with the maximum weightvtx.
  -------------------------------------------------------------------*/
float sclustMaxClusterWeightVtx(SURFCLUSTERSUM *scs, int nClusters, int thsign)
{
  int n;
  float maxw, w;

  if (nClusters == 0) return (0);

  if (thsign == 0)
    maxw = 0;
  else
    maxw = -thsign * 10e10;
  for (n = 0; n < nClusters; n++) {
    w = scs[n].weightvtx;
    if (thsign == 0 && fabs(maxw) < fabs(w)) maxw = w;
    if (thsign == +1 && maxw < w) maxw = w;
    if (thsign == -1 && maxw > w) maxw = w;
  }
  return (maxw);
}

/*---------------------------------------------------------------*/
SCS *sclustPruneByCWPval(SCS *ClusterList, int nclusters, double cwpvalthresh, int *nPruned, MRIS *surf)
{
  int n, nth, vtxno, map[10000];
  SCS *scs;

  // Construct a new SCS with pruned clusters
  nth = 0;
  for (n = 0; n < nclusters; n++)
    if (ClusterList[n].pval_clusterwise <= cwpvalthresh) nth++;
  *nPruned = nth;

  scs = (SCS *)calloc(*nPruned, sizeof(SCS));
  nth = 0;
  for (n = 0; n < nclusters; n++) {
    if (ClusterList[n].pval_clusterwise <= cwpvalthresh) {
      memmove(&scs[nth], &ClusterList[n], sizeof(SCS));
      map[n] = nth;
      nth++;
    }
  }

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    n = surf->vertices[vtxno].undefval;  // 1-based
    if (n == 0) continue;
    if (ClusterList[n - 1].pval_clusterwise > cwpvalthresh) {
      // Remove clusters/values from surface
      surf->vertices[vtxno].undefval = 0;
      surf->vertices[vtxno].val = 0;
    }
    else {
      // Re-number
      surf->vertices[vtxno].undefval = map[n - 1] + 1;
    }
  }

  return (scs);
}

/* ------------------------------------------------------------
   int sclustAnnot(MRIS *surf, int NClusters)
   Convert clusters into annotation
   ------------------------------------------------------------*/
int sclustAnnot(MRIS *surf, int NClusters)
{
  COLOR_TABLE *ct;
  int vtxno, vtx_clusterno, annot, n;

  ct = CTABalloc(NClusters + 1);
  surf->ct = ct;

  for (n = 1; n < NClusters + 1; n++)  // no cluster 0
    sprintf(surf->ct->entries[n]->name, "%s-%03d", "cluster", n);

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    vtx_clusterno = surf->vertices[vtxno].undefval;
    if (vtx_clusterno == 0 || vtx_clusterno > NClusters) {
      surf->vertices[vtxno].annotation = 0;
      continue;
    }
    CTABannotationAtIndex(surf->ct, vtx_clusterno, &annot);
    surf->vertices[vtxno].annotation = annot;
  }
  return (0);
}

/* --------------------------------------------------------------*/
/*!
  \fn int sclustGrowByDist(MRIS *surf, int seedvtxno, double dthresh,
        int shape, int vtxno, int init, int *vtxlist)
  \brief Grow a cluster from the seed vertex out to a given distance.
  The undefval of the surf will be zeroed and all vertices in
  within distance are set to 1. vtxlist is filled with the vertex
  numbers, and the return value is the number of vertices in the
  cluster. If shape is SPHERICAL_COORDS, then distance is computed
  along the arc of a sphere, otherwise 2d flat is assumed. This
  is a recursive function. Recursive calls will have the
  the vtxno >= 0, so set vtxno = -1 when you call this function.
*/
int sclustGrowByDist(MRIS *surf, int seedvtxno, double dthresh, int shape, int vtxno, int *vtxlist)
{
  static double radius = 0, radius2 = 0;
  static int nhits = 0, ncalls = 0;
  static VERTEX *v1 = NULL;
  double theta = 0, costheta = 0, d = 0;
  int nthnbr, nbrvtxno;

  if (vtxlist == NULL) {
    printf("ERROR: vtxlist is NULL\n");
    return (-1);
  }

  // This is just a check to make sure that the caller has done the right thing.
  // Note: this check does not work after the first call.
  if (ncalls == 0 && vtxno >= 0) {
    printf("ERROR: set vtxno to a negative value when calling this function!\n");
    return (-1);
  }

  // A negative vtxno indicates a non-recursive call, Init required
  if (vtxno < 0) {
    v1 = &(surf->vertices[seedvtxno]);
    ncalls = 0;
    nhits = 0;
    if (shape == SPHERICAL_COORDS) {
      radius2 = (v1->x * v1->x + v1->y * v1->y + v1->z * v1->z);
      radius = sqrt(radius2);
    }
    for (vtxno = 0; vtxno < surf->nvertices; vtxno++) surf->vertices[vtxno].undefval = 0;
    vtxno = seedvtxno;
  }

  VERTEX_TOPOLOGY const * const v2t = &surf->vertices_topology[vtxno];
  VERTEX    	        * const v2  = &surf->vertices         [vtxno];
  ncalls++;

  if (v2->undefval) return (0);

  if (shape == SPHERICAL_COORDS) {
    costheta = ((v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z)) / radius2;
    theta = acos(costheta);
    d = radius * theta;
    // printf("%g %g %g %g\n",costheta,theta,radius,radius2);
  }
  else {
    d = sqrt((v1->x - v2->x) * (v1->x - v2->x) + (v1->y - v2->y) * (v1->y - v2->y));
  }

  // Check distance against threshold
  if (d > dthresh) return (0);

  // Add to the list
  vtxlist[nhits] = vtxno;
  nhits++;

  // printf("%3d %3d %6d %6d   %g %g %g   %g %g %g   %g\n",
  // ncalls,nhits,seedvtxno,vtxno,
  // v1->x,v1->y,v1->z,v2->x,v2->y,v2->z, d);

  // Go throught the neighbors ...
  v2->undefval = 1;
  for (nthnbr = 0; nthnbr < v2t->vnum; nthnbr++) {
    nbrvtxno = v2t->v[nthnbr];
    sclustGrowByDist(surf, seedvtxno, dthresh, shape, nbrvtxno, vtxlist);
  }

  return (nhits);
}

/*!
  \fn int sclustSaveAsPointSet(char *fname, SCS *scslist, int NClusters, MRIS *surf)
  \brief Outputs a file that can be loaded into freeview with -c
 */
int sclustSaveAsPointSet(char *fname, SCS *scslist, int NClusters, MRIS *surf)
{
  int n,vtxno;
  VERTEX *v;
  FILE *fp;

  fp = fopen(fname,"w");
  for (n=0; n < NClusters; n++) {
    vtxno = scslist[n].vtxmaxval;
    v = &(surf->vertices[vtxno]);
    fprintf(fp,"%g %g %g\n",v->x,v->y,v->z);
  }
  fprintf(fp,"info\n");
  fprintf(fp,"numpoints %d\n",NClusters);
  fprintf(fp,"useRealRAS 0\n");
  fclose(fp);

  return(0);
}
