/**
 * @brief removes surface intersections
 *
 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "surfcluster.h"
#include "mrisurf_metricProperties.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int GetNotches(MRIS *surf, double projdistmm, double k1thresh, LABEL *lb, int nmin, int nmax, double stdthresh, char *psfile, char *ocnfile, MRIS *csurf=NULL);

const char *Progname ;
int FillHoles = 0;

int main(int argc, char *argv[])
{
  char         **av, *in_surf_fname, *out_fname ;
  int          ac, nargs, msec ;
  MRI_SURFACE  *mris ;
  Timer then ;


  std::string cmdline = getAllInfo(argc, argv, "mris_remove_intersection");

  nargs = handleVersionOption(argc, argv, "mris_remove_intersection");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Gdiag = DIAG_SHOW ;

  then.reset() ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit() ;
  }

  in_surf_fname = argv[1] ;
  out_fname = argv[2] ;

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;

  MRISaddCommandLine(mris, cmdline) ;

  mrisMarkIntersections(mris,FillHoles);
  int n, nintersections=0;
  for(n=0; n < mris->nvertices; n++) if(mris->vertices[n].marked) nintersections++;
  printf("Found %d intersections\n",nintersections);

  // MRISsetNeighborhoodSizeAndDist(mris, 2) ;
  MRISremoveIntersections(mris,FillHoles) ;

  printf("writing corrected surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  msec = then.milliseconds() ;
  fprintf(stderr, "intersection removal took %2.2f hours\n",
          (float)msec/(1000.0f*60.0f*60.0f));

  exit(0) ;

  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage")) print_help() ;
  else if (!stricmp(option, "-version")) print_version() ;
  else if (!stricmp(option, "fill-holes"))  FillHoles = 1;
  else if (!stricmp(option, "no-fill-holes"))  FillHoles = 0;
  else if (!stricmp(option, "map"))
  {
    MRIS *surf = MRISread(argv[2]);
    if(surf==NULL) exit(1);
    if(argc > 4) {
      double projdistmm = 0;
      sscanf(argv[4],"%lf",&projdistmm);
      printf("Projecting surface %g mm\n",projdistmm);
      for(int n=0; n < surf->nvertices; n++){
	VERTEX *v = &(surf->vertices[n]);
	v->x += projdistmm*v->nx;
	v->y += projdistmm*v->ny;
	v->z += projdistmm*v->nz;
      }
    }
    mrisMarkIntersections(surf,FillHoles);
    int n, nintersections=0;
    for(n=0; n < surf->nvertices; n++){
      if(surf->vertices[n].marked) nintersections++;
    }
    printf("Found %d intersections\n",nintersections);
    #if 0
    printf("Filling any holes\n");
    SURFCLUSTERSUM *SurfClustList;
    int nClusters;
    SurfClustList = sclustMapSurfClusters(surf,0.5,-1,1,0,&nClusters,NULL,NULL);
    printf("Found %d clusters\n",nClusters);
    for(n=0; n < surf->nvertices; n++){
      if(surf->vertices[n].undefval > 1) surf->vertices[n].marked = 1;
    }
    #endif
    MRI *mri = MRIcopyMRIS(NULL,surf,0,"marked");
    int err = MRIwrite(mri,argv[3]);
    exit(err);
  }
  else if (!stricmp(option, "notches"))
  {
    // 2=insurf 3=projdist(1) 4=k1thresh(.05) 5=label(cortex) 6=nmin(3) 7=nmax 8=stdthresh(.2) 9=pointset 10=ocn <11=pial>
    MRIS *surf = MRISread(argv[2]);
    if(surf==NULL) exit(1);
    double projdistmm, k1thresh, stdthresh;
    int nmin, nmax = -1;
    sscanf(argv[3],"%lf",&projdistmm);
    sscanf(argv[4],"%lf",&k1thresh);
    LABEL *lb = LabelRead(NULL,argv[5]);
    if(!lb) exit(1);
    sscanf(argv[6],"%d",&nmin);
    sscanf(argv[7],"%d",&nmax);
    sscanf(argv[8],"%lf",&stdthresh);
    MRIS *csurf=NULL;
    if(argc > 11){
      csurf = MRISread(argv[11]);
      if(!csurf) exit(1);
    }
    int err = GetNotches(surf, projdistmm, k1thresh, lb, nmin, nmax, stdthresh, argv[9], argv[10], csurf);
    if(err) exit(err);
    printf("\n");
    printf("fsvglrun freeview --hide-3d-slices -neuro-view -rotate-around-cursor ");
    printf("-f %s:annot=%s ",argv[2],argv[10]);
    if(argc > 11) printf("-f %s:annot=%s ",argv[11],argv[10]);
    printf("-c %s\n",argv[9]);
    printf("\n");
    exit(0);
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

#include "mris_remove_intersection.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_remove_intersection_help_xml,mris_remove_intersection_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

/*!
  \fn GetNotches()
  \brief This function finds "notches" in the white surface. A notch
   is a break in the continuous ridge of a gyrus. The notch itself is
   usually a slight error in the white surface placement. The bigger
   error happens when the pial surface is placed as the pial starts as
   the white surface. In these areas, the pial surface often cannot
   escape from the white creating inaccuracies in the pial surface
   (sometimes called "pinches"). This routine tries to localize the
   notches on the white surface. The first step is to project the
   surface along the normal by projdistmm (eg, 1mm) and then find the
   places where intersections are created. The idea is that a normal
   white surface will not have a bunch of intersections but there will
   be intersections around a notch.  Notches form along a gyrus, so
   intersections are limited to be in areas where the K1 curvature is
   less than -k1thresh (eg, .05). A label (eg, cortex.label) is used
   to exclude intersections from the medial wall. Clusters are formed
   from the intersections. Clusters are excluded if the size is less
   than nmin (eg, 3). Even with all this, there are false positives
   (eg, in the elbow where a gyrus turns), so the last step is to
   compute the stddev of the dot product of the normals with the mean
   normal in the cluster. Most of the white surface will have a pretty
   consistent normal in an area. In an elbow, there will be more
   variation, but there will be a lot in a notch, so the final step is
   to remove clusters whose dotstddev is less than a stdthresh (eg,
   0.2). The centroids of the surviving clusters are saved in a point
   set (psfile, can be NULL) and surface seg (ocnfile, can be NULL).
   If csurf is non-NULL, then the cenroids will be computed with
   coordinates from csurf (eg, to have the point set on the pial
   instead of the white).
 */
int GetNotches(MRIS *surf, double projdistmm, double k1thresh, LABEL *lb, int nmin, int nmax, double stdthresh, char *psfile, char *ocnfile, MRIS *csurf)
{
  int FillHoles = 1;

  printf("GetNotches(): %g %g %d %d %g\n",projdistmm,k1thresh,nmin,nmax,stdthresh);

  // This is needed to compute K1
  MRISsetNeighborhoodSizeAndDist(surf,2); // neighborhood size 2 
  MRIScomputeSecondFundamentalFormDiscrete(surf, 0);

  // Project surfaces so that they intersect
  MRISsaveVertexPositions(surf, TMP_VERTICES);
  printf("Projecting surface %g mm\n",projdistmm);
  for(int n=0; n < surf->nvertices; n++){
    VERTEX *v = &(surf->vertices[n]);
    v->x += projdistmm*v->nx;
    v->y += projdistmm*v->ny;
    v->z += projdistmm*v->nz;
  }

  // Mark the intersections
  mrisMarkIntersections(surf,FillHoles);
  int nintersections=0;
  for(int n=0; n < surf->nvertices; n++) if(surf->vertices[n].marked) nintersections++;
  printf("Found %d vertices with intersections\n",nintersections);
  MRISrestoreVertexPositions(surf, TMP_VERTICES);

  // If label is passed, make a binary mask from it
  MRI *labelmask = NULL;
  if(lb) labelmask = MRISlabel2Mask(surf, lb, NULL);

  // Remove marked vertices that are out of the mask or don't meet the K1 threshold
  nintersections=0;
  for(int vtxno = 0; vtxno < surf->nvertices; vtxno++){
    VERTEX *v = &(surf->vertices[vtxno]);
    v->val = 0;
    if(! v->marked) continue;
    if((labelmask && MRIgetVoxVal(labelmask,vtxno,0,0,0) < 0.5) || -v->k1 < k1thresh){
      v->marked = 0;
      continue;
    }
    v->val = 1;
    nintersections++;
  }
  printf("Found %d vertices with intersections\n",nintersections);

  //MRI *mri = MRIcopyMRIS(NULL,surf,0,"marked");
  //MRIwrite(mri,"dng.x.mgz");

  // Clusterize the remaining vertices
  SURFCLUSTERSUM *SurfClustList;
  int nClusters;
  SurfClustList = sclustMapSurfClusters(surf,0.5,-1,1,0,&nClusters,NULL,NULL);
  printf("Found %d clusters\n",nClusters);

  // Compute the mean normal and the centroid within each cluster
  double *nx = (double*)calloc(nClusters,sizeof(double));
  double *ny = (double*)calloc(nClusters,sizeof(double));
  double *nz = (double*)calloc(nClusters,sizeof(double));
  double *cx = (double*)calloc(nClusters,sizeof(double));
  double *cy = (double*)calloc(nClusters,sizeof(double));
  double *cz = (double*)calloc(nClusters,sizeof(double));
  for(int vtxno = 0; vtxno < surf->nvertices; vtxno++){
    VERTEX *v = &(surf->vertices[vtxno]);
    int cno = v->undefval-1;
    if(cno == -1) continue;
    nx[cno] += v->nx;
    ny[cno] += v->ny;
    nz[cno] += v->nz;
    if(csurf) v = &(csurf->vertices[vtxno]);
    cx[cno] += v->x;
    cy[cno] += v->y;
    cz[cno] += v->z;
  }
  for(int cno = 0; cno < nClusters; cno++){
    double mag = sqrt(nx[cno]*nx[cno] + ny[cno]*ny[cno] + nz[cno]*nz[cno]);
    nx[cno] /= mag;
    ny[cno] /= mag;
    nz[cno] /= mag;
    cx[cno] /= SurfClustList[cno].nmembers;
    cy[cno] /= SurfClustList[cno].nmembers;
    cz[cno] /= SurfClustList[cno].nmembers;
    // Save the centroid in the SurfCluster struct
    SurfClustList[cno].cx = cx[cno];
    SurfClustList[cno].cy = cy[cno];
    SurfClustList[cno].cz = cz[cno];
  }
  // Compute the mean and std of the dot product within the cluster 
  double *dotsum  = (double*) calloc(nClusters,sizeof(double));
  double *dotsum2 = (double*) calloc(nClusters,sizeof(double));
  for(int vtxno = 0; vtxno < surf->nvertices; vtxno++){
    VERTEX *v = &(surf->vertices[vtxno]);
    int cno = v->undefval-1;
    if(cno == -1) continue;
    double dot = nx[cno]*v->nx + ny[cno]*v->ny + nz[cno]*v->nz;
    dotsum[cno]  += dot;
    dotsum2[cno] += (dot*dot);
  }
  double *dotavg = (double*) calloc(nClusters,sizeof(double));
  double *dotstd = (double*) calloc(nClusters,sizeof(double));
  for(int cno = 0; cno < nClusters; cno++){
    dotavg[cno] = dotsum[cno]/SurfClustList[cno].nmembers; 
    if(SurfClustList[cno].nmembers > 2) dotstd[cno] = sum2stddev(dotsum[cno], dotsum2[cno],SurfClustList[cno].nmembers); 
    SurfClustList[cno].maxval = dotstd[cno];
    SurfClustList[cno].pval_clusterwise = 1.0/(dotstd[cno] + FLT_EPSILON); // sort by this
    //printf("%3d %4d  %6.4lf %6.4lf\n",cno,SurfClustList[cno].nmembers,dotavg[cno],dotstd[cno]);
  }

  // Now sort the clusters based on dot stddev
  SCS *SurfClustSorted = SortSurfClusterSum(SurfClustList, nClusters);

  // Create a point set of those clusters with enough vertices and whose dot stddev meets threshold
  int kClusters = 0; // number of clusters after excluding small and < stdthresh clusters
  std::vector<int> kcnolist;
  fsPointSet  ps;
  ps.vox2ras = "tkreg";
  for(int cno=0; cno < nClusters; cno++){
    SCS *scs = &(SurfClustSorted[cno]);
    if(scs->nmembers < nmin) continue;
    if(nmax > 0 && scs->nmembers > nmax) continue;
    if(scs->maxval < stdthresh) continue;
    fsPointSet::Point p = fsPointSet::Point();
    printf("%3d %3d %4d  %6.4lf %6.4lf  %6.2lf %6.2lf %6.2lf\n",
	   kClusters+1,scs->clusterno,scs->nmembers,scs->maxval,scs->pval_clusterwise,scs->cx,scs->cy,scs->cz);
    p.index = kClusters;
    p.count = scs->nmembers;
    p.value = scs->maxval;
    p.x = scs->cx;
    p.y = scs->cy;
    p.z = scs->cz;
    ps.add(p);
    kcnolist.push_back(scs->clusterno);
    kClusters++;
  }
  // Write the point set file
  if(psfile) ps.save(psfile);

  if(ocnfile){
    // Write out the output cluster number file. Each vertex/voxel has 
    // the associated cluster number
    MRI *mritmp = MRIalloc(surf->nvertices,1,1,MRI_INT);
    for(int vtxno = 0; vtxno < surf->nvertices; vtxno++){
      int cno0 = surf->vertices[vtxno].undefval;
      if(cno0 <= 0) {
	MRIsetVoxVal(mritmp,vtxno,0,0,0,0);
	continue;
      }
      int k, ok = 0; // sorted
      for(k=0; k < kClusters; k++) {
	if(kcnolist[k] == cno0){
	  ok = 1;
	  //printf("%d %d %d\n",vtxno,cno0,k);
	  break;
	}
      }
      if(ok) MRIsetVoxVal(mritmp,vtxno,0,0,0,k+1);
    }
    mritmp->ct = CTABalloc(kClusters+1);
    CTABunique(mritmp->ct, kClusters+1);
    MRIwrite(mritmp,ocnfile); // must be loaded as an annot
    MRIfree(&mritmp);
  }

  // At this point, it would probably be good to reset the undefval to
  // be the proper cluster number so that it would be available to the
  // caller.

  free(nx);
  free(ny);
  free(nz);
  free(cx);
  free(cy);
  free(cz);
  free(dotsum);
  free(dotsum2);
  free(dotavg);
  free(dotstd);
  // Free SurfClusts (or return one)

  return(0);
}
