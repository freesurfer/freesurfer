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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>
#include "error.h"
#include "diag.h"
#include "surfcluster.h"
#include "volcluster.h"
#include "romp_support.h"
#include "tfce.h"
#undef private


/*!
  \func MRI *TFCE::voxcor(MRI *statmap, std::vector<double> maxstatlist)
  \brief Computes a corrected voxel-wise p-value map given the statmap and
  a list of maximum voxelwise stats. The list is sorted. This is not 
  specific to TFCE.
*/
MRI *TFCE::voxcor(MRI *statmap, std::vector<double> maxstatlist)
{
  std::sort(maxstatlist.begin(), maxstatlist.end());
  MRI *statmapcor = MRIallocSequence(statmap->width,statmap->height,statmap->depth,MRI_FLOAT,statmap->nframes);
  MRIcopyHeader(statmap, statmapcor);
  MRIcopyPulseParameters(statmap, statmapcor);
  for(int c = 0; c < statmap->width; c++){
    for(int r = 0; r < statmap->height; r++){
      for(int s = 0; s < statmap->depth; s++){
	if(mask){
	  int m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m<0.5) continue;
	}
	double v = MRIgetVoxVal(statmap,c,r,s,0);
	int n;
	for(n=0; n < maxstatlist.size(); n++) if(v < maxstatlist[n]) break;
	double sig,p=-1;
	if(n!=maxstatlist.size()){
	  p = ((double)maxstatlist.size()-n)/maxstatlist.size();
	  sig = -log10(p);
	}
	else sig = 100;
	MRIsetVoxVal(statmapcor,c,r,s,0, sig);
	if(c == 10366) printf("vno=%d v=%g n=%d p=%g sig=%g %d\n",c,v,n,p,sig,(int)maxstatlist.size());
      }
    }
  }
  return(statmapcor);
}

/*!
  \func int TFCE::write_vector_double(char *fname, std::vector<double> vlist)
  \brief Simply writes a double vector into a text file.
  This is not specific to TFCE.
*/
int TFCE::write_vector_double(char *fname, std::vector<double> vlist)
{
  FILE *fp;
  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: cannot open %s for writing\n",fname);
    return(1);
  }
  for(int n=0; n<vlist.size(); n++) fprintf(fp,"%20.10lf\n",vlist[n]);
  fclose(fp);
  return(0);
}
/*!
  \func std::vector<double> TFCE::read_vector_double(char *fname)
  \brief Simply reads a double vector from a text file.
  This is not specific to TFCE.
*/
std::vector<double> TFCE::read_vector_double(char *fname)
{
  std::vector<double> vlist;
  FILE *fp;
  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s for reading\n",fname);
    return(vlist);
  }
  while(!feof(fp)){
    double v;
    fscanf(fp,"%lf",&v);
    vlist.push_back(v);
  }
  fclose(fp);
  return(vlist);
}

/*!
  \func double TFCE::getmax(MRI *statmap, MRI *mask)
  \brief Simply returns the maximum in the statmap within 
  the mask. This is not specific to TFCE.
*/
double TFCE::getmax(MRI *statmap, MRI *mask)
{
  double maxstat=-1;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(max:maxstat)
  #endif
  for(int c = 0; c < statmap->width; c++){
    for(int r = 0; r < statmap->height; r++){
      for(int s = 0; s < statmap->depth; s++){
        if(mask){
          int m = MRIgetVoxVal(mask,c,r,s,0);
          if(m<0.5) continue;
        }
        double v = MRIgetVoxVal(statmap,c,r,s,0);
        if(maxstat < v) maxstat = v;
      }//s
    }//r
  }//c
  return(maxstat);
}

/*!
  \func std::vector<double> TFCE::maxstatsim(MRI *temp, int niters)
  \brief Simulates TFCE assuming unsmoothed white Gaussian noise
  (ie, a z map). Returns the sorted list of maximum voxelwise stats.
*/
std::vector<double> TFCE::maxstatsim(MRI *temp, int niters)
{
  std::vector<double> maxstatlist;
  //hmin = FLT_EPSILON;  hmax = 4;  nh = 50;  hlistUniform();
  if(hlist.size()==0){
    printf("ERROR: TFCE::maxstatsim(): hlist has not been set up\n");
    return(maxstatlist);
  }
  //cant run in parallel because sclust uses surface val and undef;
  //need to copy surface for each threads/iter
  for(int n=0; n < niters; n++){
    MRI *zmap = MRIrandn(temp->width, temp->height, temp->depth, 1, 0.0, 1.0, NULL);
    MRIcopyHeader(temp, zmap); // needs to have proper voxel size for volume topo
    //MRIcopyPulseParameters(temp, zmap);
    //MRIwrite(zmap,"zmap.mgh");
    MRI *tfcemap = compute(zmap);
    double maxstat = getmax(tfcemap,mask);
    maxstatlist.push_back(maxstat);
    printf("iter=%d  maxstat %g\n",n,maxstat); fflush(stdout);
    //MRIwrite(tfcemap,"zmap.tfce.mgh");
    MRIfree(&zmap);
    MRIfree(&tfcemap);
  }// iter
  std::sort(maxstatlist.begin(), maxstatlist.end());
  return(maxstatlist);
}

/*!
  \func int TFCE::hlistUniform(void)
  \brief Create a nh-length vector of thresholds uniformly 
  distributed between hmin and hmax.
*/
int TFCE::hlistUniform(void)
{
  hlist.clear();
  double dh = (hmax-hmin)/(nh-1);
  if(debug) printf("hlistUniform(): nh=%d, hmin=%g, hmax=%g, dh=%g\n",nh,hmin,hmax,dh);
  for(int n=0; n<nh; n++){
    double h = hmin + n*dh;
    hlist.push_back(h);
  }
  return(0);
}

/*!
  \func int TFCE::hlistAuto(MRI *map)
  \brief Creates a threshold list where the number of voxels in each bin
  is approximately equal.
*/
int TFCE::hlistAuto(MRI *map)
{
  hlist.clear();
  int nhits = 0;
  std::vector<double> vlist;
  for(int c = 0; c < map->width; c++){
    for(int r = 0; r < map->height; r++){
      for(int s = 0; s < map->depth; s++){
	if(mask){
	  int m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m<0.5) continue;
	}
	nhits++;
	double v = MRIgetVoxVal(map,c,r,s,0);
	if(surf) surf->vertices[c].val = v; // don't apply a sign here
	if(thsign ==  0) v = fabs(v);
	if(thsign == -1) v = -v;
        if(v>FLT_EPSILON) vlist.push_back(v);
      }
    }
  }
  int nvlist = (int)vlist.size();
  std::sort(vlist.begin(), vlist.end());
  hmin = vlist[0]-FLT_EPSILON;
  if(hmin<=0) hmin = FLT_EPSILON;
  hmax = vlist[nvlist-1]-FLT_EPSILON;
  hlist.push_back(hmin);
  int nskip = round((double)nvlist/nh);
  printf("hlistAuto(): nh=%d, nvlist=%d, nskip=%d\n",nh,nvlist,nskip);
  if(nskip == 0){
    printf("ERROR: hlistAuto(): nskip=0\n");
    return(1);
  }
  for(int n=1; n < nh-1; n++) hlist.push_back(vlist[n*nskip]);
  hlist.push_back(hmax);
  printf("hlistAuto(): nhits=%d, thsign=%d, nh=%d hmin=%g, hmax=%g\n",thsign,nhits,nh,hmin,hmax);

  return(0);
}

/*!
  \func MRI *TFCE::compute(MRI *map)
  \brief Computes the TFCE map
*/
MRI *TFCE::compute(MRI *map)
{
  MRI *tfcemap = NULL;

  if(debug) printf("Entering TFCE::compute() hlist.size()=%d\n",(int)hlist.size());
  if(hlist.size()==0){
    printf("ERROR: TFCE::compute(): hlist has not been set up\n");
    return(NULL);
  }
  
  int nhits = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+:nhits)
  #endif
  for(int c = 0; c < map->width; c++){
    for(int r = 0; r < map->height; r++){
      for(int s = 0; s < map->depth; s++){
	if(mask){
	  int m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m<0.5) {
	    if(surf) surf->vertices[c].val = 0; // make sure=0
	    continue;
	  }
	}
	nhits++;
	double v = MRIgetVoxVal(map,c,r,s,0);
	if(surf) surf->vertices[c].val = v; // don't apply a sign here
      }
    }
  }
  if(debug) printf("E=%g, H=%g, nhits=%d, thsign=%d, nh=%d\n",E,H,thsign,nhits,nh);

  MRI *tfcemaps = MRIallocSequence(map->width,map->height,map->depth,MRI_FLOAT,nh);
  MRIcopyHeader(map, tfcemaps);
  MRIcopyPulseParameters(map, tfcemaps);

  if(debug) printf("Entering TFCE::compute() hlist.size()=%d\n",(int)hlist.size());
  // Cannot be parallized here because MapSurfClust is not thread safe
  for(int nthh=0; nthh < nh; nthh++){
    double h = hlist[nthh];
    double powhH = pow(h,H);// default: E=0.5, H=2
    int nClusters;
    if(surf){
      SCS *scs=NULL;
      scs = sclustMapSurfClusters(surf, h, -1, thsign, 0, &nClusters, NULL, NULL);
      if(nClusters==0) continue;
      int nhits = 0;
      for(int c=0; c<nClusters; c++) nhits += scs[0].nmembers;
      if(debug) printf("%2d h=%g, nc=%d nhits=%d  c0nm=%d  c0area=%6.1f\n",nthh,h,nClusters,nhits,scs[0].nmembers,scs[0].area);
      #ifdef HAVE_OPENMP
      #pragma omp parallel for 
      #endif
      for(int vno=0; vno < surf->nvertices; vno++){
        VERTEX *vtx = &surf->vertices[vno];
        int cno = vtx->undefval;
        if(cno == 0) continue;
        double v = pow(scs[cno-1].area,E)*powhH;
        if(debug && vno == vnodebug) printf("   h=%g vno=%d cno=%d nm=%d a=%6.3f v=%g %6.3f %6.3f\n",
  					  h,vno,cno,scs[cno-1].nmembers,scs[cno-1].area,v,pow(scs[cno-1].area,E),pow(h,H));
        MRIsetVoxVal(tfcemaps,vno,0,0,nthh,v);
      }
      free(scs);
    }
    else {
      VOLCLUSTER **VCList = clustGetClusters(map, 0, h,-1,thsign,0,mask, &nClusters, NULL);
      if(nClusters==0) continue;
      for(int cno=0; cno<nClusters; cno++){
        VOLCLUSTER *VC = VCList[cno];
        double csize = VC->nmembers * VC->voxsize;
        double v = pow(csize,E)*powhH;
	for(int n=0; n < VC->nmembers; n++)
          MRIsetVoxVal(tfcemaps,VC->col[n],VC->row[n],VC->slc[n],nthh,v);
      }
      clustFreeClusterList(&VCList,nClusters);
    }
  }

  tfcemap = MRIallocSequence(map->width,map->height,map->depth,MRI_FLOAT,1);
  MRIcopyHeader(map, tfcemap);
  MRIcopyPulseParameters(map, tfcemap);
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(int c = 0; c < map->width; c++){
    for(int r = 0; r < map->height; r++){
      for(int s = 0; s < map->depth; s++){
        if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
        double vsum=0;
        for(int nthh=0; nthh < nh-1; nthh++){
          double v1 = MRIgetVoxVal(tfcemaps,c,r,s,nthh);
          double v2 = MRIgetVoxVal(tfcemaps,c,r,s,nthh+1);
          double hd = hlist[nthh+1]-hlist[nthh];
          vsum += hd*(v1 + (v2-v1)/2);
        }
        MRIsetVoxVal(tfcemap,c,r,s,0, vsum);
        if(debug && c == vnodebug) printf("  nh=%d final tfce stat c=%d %g\n",nh,c,vsum);
      }
    }
  }
  MRIfree(&tfcemaps);

  return(tfcemap);
}
