/**
 * @file  dngtester.c
 * @brief dougs super special test code
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/04/26 16:15:55 $
 *    $Revision: 1.59 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>
#include "mrisurf.h"
#include "geodesics.h"
#include "timer.h"
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
  int vtxno; // surface vertex no
  int volindex; // voxel index from volume
  float dist;
}  VTXVOLINDEX;

MRI *GeoSmooth(MRI *src, double fwhm, MRIS *surf, Geodesics *geod, MRI *volindex, MRI *out);
int GeoCount(Geodesics *geod, int nvertices);
int GeoDumpVertex(char *fname, Geodesics *geod, int vtxno);
int geodesicsWriteV2(Geodesics* geo, int nvertices, char* fname) ;
Geodesics* geodesicsReadV2(char* fname, int *pnvertices) ;
double geodesicsCheckSphereDist(MRIS *sphere, Geodesics *geod);
double MRISsphereDist(MRIS *sphere, VERTEX *vtx1, VERTEX *vtx2);

int VtxVolIndexCompare(const void *a, const void *b);
int VtxVolIndexSort(VTXVOLINDEX *vvi, int nlist);
int VtxVolIndexPrint(VTXVOLINDEX *vvi, int nlist);
int VtxVolIndexSortTest(int nlist);
VTXVOLINDEX *VtxVolIndexUnique(VTXVOLINDEX *vvi, int nlist, int *nunique);
VTXVOLINDEX *VtxVolIndexPack(Geodesics *geod, int vtxno, MRI *volindex);
Geodesics *VtxVolPruneGeod(Geodesics *geod, int vtxno, MRI *volindex);

/*----------------------------------------*/
int main(int argc, char **argv) 
{
  MRIS *surf;
  int msec, nvertices; //vtxno=0;
  struct timeb  mytimer;
  MRI *mri, *mriindex, *mri2;
  Geodesics *geod;
  float maxdist;
  double d;

  // good vertex on the lateral side 108489
  omp_set_num_threads(10);

  // apply smoothing: 5 args: surf geod input index output
  surf = MRISread(argv[1]);
  printf("reading geo\n"); fflush(stdout);
  TimerStart(&mytimer) ;
  geod = geodesicsRead(argv[2], &nvertices);
  msec = TimerStop(&mytimer) ;
  printf("t = %g min\n",msec/(1000.0*60));
  mri  = MRIread(argv[3]);
  mriindex  = MRIread(argv[4]);
  printf("Smoothing\n");
  TimerStart(&mytimer) ;
  mri2 = GeoSmooth(mri, 10, surf, geod, mriindex, NULL);
  msec = TimerStop(&mytimer) ;
  printf("t = %g min\n",msec/(1000.0*60));
  fflush(stdout);
  MRIwrite(mri2,argv[5]);
  mri2 = MRIcopyMRIS(NULL,surf,0,"val2bak");
  MRIwrite(mri2,"nnbrs.mgh");

  exit(0);  //-------------------------

  // check distances: surf geod
  surf = MRISread(argv[1]);
  d = MRISsphereDist(surf, &surf->vertices[220], &surf->vertices[3140]);
  geod = geodesicsRead(argv[2], &nvertices);
  geodesicsCheckSphereDist(surf, geod);
  exit(0);  //-------------------------

  // create geod file: 3 args: surf distmax output
  surf = MRISread(argv[1]);
  sscanf(argv[2],"%f",&maxdist);
  TimerStart(&mytimer) ;
  geod = computeGeodesics(surf, maxdist);
  msec = TimerStop(&mytimer) ;
  printf("done t = %g min\n",msec/(1000.0*60));
  geodesicsWrite(geod, surf->nvertices, argv[3]);
  msec = TimerStop(&mytimer) ;
  printf("done write t = %g min\n",msec/(1000.0*60));
  exit(0); //-------------------------


  //GeoDumpVertex("uvtx.108489.dist.dat", geod, 108489);
  TimerStart(&mytimer) ;
  geod = geodesicsReadV2(argv[1], &nvertices);
  msec = TimerStop(&mytimer) ;
  printf(" read2 t = %g min\n",msec/(1000.0*60));
  exit(0); //-------------------------

  TimerStart(&mytimer) ;
  geodesicsWrite(geod, nvertices, "u1.geod");
  msec = TimerStop(&mytimer) ;
  printf(" write1 t = %g min\n",msec/(1000.0*60));

  TimerStart(&mytimer) ;
  geodesicsWriteV2(geod, nvertices, "u2.geod");
  msec = TimerStop(&mytimer) ;
  printf(" write2 t = %g min\n",msec/(1000.0*60));

  exit(0); //-------------------------



}

MRI *GeoSmooth(MRI *src, double fwhm, MRIS *surf, Geodesics *geod, MRI *volindex, MRI *out)
{
  int vtxno;
  double gvar, gstd, gf;

  if(out == NULL)  {
    out = MRIallocSequence(src->width, src->height, src->depth,
                            MRI_FLOAT, src->nframes);
    if (out==NULL){
      printf("ERROR: GeoSmooth: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(src,out);
    MRIcopyPulseParameters(src,out);
  }
  MRIScomputeMetricProperties(surf);

  gstd = fwhm/sqrt(log(256.0));
  gvar = gstd*gstd;
  // scale factor does not really matter but good check because kernel should sum to 1
  gf = sqrt(pow(2*M_PI,2)*gvar*gvar); // = 2*M_PI*gvar
  printf("fwhm %g gstd %g %g %g, nv=%d\n",fwhm,gstd,gvar,gf,surf->nvertices);

#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    int nthnbr,nbrvtxno,frame;
    double ksum, *sum, d, vkern;
    Geodesics *vtxgeod;
    if(vtxno % 10000 == 0) printf("vtxno %d\n",vtxno);
    if(surf->vertices[vtxno].ripflag) continue;

    // Remove replicate neighbors that have the same volume vertex no
    if(volindex) vtxgeod = VtxVolPruneGeod(geod, vtxno, volindex);
    else         vtxgeod = &geod[vtxno];

    // Set up init using self
    //vkern = surf->vertices[vtxno].area/gf; // scale by the area
    vkern = 1/gf;
    ksum = vkern;
    sum = (double *) calloc(sizeof(double),src->nframes);
    for(frame = 0; frame < src->nframes; frame++)
      sum[frame] = (vkern*MRIgetVoxVal(src,vtxno,0,0,frame));

    for(nthnbr = 0 ; nthnbr < vtxgeod->vnum; nthnbr++) {
      nbrvtxno = vtxgeod->v[nthnbr];
      if(surf->vertices[nbrvtxno].ripflag) continue;
      d = vtxgeod->dist[nthnbr];
      //vkern = surf->vertices[nbrvtxno].area*exp(-(d*d)/(2*gvar))/gf; // scale by the area
      vkern = exp(-(d*d)/(2*gvar))/gf; 
      ksum += vkern;
      for(frame = 0; frame < src->nframes; frame++)
	sum[frame] += (vkern*MRIgetVoxVal(src,nbrvtxno,0,0,frame));
    }

    for(frame = 0; frame < src->nframes; frame++)
      MRIsetVoxVal(out,vtxno,0,0,frame,(sum[frame]/ksum));
    surf->vertices[vtxno].valbak = ksum;
    surf->vertices[vtxno].val2bak = vtxgeod->vnum;
    free(sum);
    if(volindex) free(vtxgeod);
  } // vtxno

  return(out);
}

int GeoCount(Geodesics *geod, int nvertices)
{
  int c=0,cmax=0,vtxno;

  for(vtxno = 0; vtxno < nvertices; vtxno++){
    c += geod[vtxno].vnum;
    if(cmax < geod[vtxno].vnum) cmax = geod[vtxno].vnum;
  }
  printf(" GeoCount %d %d\n",c,cmax);
  return(c);
}

int GeoDumpVertex(char *fname, Geodesics *geod, int vtxno)
{
  FILE *fp;
  int nthnbr;
  fp = fopen(fname,"w");
  for(nthnbr = 0; nthnbr < geod[vtxno].vnum; nthnbr++){
    fprintf(fp,"%5d %g\n",geod[vtxno].v[nthnbr],geod[vtxno].dist[nthnbr]);
  }
  fclose(fp);
  return(0);
}


int geodesicsWriteV2(Geodesics* geo, int nvertices, char* fname) 
{
  int vtxno,*vnum,*vlist,nth,nthnbr,nnbrstot;
  float *dist;
  FILE *fp;

  nnbrstot = 0;
  for(vtxno = 0; vtxno < nvertices; vtxno++) nnbrstot += geo[vtxno].vnum;
  printf(" GeoCount %d\n",nnbrstot);

  fp = fopen(fname, "wb");
  fprintf(fp,"FreeSurferGeodesics-V2\n");
  fprintf(fp,"%d\n",-1);
  fprintf(fp,"%d\n",nvertices);
  fprintf(fp,"%d\n",nnbrstot);

  // Pack the number of neighbors into an array, then fwrite
  vnum = (int *) calloc(sizeof(int),nvertices);
  for(vtxno = 0; vtxno < nvertices; vtxno++) vnum[vtxno] = geo[vtxno].vnum;
  fwrite(vnum,sizeof(int), nvertices, fp);
  free(vnum);

  // Pack the neighbor vertex numbers and dist into an arrays, then fwrite
  vlist = (int *)   calloc(sizeof(int),  nnbrstot);
  dist  = (float *) calloc(sizeof(float),nnbrstot);
  nth = 0;
  for(vtxno = 0; vtxno < nvertices; vtxno++){
    for(nthnbr = 0 ; nthnbr < geo[vtxno].vnum; nthnbr++){
      vlist[nth] = geo[vtxno].v[nthnbr];
      dist[nth]  = geo[vtxno].dist[nthnbr];
      nth ++;
    }
  }
  fwrite(vlist,sizeof(int), nnbrstot, fp);
  free(vlist);
  fwrite(dist,sizeof(float), nnbrstot, fp);
  free(dist);

  fclose(fp);
  return(0);
}

Geodesics* geodesicsReadV2(char* fname, int *pnvertices) 
{
  int magic;
  char tmpstr[1000];
  FILE *fp;
  int vtxno,*vnum,*vlist,nth,nthnbr,nnbrstot;
  float *dist;
  struct timeb  mytimer;
  int msec; 

  fp = fopen(fname, "rb");
  fscanf(fp,"%s",tmpstr);
  if(strcmp(tmpstr,"FreeSurferGeodesics-V2")){
    fclose(fp);
    printf("   %s\n",tmpstr);
    printf("ERROR: %s not a geodesics file\n",fname);
    return(NULL);
  }
  fscanf(fp,"%d",&magic);
  if(magic != -1){
    fclose(fp);
    printf("ERROR: %s wrong endian\n",fname);
    return(NULL);
  }
  fscanf(fp,"%d",pnvertices);
  fscanf(fp,"%d",&nnbrstot);
  fgetc(fp); // swallow the new line
  printf("    geodesicsReadV2(): %s nvertices = %d, magic = %d\n",fname,*pnvertices,magic);fflush(stdout);

  printf("  alloc vnum \n");fflush(stdout);
  vnum = (int *) calloc(sizeof(int),*pnvertices);
  printf("  reading in vnum %d \n",*pnvertices);fflush(stdout);
  TimerStart(&mytimer) ;
  fread(vnum,sizeof(int), *pnvertices, fp);
  msec = TimerStop(&mytimer) ;
  printf("  t = %g min\n",msec/(1000.0*60));
  //printf(" setting vnum %d\n",*pnvertices);
  //for(vtxno = 0; vtxno < *pnvertices; vtxno++) 
  //  geo[vtxno].vnum = vnum[vtxno];
  //free(vnum);
  //msec = TimerStop(&mytimer) ;
  //printf("  t = %g min\n",msec/(1000.0*60));

  TimerStart(&mytimer) ;
  printf("  allocing vlist and dlist %d\n",nnbrstot);fflush(stdout);
  vlist = (int *)   calloc(sizeof(int),  nnbrstot);
  dist  = (float *) calloc(sizeof(float),nnbrstot);
  printf("  reading in vlist %d\n",nnbrstot);fflush(stdout);
  fread(vlist,sizeof(int), nnbrstot, fp);
  msec = TimerStop(&mytimer); printf("  t = %g min\n",msec/(1000.0*60));
  printf("  reading in dist %d\n",nnbrstot);fflush(stdout);
  fread(dist, sizeof(float), nnbrstot, fp);
  msec = TimerStop(&mytimer); printf("  t = %g min\n",msec/(1000.0*60));
  printf("Alloc geo\n");
  Geodesics *geo = (Geodesics*) calloc(*pnvertices, sizeof(Geodesics));
  msec = TimerStop(&mytimer); printf("  t = %g min\n",msec/(1000.0*60));
  printf("Setting\n");
  nth = 0;
  for(vtxno = 0; vtxno < *pnvertices; vtxno++){
    geo[vtxno].vnum = vnum[vtxno];
    for(nthnbr = 0 ; nthnbr < geo[vtxno].vnum; nthnbr++){
      geo[vtxno].v[nthnbr]    = vlist[nth];
      geo[vtxno].dist[nthnbr] = dist[nth];
      nth ++;
    }
  }
  free(vlist);
  free(dist);
  free(vnum);
  msec = TimerStop(&mytimer); printf("  t = %g min\n",msec/(1000.0*60));

  fclose(fp);

  return(geo);
}

double MRISsphereDist(MRIS *sphere, VERTEX *vtx1, VERTEX *vtx2)
{
  double r, dot, d;
  r = sqrt(vtx1->x * vtx1->x + vtx1->y * vtx1->y + vtx1->z * vtx1->z);
  dot = (vtx1->x * vtx2->x + vtx1->y * vtx2->y + vtx1->z * vtx2->z)/(r*r);
  d = acos(dot)*r;
  //printf("r = %g, dot = %g, d = %g\n",r,dot,d);
  return(d);
}

double geodesicsCheckSphereDist(MRIS *sphere, Geodesics *geod)
{
  int vtxno;
  double e;

  e = 0;
  for(vtxno = 0; vtxno < sphere->nvertices; vtxno++){
    int nbrvtxno,nthnbr;
    VERTEX *vtx1,*vtx2;
    double dA, dB;
    vtx1 = &sphere->vertices[vtxno];
    for(nthnbr = 0 ; nthnbr < geod[vtxno].vnum; nthnbr++) {
      nbrvtxno = geod[vtxno].v[nthnbr];
      dA = geod[vtxno].dist[nthnbr];
      vtx2 = &sphere->vertices[nbrvtxno];
      dB = MRISsphereDist(sphere, vtx1, vtx2);
      printf("%5d %5d %6.4f %6.4f  %8.6f %8.4f\n",vtxno,nbrvtxno,dA,dB,dA-dB,e);
      e += fabs(dA-dB);
    }
  }

  return(e);
}


int VtxVolIndexCompare(const void *a, const void *b)
{

  VTXVOLINDEX *vvi1;
  VTXVOLINDEX *vvi2;

  vvi1 = ((VTXVOLINDEX *) a);
  vvi2 = ((VTXVOLINDEX *) b);

  /* Sort by Volume Index */
  if(vvi1->volindex < vvi2->volindex) return(-1);
  else if(vvi1->volindex > vvi2->volindex) return(+1);
  else{
    if(vvi1->vtxno < vvi2->vtxno) return(-1);
    return(+1);
  }
}

int VtxVolIndexSort(VTXVOLINDEX *vvi, int nlist)
{
  qsort(vvi,nlist,sizeof(VTXVOLINDEX),VtxVolIndexCompare);
  return(0);
}

int VtxVolIndexPrint(VTXVOLINDEX *vvi, int nlist)
{
  int n;
  for(n=0; n < nlist; n++)
    printf("%3d %5d %7d  %6.4f\n",n,vvi[n].vtxno,vvi[n].volindex,vvi[n].dist);
  return(0);
}

int VtxVolIndexSortTest(int nlist)
{
  VTXVOLINDEX *vvi,*vvi2;
  int n,nunique;
  vvi = (VTXVOLINDEX *) calloc(sizeof(VTXVOLINDEX),nlist);
  for(n=0; n < nlist; n++){
    vvi[n].vtxno = round(100*drand48());
    vvi[n].volindex = round(100*drand48());
    vvi[n].dist = 10*drand48();
  }
  vvi[0].volindex = vvi[nlist-1].volindex;
  vvi[1].volindex = vvi[nlist-1].volindex;
  vvi[2].volindex = vvi[nlist-2].volindex;
  
  VtxVolIndexPrint(vvi, nlist);
  VtxVolIndexSort(vvi, nlist);
  printf("-----------------\n");
  VtxVolIndexPrint(vvi, nlist);

  vvi2 = VtxVolIndexUnique(vvi, nlist, &nunique);
  printf("%d -----------------\n", nunique);
  VtxVolIndexPrint(vvi2, nunique);

  return(0);
}

VTXVOLINDEX *VtxVolIndexUnique(VTXVOLINDEX *vvi, int nlist, int *nunique)
{
  int n,nthu,nhits;
  VTXVOLINDEX *uvvi;
  double distsum;

  VtxVolIndexSort(vvi, nlist);

  uvvi = (VTXVOLINDEX *) calloc(sizeof(VTXVOLINDEX),nlist);

  nthu = 0;
  n = 0;
  memcpy(&uvvi[nthu],&vvi[n],sizeof(VTXVOLINDEX));
  distsum = vvi[n].dist;
  nhits = 1;
  for(n=1; n<nlist; n++) {
    if(uvvi[nthu].volindex != vvi[n].volindex) {
      uvvi[nthu].dist = distsum/nhits;
      nthu ++;
      memcpy(&uvvi[nthu],&vvi[n],sizeof(VTXVOLINDEX));
      distsum = vvi[n].dist;
      nhits = 1;
    }
    else {
      distsum += vvi[n].dist;
      nhits ++;
    }
  }
  *nunique = nthu+1;
  return(uvvi);
}

VTXVOLINDEX *VtxVolIndexPack(Geodesics *geod, int vtxno, MRI *volindex)
{
  int nthnbr,nbrvtxno;
  VTXVOLINDEX *vvi;
  //Geodesics *geodvtx = &geod[vtxno];

  vvi = (VTXVOLINDEX *) calloc(sizeof(VTXVOLINDEX),geod[vtxno].vnum);

  for(nthnbr = 0 ; nthnbr < geod[vtxno].vnum; nthnbr++) {
    nbrvtxno = geod[vtxno].v[nthnbr];
    vvi[nthnbr].vtxno = nbrvtxno;
    vvi[nthnbr].dist = geod[vtxno].dist[nthnbr];
    vvi[nthnbr].volindex = MRIgetVoxVal(volindex,nbrvtxno,0,0,0);
  }
  //vvi[nthnbr].vtxno = vtxno;
  //vvi[nthnbr].dist = 0;


  return(vvi);
}

Geodesics *VtxVolPruneGeod(Geodesics *geod, int vtxno, MRI *volindex)
{
  VTXVOLINDEX *vvi, *uvvi; 
  Geodesics *ugeod;
  int nunique;
  int nthnbr;

  vvi = VtxVolIndexPack(geod, vtxno, volindex);
  //if(vtxno == 123093)  VtxVolIndexPrint(vvi, geod[vtxno].vnum);

  uvvi = VtxVolIndexUnique(vvi, geod[vtxno].vnum, &nunique);

  ugeod = (Geodesics *) calloc(sizeof(Geodesics),1);
  ugeod->vnum = nunique;
  for(nthnbr = 0 ; nthnbr < ugeod->vnum; nthnbr++) {
    ugeod->v[nthnbr] = uvvi[nthnbr].vtxno;
    ugeod->dist[nthnbr] = uvvi[nthnbr].dist;
  }

  free(vvi);
  free(uvvi);

  return(ugeod);
}


/*
  \fn LABEL *LabelAddVoxel(MRIS *surf, MRI *volsurf, LTA *lta, MRI *vol, int col, int row, int slice, double dmax, LABEL *slabel)
  surf is the surface, get geom from vg field
  lta - registration between volsurf and vol (NULL to align based on scanner coords)
  vol is the geometry of the voxel to add
  col, row, slice are indices of the voxel to add
  dmax - maximum allowable distance between voxel and surface to be added to label
  slabel - surface-based label.
 */

//LABEL *LabelAddVoxelToSurfLabel(MRIS *surf, LTA *lta, MRI *vol, int col, int row, int slice, double dmax, LABEL *slabel) 
//{
//}
