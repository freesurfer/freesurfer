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
  MRIS **surfs;  // list of surfaces (eg, left and right)
  int nsurfs;    // number of surfaces in list
  MRI *template; // volume template for voxel to add
  double dmax;   // max allowable distance in mm
  LTA *vol2surf; // reg bet vol and surf, direction unimportant, set to null to use header reg
  LABEL **labels;  // list of labels, one for each surface
  MHT **hashes;  // hash tables, one for each vertex
  float hashres; // eg, 16
  MRI **masks;   // label masks foreach surface
  int nsegs; // subsample each voxel into nsegs along each dim to fill holes
  MATRIX *volcrs2surfxyz; // precomputed matrix converts vol CRS to surface tkRAS
  MATRIX *surfxyz2volcrs; // inverse
  int nhopsmax; 
  int DilateWithinVoxel; // require neighboring vertices to be within voxel to include in label
  int DilateContiguous; // require neighboring vertices to be continguous with the center vertex
  int debug; // set to 1 to print out a bunch of stuff
} LABEL2SURF;
LABEL2SURF *L2Salloc(int nsurfs, char *subject);
int L2Sinit(LABEL2SURF *l2s);
int L2SaddPoint(LABEL2SURF *l2s, double col, double row, double slice, int Operation);
int L2SaddVoxel(LABEL2SURF *l2s, double col, double row, double slice, int Operation);
int L2Sfree(LABEL2SURF **pl2s);

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
  MRIS *surf, *surf2;
  int msec, nvertices; //vtxno=0;
  struct timeb  mytimer;
  MRI *mri, *mriindex, *mri2;
  Geodesics *geod;
  float maxdist;
  double d;
  int c,r,s;
  LABEL2SURF *l2s;
  LTA *lta;

  vg_isEqual_Threshold = 10e-4;

  if(0){
    mri  = MRIread(argv[1]);
    surf = MRISread(argv[2]);
    surf2 = MRISread(argv[3]);
    lta = LTAread(argv[4]);
  }else{
    mri  = MRIread("template.nii.gz");
    surf = MRISread("lh.white");
    surf2 = MRISread("rh.white");
    lta = LTAread("register.dof6.lta");
  }
  printf("\n");
  printf("alloc \n");
  l2s = L2Salloc(2, "");
  l2s->template = mri;
  l2s->surfs[0] = surf;
  l2s->surfs[1] = surf2;
  l2s->dmax = 3;
  l2s->hashres = 16;
  l2s->vol2surf = lta;
  l2s->nsegs = 1;
  l2s->debug = 0;
  //l2s->nhopsmax = 10;
  printf("init \n");
  L2Sinit(l2s);
  printf("nhopsmax %d\n",l2s->nhopsmax);
  printf("loop \n");
  s = 17;
  for(c = 0; c < mri->width; c++){
    //printf("%2d ",c); fflush(stdout);
    for(r = 0; r < mri->height; r++) L2SaddVoxel(l2s, c, r, s, 1);
    //printf("\n");
  }
  // remove, erase a few
  for(c = 10; c < mri->width-10; c++) L2SaddVoxel(l2s, c, 15, s, 0);
  L2Sfree(&l2s);
  exit(0);

  LabelWrite(l2s->labels[0],"./my.label0");
  LabelWrite(l2s->labels[1],"./my.label1");
  MRIwrite(l2s->masks[0],"mask0.mgh");
  MRIwrite(l2s->masks[1],"mask1.mgh");
  L2Sfree(&l2s);

  exit(0);


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


LABEL2SURF *L2Salloc(int nsurfs, char *subject)
{
  int n;
  LABEL2SURF *l2s;
  l2s = (LABEL2SURF *) calloc(sizeof(LABEL2SURF),1);
  l2s->surfs  = (MRIS **) calloc(sizeof(MRIS*),nsurfs);
  l2s->hashes = (MHT **)  calloc(sizeof(MHT*), nsurfs);
  l2s->nsurfs = nsurfs;
  l2s->labels = (LABEL **) calloc(sizeof(LABEL*), nsurfs);
  for(n=0; n < nsurfs; n++) l2s->labels[n] = LabelAlloc(100,subject,NULL);
  l2s->masks = (MRI **) calloc(sizeof(MRI*), nsurfs);
  l2s->nsegs = 1;
  l2s->nhopsmax = -1; // use -1 to keep track of whether it has been set or not
  l2s->DilateWithinVoxel = 1; // turn on by default
  l2s->DilateContiguous = 0; // turn off by default
  return(l2s);
}
int L2Sfree(LABEL2SURF **pl2s)
{
  int n;
  LABEL2SURF *l2s = *pl2s;

  free(l2s->surfs); // only free the pointer to the pointers
  for(n = 0; n < l2s->nsurfs; n++){
    LabelFree(&l2s->labels[n]);
    MHTfree(&l2s->hashes[n]);
    MRIfree(&l2s->masks[n]);
  }
  free(l2s->labels);
  free(l2s->hashes);
  free(l2s->masks);
  MatrixFree(&l2s->volcrs2surfxyz);
  free(l2s);
  return(0);
}

/*!
  \fn int L2Sinit(LABEL2SURF *l2s)
  \brief Initializes the L2S structure. The caller must have run
  L2Salloc() and set the surfs and vol2surf structures. If a header
  registration is good enough, then set vol2surf=NULL. Note that
  vol2surf can go in either direction; this function will determine
  which way to go from the template volume and surface vg geometries
  and handle it appropriately. This also computes the vertex coords in
  CRS space and stores them in the tx, ty, tz elements of the vertex,
  so the caller needs make sure that this is not overwriting something
  needed and to make sure not to overwrite it. Also computes a value
  for nhopsmax; the value is a reasonable guess based on the voxel
  size.
*/
int L2Sinit(LABEL2SURF *l2s)
{
  int n;

  // initialize hashes 
  for(n = 0; n < l2s->nsurfs; n++){
    l2s->hashes[n] = MHTfillVertexTableRes(l2s->surfs[n], NULL,
					   CURRENT_VERTICES,l2s->hashres);
    if(l2s->hashes[n] == NULL){
      printf("ERROR: L2Sinit(): MHTfillVertexTableRes() failed\n");
      return(-1);
    }
    l2s->masks[n] = MRIalloc(l2s->surfs[n]->nvertices,1,1,MRI_INT);
  }

  // compute the matrix that maps the template volume CRS to surface RAS
  // volcrs2surfxyz = K*inv(Vs)*R*Vv
  MATRIX *K,*Vs,*invVs,*Vv,*R;
  LTA *lta;
  K = TkrVox2RASfromVolGeom(&l2s->surfs[0]->vg); // vox2tkras of surface
  Vs = vg_i_to_r(&l2s->surfs[0]->vg); // vox2scanneras of surface
  Vv = MRIxfmCRS2XYZ(l2s->template,0); // vox2scanneras of template volume
  invVs = MatrixInverse(Vs,NULL);
  if(l2s->vol2surf == NULL) R = MatrixIdentity(4,NULL);
  else {
    // A registration LTA has been passed, make sure it is consisent
    // with the input geometries and that it points in the right
    // direction.  Note that these function use the global
    // vg_isEqual_Threshold variable which should be set to something
    // small, like 10^-4
    int DoInvert=0;
    VOL_GEOM vgvol;
    getVolGeom(l2s->template, &vgvol);
    if(!vg_isEqual(&vgvol, &(l2s->vol2surf->xforms[0].src))){
      // The src does not match the template, so try the dst
      if(!vg_isEqual(&l2s->surfs[0]->vg, &(l2s->vol2surf->xforms[0].src))){
	printf("ERROR: L2Sinit(): neither registration vgs match template %g\n",vg_isEqual_Threshold);  
	return(-1);
      }
      // Verify that the template matches the dst
      if(!vg_isEqual(&vgvol, &(l2s->vol2surf->xforms[0].dst))){
	printf("ERROR: L2Sinit(): registration does not match volume vg %g\n",vg_isEqual_Threshold);
	return(-1);
      }
      DoInvert = 1;
    }
    else{
      // The source matches, but does the target?
      if(!vg_isEqual(&l2s->surfs[0]->vg, &(l2s->vol2surf->xforms[0].dst))){
	printf("ERROR: L2Sinit(): registration does not match surface vg %g\n",vg_isEqual_Threshold);
	return(-1);
      }
    }
    // Copy the LTA
    lta = LTAcopy(l2s->vol2surf,NULL);
    // Make sure the type is RAS2RAS
    LTAchangeType(lta, LINEAR_RAS_TO_RAS);
    if(DoInvert){
      if(l2s->debug) printf("L2Sinit(): inverting reg\n");
      R = MatrixInverse(lta->xforms[0].m_L,NULL);
    }
    else R = MatrixCopy(lta->xforms[0].m_L,NULL);
    LTAfree(&lta);
  }
  // Now finally compute it
  l2s->volcrs2surfxyz = MatrixMultiplyD(K,invVs,NULL);
  MatrixMultiplyD(l2s->volcrs2surfxyz,R,l2s->volcrs2surfxyz);
  MatrixMultiplyD(l2s->volcrs2surfxyz,Vv,l2s->volcrs2surfxyz);
  l2s->surfxyz2volcrs = MatrixInverse(l2s->volcrs2surfxyz,NULL);
  if(l2s->debug){
    printf("L2Sinit(): \n");
    printf("K = [\n");
    MatrixPrint(stdout,K);
    printf("];\n");
    printf("invVs = [\n"); 
    MatrixPrint(stdout,invVs);
    printf("];\n");
    printf("R = [\n");
    MatrixPrint(stdout,R);
    printf("];\n");
    printf("Vv = [\n");
    MatrixPrint(stdout,Vv);
    printf("];\n");
    printf("volcrs2surfxyz = [\n");
    MatrixPrint(stdout,l2s->volcrs2surfxyz);
    printf("];\n");
  }
  MatrixFree(&K);
  MatrixFree(&Vs);
  MatrixFree(&Vv);
  MatrixFree(&invVs);
  MatrixFree(&R);

  // Now precompute the CRS coords of each vertex and save in the
  // t{xyz} part of the vertex structure. This makes constraining the
  // label dilation to within the voxel easier
  int vtxno;
  MATRIX *ras, *crs=NULL;
  ras = MatrixAlloc(4,1,MATRIX_REAL);
  ras->rptr[4][1] = 1;
  for(n = 0; n < l2s->nsurfs; n++){
    for(vtxno = 0; vtxno < l2s->surfs[n]->nvertices; vtxno++){
      ras->rptr[1][1] = l2s->surfs[n]->vertices[vtxno].x;
      ras->rptr[2][1] = l2s->surfs[n]->vertices[vtxno].y;
      ras->rptr[3][1] = l2s->surfs[n]->vertices[vtxno].z;
      crs = MatrixMultiplyD(l2s->surfxyz2volcrs,ras,crs);
      l2s->surfs[n]->vertices[vtxno].tx = crs->rptr[1][1];
      l2s->surfs[n]->vertices[vtxno].ty = crs->rptr[2][1];
      l2s->surfs[n]->vertices[vtxno].tz = crs->rptr[3][1];
    }
  }
  MatrixFree(&ras);
  MatrixFree(&crs);

  // Take a stab at automatically computing the maximum number of hops. The basic
  // idea is to scale the number based upon the voxel size. The 3 is just something
  // I guessed. 
  if(l2s->nhopsmax == -1)
    l2s->nhopsmax = 3*round(sqrt(pow(l2s->template->xsize,2.0)+pow(l2s->template->ysize,2)+
				 pow(l2s->template->zsize,2)));


  return(0);
}

/*!
  \fn int L2SaddVoxel(LABEL2SURF *l2s, double col, double row, double slice, int Operation)
  \brief Divides the voxel into l2s->nsegs parts and then runs L2SaddPoint(). This is an 
  alternative way to include vertices inside the voxel that are not necessarily the closest
  vertex. The problem is that it adds any vertex regardless of how far away it is along the
  surface from the center vertex. L2SaddPoint() has a built-in mechanism for solving this
  problem, and it may be better. See L2SaddPoint() for argument descriptions.
*/
int L2SaddVoxel(LABEL2SURF *l2s, double col, double row, double slice, int Operation)
{
  double c, r, s, dseg;
  int ret, kc, kr, ks;

  if(l2s->nsegs == 1){
    ret = L2SaddPoint(l2s, col, row, slice, Operation);
    return(ret);
  }

  dseg = 1.0/(l2s->nsegs-1);
  // using l2s->nsegs as the upper limit (instead of l2s->nsegs+1) means that
  // the "end" side of the voxel is not included (but the "start" side
  // is).
  for(kc=0; kc < l2s->nsegs; kc++){
    c = col + kc*dseg - 0.5;
    for(kr=0; kr < l2s->nsegs; kr++){
      r = row + kr*dseg - 0.5;
      for(ks=0; ks < l2s->nsegs; ks++){
	s = slice + ks*dseg - 0.5;
	ret = L2SaddPoint(l2s, c, r, s, Operation);
	if(ret < 0) return(ret);
      }
    }
  }
  // Make sure to add the center voxel
  ret = L2SaddPoint(l2s, col, row, slice, Operation);

  return(ret);  
}

/*!
  \fn int L2SaddPoint(LABEL2SURF *l2s, double col, double row, double slice, int Operation)
  \brief Adds (Operation==1) or removes (Operation!=1) a voxel from a label based
  on its proximity to a surface. The caller must have run L2Salloc() and L2Sinit() and
  set the appropriate variables in the L2S structure. There is a label for each 
  surface. If the CRS is within dmax to a surface, then it is assigned to the label
  of the closest surface. This allows a label to be assigned, for example, to the lh
  but not the rh if both lh and rh surfaces are in the voxel. 

*/
int L2SaddPoint(LABEL2SURF *l2s, double col, double row, double slice, int Operation)
{
  int n,nmin,vtxnominmin,vtxno,pointno;
  VERTEX v;
  static MATRIX *crs=NULL, *ras=NULL;
  float dminsurf,dminmin;
  LV *lv;
  LABEL *label;

  if(crs==NULL) {
    crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
  }
  crs->rptr[1][1] = col;
  crs->rptr[2][1] = row;
  crs->rptr[3][1] = slice;

  // Compute the TkReg RAS of this CRS in the surface volume space
  ras = MatrixMultiplyD(l2s->volcrs2surfxyz,crs,ras);

  // Load it into a vertex structure
  v.x = ras->rptr[1][1];
  v.y = ras->rptr[2][1];
  v.z = ras->rptr[3][1];
  if(l2s->debug > 1) printf("---------------\n%g %g %g   (%5.2f %5.2f %5.2f)\n",col,row,slice,v.x,v.y,v.z);

  // Go through each surface to find the surface and vertex in the surface that
  // the CRS is closest to. 
  dminmin = 10e10; // minimum distance over all surfaces
  vtxnominmin = -1; // number of closest vertex
  nmin = -1; // index of the surface with the closest vertex
  for(n = 0; n < l2s->nsurfs; n++){
    vtxno = MHTfindClosestVertexNo(l2s->hashes[n],l2s->surfs[n],&v,&dminsurf);
    if(vtxno >= 0){
      if(l2s->debug > 1) printf("%3d %6d (%5.2f %5.2f %5.2f) %g\n",n,vtxno,
	l2s->surfs[n]->vertices[vtxno].x,l2s->surfs[n]->vertices[vtxno].y,
	l2s->surfs[n]->vertices[vtxno].z,dminsurf);
    }
    // should check whether there was a hash error?
    if(dminsurf < dminmin && dminsurf < l2s->dmax){
      dminmin = dminsurf;
      vtxnominmin = vtxno;
      nmin = n;
    }
  }

  // If it does not meet the distance criteria, then return 1
  if(vtxnominmin == -1) return(1);

  // Select the label of the winning surface
  label = l2s->labels[nmin];

  // Below is a method to expand the label along the surface to
  // include other nearby vertices to the winning vertex. The hoplist
  // is a list of vertices for each ring of nearest neighbors. The 0th
  // hop is always the center vertex.
  SURFHOPLIST *shl;
  int nthnbr, nnbrs, nthhop, nthhoplast,nhits;
  shl = SetSurfHopList(vtxnominmin, l2s->surfs[nmin], l2s->nhopsmax);
  nthhoplast = 0;
  for(nthhop = 0; nthhop < l2s->nhopsmax; nthhop++){
    nnbrs = shl->nperhop[nthhop]; // number of neighbrs in the nthhop ring
    nhits = 0;

    // loop through the neighbors nthhop rings away
    for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
      vtxno = shl->vtxlist[nthhop][nthnbr];

      if(l2s->DilateWithinVoxel){
	// Require nearby vertices to be within the given voxel, where
	// "within" is defined as +/-0.5 of the voxel because the
	// voxel coords are at the center of the voxel
	v.x = l2s->surfs[nmin]->vertices[vtxno].tx;
	v.y = l2s->surfs[nmin]->vertices[vtxno].ty;
	v.z = l2s->surfs[nmin]->vertices[vtxno].tz;
	if(nthhop != 0){
	  // have to allow nthhop=0 because the closest vertex may be outside of the
	  // voxel, depending upon the distance threshold
	  if(fabs(v.x-col) > 0.501 || fabs(v.y-row) > 0.501 || fabs(v.z-slice) > 0.501)
	    continue;
	}
      }

      nthhoplast = nthhop; // keep track of the most recent hop with a hit
      nhits++; // keep track of the number of vertices hit in this hop

      // Get the 1-based label point number. Using 1-based so that 0
      // can code for no label at that point (yet)
      pointno = MRIgetVoxVal(l2s->masks[nmin],vtxno,0,0,0);
      
      if(Operation == 1){ // Add vertex to label
	if(pointno!=0) continue; //already there
	// If it gets here, then add the vertex
	if(l2s->debug) printf("Adding surf=%d hop=%d vtxno=%d  np=%5d   %g %g %g   (%5.2f %5.2f %5.2f)\n",
			      nmin,nthhop,vtxno,label->n_points,col,row,slice,v.x,v.y,v.z);
	// Set value of the mask to the label point number + 1 (so 0=nolabel)
	MRIsetVoxVal(l2s->masks[nmin],vtxno,0,0,0, label->n_points+1);
	
	// Check whether need to alloc more points in label
	if(label->n_points >= label->max_points)
	  LabelRealloc(label, nint(label->max_points*1.5)) ;
	
	// Finally, add this point
	lv = &(label->lv[label->n_points]);
	lv->vno = vtxno;
	lv->x = l2s->surfs[nmin]->vertices[vtxno].x;
	lv->y = l2s->surfs[nmin]->vertices[vtxno].y;
	lv->z = l2s->surfs[nmin]->vertices[vtxno].z;
	lv->stat = dminmin;
	
	// Incr the number of points in the label
	label->n_points++;
      }
      else { // Operation != 1, Remove vertex from label
	if(pointno==0) continue; // not already there, cant remove it
	lv = &(label->lv[pointno-1]);
	// If it gets here, then remove the vertex by copying the last
	// point into this point, then decrementing the number of points
	if(l2s->debug) printf("Removing surf=%d vtxno=%d %g %g %g   (%5.2f %5.2f %5.2f)\n",
			      nmin,vtxno,col,row,slice,v.x,v.y,v.z);
	memcpy(&label->lv[pointno-1],&(label->lv[label->n_points-1]),sizeof(LV));
	MRIsetVoxVal(l2s->masks[nmin],vtxno,0,0,0,0);
	label->n_points --;
      }
    }

    // If it has gone through a ring and not gotten a hit (nhits==0), then future hits
    // will not be contiguous with the current label. This can happen if the surface
    // folds back into the voxel. This can easily happen in fMRI-sized voxels. It is not
    // clear how this should be handled, so I put l2s->DilateContiguous into the structure
    // to allow or forbid disconinuous labels. 
    if(l2s->DilateContiguous && nhits == 0) break;
  }
  SurfHopListFree(&shl);

  // check whether it might be a good idea to increase the number of hops
  if(nthhoplast == l2s->nhopsmax-1)
    if(l2s->debug) printf("WARNING: hop saturation nthsurf=%d, cvtxno = %d\n",nmin,vtxnominmin);

  return(0);
}


