/**
 * @file  dngtester.c
 * @brief dougs super special test code
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.49 $
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
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "mrisurf.h"
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "region.h"
#include "machine.h"
#include "fio.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "fmriutils.h"
#include "gca.h"
#include "gcsa.h"
#include "fsgdf.h"
#include "icosahedron.h"
#include "gca.h"
#include "gcamorph.h"
#include "DICOMRead.h"
#include "fsenv.h"
#include "qdecutils.h"
#include "dti.h"
#include "registerio.h"
#include "timer.h"
#include "evschutils.h"
#include "matrix.h"
#include "matfile.h"
#include "randomfields.h"
#include "mri2.h"
#include "annotation.h"
#include "mrisutils.h"
#include "image.h"
#include "retinotopy.h"
#include "bfileio.h"
#include "cma.h"

// setenv SUBJECTS_DIR /space/greve/1/users/greve/subjects
// /autofs/space/greve_001/users/greve/dev/trunk/dngtester
// ./dngtester register.dat func.mgz func-in-m3z.mgh
// ./dngtester identity.dat ~/subjects/fbirn-anat-101/mri/orig.mgz orig-in-m3z.mgh
// tkregister2 --targ orig-in-m3z.mgh --mov func-in-m3z.mgh --regheader --reg tmp.reg
// Atlas: $SUBJECTS_DIR/avgtst/mri/T1MLV.mgz

MRI *MRIsetSliceNo(MRI *mri, MRI *out);

LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R);
char *Progname = "dngtester";

char *outvolfile=NULL;
char *subject = NULL;
char *regfile = NULL;
TRANSFORM *Rtransform;  //types : M3D, M3Z, LTA, FSLMAT, DAT, OCT(TA), XFM
FSENV     *fsenv;
char gcamfile[1000];
char origfile[1000];
char *sourcefile;
GCAM      *gcam;
MRI *mri, *mri2, *mask=NULL;
MATRIX *V, *W, *m_tmp;
float ipr, bpr, intensity;
MATRIX *R;
int float2int;
int err,n;
struct timeb  mytimer;
int msecFitTime;
double a,b,zthresh,pthresh;
char *SUBJECTS_DIR;
char tmpstr[2000];
MRIS *surf, *surf2;
char *hemi;

MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var);
void printrgb(void);

COLOR_TABLE *CTABaddEntry(COLOR_TABLE *ctold, const char *name);
int MRISmercator(MRIS *surf);
IMAGE *I;
MHT *lhwhite_hash;
MRIS *lhwhite,*white,*sphere;
MRIS *surfs[100], *mris;
int *XNbrVtxNo, nXNbrs;
double *XNbrDotProd;

MRIS *MRISaverageSurfaces(int nsurfaces, MRIS **surfs);
MATRIX *MatrixLoadFSL(char *fname);
int find_path ( int* vert_vno, int num_vno, int max_path_length,
                int* path, int* path_length, MRIS *mris );
MRI *seg1, *seg2;
int mygd(char *fname);
MRI *MRIaddB(MRI *mri1, MRI *mri2, MRI *mriadd);

/*----------------------------------------*/
int main(int argc, char **argv) 
{
  int err,area32p,area32v,superiorfrontal,medialorbitofrontal;
  int rostralanteriorcingulate, rostralmiddlefrontal;
  int index,k,c,r,s,nvox,a1,a2,k1,k2;
  int *nunits;
  char *parcnames[10], *annot1, *annot2;
  COLOR_TABLE *ct ;
  VERTEX *vtx,*vtx1w,*vtx2w,*vtx1s,*vtx2s;
  float dlhw,DotProdThresh;
  int  lhwvtx;
  double sumval;
  int nsegid1, *segidlist1;
  int nsegid2, *segidlist2;
  int nsegs, *segidlist,vtxno1,vtxno2;
  double *area1, *area2, *area12, *dice;
  double f,radius,radius2,DotProd,theta,d2,d3,d3Sqr;
  double fwhm,fwhmSqr,p,z;
  COLOR_TABLE *ctab = NULL;
  RFS *rfs;

  rfs = RFspecInit(53,NULL);
  rfs->name = strcpyalloc("gaussian");
  rfs->params[0] = 0;
  rfs->params[1] = 1;
  for(z=0; z < 10.0; z+=.1){
    //p = sc_cdf_flat_Q(z,0,1);
    p = RFstat2PVal(rfs, z);
    printf("%7.4f %8.4f\n",z,-log10(p));
  }
  exit(1);



  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  subject = argv[1];
  hemi = argv[2];
  fwhm = 10;
  fwhmSqr = fwhm*fwhm;

  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  white = MRISread(tmpstr);
  if(!white) exit(1);
  printf("White Total Area = %g \n",white->total_area);

  sprintf(tmpstr,"%s/%s/surf/%s.sphere",SUBJECTS_DIR,subject,hemi);
  sphere = MRISread(tmpstr);
  if(!sphere) exit(1);
  printf("Sphere Total Area = %g \n",sphere->total_area);

  f = sqrt(white->total_area/sphere->total_area);
  printf("f = %g \n",f);

  // normalize radius to 1
  vtx1s = &(sphere->vertices[100]);
  for(vtxno1 = 0; vtxno1 < sphere->nvertices; vtxno1++){
    vtx1s = &(sphere->vertices[vtxno1]);
    radius = sqrt(vtx1s->x*vtx1s->x + vtx1s->y*vtx1s->y + vtx1s->z*vtx1s->z);
    vtx1s->x *= f;
    vtx1s->y *= f;
    vtx1s->z *= f;
  }
  MRIScomputeMetricProperties(sphere);
  printf("Sphere Total Area = %g \n",sphere->total_area);

  radius = sqrt(vtx1s->x*vtx1s->x + vtx1s->y*vtx1s->y + vtx1s->z*vtx1s->z);
  printf("Radius %f\n",radius);
  radius2 = radius*radius;

  printf("Alloc\n");
  mri = MRIallocSequence(white->nvertices, 1, 1, MRI_FLOAT, 1);

  printf("loop\n");
  for(vtxno1 = 0; vtxno1 < white->nvertices-1; vtxno1++){
    if(vtxno1%100 ==0) printf("%6d \n",vtxno1);
    vtx1w = &(white->vertices[vtxno1]);
    vtx1s = &(sphere->vertices[vtxno1]);

    for(vtxno2 = vtxno1+1; vtxno2 < white->nvertices; vtxno2++){
      vtx2w = &(white->vertices[vtxno2]);
      vtx2s = &(sphere->vertices[vtxno2]);

      d3Sqr = SQR(vtx1w->x-vtx2w->x)+SQR(vtx1w->y-vtx2w->y)+SQR(vtx1w->z-vtx2w->z);
      if(d3Sqr > fwhmSqr) continue;
      d3 = sqrt(d3Sqr);

      DotProd = (vtx1s->x*vtx2s->x + vtx1s->y*vtx2s->y + vtx1s->z*vtx2s->z)/radius2;
      if(DotProd > +1) DotProd = +1;
      if(DotProd < -1) DotProd = -1;
      theta = acos(DotProd);
      d2 = f * radius * theta;
      if(d2 < fwhm) continue;

      k = MRIgetVoxVal(mri,vtxno1,0,0,0);
      MRIsetVoxVal(mri,vtxno1,0,0,0, k+1);

      //printf("%6d %6d  %4.1f %4.1f\n",vtxno1,vtxno2,d3,d2);

    }
    //if(vtxno1 > 5000 ) break;
  }

  MRIwrite(mri,"count.mgh");

  exit(0);




  k = 1;
  printf("%d %s\n",k,argv[k]);
  seg1 = MRIread(argv[k]);
  if(seg1 == NULL) exit(1);

  k = 2;
  printf("%d %s\n",k,argv[k]);
  seg2 = MRIread(argv[k]);
  if(seg2 == NULL) exit(1);

  mri = MRIadd(seg1,seg2,NULL);
  if(mri == NULL) exit(1);

  for(k=3; k < argc-1; k++){
    printf("%d %s\n",k,argv[k]);
    seg2 = MRIread(argv[k]);
    if(seg2 == NULL) exit(1);
    mri = MRIadd(mri,seg2,mri);
    if(mri == NULL) exit(1);
  }

  f = 1.0/((double)(argc-2));
  printf("Dividing by %lf\n",f);
  MRImultiplyConst(mri, f, mri);
  printf("Saving to %s\n",argv[argc-1]);
  MRIwrite(mri,argv[argc-1]);
  exit(0);


  mri = MRIsegDiff(seg1,seg2,&k);
  printf("DiffFlag %d\n",k);
  MRIwrite(mri,"aseg.diff.mgz");

  mri2 = MRIsegMergeDiff(seg1, mri);
  MRIwrite(mri2,"aseg.new.mgz");

  exit(0);

  seg1 = MRIread(argv[1]);
  seg2 = MRIread(argv[2]);
  dice = MRIsegDice(seg1, seg2, &nsegs, &segidlist);

  ctab = CTABreadASCII("/space/greve/1/users/greve/freesurfer/FreeSurferColorLUT.txt");
  if (ctab == NULL) {
    printf("ERROR: reading ctab\n");
    exit(1);
  }

  for(k=0; k < nsegs; k++){
    if(segidlist[k] >= 0){
      CTABcopyName(ctab,segidlist[k],tmpstr,sizeof(tmpstr));
      printf("%2d %4d %-30s %7.5lf\n",k,segidlist[k],tmpstr,dice[k]);
    }
  }


  return(0);


  //----------------------------------------------------
  subject = argv[1];
  hemi = argv[2];
  annot1 = argv[3];
  annot2 = argv[4];

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);
  MRIScomputeMetricProperties(surf);

  surf2 = MRISread(tmpstr);
  if(!surf2) exit(1);

  sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot1);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);

  sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot2);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf2, tmpstr);
  if(err) exit(1);

  dice = MRISannotDice(surf, surf2, &nsegs, &segidlist);

  for(k=0; k < nsegs; k++){
    if(segidlist[k] >= 0){
      printf("%2d %4d %-25s  %7.5lf\n",k,segidlist[k],
	     surf->ct->entries[segidlist[k]]->name,
	     dice[k]);
    }
  }



  return(0);
  //----------------------------------------------------

  seg1 = MRISannot2seg(surf,1000);
  seg2 = MRISannot2seg(surf2,1000);

  printf("Generating list of segmentation ids\n");
  segidlist1 = MRIsegIdList(seg1, &nsegid1,0);
  printf("Found %d\n",nsegid1);
  fflush(stdout);

  printf("Generating list of segmentation ids\n");
  segidlist2 = MRIsegIdList(seg1, &nsegid2,0);
  printf("Found %d\n",nsegid2);
  fflush(stdout);

  area1 = (double *) calloc(nsegid1,sizeof(double));
  area2 = (double *) calloc(nsegid1,sizeof(double));
  area12 = (double *) calloc(nsegid1,sizeof(double));

  for(c=0; c < seg1->width; c++){
    for(r=0; r < seg1->height; r++){
      for(s=0; s < seg1->depth; s++){
	a1 = MRIgetVoxVal(seg1,c,r,s,0);
	k1 = -1;
	for(k=0; k < nsegid1; k++) {
	  if(a1 == segidlist1[k]){
	    k1 = k; 
	    break;
	  }
	}
	a2 = MRIgetVoxVal(seg2,c,r,s,0);
	k2 = -1;
	for(k=0; k < nsegid2; k++) {
	  if(a2 == segidlist2[k]){
	    k2 = k; 
	    break;
	  }
	}
	area1[k1] += surf->vertices[c].area;
	area2[k2] += surf->vertices[c].area;
	if(a1 == a2) area12[k1] += surf->vertices[c].area;
      }
    }
  }

  for(k=0; k < nsegid1; k++) {
    printf("%2d %4d   %7.1lf %7.1lf %7.1lf   %6.4f\n",
	   k,segidlist1[k],
	   area1[k],area2[k],area12[k],
	   (float)area12[k]/((area1[k]+area2[k])/2.0));
  }


  exit(1);

  sprintf(tmpstr,"%s/%s/label/%s.aparc.annot",SUBJECTS_DIR,subject,hemi);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);


  MRISfbirnAnnot(surf);
  sprintf(tmpstr,"%s.fbirn.annot",hemi);
  MRISwriteAnnotation(surf, tmpstr);

  exit(1);
  //----------------------------------------------

  mri = MRIread(argv[1]);
  if(mri == NULL) exit(1);

  nvox = 0;
  sumval = 0;
  for (c=0; c < mri->width; c++) {
    for (r=0; r < mri->height; r++) {
      for (s=0; s < mri->depth; s++) {
	sumval += MRIgetVoxVal(mri,c,r,s,0);
	nvox ++;
      }
    }
  }
  printf("%lf\n",sumval/nvox);
  exit(0);

  if(argc <= 1){
    printf("dngtester subject hemi\n");
    exit(1);
  }
  subject = argv[1];
  hemi = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.sphere",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

  sprintf(tmpstr,"%s/%s/surf/%s.occip.patch.flat",SUBJECTS_DIR,subject,hemi);
  MRISreadPatchNoRemove(surf, tmpstr) ;

  RETlogMap(surf, 15, .7, 0, 0);
  mri = MRIcopyMRIS(NULL, surf, 0, "val");
  MRIwrite(mri,"logmap.mgh");

  exit(1);

  //-----------------------------------------------------
  printf("nsurfs %d\n",argc-1);
  for(k=1; k<argc; k++){
    printf("Loading %s\n",argv[k]);
    surf = MRISread(argv[k]);
    if(surf == NULL) exit(1);
    surfs[k-1] = surf;
  }
  surf = MRISaverageSurfaces(argc-1, surfs);
  MRISwrite(surf,"lh.avgsurf");

  exit(1);

  vtx = &(surf->vertices[12282]);
  dlhw = sqrt((vtx->x * vtx->x) + (vtx->y * vtx->y) + (vtx->z * vtx->z));
  DotProdThresh = (100*100)*cos(10/dlhw)*(1.0001);

  XNbrVtxNo   = (int *) calloc(surf->nvertices,sizeof(int));
  XNbrDotProd = (double *) calloc(surf->nvertices,sizeof(double));
  mri = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, 1);

  //-----------------------------------------------
  nXNbrs = 0;
  MRISextendedNeighbors(surf, 12282, 12282, DotProdThresh,
			XNbrVtxNo, XNbrDotProd, 
			&nXNbrs,surf->nvertices,1);
  printf("Found %d neighbors\n",nXNbrs);

  for(k = 0; k < nXNbrs; k++)
    MRIsetVoxVal(mri,XNbrVtxNo[k],0,0,0,1);
  MRIwrite(mri,"xnbr.mgh");

  //-----------------------------------------------
  err = MRISreadPatch(surf, argv[3]);

  for(lhwvtx = 0; lhwvtx < surf->nvertices; lhwvtx++){
    if(surf->vertices[lhwvtx].ripflag) continue;
    for(k = 0; k < surf->nvertices; k++) surf->vertices[k].val2bak = -1;
    nXNbrs = 0;
    MRISextendedNeighbors(surf, lhwvtx, lhwvtx, 5*5,
			  XNbrVtxNo, XNbrDotProd, 
			  &nXNbrs,surf->nvertices,2);
    if(lhwvtx%1000 == 1) printf("%5d %5d neighbors\n",lhwvtx,nXNbrs);
  }

  for(k = 0; k < nXNbrs; k++)
    MRIsetVoxVal(mri,XNbrVtxNo[k],0,0,0,1);

  MRIwrite(mri,"xnbr2.mgh");

  exit (1);

  //----------------------------------------------------
  //R = MatrixDRand48(20,20,NULL);
  R = MatrixIdentity(100,NULL);
  MatrixPrint(stdout,R);
  I = ImageFromMatrix(R, NULL);
  ImageWrite(I, "./myimage.jpg") ;
  return(1);


  //----------------------------------------------------
  subject = argv[1];
  hemi = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

  sprintf(tmpstr,"%s/%s/mri/aseg.mgz",SUBJECTS_DIR,subject);
  mri = MRIread(tmpstr);
  if(!mri) exit(1);

  //#define MIN_NONCORTEX_VERTICES 10
  //lcortex = MRIScortexLabel(surf, mri, MIN_NONCORTEX_VERTICES) ;
  //sprintf(tmpstr,"%s/%s/label/%s.%s.label",SUBJECTS_DIR, subject,hemi,"cortex");
  //printf("writing cortex label to %s...\n", tmpstr) ;
  //LabelWrite(lcortex, tmpstr) ;

  sprintf(tmpstr,"%s/%s/label/%s.aparc.annot",SUBJECTS_DIR,subject,hemi);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);

  parcnames[0] = "caudalmiddlefrontal";
  parcnames[1] = "superiorfrontal";
  parcnames[2] = "rostralmiddlefrontal";
  parcnames[3] = "parsopercularis";
  parcnames[4] = "parstriangularis";
  parcnames[5] = "parsorbitalis";
  parcnames[6] = "lateralorbitofrontal";
  parcnames[7] = "medialorbitofrontal";
  parcnames[8] = "paracentral";
  parcnames[9] = "frontalpole";
  MRISmergeAnnotations(surf, 10, parcnames, "frontal");

  parcnames[0] = "superiortemporal";
  parcnames[1] = "entorhinal";
  parcnames[2] = "temporalpole";
  parcnames[3] = "fusiform";
  parcnames[4] = "inferiortemporal";
  parcnames[5] = "middletemporal";
  parcnames[6] = "parahippocampal";
  parcnames[7] = "bankssts";
  parcnames[8] = "transversetemporal";
  MRISmergeAnnotations(surf, 9, parcnames, "temporal");

  parcnames[0] = "supramarginal";
  parcnames[1] = "inferiorparietal";
  parcnames[2] = "superiorparietal";
  parcnames[3] = "precuneus";
  MRISmergeAnnotations(surf, 4, parcnames, "parietal");

  parcnames[0] = "pericalcarine";
  parcnames[1] = "cuneus";
  parcnames[2] = "lingual";
  parcnames[3] = "lateraloccipital";
  MRISmergeAnnotations(surf, 4, parcnames, "occipital");

  parcnames[0] = "isthmuscingulate";  
  parcnames[1] = "posteriorcingulate";
  parcnames[2] = "caudalanteriorcingulate";
  parcnames[3] = "rostralanteriorcingulate";
  MRISmergeAnnotations(surf, 4, parcnames, "cingulate");

  sprintf(tmpstr,"%s.lobes.annot",hemi);
  MRISwriteAnnotation(surf, tmpstr);

  //------------------------------------------------------------
  ct = CTABaddEntry(surf->ct,"area32p");
  CTABfree(&surf->ct);
  surf->ct = ct;

  ct = CTABaddEntry(surf->ct,"area32v");
  CTABfree(&surf->ct);
  surf->ct = ct;

  CTABfindName(surf->ct, "area32p", &area32p);
  CTABfindName(surf->ct, "area32v", &area32v);
  CTABfindName(surf->ct, "superiorfrontal", &superiorfrontal);
  CTABfindName(surf->ct, "medialorbitofrontal", &medialorbitofrontal);
  CTABfindName(surf->ct, "rostralanteriorcingulate", &rostralanteriorcingulate);
  CTABfindName(surf->ct, "rostralmiddlefrontal", &rostralmiddlefrontal);

  mri = MRISfbirnMask_MOF_RACing(surf);
  MRISdilateConfined(surf, mri, medialorbitofrontal, 12, area32v);

  mri = MRISfbirnMask_SFG_Cing(surf);
  MRISdilateConfined(surf, mri, superiorfrontal, 12, area32p);


  nunits = (int *)calloc(surf->ct->nentries, sizeof(int)) ;
  nunits[rostralanteriorcingulate] = 3;
  nunits[rostralmiddlefrontal] = 3;
  nunits[superiorfrontal] = 5;
  nunits[area32p] = 2;
  //nunits[area32v] = 2;

  MRISdivideAnnotation(surf, nunits) ;

  CTABfindName(surf->ct, "area32p", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area32p_pseudo");

  CTABfindName(surf->ct, "area32p_div2", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area32a_pseudo");

  CTABfindName(surf->ct, "area32v", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area32v_pseudo");

  //CTABfindName(surf->ct, "area32v_div2", &index);
  //sprintf(surf->ct->entries[index]->name, "%s","area32v_pseudo");


  CTABfindName(surf->ct, "rostralanteriorcingulate", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area24d_pseudo");

  CTABfindName(surf->ct, "rostralanteriorcingulate_div2", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area24pg_pseudo");

  CTABfindName(surf->ct, "rostralanteriorcingulate_div3", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area24v_pseudo");

  sprintf(tmpstr,"%s.fbirn0.annot",hemi);
  MRISwriteAnnotation(surf, tmpstr);

  return(0);

  printrgb();
  return(0);

  printf("Reading\n");
  mri = MRIread(argv[1]);
  printf("chunck %d\n",mri->ischunked);

  printf("Done\n");

  exit(0);


  mri2 = MRIvote(mri, mask, NULL);
  MRIwrite(mri2,"vote.mgh");

  return(0);
  exit(0);
}

/*--------------------------------------------------------*/
MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var) {
  int c, r, s, f;
  double val,sumsqval, sumval;
  int nvox_per_row, nvox_per_slice, bytes_per_vol;
  void *p;
  val = 0;

  if (DOF < 0) DOF = fmri->nframes;

  if (var==NULL) {
    var = MRIallocSequence(fmri->width, fmri->height, fmri->depth,
                           MRI_FLOAT, 1);
    if (var==NULL) {
      printf("ERROR: fMRIvariance: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(fmri,var);
  } else {
    if (var->width  != fmri->width ||
        var->height != fmri->height ||
        var->depth  != fmri->depth) {
      printf("ERROR: fMRIvariance: output dimension mismatch\n");
      return(NULL);
    }
    if (var->type != MRI_FLOAT) {
      printf("ERROR: fMRIvariance: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  nvox_per_row = fmri->width;
  nvox_per_slice = fmri->width * fmri->height;
  bytes_per_vol = fmri->width * fmri->height * fmri->depth * fmri->bytes_per_vox;
  for (c=0; c < fmri->width; c++) {
    for (r=0; r < fmri->height; r++) {
      for (s=0; s < fmri->depth; s++) {
        sumval = 0;
        sumsqval = 0;
        p = fmri->chunk + (c + r*nvox_per_row + s*nvox_per_slice)*fmri->bytes_per_vox;
        for (f=0; f < fmri->nframes; f++) {
          if (fmri->ischunked) {
            switch (fmri->type) {
            case MRI_UCHAR:
              val = (double)(*((char *)p));
              break;
            case MRI_SHORT:
              val = (double)(*((short*)p));
              break;
            case MRI_INT:
              val = (double)(*((int  *)p));
              break;
            case MRI_LONG:
              val = (double)(*((long *)p));
              break;
            case MRI_FLOAT:
              val = (double)(*((float*)p));
              break;
            }
          } else val = MRIgetVoxVal(fmri, c, r, s, f);
          //printf("%d  %lu   %g  %g\n",f,(long int)p,val,MRIgetVoxVal(fmri,c,r,s,f));
          sumsqval += (val*val);
          if (RmMean) sumval += val;
          p += bytes_per_vol;
        }
        MRIFseq_vox(var,c,r,s,0) = sumsqval/DOF;
        if (RmMean)
          MRIFseq_vox(var,c,r,s,0) -= ((sumval/DOF)*(sumval/DOF));
      }
    }
  }

  return(var);
}



/*---------------------------------------------------------------*/
/*
LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R)
{
  LTA *lta;
  MATRIX *vox2vox; // Targ->Mov
  MATRIX *Ttarg, *Tmov, *invTmov;

  Ttarg = MRIxfmCRS2XYZtkreg(targ);
  Tmov  = MRIxfmCRS2XYZtkreg(mov);
  invTmov = MatrixInverse(Tmov,NULL);

  // vox2vox = invTmov * R * Ttarg
  vox2vox = MatrixMultiply(invTmov,R,NULL);
  MatrixMultiply(vox2vox,Ttarg,vox2vox);

  lta = LTAalloc(1,NULL);
  lta->type = LINEAR_VOX_TO_VOX;
  lta->xforms[0].type = LINEAR_VOX_TO_VOX;
  getVolGeom(targ,&lta->xforms[0].src);
  getVolGeom(mov,&lta->xforms[0].dst);
  lta->xforms[0].m_L = MatrixCopy(vox2vox,NULL);

  MatrixFree(&Ttarg);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&vox2vox);

  return(lta);
}
*/

/*---------------------------------------------------------------*/
MRI *MRIsetSliceNo(MRI *mri, MRI *out) {
  int c, r, s, f, n, ncols, nrows, nslices,nframes;
  void   *pmri=NULL, *pout=NULL;
  int sz, szout;

  ncols   = mri->width;
  nrows   = mri->height;
  nslices = mri->depth;
  nframes = mri->nframes;

  if (out==NULL) {
    out = MRIallocSequence(ncols, nrows, nslices, mri->type, nframes);
    if (out==NULL) {
      printf("ERROR: MRIsetSliceNo: could not alloc output\n");
      return(NULL);
    }
    MRIcopyHeader(mri,out); // ordinarily would need to change nframes
  }
  if (out->width != ncols   || out->height != nrows ||
      out->depth != nslices || out->nframes != nframes) {
    printf("ERROR: MRIsetSliceNo: dimension mismatch\n");
    return(NULL);
  }

  // Number of bytes in the mri data types
  sz   = MRIsizeof(mri->type);
  szout = MRIsizeof(out->type);

  n = 0;
  for (f=0; f<nframes; f++) {
    for (s=0; s<nslices; s++) {
      for (r=0; r<nrows; r++) {
        // Pointers to the start of the column
        pmri  = (void *) mri->slices[n][r];
        pout  = (void *) out->slices[n][r];
        for (c=0; c<ncols; c++) {
          MRIdbl2ptr(s, pout, out->type);
          pmri += sz;
          pout  += szout;
        } // cols
      } // rows
      n++;
    } // slices
  } // frames

  return(out);
}

/*-------------------------------------------------------*/
//  if(1){
//  MRISmercator(surf);
//  sprintf(tmpstr,"%s/%s/surf/%s.mercator",SUBJECTS_DIR,subject,hemi);
//  MRISwrite(surf,tmpstr);
//  exit(0);
//  }
int MRISmercator(MRIS *surf)
{
  int vtxno;
  VERTEX *vtx;
  double x,y,z;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  double dthresh;
  dthresh = 5;

  xmin=ymin=zmin =  1e10;
  xmax=ymax=zmax = -1e10;
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    vtx = &(surf->vertices[vtxno]);
    x = vtx->x;
    y = vtx->y;
    z = vtx->z;
    vtx->x = 0;
    vtx->y = (100*(atan2(y,-x)))/M_PI;
    // z stays the same

    x = vtx->x;
    y = vtx->y;
    z = vtx->z;
    if(xmin > x) xmin = x;
    if(xmax < x) xmax = x;
    if(ymin > y) ymin = y;
    if(ymax < y) ymax = y;
    if(zmin > z) zmin = z;
    if(zmax < z) zmax = z;
  }

  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    vtx = &(surf->vertices[vtxno]);
    if(0){
      if(fabs(vtx->x - xmax) < dthresh) vtx->y = -100;
      if(fabs(vtx->x - xmin) < dthresh) vtx->y = -100;
      if(fabs(vtx->z - zmax) < dthresh) vtx->y = -100;
      if(fabs(vtx->z - zmin) < dthresh) vtx->y = -100;
    } else {
      if(fabs(vtx->y - ymax) < dthresh) vtx->x = -100;
      if(fabs(vtx->y - ymin) < dthresh) vtx->x = -100;
      if(fabs(vtx->z - zmax) < dthresh) vtx->x = -100;
      if(fabs(vtx->z - zmin) < dthresh) vtx->x = -100;
    }
  }

  return(0);
}

MRIS *MRISaverageSurfaces(int nsurfaces, MRIS **surfs)
{
  MRIS *surf, *avgsurf;
  int n,k;
  float average_surface_area;

  surf = surfs[0];
  avgsurf = MRISalloc(surf->nvertices, surf->nfaces) ;

  // Make sure xyz is 0, copy faces and vertices
  for(k=0; k < surf->nvertices; k++){
    avgsurf->vertices[k].x = 0.0;
    avgsurf->vertices[k].y = 0.0;
    avgsurf->vertices[k].z = 0.0;
    avgsurf->vertices[k].num = surf->vertices[k].num;
    avgsurf->vertices[k].f = (int*) calloc(surf->vertices[k].num,sizeof(int));
    avgsurf->vertices[k].n = (uchar*) calloc(surf->vertices[k].num,sizeof(uchar));
    for(n=0; n < surf->vertices[k].num; n++){
      avgsurf->vertices[k].f[n] = surf->vertices[k].f[n];
      avgsurf->vertices[k].n[n] = surf->vertices[k].n[n];
    }
    avgsurf->vertices[k].vnum = surf->vertices[k].vnum;
    avgsurf->vertices[k].v = (int*) calloc(surf->vertices[k].vnum,sizeof(int));
    avgsurf->vertices[k].dist = (float*) calloc(surf->vertices[k].vnum,sizeof(float));
    for(n=0; n < surf->vertices[k].vnum; n++)
      avgsurf->vertices[k].v[n] = surf->vertices[k].v[n];
  }
  for(k=0; k < surf->nfaces; k++) {
    for(n=0; n < VERTICES_PER_FACE; n++)
      avgsurf->faces[k].v[n] = surf->faces[k].v[n];
  }

  // Loop thru all surfaces, sume xyz 
  for(n=0; n < nsurfaces; n++){
    printf("%2d  surface area %g\n",n,surf->total_area);
    surf = surfs[0];
    if(surf->nvertices != surfs[0]->nvertices){
      printf("ERROR: MRISaverageSurfaces(): dimension mismatch surface %d\n",n);
      printf(" number of vertices %d vs %d\n",surf->nvertices,surfs[0]->nvertices);
      return(NULL);
    }
    for(k=0; k < surf->nvertices; k++){
      avgsurf->vertices[k].x += surf->vertices[k].x;
      avgsurf->vertices[k].y += surf->vertices[k].y;
      avgsurf->vertices[k].z += surf->vertices[k].z;
    }
    average_surface_area += surf->total_area ;
  }

  average_surface_area /= nsurfaces;
  printf("average surface area %g\n",surf->total_area);

  // Now divide by number of surfaces
  for(k=0; k < surf->nvertices; k++){
    avgsurf->vertices[k].x /= nsurfaces;
    avgsurf->vertices[k].y /= nsurfaces;
    avgsurf->vertices[k].z /= nsurfaces;
  }

  MRIScomputeMetricProperties(avgsurf);
  printf("avg  surface area %g\n",avgsurf->total_area);
  return(avgsurf);
}
/*!
  \fn MATRIX *MatrixLoadFSL(char *fname)
  \brief Loads in an FSL matrix. Not a registration matrix but one
  stored in their design.con or design.mat files used in the GLM 
  estimation. The format is pretty simple, it just has a /Matrix 
  keyword followed by the data.
*/
MATRIX *MatrixLoadFSL(char *fname)
{
  FILE *fp;
  MATRIX *M0,*M;
  char line[2000], *s, tag[2000];
  int hit,Nc,Nc0,c,nrows,r;

  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",fname);
    return(NULL);
  }

  // Get past the "/Matrix" string
  hit = 0;
  while(1){
    s = fgets(line, 1999, fp);
    if(s == NULL){
      printf("ERROR: %s is not formatted correctly\n",fname);
      fclose(fp);
      return(NULL);
    }
    sscanf(line,"%s",tag);
    if(! strcmp(tag,"/Matrix")){
      hit = 1;
      break;
    }
  }
  if(!hit){
    printf("ERROR: %s is not formatted correctly\n",fname);
    fclose(fp);
    return(NULL);
  }

  // Now read in each line
  M = NULL;
  Nc0 = 0; Nc=0;
  nrows = 0;
  while(1){
    if(nrows > 5000){
      printf("ERROR: %s, nrows exceeds 5000\n",fname);
      fclose(fp);
      MatrixFree(&M0);
      return(NULL);
    }
    // read in line
    s = fgets(line, 1999, fp);
    if(s == NULL) break;
    // Count number of columns
    Nc = gdfCountItemsInString(line);
    if(Nc == 0){
      printf("ERROR: %s is not formatted correctly\n",fname);
      fclose(fp);
      return(NULL);
    }
    // If first pass, alloc matrix, etc
    if(nrows == 0){
      Nc0 = Nc;
      //printf("%d cols\n",Nc);
      M0 = MatrixAlloc(5000,Nc,MATRIX_REAL);
      if(M0==NULL){
	printf("ERROR: could not alloc %d cols\n",Nc);
	fclose(fp);
	return(NULL);
      }
    }
    // Make sure this row is conistent with previous rows
    if(Nc0 != Nc){
      printf("ERROR: %s is not formatted correctly\n",fname);
      fclose(fp);
      MatrixFree(&M0);
      return(NULL);
    }
    // Read in each colum for this row
    for(c=0; c < Nc0; c++){
      s = gdfGetNthItemFromString(line, c);
      sscanf(s,"%f",&M0->rptr[nrows+1][c+1]);
      free(s);
    }
    nrows ++;
  }
  fclose(fp);
  //printf("%d rows\n",nrows);

  if(nrows == 0){
    printf("ERROR: %s is not formatted correctly\n",fname);
    return(NULL);
  }

  // Pack data into a matrix of the correct size
  M = MatrixAlloc(nrows,Nc,MATRIX_REAL);
  for(r=0; r < nrows; r++)
    for(c=0; c < Nc0; c++)
      M->rptr[r+1][c+1] = M0->rptr[r+1][c+1];
  MatrixFree(&M0);

  MatrixPrint(stdout,M);
  
  return(M);
}
/* ---------------------------------------------------------------------- */

int find_path ( int* vert_vno, int num_vno, int max_path_length,
                int* path, int* path_length, MRIS *mris ) {
  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char* check;
  float* dist;
  int* pred;
  char done;
  VERTEX* v;
  VERTEX* u;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  int num_path = 0;
  int num_checked;
  float vu_x, vu_y, vu_z;
  int flag2d = 0; // for flattend surface?

  dist = (float*) calloc (mris->nvertices, sizeof(float));
  pred = (int*) calloc (mris->nvertices, sizeof(int));
  check = (char*) calloc (mris->nvertices, sizeof(char));
  num_path = 0;
  num_checked = 0;
  (*path_length) = 0;

  for (cur_vert_vno = 0; cur_vert_vno < num_vno-1; cur_vert_vno++) {
    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++) {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno+1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices ||
        dest_vno < 0 || dest_vno >= mris->nvertices)
      continue;

    if (src_vno == dest_vno)
      continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done) {

      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist) {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      v = &(mris->vertices[closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno) {
        done = TRUE;
      } else {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < v->vnum; neighbor++) {
          neighbor_vno = v->v[neighbor];
          u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          vu_x = u->x - v->x;
          vu_y = u->y - v->y;
	  if (flag2d)	    vu_z = 0;
	  else     	    vu_z = u->z - v->z;

          /* recalc the weight. */
	  if (flag2d)
	    dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
			   ((v->y - u->y) * (v->y - u->y)));
	  else
	    dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
			   ((v->y - u->y) * (v->y - u->y)) +
			   ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno]) {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0) {
        printf (".");
        fflush (stdout);
      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno &&
           (*path_length) < max_path_length ) {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }
  printf (" done\n");
  fflush (stdout);

  free (dist);
  free (pred);
  free (check);

  return (ERROR_NONE);
}

int mygd(char *fname)
{
  FILE *fp;
  char *match = "DiffusionGradientDirection\0";
  char *check;
  int matchlen, n, m;
  int c, hit;
  float f;
  char buf[4];
  
  fp = fopen(fname,"r");
  
  matchlen = strlen(match);
  check = calloc(sizeof(char),matchlen);

  n = -1;
  while(fp){
    n++;
    c = fgetc(fp);
    if(n < matchlen){
      check[n] = c;
      continue;
    }
    for(m=0; m < matchlen-1; m++) check[m] = check[m+1];
    check[m] = c;

    hit = 1;
    for(m=0; m < matchlen; m++){
      if(check[m] != match[m]) hit = 0;
    }
    if(hit){
      for(m = 0; m < 10; m++){
	fread(buf,sizeof(float),1,fp);
	byteswapbuffloat(buf,1);
	f = (float)buf[0];
	//printf("%s %7.4f %7.4f %7.4f \n",fname,f[0],f[1],f[2]);
	printf("%d %7.4f\n",m,f);
      }
      return(0);
    }
  }
  return(1);
}

MRI *MRIaddB(MRI *mri1, MRI *mri2, MRI *mriadd)
{
  int c,r,s,f;
  double v1,v2,vadd;

  if(mriadd == NULL){
    mriadd = MRIallocSequence(mri1->width, mri1->height,mri1->depth,
			      MRI_FLOAT, mri1->nframes);
    MRIcopyHeader(mri1, mriadd);
  }

  for(c=0; c < mri1->width; c++){
    for(r=0; r < mri1->height; r++){
      for(s=0; s < mri1->depth; s++){
	for(f=0; f < mri1->nframes; f++){
	  v1 = MRIgetVoxVal(mri1,c,r,s,f);
	  v2 = MRIgetVoxVal(mri2,c,r,s,f);
	  vadd = v1+v2;
	  MRIsetVoxVal(mriadd,c,r,s,f,vadd);
	}
      }
    }
  }
  return(mriadd);
}

int RunSurfMCZ(MRIS *surf, MRI *mask, int nRepetitions,
	       double *ThreshList, int nThreshList,
	       double *FWHMList, int nFWHMList,
	       int SynthSeed, char *csdbase)
{
  RFS *rfs;
  int *nSmoothsList, nSmoothsPrev, nSmoothsDelta, nthRep;
  MRI *z, *zabs, *sig=NULL, *p=NULL;
  int SignList[3] = {-1,0,1}, nthSign, nthFWHM, nthThresh;
  CSD *csdList[22][20][3], *csd;
  double sigmax, zmax, threshadj, csize;
  int nClusters, cmax,rmax,smax;
  SURFCLUSTERSUM *SurfClustList;

  rfs = RFspecInit(SynthSeed,NULL);
  rfs->name = strcpyalloc("gaussian");
  rfs->params[0] = 0;
  rfs->params[1] = 1;

  nSmoothsList = (int *) calloc(sizeof(int),nFWHMList);
  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
    nSmoothsList[nthFWHM] = MRISfwhm2niters(FWHMList[nthFWHM], surf);
    printf("%2d %5.1f  %4d\n",nFWHMList,FWHMList[nthFWHM],nSmoothsList[nthFWHM]);
  }

  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
    for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
      for(nthSign = 0; nthSign < 3; nthSign++){
	csd = CSDalloc();
	sprintf(csd->simtype,"%s","null-z");
	sprintf(csd->anattype,"%s","surface");
	sprintf(csd->subject,"%s","??");
	sprintf(csd->hemi,"%s","??");
	sprintf(csd->contrast,"%s","N/A");
	csd->seed = SynthSeed;
	csd->nreps = nRepetitions;
	csd->thresh = ThreshList[nthThresh];
	csd->threshsign = SignList[nthSign];
	csd->nullfwhm = FWHMList[nthFWHM];
	csd->varfwhm = -1;
	csd->searchspace = -1;
	CSDallocData(csd);
	csdList[nthFWHM][nthThresh][nthSign] = csd;
      }
    }
  }

  z = MRIcloneBySpace(mask,MRI_FLOAT,1);

  for(nthRep = 0; nthRep < nRepetitions; nthRep++){
    RFsynth(z,rfs,mask); // synth unsmoothed z
    nSmoothsPrev = 0;
    for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
      nSmoothsDelta = nSmoothsList[nthFWHM] - nSmoothsPrev;
      nSmoothsPrev = nSmoothsList[nthFWHM];
      MRISsmoothMRI(surf, z, nSmoothsDelta, mask, z); // smooth z
      RFrescale(z,rfs,mask,z);
      zabs = MRIabs(z,NULL);
      // Slightly tortured way to get the right p-values because
      //   RFstat2P() computes one-sided, but I handle sidedness
      //   during thresholding.
      // First, use zabs to get a two-sided pval bet 0 and 0.5
      zabs = MRIabs(z,zabs);
      p = RFstat2P(zabs,rfs,mask,0,p);
      // Next, mult pvals by 2 to get two-sided bet 0 and 1
      MRIscalarMul(p,p,2.0);
      sig = MRIlog10(p,NULL,sig,1); // sig = -log10(p)
      for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
	for(nthSign = 0; nthSign < 3; nthSign++){
	  csd = csdList[nthFWHM][nthThresh][nthSign];
	  // If test is not ABS then apply the sign
	  if(csd->threshsign != 0) MRIsetSign(sig,z,0);
	  sigmax = MRIframeMax(sig,0,mask,csd->threshsign,
			       &cmax,&rmax,&smax);
	  zmax = MRIgetVoxVal(z,cmax,rmax,smax,0);
	  if(csd->threshsign == 0) zmax = fabs(zmax);
	  if(mask) MRImask(sig,mask,sig,0.0,0.0);
	  MRIScopyMRI(surf, sig, 0, "val");
	  if(csd->threshsign == 0) threshadj = csd->thresh;
	  else threshadj = csd->thresh - log10(2.0); // one-sided test
	  SurfClustList = sclustMapSurfClusters(surf,threshadj,-1,csd->threshsign,
						0,&nClusters,NULL);
	  csize = sclustMaxClusterArea(SurfClustList, nClusters);
	  csd->nClusters[nthRep] = nClusters;
	  csd->MaxClusterSize[nthRep] = csize;
	  csd->MaxSig[nthRep] = sigmax;
	  csd->MaxStat[nthRep] = zmax;
	  
	} // Sign
      } // Thresh
    } // FWHM
  } // Iteration

  return(0);
}

