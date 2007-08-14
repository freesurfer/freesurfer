/**
 * @file  dngtester.c
 * @brief dougs super special test code
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/08/14 03:40:05 $
 *    $Revision: 1.42 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
MRIS *surf;
char *hemi;

MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var);
void printrgb(void);
MRI *GetMyMask(MRIS *surf);
MRI *GetMyMask2(MRIS *surf);

MRI *MRISdilateMask(MRIS *surf, MRI *mask, int annotidmask, int niters, int newid);
COLOR_TABLE *CTABaddEntry(COLOR_TABLE *ctold, char *name);
int MRISmercator(MRIS *surf);
IMAGE *I;
MHT *lhwhite_hash;
MRIS *lhwhite;
MRIS *surfs[100];
int *XNbrVtxNo, nXNbrs;
double *XNbrDotProd;

MRIS *MRISaverageSurfaces(int nsurfaces, MRIS **surfs);

/*----------------------------------------*/
int main(int argc, char **argv) {
  int err,area32p,area32v,superiorfrontal,medialorbitofrontal;
  int rostralanteriorcingulate, rostralmiddlefrontal;
  int index,k;
  int *nunits;
  char *parcnames[10];
  COLOR_TABLE *ct ;
  LABEL *lcortex;
  VERTEX *vtx;
  float dlhw,DotProdThresh;
  int  lhwvtx;

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

  subject = argv[1];
  hemi = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.sphere",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

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
  //---------------------------------------------------------
#if 0
  printf("nvertices %d\n",surf->nvertices);
  printf("patch %d\n",surf->patch);
  printf("radius %g\n",dlhw);
  printf("patch %d\n",surf->patch);
  printf("rip %d\n",vtx->ripflag);
#endif



  printf("Building hash of lh white\n");
  lhwhite_hash = MHTfillVertexTableRes(lhwhite, NULL,CURRENT_VERTICES,16);

  vtx = &(lhwhite->vertices[1000]);

  // Find closest when not ripped
  vtx->ripflag = 0;
  lhwvtx = MHTfindClosestVertexNo(lhwhite_hash,lhwhite,vtx,&dlhw);
  printf("%d\n",lhwvtx);

  // Find closest when ripped
  vtx->ripflag = 1;
  lhwvtx = MHTfindClosestVertexNo(lhwhite_hash,lhwhite,vtx,&dlhw);
  printf("%d\n",lhwvtx);

  // Brute-force find closest when ripped
  lhwvtx = MRISfindClosestVertex(lhwhite,vtx->x,vtx->y,vtx->z,&dlhw);
  printf("%d\n",lhwvtx);



  return(1);



  //R = MatrixDRand48(20,20,NULL);
  R = MatrixIdentity(100,NULL);
  MatrixPrint(stdout,R);
  I = ImageFromMatrix(R, NULL);
  ImageWrite(I, "./myimage.jpg") ;
  return(1);


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

#define MIN_NONCORTEX_VERTICES 10

  lcortex = MRIScortexLabel(surf, mri, MIN_NONCORTEX_VERTICES) ;
  sprintf(tmpstr,"%s/%s/label/%s.%s.label",SUBJECTS_DIR, subject,hemi,"cortex");
  printf("writing cortex label to %s...\n", tmpstr) ;
  LabelWrite(lcortex, tmpstr) ;

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


  return(1);

  if(1){
  MRISmercator(surf);
  sprintf(tmpstr,"%s/%s/surf/%s.mercator",SUBJECTS_DIR,subject,hemi);
  MRISwrite(surf,tmpstr);
  exit(0);
  }


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

  mri = GetMyMask2(surf);
  MRISdilateMask(surf, mri, medialorbitofrontal, 12, area32v);

  mri = GetMyMask(surf);
  MRISdilateMask(surf, mri, superiorfrontal, 12, area32p);


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

  sprintf(tmpstr,"%s.fbirn.annot",hemi);
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
/*-------------------------------------------------------------
  MRI *GetMyMask(MRIS *surf) - creates a mask at the intersection
  of SFG and CAcing, PAcing, and RAcing.
  -------------------------------------------------------------*/
MRI *GetMyMask(MRIS *surf)
{
  int vtxno, annot, annotid, nnbrs, nbrvtx, nthnbr, nbrannotid; 
  VERTEX *vtx;
  int superiorfrontal, posteriorcingulate;
  int caudalanteriorcingulate, rostralanteriorcingulate;
  MRI *mri;

  superiorfrontal = 28;
  posteriorcingulate = 23;
  caudalanteriorcingulate = 2;
  rostralanteriorcingulate = 26;

  mri = MRIalloc(surf->nvertices, 1, 1, MRI_INT) ;

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    MRIsetVoxVal(mri,vtxno,0,0,0, 0);

    vtx = &(surf->vertices[vtxno]);
    annot = surf->vertices[vtxno].annotation;
    CTABfindAnnotation(surf->ct, annot, &annotid);

    if(annotid != superiorfrontal && 
       annotid != posteriorcingulate &&
       annotid != caudalanteriorcingulate && 
       annotid != rostralanteriorcingulate) continue;

    nnbrs = surf->vertices[vtxno].vnum;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++){
      nbrvtx = surf->vertices[vtxno].v[nthnbr];
      if (surf->vertices[nbrvtx].ripflag) continue; //skip ripped vtxs
      annot = surf->vertices[nbrvtx].annotation;
      CTABfindAnnotation(surf->ct, annot, &nbrannotid);
      if(annotid == superiorfrontal && 
	 (nbrannotid == posteriorcingulate ||
	  nbrannotid == caudalanteriorcingulate ||
	  nbrannotid == rostralanteriorcingulate) ){
	MRIsetVoxVal(mri,vtxno,0,0,0,1);
      }
    }
  }

  //MRIwrite(mri,"mymask.mgh");

  return(mri);
}

/*-------------------------------------------------------------
  MRI *GetMyMask(MRIS *surf) - creates a mask at the intersection
  of SFG and CAcing, PAcing, and RAcing.
  -------------------------------------------------------------*/
MRI *GetMyMask2(MRIS *surf)
{
  int vtxno, annot, annotid, nnbrs, nbrvtx, nthnbr, nbrannotid; 
  VERTEX *vtx;
  int medialorbitofrontal, rostralanteriorcingulate;
  MRI *mri;

  CTABfindName(surf->ct, "medialorbitofrontal", &medialorbitofrontal);
  CTABfindName(surf->ct, "rostralanteriorcingulate", &rostralanteriorcingulate);

  printf("medialorbitofrontal %d\n",medialorbitofrontal);
  printf("rostralanteriorcingulate %d\n",rostralanteriorcingulate);

  mri = MRIalloc(surf->nvertices, 1, 1, MRI_INT) ;

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    MRIsetVoxVal(mri,vtxno,0,0,0, 0);

    vtx = &(surf->vertices[vtxno]);
    annot = surf->vertices[vtxno].annotation;
    CTABfindAnnotation(surf->ct, annot, &annotid);

    if(annotid != medialorbitofrontal) continue;

    nnbrs = surf->vertices[vtxno].vnum;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++){

      nbrvtx = surf->vertices[vtxno].v[nthnbr];
      if (surf->vertices[nbrvtx].ripflag) continue; //skip ripped vtxs

      annot = surf->vertices[nbrvtx].annotation;
      CTABfindAnnotation(surf->ct, annot, &nbrannotid);
      if(nbrannotid == rostralanteriorcingulate)
	MRIsetVoxVal(mri,vtxno,0,0,0,1);

    }
  }

  //MRIwrite(mri,"area32v.mgh");

  return(mri);
}
/*---------------------------------------------------------------------*/
MRI *MRISdilateMask(MRIS *surf, MRI *mask, int annotidmask, int niters, int newid)
{
  int vtxno, annot, annotid, nnbrs, nbrvtxno, nthnbr, nthiter,new_annot ;
  VERTEX *vtx;
  MRI *mri1, *mri2;
  float val;

  mri1 = MRIcopy(mask,NULL);
  mri2 = MRIcopy(mask,NULL);

  for(nthiter = 0; nthiter < niters; nthiter++){
    printf("iter %d\n",nthiter);

    for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
      vtx = &(surf->vertices[vtxno]);

      if(annotidmask > -1){
	annot = surf->vertices[vtxno].annotation;
	CTABfindAnnotation(surf->ct, annot, &annotid);
	if(annotid != annotidmask){
	  MRIsetVoxVal(mri2,vtxno,0,0,0,0);
	  continue;
	}
      }

      val = MRIgetVoxVal(mri1,vtxno,0,0,0);

      if(val){
	MRIsetVoxVal(mri2,vtxno,0,0,0,1);
	continue;
      }
      // If it gets here, the vtx is in the annot and has not been set 
      nnbrs = surf->vertices[vtxno].vnum;
      for (nthnbr = 0; nthnbr < nnbrs; nthnbr++){
	nbrvtxno = surf->vertices[vtxno].v[nthnbr];
	if (surf->vertices[nbrvtxno].ripflag) continue; //skip ripped vtxs
	val = MRIgetVoxVal(mri1,nbrvtxno,0,0,0);
	if(val){
	  MRIsetVoxVal(mri2,vtxno,0,0,0,1);
	  continue;
	}
      }
    }
    MRIcopy(mri2,mri1);
  
  }

  MRIfree(&mri2);

  //MRIwrite(mri1,"mymask2.mgh");
  
  CTABannotationAtIndex(surf->ct, newid, &new_annot) ;

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    if(MRIgetVoxVal(mri1,vtxno,0,0,0))
      surf->vertices[vtxno].annotation = new_annot;
  }

  return(mri1);
}

/*-------------------------------------------------------*/
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
