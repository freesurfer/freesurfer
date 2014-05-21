/**
 * @file  gtm.c
 * @brief Routines to create and analyze the Geometric Transfer Matrix (GTM)
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/05/21 17:56:15 $
 *    $Revision: 1.9 $
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
#include <sys/utsname.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "gtm.h"
#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "timer.h"
#include "fmriutils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "cma.h"
#include "matfile.h"
#include "region.h"
#include "mrimorph.h"
#include "numerics.h"
#include "resample.h"

/*------------------------------------------------------------------------------------*/
int GTMSEGprint(GTMSEG *gtmseg, FILE *fp)
{
  fprintf(fp,"subject %s\n",gtmseg->subject);
  fprintf(fp,"USF %d\n",gtmseg->USF);
  fprintf(fp,"OutputUSF %d\n",gtmseg->OutputUSF);
  fprintf(fp,"apasfile %s\n",gtmseg->apasfile);
  if(gtmseg->wmannotfile != NULL){
    fprintf(fp,"wmannotfile %s\n",gtmseg->wmannotfile);
    fprintf(fp,"wmlhbase %d\n",gtmseg->wmlhbase);
    fprintf(fp,"wmrhbase %d\n",gtmseg->wmrhbase);
  }
  else fprintf(fp,"wmannotfile NULL\n");
  fprintf(fp,"ctxannotfile %s\n",gtmseg->ctxannotfile);
  fprintf(fp,"ctxlhbase %d\n",gtmseg->ctxlhbase);
  fprintf(fp,"ctxrhbase %d\n",gtmseg->ctxrhbase);
  fprintf(fp,"SubSegWM %3d\n",gtmseg->SubSegWM);
  fprintf(fp,"KeepHypo %3d\n",gtmseg->KeepHypo);
  fprintf(fp,"KeepCC %3d\n",gtmseg->KeepCC);
  fprintf(fp,"dmax %f\n",gtmseg->dmax);
  fprintf(fp,"nlist %3d\n",gtmseg->nlist);
  fprintf(fp,"lhmin %5d\n",gtmseg->lhmin);
  fprintf(fp,"lhmax %5d\n",gtmseg->lhmax);
  fprintf(fp,"rhmin %5d\n",gtmseg->rhmin);
  fprintf(fp,"rhmax %5d\n",gtmseg->rhmax);
  fflush(fp);
  return(0);
}

/*------------------------------------------------------------------------------------*/
int MRIgtmSeg(GTMSEG *gtmseg)
{
  int err,*segidlist,nsegs,n;
  char *SUBJECTS_DIR, tmpstr[5000];
  MRI *apas, *aseg, *hrseg, *ctxseg;
  struct timeb timer;
  TimerStart(&timer);

  printf("Starting MRIgtmSeg() USF=%d\n",gtmseg->USF);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->apasfile);
  printf("Loading %s\n",tmpstr);
  apas = MRIread(tmpstr);
  if(apas==NULL) return(1);

  gtmseg->anat = MRIcopyHeader(apas,NULL); // keep anat header around

  aseg = MRIunsegmentCortex(apas, gtmseg->lhmin, gtmseg->lhmax, gtmseg->rhmin, gtmseg->rhmax, NULL); // aseg+head
  if(aseg == NULL) return(1);
  MRIfree(&apas);

  printf("Loading surfaces ");
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,gtmseg->subject);
  gtmseg->lhw = MRISread(tmpstr); if(gtmseg->lhw==NULL) return(1);

  sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,gtmseg->subject);
  gtmseg->lhp = MRISread(tmpstr);  if(gtmseg->lhp==NULL) return(1);

  sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,gtmseg->subject);
  gtmseg->rhw = MRISread(tmpstr);  if(gtmseg->rhw==NULL) return(1);

  sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,gtmseg->subject);
  gtmseg->rhp = MRISread(tmpstr);  if(gtmseg->rhp==NULL) return(1);

  printf("Loading annotations ");
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  if(gtmseg->wmannotfile != NULL){
    sprintf(tmpstr,"%s/%s/label/lh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->wmannotfile);
    err = MRISreadAnnotation(gtmseg->lhw,tmpstr);  
    if(err) {
      printf("Try running mri_annotation2label --s %s --hemi lh --lobesStrict lobefilename\n",gtmseg->subject);
      return(1);
    }
    sprintf(tmpstr,"%s/%s/label/rh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->wmannotfile);
    err = MRISreadAnnotation(gtmseg->rhw,tmpstr);
    if(err) {
      printf("Try running mri_annotation2label --s %s --hemi rh --lobesStrict lobefilename\n",gtmseg->subject);
      return(1);
    }
    MRISripUnknown(gtmseg->lhw);
    MRISripUnknown(gtmseg->rhw);
    gtmseg->lhw->ct->idbase = gtmseg->wmlhbase; 
    gtmseg->rhw->ct->idbase = gtmseg->wmrhbase; 
  }
  else printf("Not segmenting WM\n");

  sprintf(tmpstr,"%s/%s/label/lh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->ctxannotfile);
  err = MRISreadAnnotation(gtmseg->lhp,tmpstr);  if(err) return(1);
  sprintf(tmpstr,"%s/%s/label/rh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->ctxannotfile);
  err = MRISreadAnnotation(gtmseg->rhp,tmpstr);  if(err) return(1);

  gtmseg->lhp->ct->idbase = gtmseg->ctxlhbase;
  gtmseg->rhp->ct->idbase = gtmseg->ctxrhbase;

  MRISsetPialUnknownToWhite(gtmseg->lhw, gtmseg->lhp);
  MRISsetPialUnknownToWhite(gtmseg->rhw, gtmseg->rhp);
  MRISripUnknown(gtmseg->lhp);
  MRISripUnknown(gtmseg->rhp);

  // relabel WM_hypointensities as {Left,Right}_WM_hypointensities as 
  printf(" Relabeling any unlateralized hypointenities as lateralized hypointenities\n");
  MRIrelabelHypoHemi(aseg, gtmseg->lhw, gtmseg->rhw, NULL, aseg);

  int nlist, list[100];
  nlist = 0;
  if(!gtmseg->KeepCC){
    printf(" Relabeling CC as WM\n");
    list[nlist] = 192; nlist++;
    list[nlist] = 251; nlist++;
    list[nlist] = 252; nlist++;
    list[nlist] = 253; nlist++;
    list[nlist] = 254; nlist++;
    list[nlist] = 255; nlist++;
  }
  if(!gtmseg->KeepHypo){
    printf(" Relabeling any hypointensities as WM\n");
    list[nlist] = 77; nlist++;
    list[nlist] = 78; nlist++;
    list[nlist] = 79; nlist++;
  }
  if(nlist > 0) MRIunsegmentWM(aseg, gtmseg->lhw, gtmseg->rhw, list, nlist, NULL, aseg);

  // Upsample the segmentation
  printf("Upsampling segmentation USF = %d",gtmseg->USF);fflush(stdout);
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  hrseg = MRIhiresSeg(aseg, gtmseg->lhw, gtmseg->lhp, gtmseg->rhw, gtmseg->rhp, gtmseg->USF, &gtmseg->anat2seg);
  if(hrseg == NULL) return(1);
  MRIfree(&aseg);

  // Label cortex (like aparc+aseg)
  printf("Beginning cortical segmentation using %s",gtmseg->ctxannotfile);fflush(stdout);
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  ctxseg = MRIannot2CorticalSeg(hrseg, gtmseg->lhw, gtmseg->lhp, gtmseg->rhw, gtmseg->rhp, NULL, NULL);
  MRIfree(&hrseg);

  // Label wm (like wmaparc)
  if(gtmseg->wmannotfile != NULL){
    printf("Beginning WM segmentation using %s",gtmseg->wmannotfile);fflush(stdout);
    printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
    ctxseg = MRIannot2CerebralWMSeg(ctxseg, gtmseg->lhw, gtmseg->rhw, gtmseg->dmax, NULL, ctxseg);
  } else     printf("Not subsegmenting WM\n");

  gtmseg->seg = MRIreplaceList(ctxseg, gtmseg->srclist, gtmseg->targlist, gtmseg->nlist, NULL);
  if(gtmseg == NULL) return(1);
  MRIfree(&ctxseg);

  segidlist = MRIsegIdListNot0(gtmseg->seg, &nsegs, 0);
  printf("Found %d segs in the final list\n",nsegs);
  for(n=0; n < nsegs; n++){
    if(segidlist[n] == Left_Cerebral_Cortex){
      printf("ERROR: MRIgtmSeg() found left cortical label\n");
      err = 1;
    }
    if(segidlist[n] == Right_Cerebral_Cortex){
      printf("ERROR: MRIgtmSeg() found right cortical label\n");
      err = 1;
    }
    if(gtmseg->SubSegWM){
      if(segidlist[n] == Left_Cerebral_White_Matter){
	printf("ERROR: MRIgtmSeg() found left cerebral WM label\n");
	err = 1;
      }
      if(segidlist[n] == Right_Cerebral_White_Matter){
	printf("ERROR: MRIgtmSeg() found right cerebral WM label\n");
	err = 1;
      }
    }
    if(segidlist[n]==WM_hypointensities){
      printf("ERROR: MRIgtmSeg() found unlateralized WM hypointensity label\n");
      err = 1;
    }
    if(err) return(1);
  }
  gtmseg->segidlist = segidlist;
  gtmseg->nsegs = nsegs;

  printf("MRIgtmSeg() done, t = %6.4f\n",TimerStop(&timer)/1000.0);
  fflush(stdout);
  return(0);
}

/*-----------------------------------------------------------------------------*/
/* 
\fn int GTMdefaultSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat)
\brief Creates a list of segids to replace with other seg ids. This is
mainly to merge very small segs (like temporal pole) with bigger segs
so that the GTM does not become too ill-conditioned. All ventricular
CSF segs are merged. CC subsegs are merged into a single label 192. It
is assumed that ReplaceThis and WithThat are arrays that have already
been allocated. It is also assumed that nReplace has been initialized.
The result is that items are added to the list.
 */
int GTMdefaultSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat)
{
  int nlist;

  nlist = *nReplace;
  ReplaceThis[nlist] = 1033; WithThat[nlist] = 1030; nlist++; // temppole=stg
  ReplaceThis[nlist] = 2033; WithThat[nlist] = 2030; nlist++; // temppole=stg
  ReplaceThis[nlist] = 1034; WithThat[nlist] = 1030; nlist++; // transtemp=stg
  ReplaceThis[nlist] = 2034; WithThat[nlist] = 1030; nlist++; // transtemp=stg
  ReplaceThis[nlist] = 1001; WithThat[nlist] = 1015; nlist++; // bankssts=mtg
  ReplaceThis[nlist] = 2001; WithThat[nlist] = 2015; nlist++; // bankssts=mtg
  ReplaceThis[nlist] = 1032; WithThat[nlist] = 1027; nlist++; // frontpole=rmf
  ReplaceThis[nlist] = 2032; WithThat[nlist] = 2027; nlist++; // frontpole=rmf
  //ReplaceThis[nlist] = 1016; WithThat[nlist] = 1006; nlist++; // parahip=entorhinal ?
  //ReplaceThis[nlist] = 2016; WithThat[nlist] = 2006; nlist++; // parahip=entorhinal ?

  // There should not be any cortex unknown after MRIannot2CorticalSeg()
  ReplaceThis[nlist] = 1000; WithThat[nlist] =    0; nlist++; // cortex unknown
  ReplaceThis[nlist] = 2000; WithThat[nlist] =    0; nlist++; // cortex unknown

  ReplaceThis[nlist] =   85; WithThat[nlist] =    0; nlist++; // optic chiasm

  // Merge ventricular CSF into one label
  ReplaceThis[nlist] =    4; WithThat[nlist] =   24; nlist++; // LLatVent
  ReplaceThis[nlist] =    5; WithThat[nlist] =   24; nlist++; // LInfLatVent
  ReplaceThis[nlist] =   14; WithThat[nlist] =   24; nlist++; // 3rd
  ReplaceThis[nlist] =   15; WithThat[nlist] =   24; nlist++; // 4th
  ReplaceThis[nlist] =   72; WithThat[nlist] =   24; nlist++; // 5th
  ReplaceThis[nlist] =   43; WithThat[nlist] =   24; nlist++; // RLatVent
  ReplaceThis[nlist] =   44; WithThat[nlist] =   24; nlist++; // RInfLatVent

  // Not sure about choriod plexi. Make part of CSF?
  ReplaceThis[nlist] =   31; WithThat[nlist] =   24; nlist++; // LChoroidP ?
  ReplaceThis[nlist] =   63; WithThat[nlist] =   24; nlist++; // RChoroidP ?

  // And these?
  ReplaceThis[nlist] =   30; WithThat[nlist] =   24; nlist++; // LVessel ?
  ReplaceThis[nlist] =   62; WithThat[nlist] =   24; nlist++; // RVessel ?
  ReplaceThis[nlist] =   80; WithThat[nlist] =   24; nlist++; // non-WM-hypo ?

  /* Merge multiple CC subsegments into one CC */
  ReplaceThis[nlist] =  251; WithThat[nlist] =  192; nlist++; 
  ReplaceThis[nlist] =  252; WithThat[nlist] =  192; nlist++; 
  ReplaceThis[nlist] =  253; WithThat[nlist] =  192; nlist++; 
  ReplaceThis[nlist] =  254; WithThat[nlist] =  192; nlist++; 
  ReplaceThis[nlist] =  255; WithThat[nlist] =  192; nlist++; 

  *nReplace += nlist;
  return(0);
}

/*
\fn COLOR_TABLE *GTMSEGctab(GTMSEG *gtmseg, COLOR_TABLE *ctSubCort)
\brief Creates a color table to represent the gtm segmentation as
created by MRIgtmSeg(). Cortical segmentations are determined from the
pial surf annot. If WM is subsegmented, then the subsegments are
determined from the white surf annot. Entries for subcortical
structures are determined from ctSubCort (which should have tissue
type set for all relevant structures).
*/
COLOR_TABLE *GTMSEGctab(GTMSEG *gtmseg, COLOR_TABLE *ctSubCort)
{
  int nsegs, *segidlist,*segidhit,n,m,hit,segid,err;
  COLOR_TABLE *ct, *ct0;
  CTE *cte,*cte0;

  segidlist = MRIsegIdListNot0(gtmseg->seg, &nsegs, 0);
  segidhit  = (int *)calloc(nsegs,sizeof(int));

  //FILE *fp = fopen("segidlist.dat","w");
  //for(m=0; m < nsegs; m++) fprintf(fp,"%4d\n",segidlist[m]);
  //fclose(fp);
  
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = segidlist[nsegs-1] + 1;
  ct->entries = (COLOR_TABLE_ENTRY**) calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY*));

  ct->ctabTissueType = CTABdeepCopy(ctSubCort->ctabTissueType);

  // Add an entry for unknown
  segid = 0;
  ct->entries[segid] = (CTE*) calloc(1, sizeof(COLOR_TABLE_ENTRY));
  cte = ct->entries[segid];
  cte->ri = 0; cte->gi = 0; cte->bi = 0; cte->ai = 255;
  cte->rf = 0; cte->gf = 0; cte->bf = 0; cte->af = 255;
  cte->TissueType = 0;
  sprintf(cte->name,"Unknown");

  fflush(stdout);
  ct0 = gtmseg->lhp->ct;
  // Go thru each entry in the left cortical annotation
  for(n=0; n < ct0->nentries; n++){
    cte0 = ct0->entries[n];
    if(cte0 == NULL) continue;
    // compute the segid
    segid = gtmseg->ctxlhbase + n;
    // check whether it is in the segidlist
    hit = 0;
    for(m=0; m < nsegs; m++) {
      if(segid == segidlist[m]) {
	hit = 1;
	break;
      }
    }
    if(hit==0){
      /* It does not have to be there. Eg, CC will be removed. Also
	 there are a bunch of entries that are non-null but not
	 actually represented in the annot (they have the name
	 clusterXX. */
      if(Gdiag_no > 1) printf("INFO: GTMSEGctab(): lhp n=%d, segid=%d, %s not in segidlist\n",n,segid,cte0->name);
      continue;
    }
    segidhit[m]++; // keep track of number of times this seg id represented
    ct->entries[segid] = (CTE*) calloc(1, sizeof(COLOR_TABLE_ENTRY));
    cte = ct->entries[segid];
    memcpy(cte,cte0,sizeof(CTE));
    sprintf(cte->name,"ctx-lh-%s",cte0->name); // new name reflects cortex and hemi
    cte->TissueType = 1;
  }

  // Do the same thing for the right hemi
  ct0 = gtmseg->rhp->ct;
  for(n=0; n < ct0->nentries; n++){
    cte0 = ct0->entries[n];
    if(cte0 == NULL) continue;
    segid = gtmseg->ctxrhbase + n;
    hit = 0;
    for(m=0; m < nsegs; m++) {
      if(segid == segidlist[m]) {
	hit = 1;
	break;
      }
    }
    if(hit==0){
      if(Gdiag_no > 1) printf("INFO: GTMSEGctab(): rhp n=%d, segid=%d, %s not in segidlist\n",n,segid,cte0->name);
      continue;
    }
    segidhit[m]++;
    ct->entries[segid] = (CTE*) calloc(1, sizeof(COLOR_TABLE_ENTRY));
    cte = ct->entries[segid];
    memcpy(cte,cte0,sizeof(CTE));
    sprintf(cte->name,"ctx-rh-%s",cte0->name);
    cte->TissueType = 1;
  }

  // Do the same thing if subsegmenting WM based on proximity to cortex
  if(gtmseg->SubSegWM){
    ct0 = gtmseg->lhw->ct;
    for(n=0; n < ct0->nentries; n++){
      cte0 = ct0->entries[n];
      if(cte0 == NULL) continue;
      segid = gtmseg->wmlhbase + n;
      hit = 0;
      for(m=0; m < nsegs; m++) {
	if(segid == segidlist[m]) {
	  hit = 1;
	  break;
	}
      }
      if(hit==0){
	if(Gdiag_no > 1) printf("INFO: GTMSEGctab(): lhw n=%d, segid=%d, %s not in segidlist\n",n,segid,cte0->name);
	continue;
      }
      segidhit[m]++;
      ct->entries[segid] = (CTE*) calloc(1, sizeof(COLOR_TABLE_ENTRY));
      cte = ct->entries[segid];
      memcpy(cte,cte0,sizeof(CTE));
      sprintf(cte->name,"wm-lh-%s",cte0->name);
      cte->TissueType = 3;
    }
    ct0 = gtmseg->rhw->ct;
    for(n=0; n < ct0->nentries; n++){
      cte0 = ct0->entries[n];
      if(cte0 == NULL) continue;
      segid = gtmseg->wmrhbase + n;
      hit = 0;
      for(m=0; m < nsegs; m++) {
	if(segid == segidlist[m]) {
	  hit = 1;
	  break;
	}
      }
      if(hit==0){
	if(Gdiag_no > 1) printf("INFO: GTMSEGctab(): lhp n=%d, segid=%d, %s not in segidlist\n",n,segid,cte0->name);
	continue;
      }
      segidhit[m]++;
      ct->entries[segid] = (CTE*) calloc(1, sizeof(COLOR_TABLE_ENTRY));
      cte = ct->entries[segid];
      memcpy(cte,cte0,sizeof(CTE));
      sprintf(cte->name,"wm-rh-%s",cte0->name);
      cte->TissueType = 3;
    }
  }

  // Add entries for subcortical regions
  ct0 = ctSubCort;
  for(m=0; m < nsegs; m++){
    if(segidhit[m]!=0) continue; // skip of already hit
    segid = segidlist[m];
    cte0 = ct0->entries[segid];
    if(cte0 == NULL) {
      printf("ERROR: cannot find match for subcortical segid %d\n",segid);
      return(NULL);
    }
    ct->entries[segid] = (CTE*) calloc(1, sizeof(COLOR_TABLE_ENTRY));
    cte = ct->entries[segid];
    memcpy(cte,cte0,sizeof(CTE));
    segidhit[m]++;
  }

  // Check the final ctab
  err = 0;
  for(m=0; m < nsegs; m++){
    if(segidhit[m] > 1) {
      printf("ERROR: segid %4d is represented multiple (%d) times \n",segidhit[m],segidlist[m]);
      err++;
    }
  }
  if(err) {
    printf("ERROR: found %d segmentations with multiple representations\n",err);
    return(NULL);
  }
  for(m=0; m < nsegs; m++){
    segid = segidlist[m];
    cte = ct->entries[segid];
    if(cte->TissueType == -1)
      printf("WARNING: segid %4d %s tissue type is not set\n",segid,cte->name);
  }

  free(segidlist);
  free(segidhit);

  return(ct);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax)
  \brief Computes the variance reduction factor (VRF) for each seg, the number
  of voxels for the seg, and betavar for each seg and frame. Also computes
  the max and min VRF.
*/
int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax)
{
  int n,nvox,segid,f;
  double vrf;

  *vrfmean = 0;
  *vrfmax = 0;
  *vrfmin = 0;
  for(n=0; n < gtm->iXtX->rows; n++){
    vrf = (double) 1.0/gtm->iXtX->rptr[n+1][n+1];
    if(n==0){
      *vrfmax = vrf;
      *vrfmin = vrf;
    }
    if(*vrfmax < vrf) *vrfmax = vrf;
    if(*vrfmin > vrf) *vrfmin = vrf;
    *vrfmean += vrf;
  }
  *vrfmean /= gtm->iXtX->rows;

  if(gtm->betavar) MatrixFree(&gtm->betavar);
  gtm->betavar = MatrixAlloc(gtm->iXtX->rows,gtm->beta->cols,MATRIX_REAL);

  if(gtm->vrf) MatrixFree(&gtm->vrf);
  gtm->vrf = MatrixAlloc(gtm->iXtX->rows,1,MATRIX_REAL);

  if(gtm->nvox) MatrixFree(&gtm->nvox);
  gtm->nvox = MatrixAlloc(gtm->iXtX->rows,1,MATRIX_REAL);

  for(n=0; n < gtm->iXtX->rows; n++){
    segid = gtm->segidlist[n];
    nvox = MRIcountMatches(gtm->gtmseg, segid, 0, gtm->mask);
    vrf = (double) 1.0/gtm->iXtX->rptr[n+1][n+1];
    gtm->vrf->rptr[n+1][1]  = vrf;
    gtm->nvox->rptr[n+1][1] = nvox;
    for(f=0; f < gtm->beta->cols; f++)
      gtm->betavar->rptr[n+1][f+1] = gtm->rvar->rptr[1][f+1]/vrf;
  }

  return(0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int WriteVRFStats(char *fname, GTM *gtm)
  \brief Creates the vrf.dat file in the output folder. This is a text file
  that reports several statistics including the variance reduction factor (VRF).
*/
int WriteVRFStats(char *fname, GTM *gtm)
{
  int n, segid, nvox;
  double vrf;
  FILE *fp;
  CTE *cte;
  CT *ttctab;

  ttctab = gtm->ctGTMSeg->ctabTissueType;

  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",fname);
    return(1);
  }

  for(n=0; n < gtm->iXtX->rows; n++){
    segid = gtm->segidlist[n];
    vrf = gtm->vrf->rptr[n+1][1];
    nvox = gtm->nvox->rptr[n+1][1];
    cte = gtm->ctGTMSeg->entries[segid];
    fprintf(fp,"%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
	    ttctab->entries[cte->TissueType]->name,nvox,vrf);
    //printf("%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
    //ttctab->entries[cte->TissueType]->name,nvox,vrf);
    if(gtm->beta)    fprintf(fp,"   %10.3f",gtm->beta->rptr[n+1][1]);
    if(gtm->segrvar) fprintf(fp,"   %10.4f",sqrt(gtm->segrvar->rptr[n+1][1]));
    fprintf(fp,"\n");
  }
  fclose(fp);
  fflush(stdout);

  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMfree(GTM **pGTM)
  \brief Frees a lot of stuff, but not everything.
*/
int GTMfree(GTM **pGTM)
{
  GTM *gtm = *pGTM;

  MRIfree(&gtm->yvol);
  //MRIfree(&gtm->gtmseg);
  MRIfree(&gtm->mask);
  MatrixFree(&gtm->X);
  MatrixFree(&gtm->y);
  MatrixFree(&gtm->XtX);
  MatrixFree(&gtm->iXtX);
  MatrixFree(&gtm->Xty);
  MatrixFree(&gtm->beta);
  MatrixFree(&gtm->res);
  MatrixFree(&gtm->yhat);
  MRIfree(&gtm->ysynth);
  free(gtm);
  *pGTM=NULL;
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMmatrixY(GTM *gtm)
  \brief Converts the input gtm->yvol to a matrix using GTMvol2mat().
  It is important that this be done conistently with X, etc.
*/
int GTMmatrixY(GTM *gtm)
{
  gtm->y = GTMvol2mat(gtm, gtm->yvol, NULL);
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsetNMask(GTM *gtm)
  \brief Computes the number of voxels in the mask. If the mask is
  NULL, then just computes the number of voxels in the input.
*/
int GTMsetNMask(GTM *gtm)
{
  if(gtm->mask) gtm->nmask = MRIcountAboveThreshold(gtm->mask, 0.5);
  else gtm->nmask = gtm->yvol->width*gtm->yvol->height*gtm->yvol->depth;
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMpsfStd(GTM *gtm)
  \brief Convert the PSF {crs}FWHM to a standard deviation.
*/
int GTMpsfStd(GTM *gtm)
{
  gtm->cStd = gtm->cFWHM/sqrt(log(256.0));
  gtm->rStd = gtm->rFWHM/sqrt(log(256.0));
  gtm->sStd = gtm->sFWHM/sqrt(log(256.0));
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsegidlist(GTM *gtm)
  \brief Compute a sorted list of segmentation IDs from the segmentation
  itself (excludes 0). Just runs MRIsegIdListNot0().
*/
int GTMsegidlist(GTM *gtm)
{
  gtm->segidlist = MRIsegIdListNot0(gtm->anatseg, &gtm->nsegs, 0);
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn GTM *GTMalloc()
  \brief Allocates the GTM structure but nothing in the structure.
   sets PadThresh = .0001;
*/
GTM *GTMalloc()
{
  GTM *gtm;
  gtm = (GTM *) calloc(sizeof(GTM),1);
  gtm->PadThresh = .0001;
  return(gtm);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMnPad(GTM *gtm)
  \brief Computes the number of voxels used to pad the tightest
  fitting bounding box around the nonzero voxels of the seg. There
  must be enough padding to account for spill-over created by
  smoothing the input by the PSF. It works by determining the distance
  from the center of a Gaussian such that the kernel equals PadThresh
  (a value between 0 and 1). The FWHM of the Guassian is the maximum
  of the {crs}FWHM. The way this should be interpreted is that any
  spill-over signal outside of the brain that is excluded by the
  bounding box will be no larger than PadThresh times the unsmoothed
  signal. Eg, if PadThresh is .001, then nPad will be computed such
  that the spill-over signal will be 0.1% less than the original
  signal. Using the bounding box can greatly reduce the size of the
  input (esp for PET). This is only used to create a bounding box for
  an individual seg in GTMbuildX(). A similar method is used by
  GTMautoMask() to create a bounding box around the head mask.
*/
int GTMnPad(GTM *gtm)
{
  double maxFWHM, maxStd;
  maxFWHM = MAX(gtm->cFWHM/gtm->yvol->xsize,
		MAX(gtm->rFWHM/gtm->yvol->ysize,gtm->sFWHM/gtm->yvol->zsize));
  if(maxFWHM > 0){
    // The case where psf=0, just correcting for volume fraction
    maxStd = maxFWHM*sqrt(log(256.0));
    gtm->nPad = ceil(sqrt(-log(gtm->PadThresh*maxStd*sqrt(2*M_PI))*2*maxStd));
    printf("maxFWHM = %g (voxels), PadThresh=%g, nPad=%d\n",maxFWHM,gtm->PadThresh,gtm->nPad);
  }
  else gtm->nPad = 1;
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsolve(GTM *gtm)
  \brief Solves the GTM using a GLM. X must already have been created.
  Computes Xt, XtX, iXtX, beta, yhat, res, dof, rvar, kurtosis, and skew.
  Also will rescale if rescaling. Returns 1 and computes condition
  number if matrix cannot be inverted. Otherwise returns 0.
*/
int GTMsolve(GTM *gtm)
{
  struct timeb timer;
  int n,f;
  double sum;

  if(gtm->X == NULL){
    printf("ERROR: GTMsolve(): must build design matrix first\n");
    exit(1);
  }

  printf("Computing  XtX ... ");fflush(stdout);
  TimerStart(&timer);
  gtm->XtX = MatrixMtM(gtm->X,gtm->XtX);
  printf(" %4.1f sec\n",TimerStop(&timer)/1000.0);fflush(stdout);

  gtm->iXtX = MatrixInverse(gtm->XtX,gtm->iXtX);
  if(gtm->iXtX==NULL){
    gtm->XtXcond = MatrixConditionNumber(gtm->XtX);
    printf("ERROR: matrix cannot be inverted, cond=%g\n",gtm->XtXcond);
    return(1);
  }
  gtm->Xty  = MatrixAtB(gtm->X,gtm->y,gtm->Xty);
  gtm->beta = MatrixMultiplyD(gtm->iXtX,gtm->Xty,gtm->beta);
  if(gtm->rescale) GTMrescale(gtm);
  gtm->yhat = MatrixMultiplyD(gtm->X,gtm->beta,gtm->yhat);
  gtm->res  = MatrixSubtract(gtm->y,gtm->yhat,gtm->res);
  gtm->dof = gtm->X->rows - gtm->X->cols;
  if(gtm->rvar==NULL) gtm->rvar = MatrixAlloc(1,gtm->res->cols,MATRIX_REAL);
  for(f=0; f < gtm->res->cols; f++){
    sum = 0;
    for(n=0; n < gtm->res->rows; n++) sum += ((double)gtm->res->rptr[n+1][f+1]*gtm->res->rptr[n+1][f+1]);
    gtm->rvar->rptr[1][f+1] = sum/gtm->dof;
  }
  gtm->kurtosis = MatrixKurtosis(gtm->res,gtm->kurtosis);
  gtm->skew     = MatrixSkew(gtm->res,gtm->skew);
  return(0);
}
/*-----------------------------------------------------------------*/
/*
  \fn MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol)
  \brief Converts a matrix an MRI volume with rows=frames
  and columns=spatial dims. It is done row fastest, then col, then
  slice which makes it consistent with matlab. Any place that operates
  on the matrix data must be consistent when going from vol to matrix
  and back. See also GTMbuildX() and GTMvol2mat().
*/
MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol)
{
  int k,c,r,s,f;

  if(vol == NULL){
    vol = MRIallocSequence(gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth,MRI_FLOAT,m->cols);
    if(vol == NULL) return(NULL);
    MRIcopyHeader(gtm->yvol,vol);
    MRIcopyPulseParameters(gtm->yvol,vol);
  }
  if(MRIdimMismatch(gtm->yvol,vol,0)){
    printf("ERROR: GTMmat2vol() dim mismatch\n");
    return(NULL);
  }

  k = 0;
  for(s=0; s < gtm->yvol->depth; s++){
    for(c=0; c < gtm->yvol->width; c++){
      for(r=0; r < gtm->yvol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < gtm->yvol->nframes; f++)
	  MRIsetVoxVal(vol,c,r,s,f,m->rptr[k+1][f+1]);
	k++;
      }
    }
  }
  return(vol);
}
/*-----------------------------------------------------------------*/
/*
  \fn MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m)
  \brief Converts an MRI volume into a matrix with rows=frames
  and columns=spatial dims. It is done row fastest, then col, then
  slice which makes it consistent with matlab. Any place that operates
  on the matrix data must be consistent when going from vol to matrix
  and back. See also GTMbuildX() and GTMmat2vol(). 
*/
MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m)
{
  int k,c,r,s,f;
  
  if(m==NULL){
    m = MatrixAlloc(gtm->nmask,vol->nframes,MATRIX_REAL);
    if(m==NULL){
      printf("ERROR: GTMvol2mat(): could not alloc matrix %d %d\n",gtm->nmask,vol->nframes);
      return(NULL);
    }
  }
  if(m->rows != gtm->nmask){
    printf("ERROR: GTMvol2mat(): row mismatch %d %d\n",m->rows,gtm->nmask);
    return(NULL);
  }
  if(m->cols != vol->nframes){
    printf("ERROR: GTMvol2mat(): col mismatch %d %d\n",m->cols,vol->nframes);
    return(NULL);
  }

  k = 0;
  for(s=0; s < vol->depth; s++){ // crs order is important here!
    for(c=0; c < vol->width; c++){
      for(r=0; r < vol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < vol->nframes; f++)
	  m->rptr[k+1][f+1] = MRIgetVoxVal(vol,c,r,s,f);
	k++;
      }
    }
  }
  return(m);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMsegrvar(GTM *gtm)
  \brief Computes the residual variance in each segmentation. Not perfect
  because the resdiual has spill-over. Hopefully it is meaningful for something.
*/
int GTMsegrvar(GTM *gtm)
{
  int c,r,s,f,k;
  int nthseg=0,segid;
  double v;

  gtm->segrvar = MatrixAlloc(gtm->nsegs,gtm->beta->cols,MATRIX_REAL);
  gtm->nperseg = (int *)calloc(sizeof(int),gtm->nsegs);

  k=0;
  for(s=0; s < gtm->yvol->depth; s++){
    for(c=0; c < gtm->yvol->width; c++){
      for(r=0; r < gtm->yvol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	segid = MRIgetVoxVal(gtm->gtmseg,c,r,s,0);
	if(segid != 0){
	  for(nthseg=0; nthseg < gtm->nsegs; nthseg++) if(gtm->segidlist[nthseg]==segid) break;
	  gtm->nperseg[nthseg] ++;
	}
	for(f=0; f < gtm->beta->cols; f++){
	  v = gtm->res->rptr[k+1][f+1];
	  if(segid != 0) gtm->segrvar->rptr[nthseg+1][f+1] += v*v;
	}
	k++;
      }// r 
    }// c
  }// s

  for(f=0; f < gtm->beta->cols; f++) {
    for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
      v = gtm->segrvar->rptr[nthseg+1][f+1];
      gtm->segrvar->rptr[nthseg+1][f+1] = v/gtm->nperseg[nthseg];
    }
  }
  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrbv(GTM *gtm)
  \brief Performs Region-based Voxel-wise PVF correction. Can reduce the FoV
  of the output to a bounding box that tightly fits around the brain to
  reduce memory and disk space. The output shares a RAS space with the 
  anatseg (which shares a RAS space with the conformed anat) so a new
  registration is not needed. This is a re-write of GTMrbv0(), which
  did not manage memory very well. The RBV volume output is identical.
  Note: gtm->rbvsegmean is the mean input inside each seg for the RBV.
  It is a QA metric to compare against the GTM. If masking, it will 
  not be accurate for segs outside the brain.
 */
int GTMrbv(GTM *gtm)
{
  int c,r,s,f,nthseg,segid;
  double val,v,vhat0,vhat,v2;
  LTA *lta;
  struct timeb mytimer;
  MATRIX *nhits;
  MRI *yframe=NULL;
  MRI_REGION *region=NULL;
  MRI *yseg=NULL; // source volume trilin resampled to seg space (used with RBV)
  MRI *yhat0seg=NULL; // unsmoothed yhat created in seg space (used with RBV)
  MRI *yhatseg=NULL;  // smoothed yhat in seg space (used with RBV)

  if(gtm->rbv)      MRIfree(&gtm->rbv);

  TimerStart(&mytimer) ; 
  PrintMemUsage(stdout);

  if(gtm->mask_rbv_to_brain){
    // Reduce RBV to a tight bounding box around the brain by excluding anything
    // with "head" tissue type (hardcoded=5). This can greatly reduce 
    // memory requirements. The RAS space is still that of the rbvseg
    // (and so also that of the conformed anat) so no new registration is necessary
    printf("   masking RBV to brain\n");
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    int n, nReplace, ReplaceThis[1000], WithThat[1000];
    nReplace = 0;
    for(n=0; n < gtm->ctGTMSeg->nentries; n++){
      if(gtm->ctGTMSeg->entries[n] == NULL)  continue;
      if(gtm->ctGTMSeg->entries[n]->TissueType != 5) continue; // should not hard-code
      ReplaceThis[nReplace] = n;
      WithThat[nReplace] = 0;
      nReplace++;
    }
    printf("  replacing head voxels with 0\n");
    gtm->rbvsegmasked = MRIreplaceList(gtm->rbvseg, ReplaceThis, WithThat, nReplace, NULL);
    printf("  computing bounding box  %d %d %d ",gtm->rbvseg->width, gtm->rbvseg->height, gtm->rbvseg->depth);
    region = REGIONgetBoundingBox(gtm->rbvsegmasked,10);
    REGIONprint(stdout, region);
    MRI *tmp = MRIextractRegion(gtm->rbvsegmasked, NULL, region);
    MRIfree(&gtm->rbvsegmasked);
    gtm->rbvsegmasked = tmp;
  }
  else gtm->rbvsegmasked = gtm->rbvseg;

  //printf("writing gtm->rbvsegmasked\n");
  //MRIwrite(gtm->rbvsegmasked,"segbrain.mgh");

  printf("   Allocating RBV nvox=%d\n",gtm->rbvsegmasked->width*gtm->rbvsegmasked->height*gtm->rbvsegmasked->depth*gtm->nframes);
  fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  gtm->rbv = MRIallocSequence(gtm->rbvsegmasked->width, gtm->rbvsegmasked->height, gtm->rbvsegmasked->depth,
				MRI_FLOAT, gtm->nframes);
  if(gtm->rbv == NULL){
    printf("ERROR: GTMrbv() could not alloc rbv\n");
    return(1);
  }
  MRIcopyHeader(gtm->rbvsegmasked,gtm->rbv);
  MRIcopyPulseParameters(gtm->yvol,gtm->rbv);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  //if(gtm->mask_rbv_to_brain) MRIfree(&segbrain);

  // must be in rbvseg space
  yseg = MRIallocSequence(gtm->rbvseg->width, gtm->rbvseg->height, gtm->rbvseg->depth, MRI_FLOAT, 1);
  if(yseg == NULL){
    printf("ERROR: GTMrbv() could not alloc yseg\n");
    return(1);
  }
  MRIcopyHeader(gtm->rbvseg,yseg);
  MRIcopyPulseParameters(gtm->yvol,yseg);

  // Keep track of segmeans in RBV for QA
  gtm->rbvsegmean = MRIallocSequence(gtm->nsegs,1,1,MRI_FLOAT,gtm->nframes);

  printf("RBV looping over %d frames, t = %4.2f min \n",gtm->nframes,TimerStop(&mytimer)/60000.0);fflush(stdout);
  for(f=0; f < gtm->nframes; f++){
    printf("   f=%d t = %4.2f\n",f,TimerStop(&mytimer)/60000.0);fflush(stdout);
    yframe = fMRIframe(gtm->yvol,f,yframe); 

    printf("   Synthesizing unsmoothed input in seg space %4.2f \n",TimerStop(&mytimer)/60000.0);fflush(stdout);
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    yhat0seg = GTMsegSynth(gtm,f,yhat0seg);
    if(yhat0seg == NULL){
      printf("ERROR: GTMrbv() could not synthesize yhat0seg\n");
      return(1);
    }

    printf("   Smoothing synthesized in seg space %4.2f \n",TimerStop(&mytimer)/60000.0);fflush(stdout);
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    yhatseg = MRIgaussianSmoothNI(yhat0seg, gtm->cStd, gtm->rStd, gtm->sStd, yhatseg);
    if(yhatseg == NULL){
      printf("ERROR: GTMrbv() could not smooth yhatseg\n");
      return(1);
    }

    printf("   Sampling input to seg space with trilin %4.2f \n",TimerStop(&mytimer)/60000.0);fflush(stdout);
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    lta = LTAcopy(gtm->rbvseg2pet,NULL);
    LTAchangeType(lta,LINEAR_VOX_TO_VOX);
    MRIvol2Vol(yframe,yseg,(lta->xforms[0].m_L),SAMPLE_TRILINEAR, 0.0);
    LTAfree(&lta);

    printf("   Computing RBV %4.2f \n",TimerStop(&mytimer)/60000.0);fflush(stdout);
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    nhits = MatrixAlloc(gtm->beta->rows,1,MATRIX_REAL);
    for(c=0; c < gtm->rbvseg->width; c++){ // crs order not important
      for(r=0; r < gtm->rbvseg->height; r++){
	for(s=0; s < gtm->rbvseg->depth; s++){
	  segid = MRIgetVoxVal(gtm->rbvseg,c,r,s,0);
	  if(segid < 0.5) continue;

	  if(gtm->mask_rbv_to_brain){
	    if(c < region->x || c >= region->x+region->dx)  continue;
	    if(r < region->y || r >= region->y+region->dy)  continue;
	    if(s < region->z || s >= region->z+region->dz)  continue;
	  }

	  for(nthseg=0; nthseg < gtm->nsegs; nthseg++) if(gtm->segidlist[nthseg] == segid) break;
	  if(f==0) nhits->rptr[nthseg+1][1] ++;

	  v     = MRIgetVoxVal(yseg,c,r,s,0);
	  vhat0 = MRIgetVoxVal(yhat0seg,c,r,s,0);
	  vhat  = MRIgetVoxVal(yhatseg,c,r,s,0);
	  val = v*vhat0/(vhat+FLT_EPSILON); // RBV equation
	  if(gtm->mask_rbv_to_brain)
	    MRIsetVoxVal(gtm->rbv,c-region->x,r-region->y,s-region->z,f,val);
	  else
	    MRIsetVoxVal(gtm->rbv,c,r,s,f,val);

	  // track seg means for QA. Head Segs won't reflect QA if masking
	  v2 = MRIgetVoxVal(gtm->rbvsegmean,nthseg,0,0,f); 
	  MRIsetVoxVal(gtm->rbvsegmean,nthseg,0,0,f,v2+val);

	}
      }
    }
  }
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);
  MRIfree(&yseg);
  MRIfree(&yhat0seg);
  MRIfree(&yhatseg);
  MRIfree(&yframe);
  if(gtm->mask_rbv_to_brain) free(region); 

  // track seg means for QA
  for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
    for(f=0; f < gtm->nframes; f++){
      val = MRIgetVoxVal(gtm->rbvsegmean,nthseg,0,0,f)/nhits->rptr[nthseg+1][1];
      MRIsetVoxVal(gtm->rbvsegmean,nthseg,0,0,f,val);
    }
  }
  MatrixFree(&nhits);


  PrintMemUsage(stdout);
  printf("  RBV took %4.2f min\n",TimerStop(&mytimer)/60000.0);

  return(0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrbvseg(GTM *gtm)
  \brief Create the segmentation for the RBV by (possibly) increasing the
  voxel size of the anatseg to gtm->rbvsegres. If there is no res
  reduction, then rbvseg=anatseg and rbvseg2pet=seg2pet. This does not
  include reducing the FoV of the rbv seg to a tight box around the
  brain. This is done in GTMrbv().
 */
int GTMrbvseg(GTM *gtm)
{
  double res;
  LTA *new2seg;// *seg2new; // 

  if(gtm->rbvsegres <= 0 || gtm->rbvsegres == gtm->anatseg->xsize){
    printf("Not changing res of seg for RBV\n");
    gtm->rbvseg = gtm->anatseg;
    gtm->rbvseg2pet = gtm->seg2pet;
    return(0);
  }

  res = gtm->rbvsegres;
  printf("Changing res of seg for RBV to %g\n",res);
  gtm->rbvseg = MRIchangeSegRes(gtm->anatseg, res,res,res, gtm->ctGTMSeg, &new2seg);
  if(gtm->rbvseg == NULL) return(1);;
  gtm->rbvseg2pet = LTAconcat2(new2seg,gtm->seg2pet,1);
  if(gtm->rbvseg2pet == NULL) exit(1);
  LTAfree(&new2seg);

  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrbv0(GTM *gtm)
  \brief Performs Region-based Voxel-wise PVF correction. Can reduce the FoV
  of the output to a bounding box that tightly fits around the brain to
  reduce memory and disk space. The output shares a RAS space with the 
  anatseg (which shares a RAS space with the conformed anat) so a new
  registration is not needed. Use GTMrbv() instead since it manages
  memory much better and provides the same output.
 */
int GTMrbv0(GTM *gtm)
{
  int c,r,s,f,nthseg,segid;
  double val,v,vhat0,vhat,v2;
  LTA *lta;
  struct timeb mytimer;
  MATRIX *nhits;
  MRI *yseg; // source volume trilin resampled to seg space (used with RBV)
  MRI *yhat0seg; // unsmoothed yhat created in seg space (used with RBV)
  MRI *yhatseg;  // smoothed yhat in seg space (used with RBV)

  TimerStart(&mytimer) ; 

  printf("   Synthesizing unsmoothed input in seg space... ");fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  yhat0seg = GTMsegSynth(gtm,-1,NULL);
  if(yhat0seg == NULL){
    printf("ERROR: GTMrbv0() could not synthesize yhat0seg\n");
    return(1);
  }
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);


  printf("   Smoothing synthesized in seg space... ");fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  yhatseg = MRIgaussianSmoothNI(yhat0seg, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  if(yhatseg == NULL){
    printf("ERROR: GTMrbv0() could not smooth yhatseg\n");
    return(1);
  }
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);

  printf("   Sampling input to seg space with trilin... ");fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  yseg = MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth,
			      MRI_FLOAT, gtm->yvol->nframes);
  if(yseg == NULL){
    printf("ERROR: GTMrbv0() could not alloc yseg\n");
    return(1);
  }
  MRIcopyHeader(gtm->anatseg,yseg);
  MRIcopyPulseParameters(gtm->yvol,yseg);
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);

  lta = LTAcopy(gtm->rbvseg2pet,NULL);
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  MRIvol2Vol(gtm->yvol,yseg,(lta->xforms[0].m_L),SAMPLE_TRILINEAR, 0.0);
  LTAfree(&lta);

  printf("   Computing RBV ... ");fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  gtm->rbv = MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth,
			      MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->rbv == NULL){
    printf("ERROR: GTMrbv0() could not alloc rbv\n");
    return(1);
  }
  MRIcopyHeader(gtm->anatseg,gtm->rbv);
  MRIcopyPulseParameters(gtm->yvol,gtm->rbv);

  if(Gdiag_no > 0) PrintMemUsage(stdout);
  gtm->rbvsegmean = MRIallocSequence(gtm->nsegs,1,1,MRI_FLOAT,gtm->yvol->nframes);
  nhits = MatrixAlloc(gtm->beta->rows,1,MATRIX_REAL);
  for(s=0; s < gtm->anatseg->depth; s++){ // crs order not important
    for(c=0; c < gtm->anatseg->width; c++){
      for(r=0; r < gtm->anatseg->height; r++){
	segid = MRIgetVoxVal(gtm->anatseg,c,r,s,0);
	if(segid < 0.5) continue;
	for(nthseg=0; nthseg < gtm->nsegs; nthseg++) if(gtm->segidlist[nthseg] == segid) break;
	nhits->rptr[nthseg+1][1] ++;
	for(f=0; f < gtm->yvol->nframes; f++){
	  v     = MRIgetVoxVal(yseg,c,r,s,f);
	  vhat0 = MRIgetVoxVal(yhat0seg,c,r,s,f);
	  vhat  = MRIgetVoxVal(yhatseg,c,r,s,f);
	  val = v*vhat0/(vhat+FLT_EPSILON); // RBV equation
	  MRIsetVoxVal(gtm->rbv,c,r,s,f,val);
	  v2 = MRIgetVoxVal(gtm->rbvsegmean,nthseg,0,0,f); // track seg means for QA
	  MRIsetVoxVal(gtm->rbvsegmean,nthseg,0,0,f,v2+val);
	}
      }
    }
  }
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  MRIfree(&yseg);
  MRIfree(&yhat0seg);
  MRIfree(&yhatseg);
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);

  // track seg means for QA
  for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
    for(f=0; f < gtm->yvol->nframes; f++){
      val = MRIgetVoxVal(gtm->rbvsegmean,nthseg,0,0,f)/nhits->rptr[nthseg+1][1];
      MRIsetVoxVal(gtm->rbvsegmean,nthseg,0,0,f,val);
    }
  }
  MatrixFree(&nhits);

  if(gtm->mask_rbv_to_brain){
    // Reduce RBV to a tight mask around the brain. This can greatly reduce 
    // memory requirements. The RAS space is still that of the anatseg
    // (and so also that of the conformed anat) so no new registration is necessary
    printf("   masking RBV to brain\n");
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    int n, nReplace, ReplaceThis[1000], WithThat[1000];
    MRI *segtmp,*rbvtmp;
    MRI_REGION *region;
    nReplace = 0;
    for(n=0; n < gtm->ctGTMSeg->nentries; n++){
      if(gtm->ctGTMSeg->entries[n] == NULL)  continue;
      if(gtm->ctGTMSeg->entries[n]->TissueType != 5) continue; // should not hard-code
      ReplaceThis[nReplace] = n;
      WithThat[nReplace] = 0;
      nReplace++;
    }
    printf("  replacing head voxels with 0\n");
    segtmp = MRIreplaceList(gtm->anatseg, ReplaceThis, WithThat, nReplace, NULL);
    printf("  computing bounding box  ");
    region = REGIONgetBoundingBox(segtmp,10);
    REGIONprint(stdout, region);
    printf("  extracting bounding box\n");
    rbvtmp = MRIextractRegion(gtm->rbv, NULL, region);
    free(region); region=NULL;
    MRIfree(&segtmp);
    MRIfree(&gtm->rbv);
    gtm->rbv = rbvtmp;
  }
  PrintMemUsage(stdout);
  printf("  RBV took %4.2f min\n",TimerStop(&mytimer)/60000.0);

  return(0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMmgpvc(GTM *gtm)
  \brief Performs Muller-Gartner PVC. Hardcodes tissue type IDs to 
  be 0=cortex, 1=subcortexgm, 2=WM.
 */
int GTMmgpvc(GTM *gtm)
{
  int c,r,s,f,n,nthseg,segid,found,nhits;
  double sum,vgmpsf,vwmpsf,vwmtac,vtac,vmgtac;
  MRI *ctxpvf, *subctxpvf, *wmpvf, *gmpvf,*gmpvfpsf,*wmpvfpsf;

  // Compute the MG reference TAC
  gtm->mg_reftac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);
  for(f=0; f < gtm->yvol->nframes; f++){
    nhits = 0;
    sum = 0;
    for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      found = 0;
      for(n=0; n < gtm->n_mg_refids; n++) {
	if(segid == gtm->mg_refids[n]) {
	  found = 1;
	  nhits++;
	  break;
	}
      }
      if(!found) continue;
      sum += gtm->beta->rptr[nthseg+1][f+1];
      printf("   n=%d, nthseg=%d %g\n",n,nthseg,gtm->beta->rptr[nthseg+1][f+1]);
    }
    gtm->mg_reftac->rptr[f+1][1] = sum/nhits;
    printf("   wm tac %2d %2d %g\n",f,nhits,gtm->mg_reftac->rptr[f+1][1]);
  }

  if(gtm->mg) MRIfree(&gtm->mg);
  gtm->mg = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth,
			   MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->mg == NULL) return(1);
  MRIcopyHeader(gtm->yvol,gtm->mg);

  // Compute gray matter PVF with smoothing
  ctxpvf    = fMRIframe(gtm->ttpvf,0,NULL); // cortex PVF
  subctxpvf = fMRIframe(gtm->ttpvf,1,NULL); // subcortex GM PVF
  gmpvf = MRIadd(ctxpvf,subctxpvf,NULL); // All GM PVF
  // Smooth GM PVF by PSF
  gmpvfpsf = MRIgaussianSmoothNI(gmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);

  // WM PVF
  wmpvf = fMRIframe(gtm->ttpvf,2,NULL); 
  // Smooth GM PVF by PSF
  wmpvfpsf = MRIgaussianSmoothNI(wmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);

  // Finally, do the actual MG correction
  for(c=0; c < gtm->yvol->width; c++){ // crs order not important
    for(r=0; r < gtm->yvol->height; r++){
      for(s=0; s < gtm->yvol->depth; s++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue; 
	vgmpsf = MRIgetVoxVal(gmpvfpsf,c,r,s,0);
	if(vgmpsf < gtm->mg_gmthresh) continue; 
	vwmpsf = MRIgetVoxVal(wmpvfpsf,c,r,s,0);
	for(f=0; f < gtm->yvol->nframes; f++){
	  vwmtac = gtm->mg_reftac->rptr[f+1][1];
	  vtac = MRIgetVoxVal(gtm->yvol,c,r,s,f);
	  vmgtac = (vtac - vwmpsf*vwmtac)/vgmpsf;
	  MRIsetVoxVal(gtm->mg,c,r,s,f, vmgtac);
	}
      }
    }
  }
  MRIfree(&ctxpvf);
  MRIfree(&subctxpvf);
  MRIfree(&wmpvf);
  MRIfree(&gmpvf);
  MRIfree(&gmpvfpsf);
  MRIfree(&wmpvfpsf);
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsynth(GTM *gtm)
  \brief Synthesizes the unsmoothed PET image by computing 
   ysynth = X0*beta and then re-packing the result into a volume.
   This is then smoothed with GTMsmoothSynth() to give a full
   synthesis of the input (same as X*beta).
 */
int GTMsynth(GTM *gtm)
{
  MATRIX *yhat;

  if(gtm->ysynth) MRIfree(&gtm->ysynth);
  gtm->ysynth = MRIallocSequence(gtm->gtmseg->width,gtm->gtmseg->height,gtm->gtmseg->depth,MRI_FLOAT,
				 gtm->beta->cols);
  if(gtm->yvol) {
    MRIcopyHeader(gtm->yvol,gtm->ysynth);
    MRIcopyPulseParameters(gtm->yvol,gtm->ysynth);
  }
  yhat = MatrixMultiply(gtm->X0,gtm->beta,NULL);
  GTMmat2vol(gtm, yhat, gtm->ysynth);
  MatrixFree(&yhat);
    
  return(0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMsmoothSynth(GTM *gtm)
  \brief Smooths the synthesized volume (ysynth) to create the ysynthsm vol.
  Should give the same result as X*beta.
 */
int GTMsmoothSynth(GTM *gtm)
{
  if(gtm->ysynth == NULL) GTMsynth(gtm);
  gtm->ysynthsm = MRIgaussianSmoothNI(gtm->ysynth, gtm->cStd, gtm->rStd, gtm->sStd, gtm->ysynthsm);
  return(0);
}

/*------------------------------------------------------------------*/
/*
  \fn int GTMcheckX(MATRIX *X)
  \brief Checks that all colums sum to 1
 */
int GTMcheckX(MATRIX *X)
{
  int r,c,count;
  double sum,d,dmax;

  count=0;
  dmax = -1;
  for(r=0; r < X->rows; r++){
    sum = 0;
    for(c=0; c < X->cols; c++){
      sum += X->rptr[r+1][c+1];
    }
    d = abs(sum-1);
    if(d>.00001) count++;
    if(dmax < d) dmax = d;
  }

  printf("GTMcheckX: count=%d, dmax=%g\n",count,dmax);
  return(count);
}
/*------------------------------------------------------------------------------*/
/*
  \fn int GTMbuildX(GTM *gtm)
  \brief Builds the GTM design matrix both with (X) and without (X0) PSF.  If 
  gtm->DoVoxFracCor=1 then corrects for volume fraction effect.
*/
int GTMbuildX(GTM *gtm)
{
  int nthseg;
  struct timeb timer;

  if(gtm->X==NULL || gtm->X->rows != gtm->nmask || gtm->X->cols != gtm->nsegs){
    // Alloc or realloc X
    if(gtm->X) MatrixFree(&gtm->X);
    gtm->X = MatrixAlloc(gtm->nmask,gtm->nsegs,MATRIX_REAL);
    if(gtm->X == NULL){
      printf("ERROR: GTMbuildX(): could not alloc X %d %d\n",gtm->nmask,gtm->nsegs);
      return(1);
    }
  }
  if(gtm->X0==NULL || gtm->X0->rows != gtm->nmask || gtm->X0->cols != gtm->nsegs){
    if(gtm->X0) MatrixFree(&gtm->X0);
    gtm->X0 = MatrixAlloc(gtm->nmask,gtm->nsegs,MATRIX_REAL);
    if(gtm->X0 == NULL){
      printf("ERROR: GTMbuildX(): could not alloc X0 %d %d\n",gtm->nmask,gtm->nsegs);
      return(1);
    }
  }
  gtm->dof = gtm->X->rows - gtm->X->cols;

  TimerStart(&timer);

  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for(nthseg = 0; nthseg < gtm->nsegs; nthseg++){
    int segid,k,c,r,s;
    MRI *nthsegpvf=NULL,*nthsegpvfbb=NULL,*nthsegpvfbbsm=NULL;
    MRI_REGION *region;
    segid = gtm->segidlist[nthseg];
    //printf("nthseg = %d, %d %6.4f\n",nthseg,segid,TimerStop(&timer)/1000.0);fflush(stdout);
    if(gtm->DoVoxFracCor)
      nthsegpvf = fMRIframe(gtm->segpvf,nthseg,NULL); // extract PVF for this seg
    else
      nthsegpvf = MRIbinarizeMatch(gtm->gtmseg,&segid,1,0,NULL); // or get binary mask
    // Extract a region for the seg. This speeds up smoothing.
    region      = REGIONgetBoundingBox(nthsegpvf,gtm->nPad); // tight+pad bounding box 
    nthsegpvfbb = MRIextractRegion(nthsegpvf, NULL, region) ; // extract BB
    nthsegpvfbbsm = MRIgaussianSmoothNI(nthsegpvfbb, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
    // Fill X, creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way. See GTMvol2mat()
    k = 0;
    for(s=0; s < gtm->yvol->depth; s++){
      for(c=0; c < gtm->yvol->width; c++){
	for(r=0; r < gtm->yvol->height; r++){
	  if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	  k ++; // have to incr here in case continue below
	  if(c < region->x || c >= region->x+region->dx)  continue;
	  if(r < region->y || r >= region->y+region->dy)  continue;
	  if(s < region->z || s >= region->z+region->dz)  continue;
	  // do not use k+1 here because it has already been incr above
	  gtm->X->rptr[k][nthseg+1] = 
	    MRIgetVoxVal(nthsegpvfbbsm,c-region->x,r-region->y,s-region->z,0);
	  gtm->X0->rptr[k][nthseg+1] = 
	    MRIgetVoxVal(nthsegpvfbb,c-region->x,r-region->y,s-region->z,0);
	}
      }
    }
    MRIfree(&nthsegpvf);
    MRIfree(&nthsegpvfbb);
    MRIfree(&nthsegpvfbbsm);
  }
  printf(" Build time %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);

  return(0);

}

/*--------------------------------------------------------------------------*/
/*
  \fn MRI *GTMsegSynth(GTM *gtm, int frame, MRI *synth)
  \brief Creates a volume that is the same size as the anatseg in which the value
  at a voxel equals the beta of the seg ID at that voxel. This is used for RBV.
  If frame < 0, then all frames are are computed. If frame >= 0, then synth
  will have only one frame reprsenting the passed frame.
*/
MRI *GTMsegSynth(GTM *gtm, int frame, MRI *synth)
{
  int c,r,s,f,segid,segno,nframes;

  if(frame < 0) nframes = gtm->nframes;
  else          nframes = 1;

  if(synth == NULL){
    synth = MRIallocSequence(gtm->rbvseg->width,gtm->rbvseg->height,gtm->rbvseg->depth,MRI_FLOAT,nframes);
    if(synth == NULL){
      printf("ERROR: GTMsegSynth(): could not alloc %d frames\n",nframes);
      return(NULL);
    }
    MRIcopyHeader(gtm->rbvseg,synth);
    MRIcopyPulseParameters(gtm->yvol,synth);
  }

  for(c=0; c < gtm->rbvseg->width; c++){ // crs order does not matter here
    for(r=0; r < gtm->rbvseg->height; r++){
      for(s=0; s < gtm->rbvseg->depth; s++){
	segid = MRIgetVoxVal(gtm->rbvseg,c,r,s,0);
	if(segid == 0) continue;
	for(segno=0; segno < gtm->nsegs; segno++)
	  if(segid == gtm->segidlist[segno]) break;
	if(segno == gtm->nsegs){
	  printf("ERROR: GTMsegSynth(): could not find a match for segid=%d\n",segid);
	  for(segno=0; segno < gtm->nsegs; segno++) printf("%3d %5d\n",segno,gtm->segidlist[segno]);
	  return(NULL);
	}
	if(frame < 0){
	  for(f=0; f < gtm->beta->cols; f++)
	    MRIsetVoxVal(synth,c,r,s,f,gtm->beta->rptr[segno+1][f+1]);
	}
	else MRIsetVoxVal(synth,c,r,s,0,gtm->beta->rptr[segno+1][frame+1]);
      }
    }
  }

  return(synth);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
  \brief Prints the Muller-Gartner WM reference TAC to the given file pointer
*/
int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
{
  int f;
  for(f=0; f < gtm->yvol->nframes; f++)
    fprintf(fp,"%3d %10.5f\n",f,gtm->mg_reftac->rptr[f+1][1]);

  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMwriteMGRefTAC(GTM *gtm, char *filename)
  \brief Writes the Muller-Gartner WM reference TAC to the given file
*/
int GTMwriteMGRefTAC(GTM *gtm, char *filename)
{
  FILE *fp;
  fp = fopen(filename,"w");
  GTMprintMGRefTAC(gtm, fp);
  fclose(fp);
  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrescale(GTM *gtm)
  \brief Computes global rescaling factor and applies it to yvol, beta, and y.
  The factor = 100/mean(beta_i) where i is the list of scale seg IDs (scale_ref_ids)
*/
int GTMrescale(GTM *gtm)
{
  int f,n,nthseg,segid,found,nhits;
  double sum;

  // global rescaling
  nhits = 0;
  sum = 0;
  for(f=0; f < gtm->yvol->nframes; f++){
    for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      found = 0;
      for(n=0; n < gtm->n_scale_refids; n++) {
	if(segid == gtm->scale_refids[n]) {
	  found = 1;
	  nhits++;
	  break;
	}
      }
      if(!found) continue;
      sum += gtm->beta->rptr[nthseg+1][f+1];
    }
  }
  gtm->scale = 100/(sum/nhits);
  printf("gtm multiplicative scale %g\n",gtm->scale);

  MRImultiplyConst(gtm->yvol,gtm->scale,gtm->yvol);
  MatrixScalarMul(gtm->beta, gtm->scale,gtm->beta);
  MatrixScalarMul(gtm->y,    gtm->scale,gtm->y);

  return(0);
}

/*-------------------------------------------------------------------------------*/
/*
  \fn int GTMttest(GTM *gtm)
  \brief Computes a t-test for each contrast. This includes computing gamma, 
  gammavar, t, and p.
*/
int GTMttest(GTM *gtm)
{
  MATRIX *Ct,*CiXtX,*CiXtXCt;
  int n,f,nframes;
  GTMCON *gtmc;
  double F;

  nframes = gtm->beta->cols;

  for(n=0; n < gtm->nContrasts; n++){
    gtmc = gtm->contrasts[n];
    gtmc->gamma = MatrixMultiply(gtmc->C,gtm->beta,NULL);
    gtmc->gammavar = MatrixAlloc(1,nframes,MATRIX_REAL);
    gtmc->t = MatrixAlloc(1,nframes,MATRIX_REAL);
    gtmc->p = MatrixAlloc(1,nframes,MATRIX_REAL);
    Ct = MatrixTranspose(gtmc->C,NULL);
    CiXtX = MatrixMultiply(gtmc->C,gtm->iXtX,NULL);
    CiXtXCt = MatrixMultiply(CiXtX,Ct,NULL);
    for(f=0; f < nframes; f++){
      gtmc->gammavar->rptr[1][f+1] = CiXtXCt->rptr[1][1]*gtm->rvar->rptr[1][f+1];
      gtmc->t->rptr[1][f+1] = gtmc->gamma->rptr[1][f+1]/sqrt(gtmc->gammavar->rptr[1][f+1]);
      F = gtmc->t->rptr[1][f+1]*gtmc->t->rptr[1][f+1];
      gtmc->p->rptr[1][f+1] = sc_cdf_fdist_Q(F,gtmc->C->rows,gtm->dof);
    }
    MatrixFree(&Ct);
    MatrixFree(&CiXtX);
    MatrixFree(&CiXtXCt);
  }
  return(0);
}
/*------------------------------------------------------------------------------*/
/*
  \fn int GTMwriteContrasts(GTM *GTM)
  \brief Creates an ASCII file in the output folder for each contrast with gamma, 
  gammavar, t, and p
 */
int GTMwriteContrasts(GTM *gtm)
{
  int n,nframes,f;
  GTMCON *gtmc;
  char tmpstr[5000];
  FILE *fp;

  nframes = gtm->beta->cols;

  for(n=0; n < gtm->nContrasts; n++){
    gtmc = gtm->contrasts[n];
    sprintf(tmpstr,"%s/%s.mat",gtm->OutDir,gtmc->name);
    MatlabWrite(gtmc->C, tmpstr,"C");
    sprintf(tmpstr,"%s/%s.dat",gtm->OutDir,gtmc->name);
    fp = fopen(tmpstr,"w");
    for(f=0; f < nframes; f++){
      fprintf(fp,"%2d %20.15f %20.15f %20.15f %20.15f\n",f,
	      gtmc->gamma->rptr[1][1],gtmc->gammavar->rptr[1][1],
	      gtmc->t->rptr[1][1],gtmc->p->rptr[1][1]);
    }
    fclose(fp);
  }
  return(0);
}

/*-------------------------------------------------------------------------*/
/*
  \fn int GTMcheckRefIds(GTM *gtm)
  \brief Checks the segmentation IDs used for rescaling, MG, KM Ref, and KM HB
  to make sure that they are in the segmentation.
 */
int GTMcheckRefIds(GTM *gtm)
{
  int n,m,ok;

  if(gtm->rescale){
    for(n=0; n < gtm->n_scale_refids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->scale_refids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: scale reference id %d cannot be found in seg id list\n",gtm->scale_refids[n]);
	fprintf(gtm->logfp,"ERROR: scale reference id %d cannot be found in seg id list\n",gtm->scale_refids[n]);
	return(1);
      }
    }
  }

  if(gtm->DoMGPVC){
    for(n=0; n < gtm->n_mg_refids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->mg_refids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: MG reference id %d cannot be found in seg id list\n",gtm->mg_refids[n]);
	fprintf(gtm->logfp,"ERROR: scale reference id %d cannot be found in seg id list\n",gtm->mg_refids[n]);
	return(1);
      }
    }
  }

  if(gtm->DoKMRef){
    for(n=0; n < gtm->n_km_refids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->km_refids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: KM reference id %d cannot be found in seg id list\n",gtm->km_refids[n]);
	fprintf(gtm->logfp,"ERROR: scale reference id %d cannot be found in seg id list\n",gtm->km_refids[n]);
	return(1);
      }
    }
  }

  if(gtm->DoKMHB){
    for(n=0; n < gtm->n_km_hbids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->km_hbids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: KM high binding id %d cannot be found in seg id list\n",gtm->km_hbids[n]);
	fprintf(gtm->logfp,"ERROR: scale high binding id %d cannot be found in seg id list\n",gtm->km_hbids[n]);
	return(1);
      }
    }
  }

  return(0);
}

/*
  \fn int GTMprintRefIds(GTM *gtm, FILE *fp)
  \brief Prints the segmentation IDs used for rescaling, MG, KM Ref, and KM HB
  to the given FILE pointer.
 */
int GTMprintRefIds(GTM *gtm, FILE *fp)
{
  int n;

  if(gtm->rescale){
    fprintf(fp,"Segmentations used for rescaling\n");
    for(n=0; n < gtm->n_scale_refids; n++) 
      fprintf(fp,"%4d %s\n",gtm->scale_refids[n],gtm->ctGTMSeg->entries[gtm->scale_refids[n]]->name);
  }

  if(gtm->DoMGPVC){
    fprintf(fp,"Segmentations used for MG PVC WM reference\n");
    for(n=0; n < gtm->n_mg_refids; n++) 
      fprintf(fp,"%4d %s\n",gtm->mg_refids[n],gtm->ctGTMSeg->entries[gtm->mg_refids[n]]->name);
  }

  if(gtm->DoKMRef){
    fprintf(fp,"Segmentations used for KM reference TAC\n");
    for(n=0; n < gtm->n_km_refids; n++)
      fprintf(fp,"%4d %s\n",gtm->km_refids[n],gtm->ctGTMSeg->entries[gtm->km_refids[n]]->name);
  }

  if(gtm->DoKMHB){
    fprintf(fp,"Segmentations used for KM High Binding TAC\n");
    for(n=0; n < gtm->n_km_hbids; n++)
      fprintf(fp,"%4d %s\n",gtm->km_hbids[n],gtm->ctGTMSeg->entries[gtm->km_hbids[n]]->name);
  }
  fflush(fp);
  return(0);
}
/*
  \fn int GTMrefTAC(GTM *gtm)
  \brief Computes KM reference and hibinding TACs. Also writes them
  to the output folder. The HB TAC is written as both an ascii file
  and a nii.gz (the later needed for KM analysis with mri_glmfit)
 */
int GTMrefTAC(GTM *gtm)
{
  int f,n,nthseg,segid;
  double sum;
  char tmpstr[5000];
  FILE *fp;

  if(gtm->DoKMRef){
    sprintf(tmpstr,"%s/km.ref.tac.dat",gtm->OutDir);
    fp = fopen(tmpstr,"w");
    if(fp==NULL) return(1);
    gtm->km_reftac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);
    for(f=0; f < gtm->yvol->nframes; f++){
      sum = 0;
      for(n=0; n < gtm->n_km_refids; n++) {
	for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
	  segid = gtm->segidlist[nthseg];
	  if(segid == gtm->km_refids[n]) break;
	}
	sum += gtm->beta->rptr[nthseg+1][f+1];
      }
      gtm->km_reftac->rptr[f+1][1] = sum/gtm->n_km_refids;
      fprintf(fp,"%20.15lf\n",sum/gtm->n_km_refids);
    }
    fclose(fp);
  }

  if(gtm->DoKMHB){
    sprintf(tmpstr,"%s/km.hb.tac.dat",gtm->OutDir);
    fp = fopen(tmpstr,"w");
    if(fp==NULL) return(1);
    gtm->km_hbtac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);
    for(f=0; f < gtm->yvol->nframes; f++){
      sum = 0;
      for(n=0; n < gtm->n_km_hbids; n++) {
	for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
	  segid = gtm->segidlist[nthseg];
	  if(segid == gtm->km_hbids[n]) break;
	}
	sum += gtm->beta->rptr[nthseg+1][f+1];
      }
      gtm->km_hbtac->rptr[f+1][1] = sum/gtm->n_km_hbids;
      fprintf(fp,"%20.15lf\n",gtm->km_hbtac->rptr[f+1][1]);
    }
    fclose(fp);

    MRI *mritmp = MRIallocSequence(1, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
    for(f=0; f < gtm->yvol->nframes; f++)
      MRIsetVoxVal(mritmp,0,0,0,f, gtm->km_hbtac->rptr[f+1][1]);
    sprintf(tmpstr,"%s/km.hb.tac.nii.gz",gtm->OutDir);
    MRIwrite(mritmp,tmpstr);
    MRIfree(&mritmp);
  }

  return(0);
}

/*
  \fn int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat)
  \brief Prints the replacement list to the FILE pointer.
 */
int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat)
{
  int n;
  for(n=0; n < nReplace; n++)
    fprintf(fp,"%5d %5d\n",ReplaceThis[n],WithThat[n]);
  return(0);
}

/*
  \fn int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat)
  \brief Checks replacement list to make sure that no item in ReplaceThis list appears in
  the WithThat list.
 */
int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat)
{
  int n,m;
  for(n=0; n < nReplace; n++){
    for(m=0; m < nReplace; m++){
      if(ReplaceThis[n] == WithThat[m]){
	printf("ERROR: item %d appears as both source and target seg id in replacement list\n",ReplaceThis[n]);
	return(1);
      }
    }
  }
  return(0);
}

/*
  \fn int GTMloadReplacmentList(const char *fname, int *nReplace, int *ReplaceThis, int *WithThat)
  \brief Loads in data from a file. The file should have two columns of numbers. The first
  is the segmentation ID to replace the second is the segmentation ID to replace it with.
  It will ignore any line that begins with a #. 
 */
int GTMloadReplacmentList(const char *fname, int *nReplace, int *ReplaceThis, int *WithThat)
{
  FILE *fp;
  int nlist,r,nth;
  char tmpstr[1001], *s;

  fp = fopen(fname,"r");
  if(fp==NULL){
    int e=errno;
    printf("ERROR: could not open %s %d\n",fname,e);
    printf("%s\n",strerror(e));
    return(1);
  }
  nlist = *nReplace;

  nth = -1;
  while(1){
    nth++;
    s = fgets(tmpstr,1000,fp);
    if(s==NULL) break;
    if(tmpstr[0] == '#') continue;
    r = CountItemsInString(tmpstr);
    if(r != 2){
      printf("ERROR: line %d in %s has the wrong format (has %d elements, expecting 2)\n",nth,fname,r);
      return(1);
    }
    sscanf(tmpstr,"%d %d",&ReplaceThis[nlist],&WithThat[nlist]); 
    nlist++;
  }
  if(nlist == 0){
    printf("ERROR: could not find any replacement items in %s\n",fname);
    return(1);
  }
  printf("Read in %d replacement items from %s\n",nlist,fname);
  *nReplace += nlist;
  return(0);
}
/*------------------------------------------------------------*/
/*
  \fn int GTMautoMask(GTM *gtm)
  \brief Computes a mask in PET space based on the segmentation and
  smoothing level. Optionally, it reduces the FoV of PET to be tight
  with the mask. If this is done, it recomputes the registration
  matrices.  The header of the full FoV PET is kept in
  gtm->yvol_full_fov but the pixel data are freed to reduce storage.
*/
int GTMautoMask(GTM *gtm)
{
  LTA *lta,*pet2bbpet,*seg2bbpet,*anat2bbpet;
  double std;
  MRI *masktmp,*yvoltmp;

  gtm->mask = MRIalloc(gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth,MRI_FLOAT);
  MRIcopyHeader(gtm->yvol,gtm->mask);
  if(LTAmriIsSource(gtm->seg2pet,gtm->yvol)) lta = LTAcopy(gtm->seg2pet,NULL);
  else                                       lta = LTAinvert(gtm->seg2pet,NULL);        
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);

  // Sample the seg into the pet space
  MRIvol2Vol(gtm->anatseg,gtm->mask,lta->xforms[0].m_L,SAMPLE_NEAREST,0);
  LTAfree(&lta);
  // Threshold it at 0.5 to give the mask
  MRIbinarize(gtm->mask,gtm->mask,0.5,0,1);
  // Smooth binary mask
  std = gtm->automask_fwhm/sqrt(log(256.0)); // convert fwhm to std
  MRIgaussianSmoothNI(gtm->mask, std,std,std, gtm->mask);
  // Binarize again to get final mask
  MRIbinarize(gtm->mask,gtm->mask,gtm->automask_thresh,0,1);

  if(gtm->automask_reduce_fov){
    printf("Automask, reducing FOV\n");
    gtm->automaskRegion = REGIONgetBoundingBox(gtm->mask,1);
    printf("region %d %d %d reduced to ",gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth);
    REGIONprint(stdout, gtm->automaskRegion);
    fflush(stdout);
    masktmp = MRIextractRegion(gtm->mask, NULL, gtm->automaskRegion);
    yvoltmp = MRIextractRegion(gtm->yvol, NULL, gtm->automaskRegion);
    pet2bbpet = TransformRegDat2LTA(gtm->yvol, yvoltmp, NULL);
    seg2bbpet = LTAconcat2(gtm->seg2pet,pet2bbpet, 1);
    if(LTAmriIsSource(gtm->anat2pet,gtm->yvol)) lta = LTAinvert(gtm->anat2pet,NULL);        
    else                                        lta = LTAcopy(gtm->anat2pet,NULL);
    anat2bbpet = LTAconcat2(lta,pet2bbpet,1);
    LTAfree(&lta);
    gtm->yvol_full_fov = MRIallocHeader(gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth,MRI_FLOAT,1);
    MRIcopyHeader(gtm->yvol,gtm->yvol_full_fov);
    MRIcopyPulseParameters(gtm->yvol,gtm->yvol_full_fov);
    MRIfree(&gtm->yvol);
    gtm->yvol = yvoltmp;
    MRIfree(&gtm->mask);
    gtm->mask = masktmp;
    LTAfree(&gtm->seg2pet);
    gtm->seg2pet = seg2bbpet;
    LTAfree(&gtm->anat2pet);
    gtm->anat2pet = anat2bbpet;
  }

  return(0);
}
/*
  \fn int GTMrvarGM(GTM *gtm)
  \brief Computes residual variance only in GM structures. The measure
  is not perfect because there will be spill-over from adjacent
  structures. But it will not include much from outside of the head (as
  the standard rvar measure will). The result is stored in
  gtm->rvargm. This hard-codes GM as 1 and 2 in the tissue type CTAB.
  Also writes rvar.gm.dat in the output folder. The value is computed
  as the sum of res.^2 in GM divided by the number of GM voxels.x  
*/
int GTMrvarGM(GTM *gtm)
{
  COLOR_TABLE *ct;
  int f,s,c,r,n,nhits,segid,tt;
  double sum;
  FILE *fp;
  char tmpstr[2000];

  if(gtm->rvargm==NULL) gtm->rvargm = MatrixAlloc(1,gtm->res->cols,MATRIX_REAL);
  ct = gtm->ctGTMSeg;

  for(f=0; f < gtm->res->cols; f++){
    sum = 0;
    n = -1;
    nhits = 0;
    // slice, col, row order is important here as is skipping the mask
    for(s=0; s < gtm->yvol->depth; s++){
      for(c=0; c < gtm->yvol->width; c++){
	for(r=0; r < gtm->yvol->height; r++){
	  if(MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	  n++;
	  segid = MRIgetVoxVal(gtm->gtmseg,c,r,s,0);
	  tt = ct->entries[segid]->TissueType;
	  if(tt != 1 && tt != 2) continue; // not GM (hardcoded)
	  sum += ((double)gtm->res->rptr[n+1][f+1]*gtm->res->rptr[n+1][f+1]);
	  nhits ++;
	}
      }
    }
    gtm->rvargm->rptr[1][f+1] = sum/nhits;
    if(f==0) printf("rvargm %2d %6.4f\n",f,gtm->rvargm->rptr[1][f+1]);
  }

  sprintf(tmpstr,"%s/rvar.gm.dat",gtm->AuxDir);
  fp = fopen(tmpstr,"w");
  for(f=0; f < gtm->res->cols; f++) fprintf(fp,"%30.20f\n",gtm->rvargm->rptr[1][f+1]);
  fclose(fp);

  return(0);
}
