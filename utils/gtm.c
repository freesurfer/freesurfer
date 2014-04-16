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
 *    $Date: 2014/04/16 19:22:44 $
 *    $Revision: 1.4 $
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

/*------------------------------------------------------------------------------------*/
int GTMSEGprint(GTMSEG *gtmseg, FILE *fp)
{
  fprintf(fp,"subject %s\n",gtmseg->subject);
  fprintf(fp,"USF %d\n",gtmseg->USF);
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
  fflush(fp);
  return(0);
}

/*------------------------------------------------------------------------------------*/
int MRIgtmSeg(GTMSEG *gtmseg)
{
  int err,*segidlist,nsegs,n;
  char *SUBJECTS_DIR, tmpstr[5000];
  MRI *apas, *ribbon, *aseg, *hrseg, *ctxseg;
  struct timeb timer;
  TimerStart(&timer);

  printf("Starting MRIgtmSeg() USF=%d\n",gtmseg->USF);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->apasfile);
  printf("Loading %s\n",tmpstr);
  apas = MRIread(tmpstr);
  if(apas==NULL) return(1);

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,gtmseg->subject,"ribbon.mgz");
  printf("Loading %s\n",tmpstr);
  ribbon = MRIread(tmpstr);
  if(ribbon==NULL) return(1);

  aseg = MRIunsegmentCortex(apas, ribbon, NULL); // aseg+head
  if(aseg == NULL) return(1);
  MRIfree(&apas);
  MRIfree(&ribbon);

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
  free(segidlist);

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

