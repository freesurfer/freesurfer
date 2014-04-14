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
 *    $Date: 2014/04/14 22:13:46 $
 *    $Revision: 1.2 $
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
  MRIS *lhw, *lhp, *rhw, *rhp;
  struct timeb timer;
  TimerStart(&timer);

  printf("Starting MRIgtmSeg() USF=%d\n",gtmseg->USF);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  SUBJECTS_DIR = "/autofs/cluster/con_009/users/greve/fdg-pvc/FSMR";

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
  lhw = MRISread(tmpstr); if(lhw==NULL) return(1);

  sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,gtmseg->subject);
  lhp = MRISread(tmpstr);  if(lhp==NULL) return(1);

  sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,gtmseg->subject);
  rhw = MRISread(tmpstr);  if(rhw==NULL) return(1);

  sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,gtmseg->subject);
  rhp = MRISread(tmpstr);  if(rhp==NULL) return(1);

  printf("Loading annotations ");
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  if(gtmseg->wmannotfile != NULL){
    sprintf(tmpstr,"%s/%s/label/lh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->wmannotfile);
    err = MRISreadAnnotation(lhw,tmpstr);  
    if(err) {
      printf("Try running mri_annotation2label --s %s --hemi lh --lobesStrict lobefilename\n",gtmseg->subject);
      return(1);
    }
    sprintf(tmpstr,"%s/%s/label/rh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->wmannotfile);
    err = MRISreadAnnotation(rhw,tmpstr);
    if(err) {
      printf("Try running mri_annotation2label --s %s --hemi rh --lobesStrict lobefilename\n",gtmseg->subject);
      return(1);
    }
    MRISripUnknown(lhw);
    MRISripUnknown(rhw);
    lhw->ct->idbase = gtmseg->wmlhbase; 
    rhw->ct->idbase = gtmseg->wmrhbase; 
  }
  else printf("Not segmenting WM\n");

  sprintf(tmpstr,"%s/%s/label/lh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->ctxannotfile);
  err = MRISreadAnnotation(lhp,tmpstr);  if(err) return(1);
  sprintf(tmpstr,"%s/%s/label/rh.%s",SUBJECTS_DIR,gtmseg->subject,gtmseg->ctxannotfile);
  err = MRISreadAnnotation(rhp,tmpstr);  if(err) return(1);

  lhp->ct->idbase = gtmseg->ctxlhbase;
  rhp->ct->idbase = gtmseg->ctxrhbase;

  MRISsetPialUnknownToWhite(lhw, lhp);
  MRISsetPialUnknownToWhite(rhw, rhp);
  MRISripUnknown(lhp);
  MRISripUnknown(rhp);

  // relabel WM_hypointensities as {Left,Right}_WM_hypointensities as 
  printf(" Relabeling any unlateralized hypointenities as lateralized hypointenities\n");
  MRIrelabelHypoHemi(aseg, lhw, rhw, NULL, aseg);

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
  if(nlist > 0) MRIunsegmentWM(aseg, lhw, rhw, list, nlist, NULL, aseg);

  // Upsample the segmentation
  printf("Upsampling segmentation USF = %d",gtmseg->USF);fflush(stdout);
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  hrseg = MRIhiresSeg(aseg, lhw, lhp, rhw, rhp, gtmseg->USF, &gtmseg->anat2seg);
  if(hrseg == NULL) return(1);
  MRIfree(&aseg);

  // Label cortex (like aparc+aseg)
  printf("Beginning cortical segmentation using %s",gtmseg->ctxannotfile);fflush(stdout);
  printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
  ctxseg = MRIannot2CorticalSeg(hrseg, lhw, lhp, rhw, rhp, NULL, NULL);
  MRIfree(&hrseg);

  // Label wm (like wmaparc)
  if(gtmseg->wmannotfile != NULL){
    printf("Beginning WM segmentation using %s",gtmseg->wmannotfile);fflush(stdout);
    printf(" t = %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);
    ctxseg = MRIannot2CerebralWMSeg(ctxseg, lhw, rhw, gtmseg->dmax, NULL, ctxseg);
  } else     printf("Not subsegmenting WM\n");

  MRISfree(&lhw);
  MRISfree(&lhp);
  MRISfree(&rhw);
  MRISfree(&rhp);

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
int GTMdefaultSegReplacmentList(GTMSEG *gtmseg)
{

  int nlist, *srclist, *targlist;

  srclist  = &(gtmseg->srclist[0]);
  targlist = &(gtmseg->targlist[0]);

  nlist = 0;
  srclist[nlist] = 1033; targlist[nlist] = 1030; nlist++; // temppole=stg
  srclist[nlist] = 2033; targlist[nlist] = 2030; nlist++; // temppole=stg
  srclist[nlist] = 1034; targlist[nlist] = 1030; nlist++; // transtemp=stg
  srclist[nlist] = 2034; targlist[nlist] = 1030; nlist++; // transtemp=stg
  srclist[nlist] = 1001; targlist[nlist] = 1015; nlist++; // bankssts=mtg
  srclist[nlist] = 2001; targlist[nlist] = 2015; nlist++; // bankssts=mtg
  srclist[nlist] = 1032; targlist[nlist] = 1027; nlist++; // frontpole=rmf
  srclist[nlist] = 2032; targlist[nlist] = 2027; nlist++; // frontpole=rmf
  //srclist[nlist] = 1016; targlist[nlist] = 1006; nlist++; // parahip=entorhinal ?
  //srclist[nlist] = 2016; targlist[nlist] = 2006; nlist++; // parahip=entorhinal ?

  // There should not be any cortex unknown after MRIannot2CorticalSeg()
  srclist[nlist] = 1000; targlist[nlist] =    0; nlist++; // cortex unknown
  srclist[nlist] = 2000; targlist[nlist] =    0; nlist++; // cortex unknown

  // Should I replace subcorts before hires seg?
  srclist[nlist] =   85; targlist[nlist] =    0; nlist++; // optic chiasm
  srclist[nlist] =    4; targlist[nlist] =   24; nlist++; // LLatVent
  srclist[nlist] =    5; targlist[nlist] =   24; nlist++; // LInfLatVent
  srclist[nlist] =   14; targlist[nlist] =   24; nlist++; // 3rd
  srclist[nlist] =   15; targlist[nlist] =   24; nlist++; // 4th
  srclist[nlist] =   72; targlist[nlist] =   24; nlist++; // 5th
  srclist[nlist] =   31; targlist[nlist] =   24; nlist++; // LChoroidP ?
  srclist[nlist] =   43; targlist[nlist] =   24; nlist++; // RLatVent
  srclist[nlist] =   44; targlist[nlist] =   24; nlist++; // RInfLatVent
  srclist[nlist] =   63; targlist[nlist] =   24; nlist++; // RChoroidP ?
  srclist[nlist] =   30; targlist[nlist] =   24; nlist++; // LVessel ?
  srclist[nlist] =   62; targlist[nlist] =   24; nlist++; // RVessel ?
  srclist[nlist] =   80; targlist[nlist] =   24; nlist++; // non-WM-hypo ?

  /* Repace CC segments with one CC if not unsegmenting CC */
  if(gtmseg->KeepCC){
    srclist[nlist] =  251; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  252; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  253; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  254; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  255; targlist[nlist] =  192; nlist++; 
  }

  gtmseg->nlist = nlist;
  return(0);
}
