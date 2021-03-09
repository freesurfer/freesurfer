/**
 * @brief Routines to create and analyze the Geometric Transfer Matrix (GTM)
 *
 */
/*
 * Original Author: Douglas N. Greve
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

#include "gtm.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "cma.h"
#include "cmdargs.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "fmriutils.h"
#include "macros.h"
#include "matfile.h"
#include "mri.h"
#include "mri2.h"
#include "mrimorph.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "numerics.h"
#include "region.h"
#include "resample.h"
#include "timer.h"
#include "utils.h"
#include "version.h"

#include "romp_support.h"


/*------------------------------------------------------------------------------------*/
int GTMSEGprint(GTMSEG *gtmseg, FILE *fp)
{
  fprintf(fp, "subject %s\n", gtmseg->subject);
  fprintf(fp, "USF %d\n", gtmseg->USF);
  fprintf(fp, "OutputUSF %d\n", gtmseg->OutputUSF);
  fprintf(fp, "apasfile %s\n", gtmseg->apasfile);
  if (gtmseg->wmannotfile != NULL) {
    fprintf(fp, "wmannotfile %s\n", gtmseg->wmannotfile);
    fprintf(fp, "wmlhbase %d\n", gtmseg->wmlhbase);
    fprintf(fp, "wmrhbase %d\n", gtmseg->wmrhbase);
  }
  else
    fprintf(fp, "wmannotfile NULL\n");
  fprintf(fp, "ctxannotfile %s\n", gtmseg->ctxannotfile);
  fprintf(fp, "ctxlhbase %d\n", gtmseg->ctxlhbase);
  fprintf(fp, "ctxrhbase %d\n", gtmseg->ctxrhbase);
  fprintf(fp, "SubSegWM %3d\n", gtmseg->SubSegWM);
  fprintf(fp, "KeepHypo %3d\n", gtmseg->KeepHypo);
  fprintf(fp, "KeepCC %3d\n", gtmseg->KeepCC);
  fprintf(fp, "dmax %f\n", gtmseg->dmax);
  fprintf(fp, "nlist %3d\n", gtmseg->nlist);
  fprintf(fp, "lhmin %5d\n", gtmseg->lhmin);
  fprintf(fp, "lhmax %5d\n", gtmseg->lhmax);
  fprintf(fp, "rhmin %5d\n", gtmseg->rhmin);
  fprintf(fp, "rhmax %5d\n", gtmseg->rhmax);
  fflush(fp);
  return (0);
}

/*------------------------------------------------------------------------------------*/
int MRIgtmSeg(GTMSEG *gtmseg)
{
  int err, *segidlist, nsegs, n;
  char *SUBJECTS_DIR, tmpstr[5000];
  MRI *apas, *aseg, *hrseg, *ctxseg;
  Timer timer;

  printf("Starting MRIgtmSeg() USF=%d\n", gtmseg->USF);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr, "%s/%s/mri/%s", SUBJECTS_DIR, gtmseg->subject, gtmseg->apasfile);
  printf("Loading %s\n", tmpstr);
  apas = MRIread(tmpstr);
  if (apas == NULL) return (1);

  gtmseg->anat = MRIcopyHeader(apas, NULL);  // keep anat header around

  aseg = MRIunsegmentCortex(apas, gtmseg->lhmin, gtmseg->lhmax, gtmseg->rhmin, gtmseg->rhmax, NULL);  // aseg+head
  if (aseg == NULL) return (1);
  MRIfree(&apas);

  printf("Loading surfaces ");
  printf(" t = %6.4f\n", timer.seconds());
  fflush(stdout);
  sprintf(tmpstr, "%s/%s/surf/lh.white", SUBJECTS_DIR, gtmseg->subject);
  gtmseg->lhw = MRISread(tmpstr);
  if (gtmseg->lhw == NULL) return (1);

  sprintf(tmpstr, "%s/%s/surf/lh.pial", SUBJECTS_DIR, gtmseg->subject);
  gtmseg->lhp = MRISread(tmpstr);
  if (gtmseg->lhp == NULL) return (1);

  sprintf(tmpstr, "%s/%s/surf/rh.white", SUBJECTS_DIR, gtmseg->subject);
  gtmseg->rhw = MRISread(tmpstr);
  if (gtmseg->rhw == NULL) return (1);

  sprintf(tmpstr, "%s/%s/surf/rh.pial", SUBJECTS_DIR, gtmseg->subject);
  gtmseg->rhp = MRISread(tmpstr);
  if (gtmseg->rhp == NULL) return (1);

  printf("Loading annotations ");
  printf(" t = %6.4f\n", timer.seconds());
  fflush(stdout);
  if (gtmseg->wmannotfile != NULL) {
    sprintf(tmpstr, "%s/%s/label/lh.%s", SUBJECTS_DIR, gtmseg->subject, gtmseg->wmannotfile);
    err = MRISreadAnnotation(gtmseg->lhw, tmpstr);
    if (err) {
      printf("Try running mri_annotation2label --s %s --hemi lh --lobesStrict lobefilename\n", gtmseg->subject);
      return (1);
    }
    sprintf(tmpstr, "%s/%s/label/rh.%s", SUBJECTS_DIR, gtmseg->subject, gtmseg->wmannotfile);
    err = MRISreadAnnotation(gtmseg->rhw, tmpstr);
    if (err) {
      printf("Try running mri_annotation2label --s %s --hemi rh --lobesStrict lobefilename\n", gtmseg->subject);
      return (1);
    }
    MRISripUnknown(gtmseg->lhw);
    MRISripUnknown(gtmseg->rhw);
    gtmseg->lhw->ct->idbase = gtmseg->wmlhbase;
    gtmseg->rhw->ct->idbase = gtmseg->wmrhbase;
  }
  else
    printf("Not segmenting WM\n");

  sprintf(tmpstr, "%s/%s/label/lh.%s", SUBJECTS_DIR, gtmseg->subject, gtmseg->ctxannotfile);
  err = MRISreadAnnotation(gtmseg->lhp, tmpstr);
  if (err) return (1);
  sprintf(tmpstr, "%s/%s/label/rh.%s", SUBJECTS_DIR, gtmseg->subject, gtmseg->ctxannotfile);
  err = MRISreadAnnotation(gtmseg->rhp, tmpstr);
  if (err) return (1);

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
  if (!gtmseg->KeepCC) {
    printf(" Relabeling CC as WM\n");
    list[nlist] = 192;
    nlist++;
    list[nlist] = 251;
    nlist++;
    list[nlist] = 252;
    nlist++;
    list[nlist] = 253;
    nlist++;
    list[nlist] = 254;
    nlist++;
    list[nlist] = 255;
    nlist++;
  }
  if (!gtmseg->KeepHypo) {
    printf(" Relabeling any hypointensities as WM\n");
    list[nlist] = 77;
    nlist++;
    list[nlist] = 78;
    nlist++;
    list[nlist] = 79;
    nlist++;
  }
  if (nlist > 0) MRIunsegmentWM(aseg, gtmseg->lhw, gtmseg->rhw, list, nlist, NULL, aseg);

  // Upsample the segmentation
  printf("Upsampling segmentation USF = %d", gtmseg->USF);
  fflush(stdout);
  printf(" t = %6.4f\n", timer.seconds());
  fflush(stdout);
  hrseg = MRIhiresSeg(aseg, gtmseg->lhw, gtmseg->lhp, gtmseg->rhw, gtmseg->rhp, gtmseg->USF, &gtmseg->anat2seg);
  if (hrseg == NULL) return (1);
  strcpy(gtmseg->anat2seg->subject, gtmseg->subject);
  MRIfree(&aseg);

  // Label cortex (like aparc+aseg)
  printf("Beginning cortical segmentation using %s", gtmseg->ctxannotfile);
  fflush(stdout);
  printf(" t = %6.4f\n", timer.seconds());
  fflush(stdout);
  ctxseg = MRIannot2CorticalSeg(hrseg, gtmseg->lhw, gtmseg->lhp, gtmseg->rhw, gtmseg->rhp, NULL, NULL);
  MRIfree(&hrseg);

  // Label wm (like wmaparc)
  if (gtmseg->wmannotfile != NULL) {
    printf("Beginning WM segmentation using %s", gtmseg->wmannotfile);
    fflush(stdout);
    printf(" t = %6.4f\n", timer.seconds());
    fflush(stdout);
    ctxseg = MRIannot2CerebralWMSeg(ctxseg, gtmseg->lhw, gtmseg->rhw, gtmseg->dmax, NULL, ctxseg);
  }
  else
    printf("Not subsegmenting WM\n");

  gtmseg->seg = MRIreplaceList(ctxseg, gtmseg->srclist, gtmseg->targlist, gtmseg->nlist, NULL, NULL);
  if (gtmseg == NULL) return (1);
  MRIfree(&ctxseg);

  segidlist = MRIsegIdListNot0(gtmseg->seg, &nsegs, 0);
  printf("Found %d segs in the final list\n", nsegs);
  for (n = 0; n < nsegs; n++) {
    if (segidlist[n] == Left_Cerebral_Cortex) {
      printf("ERROR: MRIgtmSeg() found left cortical label\n");
      err = 1;
    }
    if (segidlist[n] == Right_Cerebral_Cortex) {
      printf("ERROR: MRIgtmSeg() found right cortical label\n");
      err = 1;
    }
    if (gtmseg->SubSegWM) {
      if (segidlist[n] == Left_Cerebral_White_Matter) {
        printf("ERROR: MRIgtmSeg() found left cerebral WM label\n");
        err = 1;
      }
      if (segidlist[n] == Right_Cerebral_White_Matter) {
        printf("ERROR: MRIgtmSeg() found right cerebral WM label\n");
        err = 1;
      }
    }
    if (segidlist[n] == WM_hypointensities) {
      printf("ERROR: MRIgtmSeg() found unlateralized WM hypointensity label\n");
      err = 1;
    }
    if (err) return (1);
  }
  gtmseg->segidlist = segidlist;
  gtmseg->nsegs = nsegs;

  printf("MRIgtmSeg() done, t = %6.4f\n", timer.seconds());
  fflush(stdout);
  return (0);
}

/*-----------------------------------------------------------------------------*/
/*
\fn int GTMdefaultSegReplacmentList(int *nReplace, int *ReplaceThis,int *WithThat)
\brief Creates a list of segids to replace with other seg ids. All
ventricular CSF segs are merged. Various ROIs are replaced with 0 (eg,
vessel).  CC subsegs are merged into a single label 192.  It is
assumed that ReplaceThis and WithThat are arrays that have already
been allocated. It is also assumed that nReplace has been initialized.
The result is that items are added to the list. */
int GTMdefaultSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat)
{
  int nlist;

  nlist = *nReplace;
  // ReplaceThis[nlist] = 1033; WithThat[nlist] = 1030; nlist++; // temppole=stg
  // ReplaceThis[nlist] = 2033; WithThat[nlist] = 2030; nlist++; // temppole=stg
  // ReplaceThis[nlist] = 1034; WithThat[nlist] = 1030; nlist++; // transtemp=stg
  // ReplaceThis[nlist] = 2034; WithThat[nlist] = 1030; nlist++; // transtemp=stg
  // ReplaceThis[nlist] = 1001; WithThat[nlist] = 1015; nlist++; // bankssts=mtg
  // ReplaceThis[nlist] = 2001; WithThat[nlist] = 2015; nlist++; // bankssts=mtg
  // ReplaceThis[nlist] = 1032; WithThat[nlist] = 1027; nlist++; // frontpole=rmf
  // ReplaceThis[nlist] = 2032; WithThat[nlist] = 2027; nlist++; // frontpole=rmf
  // ReplaceThis[nlist] = 1016; WithThat[nlist] = 1006; nlist++; // parahip=entorhinal ?
  // ReplaceThis[nlist] = 2016; WithThat[nlist] = 2006; nlist++; // parahip=entorhinal ?

  // Merge ventricular CSF into one label
  ReplaceThis[nlist] = 4;
  WithThat[nlist] = 24;
  nlist++;  // LLatVent
  ReplaceThis[nlist] = 5;
  WithThat[nlist] = 24;
  nlist++;  // LInfLatVent
  ReplaceThis[nlist] = 14;
  WithThat[nlist] = 24;
  nlist++;  // 3rd
  ReplaceThis[nlist] = 15;
  WithThat[nlist] = 24;
  nlist++;  // 4th
  ReplaceThis[nlist] = 72;
  WithThat[nlist] = 24;
  nlist++;  // 5th
  ReplaceThis[nlist] = 43;
  WithThat[nlist] = 24;
  nlist++;  // RLatVent
  ReplaceThis[nlist] = 44;
  WithThat[nlist] = 24;
  nlist++;  // RInfLatVent

  /* Merge multiple CC subsegments into one CC */
  ReplaceThis[nlist] = 251;
  WithThat[nlist] = 192;
  nlist++;
  ReplaceThis[nlist] = 252;
  WithThat[nlist] = 192;
  nlist++;
  ReplaceThis[nlist] = 253;
  WithThat[nlist] = 192;
  nlist++;
  ReplaceThis[nlist] = 254;
  WithThat[nlist] = 192;
  nlist++;
  ReplaceThis[nlist] = 255;
  WithThat[nlist] = 192;
  nlist++;

  // There should not be any cortex unknown after MRIannot2CorticalSeg()
  ReplaceThis[nlist] = 1000;
  WithThat[nlist] = 0;
  nlist++;  // cortex unknown
  ReplaceThis[nlist] = 2000;
  WithThat[nlist] = 0;
  nlist++;  // cortex unknown
  ReplaceThis[nlist] = 85;
  WithThat[nlist] = 0;
  nlist++;  // optic chiasm

  // And these?
  ReplaceThis[nlist] = 30;
  WithThat[nlist] = 0;
  nlist++;  // LVessel ?
  ReplaceThis[nlist] = 62;
  WithThat[nlist] = 0;
  nlist++;  // RVessel ?
  ReplaceThis[nlist] = 80;
  WithThat[nlist] = 0;
  nlist++;  // non-WM-hypo ?

  // Not sure about choriod plexus. Make part of CSF?
  // Note: the location as last two items makes it so that --default-seg-merge-choroid
  // in mri_gtmpvc works.
  // ReplaceThis[nlist] =   31; WithThat[nlist] =   24; nlist++; // LChoroidP ?
  // ReplaceThis[nlist] =   63; WithThat[nlist] =   24; nlist++; // RChoroidP ?

  *nReplace = nlist;
  return (0);
}
/*!
  \fn int GTMoptSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat)
  \brief Reduces the number of ROIs by merging left and right and
  combining cortical segs into a much fewer number (lobes,
  more-or-less); also combines some subcort. The idea here is to
  reduce the number of ROIs to make optimization faster without
  loosing a lot of anatomical detail. Calls GTMdefaultSegReplacmentList().
 */
int GTMoptSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat)
{
  int nlist;
  GTMdefaultSegReplacmentList(nReplace, ReplaceThis, WithThat);

  nlist = *nReplace;
  // Merge all occip (left and right) into a single seg (1021, lh calc)
  ReplaceThis[nlist] = 2021; WithThat[nlist] = 1021; nlist++;  // rh calc
  ReplaceThis[nlist] = 1011; WithThat[nlist] = 1021; nlist++;  // lh lat occip
  ReplaceThis[nlist] = 2011; WithThat[nlist] = 1021; nlist++;  // rh lat occip
  ReplaceThis[nlist] = 1005; WithThat[nlist] = 1021; nlist++;  // lh cuneus
  ReplaceThis[nlist] = 2005; WithThat[nlist] = 1021; nlist++;  // rh cuneus
  ReplaceThis[nlist] = 1013; WithThat[nlist] = 1021; nlist++;  // lh lingual
  ReplaceThis[nlist] = 2013; WithThat[nlist] = 1021; nlist++;  // rh lingual

  // Merge all temporal (left and right) into a single seg (1015, lh med temp)
  ReplaceThis[nlist] = 2015; WithThat[nlist] = 1015; nlist++;  // rh medtemp
  ReplaceThis[nlist] = 1009; WithThat[nlist] = 1015; nlist++;  // lh inftemp
  ReplaceThis[nlist] = 2009; WithThat[nlist] = 1015; nlist++;  // rh inftemp
  ReplaceThis[nlist] = 1030; WithThat[nlist] = 1015; nlist++;  // lh suptemp
  ReplaceThis[nlist] = 2030; WithThat[nlist] = 1015; nlist++;  // rh suptemp
  ReplaceThis[nlist] = 1033; WithThat[nlist] = 1015; nlist++;  // lh temppole
  ReplaceThis[nlist] = 2033; WithThat[nlist] = 1015; nlist++;  // rh temppole
  ReplaceThis[nlist] = 1034; WithThat[nlist] = 1015; nlist++;  // lh transtemp
  ReplaceThis[nlist] = 2034; WithThat[nlist] = 1015; nlist++;  // rh transtemp
  ReplaceThis[nlist] = 1007; WithThat[nlist] = 1015; nlist++;  // lh fusi
  ReplaceThis[nlist] = 2007; WithThat[nlist] = 1015; nlist++;  // rh fusi
  ReplaceThis[nlist] = 1001; WithThat[nlist] = 1015; nlist++;  // lh banks sts
  ReplaceThis[nlist] = 2001; WithThat[nlist] = 1015; nlist++;  // rh banks sts
  ReplaceThis[nlist] = 1016; WithThat[nlist] = 1015; nlist++;  // lh parahip
  ReplaceThis[nlist] = 2016; WithThat[nlist] = 1015; nlist++;  // rh parahip
  ReplaceThis[nlist] = 1006; WithThat[nlist] = 1015; nlist++;  // lh ento
  ReplaceThis[nlist] = 2006; WithThat[nlist] = 1015; nlist++;  // rh ento

  // Merge all parietal (left and right) into a single seg (1029, lh sup parietal)
  ReplaceThis[nlist] = 2029; WithThat[nlist] = 1029; nlist++;  // rh sup par
  ReplaceThis[nlist] = 1008; WithThat[nlist] = 1029; nlist++;  // lh inf par
  ReplaceThis[nlist] = 2008; WithThat[nlist] = 1029; nlist++;  // rh inf par
  ReplaceThis[nlist] = 1025; WithThat[nlist] = 1029; nlist++;  // lh precuen
  ReplaceThis[nlist] = 2025; WithThat[nlist] = 1029; nlist++;  // rh precuen
  ReplaceThis[nlist] = 1031; WithThat[nlist] = 1029; nlist++;  // lh sup marg
  ReplaceThis[nlist] = 2031; WithThat[nlist] = 1029; nlist++;  // rh sup marg

  // Merge all cingulate (left and right) into a single seg (1002, lh caud ant cing)
  ReplaceThis[nlist] = 2002; WithThat[nlist] = 1002; nlist++;  // rh caudantcing
  ReplaceThis[nlist] = 1026; WithThat[nlist] = 1002; nlist++;  // lh rostantcing
  ReplaceThis[nlist] = 2026; WithThat[nlist] = 1002; nlist++;  // rh rostantcing
  ReplaceThis[nlist] = 1023; WithThat[nlist] = 1002; nlist++;  // lh postcing
  ReplaceThis[nlist] = 2023; WithThat[nlist] = 1002; nlist++;  // rh postcing
  ReplaceThis[nlist] = 1010; WithThat[nlist] = 1002; nlist++;  // lh isthmuscing
  ReplaceThis[nlist] = 2010; WithThat[nlist] = 1002; nlist++;  // rh isthmuscing

  // Merge all central (left and right) into a single seg (1022, lh post cent)
  ReplaceThis[nlist] = 2022; WithThat[nlist] = 1022; nlist++;  // rh post cent
  ReplaceThis[nlist] = 1024; WithThat[nlist] = 1022; nlist++;  // lh pre cent
  ReplaceThis[nlist] = 2024; WithThat[nlist] = 1022; nlist++;  // rh pre cent
  ReplaceThis[nlist] = 1017; WithThat[nlist] = 1022; nlist++;  // lh para cent
  ReplaceThis[nlist] = 2017; WithThat[nlist] = 1022; nlist++;  // rh para cent

  // Merge all supfront (left and right) into a single seg (1028, lh sup front)
  ReplaceThis[nlist] = 2028; WithThat[nlist] = 1028; nlist++;  // rh sup front

  // Merge all midfront (left and right) into a single seg (1003, lh caudmid front)
  ReplaceThis[nlist] = 2003; WithThat[nlist] = 1003; nlist++;  // rh caudmid front
  ReplaceThis[nlist] = 1027; WithThat[nlist] = 1003; nlist++;  // lh rostmid front
  ReplaceThis[nlist] = 2027; WithThat[nlist] = 1003; nlist++;  // rh rostmid front

  // Merge all pars (left and right) into a single seg (1020, lh pars tri)
  ReplaceThis[nlist] = 2020; WithThat[nlist] = 1020; nlist++;  // rh pars tri
  ReplaceThis[nlist] = 1018; WithThat[nlist] = 1020; nlist++;  // lh pars oper
  ReplaceThis[nlist] = 2018; WithThat[nlist] = 1020; nlist++;  // rh pars oper
  ReplaceThis[nlist] = 1019; WithThat[nlist] = 1020; nlist++;  // lh pars orbit
  ReplaceThis[nlist] = 2019; WithThat[nlist] = 1020; nlist++;  // rh pars orbit

  // Merge all orb front (left and right) into a single seg (1012, lh lat orb front)
  ReplaceThis[nlist] = 2012; WithThat[nlist] = 1012; nlist++;  // rh lat orb front
  ReplaceThis[nlist] = 1014; WithThat[nlist] = 1012; nlist++;  // lh med orb front
  ReplaceThis[nlist] = 2014; WithThat[nlist] = 1012; nlist++;  // rh med orb front
  ReplaceThis[nlist] = 1032; WithThat[nlist] = 1012; nlist++;  // lh frontal pole
  ReplaceThis[nlist] = 2032; WithThat[nlist] = 1012; nlist++;  // rh frontal pole

  // Merge left and right insula into left ins 1035
  ReplaceThis[nlist] = 2035; WithThat[nlist] = 1035; nlist++;  // rh ins

  // Merge CP with CSF
  ReplaceThis[nlist] =   31; WithThat[nlist] =   24; nlist++; // LChoroidP
  ReplaceThis[nlist] =   63; WithThat[nlist] =   24; nlist++; // RChoroidP

  // Merge Cerebellar GM into left cblum 8
  ReplaceThis[nlist] =   47; WithThat[nlist] =   8; nlist++; // L Cblum GM
  ReplaceThis[nlist] =  172; WithThat[nlist] =   8; nlist++; // Vermis

  // Merge L and R Cerebellar WGM into left cblum wm 7
  ReplaceThis[nlist] =   46; WithThat[nlist] =   7; nlist++; // R Cblum WM

  // Merge Hippo and Amyg into LHip 17
  ReplaceThis[nlist] =   53; WithThat[nlist] =   17; nlist++; // R Hippo
  ReplaceThis[nlist] =   18; WithThat[nlist] =   17; nlist++; // L Amyg
  ReplaceThis[nlist] =   54; WithThat[nlist] =   17; nlist++; // R Amyg

  // Merge Putamen and Pallidum into L Put 12
  ReplaceThis[nlist] =   51; WithThat[nlist] =   12; nlist++; // R Put
  ReplaceThis[nlist] =   13; WithThat[nlist] =   12; nlist++; // L Pal
  ReplaceThis[nlist] =   52; WithThat[nlist] =   12; nlist++; // R Pal

  // Merge Caudate and NucAcc into L Caud 11
  ReplaceThis[nlist] =   50; WithThat[nlist] =   11; nlist++; // R Put
  ReplaceThis[nlist] =   26; WithThat[nlist] =   11; nlist++; // L NAcc
  ReplaceThis[nlist] =   58; WithThat[nlist] =   11; nlist++; // R NAcc

  // Merge L and R Thal into L Thal 10
  ReplaceThis[nlist] =   49; WithThat[nlist] =   10; nlist++; // R Thal

  // Merge L and R VDC into L VDC 28
  ReplaceThis[nlist] =   60; WithThat[nlist] =   28; nlist++; // R VDC

  // Merge L and R Cerebral WM into L WM 2 (and WMSAs too)
  ReplaceThis[nlist] =   41; WithThat[nlist] =   2; nlist++; // R WM
  ReplaceThis[nlist] =   77; WithThat[nlist] =   2; nlist++; // WMSA
  ReplaceThis[nlist] =   78; WithThat[nlist] =   2; nlist++; // L WMSA
  ReplaceThis[nlist] =   79; WithThat[nlist] =   2; nlist++; // R WMSA

  // Merge ExtraCerebral CSF with CSF 24
  ReplaceThis[nlist] =   257; WithThat[nlist] =   24; nlist++; // XCSF

  // Merge Head segs
  ReplaceThis[nlist] =   130; WithThat[nlist] =   258; nlist++; // Air
  ReplaceThis[nlist] =   165; WithThat[nlist] =   258; nlist++; // Skull

  *nReplace = nlist;
  return (0);
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
  int nsegs, *segidlist, *segidhit, n, m, hit, segid, err;
  COLOR_TABLE *ct, *ct0;
  CTE *cte, *cte0;

  segidlist = MRIsegIdListNot0(gtmseg->seg, &nsegs, 0);
  segidhit = (int *)calloc(nsegs, sizeof(int));

  // FILE *fp = fopen("segidlist.dat","w");
  // for(m=0; m < nsegs; m++) fprintf(fp,"%4d\n",segidlist[m]);
  // fclose(fp);

  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = segidlist[nsegs - 1] + 1;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));

  ct->ctabTissueType = CTABdeepCopy(ctSubCort->ctabTissueType);
  strcpy(ct->TissueTypeSchema,ctSubCort->TissueTypeSchema);

  // Add an entry for unknown
  segid = 0;
  ct->entries[segid] = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
  cte = ct->entries[segid];
  cte->ri = 0;
  cte->gi = 0;
  cte->bi = 0;
  cte->ai = 255;
  cte->rf = 0;
  cte->gf = 0;
  cte->bf = 0;
  cte->af = 255;
  cte->TissueType = 0;
  sprintf(cte->name, "Unknown");

  fflush(stdout);
  ct0 = gtmseg->lhp->ct;
  // Go thru each entry in the left cortical annotation
  for (n = 0; n < ct0->nentries; n++) {
    cte0 = ct0->entries[n];
    if (cte0 == NULL) continue;
    // compute the segid
    segid = gtmseg->ctxlhbase + n;
    // check whether it is in the segidlist
    hit = 0;
    for (m = 0; m < nsegs; m++) {
      if (segid == segidlist[m]) {
        hit = 1;
        break;
      }
    }
    if (hit == 0) {
      /* It does not have to be there. Eg, CC will be removed. Also
         there are a bunch of entries that are non-null but not
         actually represented in the annot (they have the name
         clusterXX. */
      if (Gdiag_no > 1) printf("INFO: GTMSEGctab(): lhp n=%d, segid=%d, %s not in segidlist\n", n, segid, cte0->name);
      continue;
    }
    segidhit[m]++;  // keep track of number of times this seg id represented
    ct->entries[segid] = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
    cte = ct->entries[segid];
    memcpy(cte, cte0, sizeof(CTE));
    // The memcpy on the previous line makes converting to std::string rather tricky
    int err = snprintf(cte->name, STRLEN-1, "ctx-lh-%s", cte0->name);  // new name reflects cortex and hemi
    if( err >= STRLEN-1 ) {
      std::cerr << __FUNCTION__ << ": Truncation prepending ctx-lh-" << std::endl;
    }
    cte->TissueType = 1;
  }

  // Do the same thing for the right hemi
  ct0 = gtmseg->rhp->ct;
  for (n = 0; n < ct0->nentries; n++) {
    cte0 = ct0->entries[n];
    if (cte0 == NULL) continue;
    segid = gtmseg->ctxrhbase + n;
    hit = 0;
    for (m = 0; m < nsegs; m++) {
      if (segid == segidlist[m]) {
        hit = 1;
        break;
      }
    }
    if (hit == 0) {
      if (Gdiag_no > 1) printf("INFO: GTMSEGctab(): rhp n=%d, segid=%d, %s not in segidlist\n", n, segid, cte0->name);
      continue;
    }
    segidhit[m]++;
    ct->entries[segid] = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
    cte = ct->entries[segid];
    memcpy(cte, cte0, sizeof(CTE));
    // The preceding memcpy is an issue with trying to change to std::string
    int err = snprintf(cte->name, STRLEN-1, "ctx-rh-%s", cte0->name);
    if( err >= STRLEN-1 ) {
      std::cerr << __FUNCTION__ << ": Truncation prepending ctx-rh-" << std::endl;
    }
    cte->TissueType = 1;
  }

  // Do the same thing if subsegmenting WM based on proximity to cortex
  if (gtmseg->SubSegWM) {
    ct0 = gtmseg->lhw->ct;
    for (n = 0; n < ct0->nentries; n++) {
      cte0 = ct0->entries[n];
      if (cte0 == NULL) continue;
      segid = gtmseg->wmlhbase + n;
      hit = 0;
      for (m = 0; m < nsegs; m++) {
        if (segid == segidlist[m]) {
          hit = 1;
          break;
        }
      }
      if (hit == 0) {
        if (Gdiag_no > 1) printf("INFO: GTMSEGctab(): lhw n=%d, segid=%d, %s not in segidlist\n", n, segid, cte0->name);
        continue;
      }
      segidhit[m]++;
      ct->entries[segid] = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
      cte = ct->entries[segid];
      memcpy(cte, cte0, sizeof(CTE));
      // The preceding memcpy is an issue with trying to change to std::string
      int err = snprintf(cte->name, STRLEN-1, "wm-lh-%s", cte0->name);
      if( err >= STRLEN-1 ) {
	std::cerr << __FUNCTION__ << ": Truncation prepending wm-lh-" << std::endl;
      }
      cte->TissueType = 3;
    }
    ct0 = gtmseg->rhw->ct;
    for (n = 0; n < ct0->nentries; n++) {
      cte0 = ct0->entries[n];
      if (cte0 == NULL) continue;
      segid = gtmseg->wmrhbase + n;
      hit = 0;
      for (m = 0; m < nsegs; m++) {
        if (segid == segidlist[m]) {
          hit = 1;
          break;
        }
      }
      if (hit == 0) {
        if (Gdiag_no > 1) printf("INFO: GTMSEGctab(): lhp n=%d, segid=%d, %s not in segidlist\n", n, segid, cte0->name);
        continue;
      }
      segidhit[m]++;
      ct->entries[segid] = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
      cte = ct->entries[segid];
      memcpy(cte, cte0, sizeof(CTE));
      // The preceding memcpy is an issue with trying to change to std::string
      int err = snprintf(cte->name, STRLEN-1, "wm-rh-%s", cte0->name);
      if( err >= STRLEN-1 ) {
	std::cerr << __FUNCTION__ << ": Truncation prepending wm-rh-" << std::endl;
      }
      cte->TissueType = 3;
    }
  }

  // Add entries for subcortical regions
  ct0 = ctSubCort;
  for (m = 0; m < nsegs; m++) {
    if (segidhit[m] != 0) continue;  // skip of already hit
    segid = segidlist[m];
    cte0 = ct0->entries[segid];
    if (cte0 == NULL) {
      printf("ERROR: cannot find match for subcortical segid %d\n", segid);
      return (NULL);
    }
    ct->entries[segid] = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
    cte = ct->entries[segid];
    memcpy(cte, cte0, sizeof(CTE));
    segidhit[m]++;
  }

  // Check the final ctab
  err = 0;
  for (m = 0; m < nsegs; m++) {
    if (segidhit[m] > 1) {
      printf("ERROR: segid %4d is represented multiple (%d) times \n", segidhit[m], segidlist[m]);
      err++;
    }
  }
  if (err) {
    printf("ERROR: found %d segmentations with multiple representations\n", err);
    return (NULL);
  }
  for (m = 0; m < nsegs; m++) {
    segid = segidlist[m];
    cte = ct->entries[segid];
    if (cte->TissueType == -1) printf("WARNING: segid %4d %s tissue type is not set\n", segid, cte->name);
  }

  free(segidlist);
  free(segidhit);

  return (ct);
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
  int n, nvox, segid, f;
  double vrf;

  *vrfmean = 0;
  *vrfmax = 0;
  *vrfmin = 0;
  for (n = 0; n < gtm->iXtX->rows; n++) {
    vrf = (double)1.0 / gtm->iXtX->rptr[n + 1][n + 1];
    if (n == 0) {
      *vrfmax = vrf;
      *vrfmin = vrf;
    }
    if (*vrfmax < vrf) *vrfmax = vrf;
    if (*vrfmin > vrf) *vrfmin = vrf;
    *vrfmean += vrf;
  }
  *vrfmean /= gtm->iXtX->rows;

  if (gtm->betavar) MatrixFree(&gtm->betavar);
  gtm->betavar = MatrixAlloc(gtm->iXtX->rows, gtm->beta->cols, MATRIX_REAL);

  if (gtm->vrf) MatrixFree(&gtm->vrf);
  gtm->vrf = MatrixAlloc(gtm->iXtX->rows, 1, MATRIX_REAL);

  if (gtm->nvox) MatrixFree(&gtm->nvox);
  gtm->nvox = MatrixAlloc(gtm->iXtX->rows, 1, MATRIX_REAL);

  for (n = 0; n < gtm->iXtX->rows; n++) {
    segid = gtm->segidlist[n];
    nvox = MRIcountMatches(gtm->gtmseg, segid, 0, gtm->mask);
    vrf = (double)1.0 / gtm->iXtX->rptr[n + 1][n + 1];
    gtm->vrf->rptr[n + 1][1] = vrf;
    gtm->nvox->rptr[n + 1][1] = nvox;
    for (f = 0; f < gtm->beta->cols; f++) gtm->betavar->rptr[n + 1][f + 1] = gtm->rvar->rptr[1][f + 1] / vrf;
  }

  return (0);
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

  fp = fopen(fname, "w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", fname);
    return (1);
  }

  for (n = 0; n < gtm->iXtX->rows; n++) {
    segid = gtm->segidlist[n];
    vrf = gtm->vrf->rptr[n + 1][1];
    nvox = gtm->nvox->rptr[n + 1][1];
    cte = gtm->ctGTMSeg->entries[segid];
    fprintf(fp,
            "%3d %4d %-31s %-13s %6d %8.3f",
            n + 1,
            segid,
            cte->name,
            ttctab->entries[cte->TissueType]->name,
            nvox,
            vrf);
    // printf("%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
    // ttctab->entries[cte->TissueType]->name,nvox,vrf);
    if (gtm->beta) fprintf(fp, "   %10.3f", gtm->beta->rptr[n + 1][1]);
    if (gtm->segrvar) fprintf(fp, "   %10.4f", sqrt(gtm->segrvar->rptr[n + 1][1]));
    fprintf(fp, "\n");
  }
  fclose(fp);
  fflush(stdout);

  return (0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMglobalStats(GTM *gtm)
  \brief Computes global mean for GM, GM+WM, and GM+WM+CSF, each
  weighted according to number of voxels in each ROI.
*/
int GTMglobalStats(GTM *gtm)
{
  int nthseg, segid, tt, f, ngm, ngmwm, ngmwmcsf;
  double v;

  gtm->glob_gm = MatrixAlloc(gtm->nframes, 1, MATRIX_REAL);
  gtm->glob_gmwm = MatrixAlloc(gtm->nframes, 1, MATRIX_REAL);
  gtm->glob_gmwmcsf = MatrixAlloc(gtm->nframes, 1, MATRIX_REAL);
  ngm = 0;
  ngmwm = 0;
  ngmwmcsf = 0;

  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    segid = gtm->segidlist[nthseg];
    tt = gtm->ctGTMSeg->entries[segid]->TissueType;
    if (tt == 1 || tt == 2) {
      ngm += gtm->nvox->rptr[nthseg + 1][1];
      for (f = 0; f < gtm->nframes; f++) {
        v = gtm->beta->rptr[nthseg + 1][f + 1] * gtm->nvox->rptr[nthseg + 1][1];
        gtm->glob_gm->rptr[f + 1][1] += v;
      }
    }
    if (tt == 1 || tt == 2 || tt == 3) {
      ngmwm += gtm->nvox->rptr[nthseg + 1][1];
      for (f = 0; f < gtm->nframes; f++) {
        v = gtm->beta->rptr[nthseg + 1][f + 1] * gtm->nvox->rptr[nthseg + 1][1];
        gtm->glob_gmwm->rptr[f + 1][1] += v;
      }
    }
    if (tt == 1 || tt == 2 || tt == 3 || tt == 4) {
      ngmwmcsf += gtm->nvox->rptr[nthseg + 1][1];
      for (f = 0; f < gtm->nframes; f++) {
        v = gtm->beta->rptr[nthseg + 1][f + 1] * gtm->nvox->rptr[nthseg + 1][1];
        gtm->glob_gmwmcsf->rptr[f + 1][1] += v;
      }
    }
  }

  printf("Global S: %f %f %f\n", gtm->glob_gm->rptr[1][1], gtm->glob_gmwm->rptr[1][1], gtm->glob_gmwmcsf->rptr[1][1]);
  printf("Global N: %d %d %d\n", ngm, ngmwm, ngmwmcsf);

  for (f = 0; f < gtm->nframes; f++) {
    gtm->glob_gm->rptr[f + 1][1] /= ngm;
    gtm->glob_gmwm->rptr[f + 1][1] /= ngmwm;
    gtm->glob_gmwmcsf->rptr[f + 1][1] /= ngmwmcsf;
  }

  return (0);
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
  // MRIfree(&gtm->gtmseg);
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
  *pGTM = NULL;
  return (0);
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
  return (0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsetNMask(GTM *gtm)
  \brief Computes the number of voxels in the mask. If the mask is
  NULL, then just computes the number of voxels in the input.
*/
int GTMsetNMask(GTM *gtm)
{
  if (gtm->mask)
    gtm->nmask = MRIcountAboveThreshold(gtm->mask, 0.5);
  else
    gtm->nmask = gtm->yvol->width * gtm->yvol->height * gtm->yvol->depth;
  return (0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMpsfStd(GTM *gtm)
  \brief Convert the PSF {crs}FWHM to a standard deviation.
*/
int GTMpsfStd(GTM *gtm)
{
  gtm->cStd = gtm->cFWHM / sqrt(log(256.0));
  gtm->rStd = gtm->rFWHM / sqrt(log(256.0));
  gtm->sStd = gtm->sFWHM / sqrt(log(256.0));
  return (0);
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
  return (0);
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
  gtm = (GTM *)calloc(sizeof(GTM), 1);
  gtm->PadThresh = .0001;
  return (gtm);
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
  maxFWHM = MAX(gtm->cFWHM / gtm->yvol->xsize, MAX(gtm->rFWHM / gtm->yvol->ysize, gtm->sFWHM / gtm->yvol->zsize));
  if (maxFWHM > 0) {
    // The case where psf=0, just correcting for volume fraction
    maxStd = maxFWHM * sqrt(log(256.0));
    gtm->nPad = ceil(sqrt(-log(gtm->PadThresh * maxStd * sqrt(2 * M_PI)) * 2 * maxStd));
    printf("maxFWHM = %g (voxels), PadThresh=%g, nPad=%d\n", maxFWHM, gtm->PadThresh, gtm->nPad);
  }
  else
    gtm->nPad = 1;
  return (0);
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
  int n, f;
  double sum;

  if (gtm->X == NULL) {
    printf("ERROR: GTMsolve(): must build design matrix first\n");
    exit(1);
  }

  if (!gtm->Optimizing) printf("Computing  XtX ... ");
  fflush(stdout);
  Timer timer;
  gtm->XtX = MatrixMtM(gtm->X, gtm->XtX);
  if (!gtm->Optimizing) printf(" %4.1f sec\n", timer.seconds());
  fflush(stdout);

  gtm->iXtX = MatrixInverse(gtm->XtX, gtm->iXtX);
  if (gtm->iXtX == NULL) {
    if (gtm->Optimizing) return (1);
    gtm->XtXcond = MatrixConditionNumber(gtm->XtX);
    printf("ERROR: matrix cannot be inverted, cond=%g\n", gtm->XtXcond);
    return (1);
  }
  gtm->Xty = MatrixAtB(gtm->X, gtm->y, gtm->Xty);
  gtm->beta = MatrixMultiplyD(gtm->iXtX, gtm->Xty, gtm->beta);
  if (gtm->rescale) GTMrescale(gtm);
  GTMrefTAC(gtm);
  if (gtm->DoSteadyState) GTMsteadyState(gtm);

  gtm->yhat = MatrixMultiplyD(gtm->X, gtm->beta, gtm->yhat);
  gtm->res = MatrixSubtract(gtm->y, gtm->yhat, gtm->res);
  gtm->dof = gtm->X->rows - gtm->X->cols;
  if (gtm->rvar == NULL) gtm->rvar = MatrixAlloc(1, gtm->res->cols, MATRIX_REAL);
  if (gtm->rvarUnscaled == NULL) gtm->rvarUnscaled = MatrixAlloc(1, gtm->res->cols, MATRIX_REAL);
  for (f = 0; f < gtm->res->cols; f++) {
    sum = 0;
    for (n = 0; n < gtm->res->rows; n++) sum += ((double)gtm->res->rptr[n + 1][f + 1] * gtm->res->rptr[n + 1][f + 1]);
    gtm->rvar->rptr[1][f + 1] = sum / gtm->dof;
    if (gtm->rescale)
      gtm->rvarUnscaled->rptr[1][f + 1] = gtm->rvar->rptr[1][f + 1] / (gtm->scale * gtm->scale);
    else
      gtm->rvarUnscaled->rptr[1][f + 1] = gtm->rvar->rptr[1][f + 1];
  }
  gtm->kurtosis = MatrixKurtosis(gtm->res, gtm->kurtosis);
  gtm->skew = MatrixSkew(gtm->res, gtm->skew);

  return (0);
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
  int k, c, r, s, f;

  if (vol == NULL) {
    vol = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT, m->cols);
    if (vol == NULL) return (NULL);
    MRIcopyHeader(gtm->yvol, vol);
    MRIcopyPulseParameters(gtm->yvol, vol);
  }
  if (MRIdimMismatch(gtm->yvol, vol, 0)) {
    printf("ERROR: GTMmat2vol() dim mismatch\n");
    return (NULL);
  }

  k = 0;
  for (s = 0; s < gtm->yvol->depth; s++) {
    for (c = 0; c < gtm->yvol->width; c++) {
      for (r = 0; r < gtm->yvol->height; r++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        for (f = 0; f < m->cols; f++) MRIsetVoxVal(vol, c, r, s, f, m->rptr[k + 1][f + 1]);
        k++;
      }
    }
  }
  return (vol);
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
  int k, c, r, s, f;

  if (m == NULL) {
    m = MatrixAlloc(gtm->nmask, vol->nframes, MATRIX_REAL);
    if (m == NULL) {
      printf("ERROR: GTMvol2mat(): could not alloc matrix %d %d\n", gtm->nmask, vol->nframes);
      return (NULL);
    }
  }
  if (m->rows != gtm->nmask) {
    printf("ERROR: GTMvol2mat(): row mismatch %d %d\n", m->rows, gtm->nmask);
    return (NULL);
  }
  if (m->cols != vol->nframes) {
    printf("ERROR: GTMvol2mat(): col mismatch %d %d\n", m->cols, vol->nframes);
    return (NULL);
  }

  k = 0;
  for (s = 0; s < vol->depth; s++) {  // crs order is important here!
    for (c = 0; c < vol->width; c++) {
      for (r = 0; r < vol->height; r++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        for (f = 0; f < vol->nframes; f++) m->rptr[k + 1][f + 1] = MRIgetVoxVal(vol, c, r, s, f);
        k++;
      }
    }
  }
  return (m);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMsegrvar(GTM *gtm)
  \brief Computes the residual variance in each segmentation. Not perfect
  because the resdiual has spill-over. Hopefully it is meaningful for something.
*/
int GTMsegrvar(GTM *gtm)
{
  int c, r, s, f, k;
  int nthseg = 0, segid;
  double v;

  gtm->segrvar = MatrixAlloc(gtm->nsegs, gtm->beta->cols, MATRIX_REAL);
  gtm->nperseg = (int *)calloc(sizeof(int), gtm->nsegs);

  k = 0;
  for (s = 0; s < gtm->yvol->depth; s++) {
    for (c = 0; c < gtm->yvol->width; c++) {
      for (r = 0; r < gtm->yvol->height; r++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        segid = MRIgetVoxVal(gtm->gtmseg, c, r, s, 0);
        if (segid != 0) {
          for (nthseg = 0; nthseg < gtm->nsegs; nthseg++)
            if (gtm->segidlist[nthseg] == segid) break;
          gtm->nperseg[nthseg]++;
        }
        for (f = 0; f < gtm->beta->cols; f++) {
          v = gtm->res->rptr[k + 1][f + 1];
          if (segid != 0) gtm->segrvar->rptr[nthseg + 1][f + 1] += v * v;
        }
        k++;
      }  // r
    }    // c
  }      // s

  for (f = 0; f < gtm->beta->cols; f++) {
    for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      v = gtm->segrvar->rptr[nthseg + 1][f + 1];
      gtm->segrvar->rptr[nthseg + 1][f + 1] = v / gtm->nperseg[nthseg];
    }
  }
  return (0);
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
  int c, r, s, f, nthseg, segid;
  double val, v, vhat0, vhat, v2;
  LTA *lta;
  MATRIX *nhits;
  MRI *yframe = NULL;
  MRI_REGION *region = NULL;
  MRI *yseg = NULL;      // source volume trilin resampled to seg space (used with RBV)
  MRI *yhat0seg = NULL;  // unsmoothed yhat created in seg space (used with RBV)
  MRI *yhatseg = NULL;   // smoothed yhat in seg space (used with RBV)

  if (gtm->rbv) MRIfree(&gtm->rbv);

  Timer mytimer;
  PrintMemUsage(stdout);

  if (gtm->mask_rbv_to_brain) {
    // Reduce RBV to a tight bounding box around the brain by excluding anything
    // with "head" tissue type (hardcoded=5). This can greatly reduce
    // memory requirements. The RAS space is still that of the rbvseg
    // (and so also that of the conformed anat) so no new registration is necessary
    printf("   masking RBV to brain\n");
    if (Gdiag_no > 0) PrintMemUsage(stdout);
    int n, nReplace, ReplaceThis[1000], WithThat[1000];
    nReplace = 0;
    for (n = 0; n < gtm->ctGTMSeg->nentries; n++) {
      if (gtm->ctGTMSeg->entries[n] == NULL) continue;
      if (gtm->ctGTMSeg->entries[n]->TissueType != 5) continue;  // should not hard-code
      ReplaceThis[nReplace] = n;
      WithThat[nReplace] = 0;
      nReplace++;
    }
    printf("  replacing head voxels with 0\n");
    gtm->rbvsegmasked = MRIreplaceList(gtm->rbvseg, ReplaceThis, WithThat, nReplace, NULL, NULL);
    printf("  computing bounding box  %d %d %d ", gtm->rbvseg->width, gtm->rbvseg->height, gtm->rbvseg->depth);
    region = REGIONgetBoundingBox(gtm->rbvsegmasked, 10);
    REGIONprint(stdout, region);
    MRI *tmp = MRIextractRegion(gtm->rbvsegmasked, NULL, region);
    MRIfree(&gtm->rbvsegmasked);
    gtm->rbvsegmasked = tmp;
  }
  else
    gtm->rbvsegmasked = gtm->rbvseg;
  gtm->anat2rbv = TransformRegDat2LTA(gtm->anatconf, gtm->rbvsegmasked, NULL);
  strcpy(gtm->anat2rbv->subject, gtm->anat2pet->subject);

  // printf("writing gtm->rbvsegmasked\n");
  // MRIwrite(gtm->rbvsegmasked,"segbrain.mgh");

  printf("   Allocating RBV nvox=%d\n",
         gtm->rbvsegmasked->width * gtm->rbvsegmasked->height * gtm->rbvsegmasked->depth * gtm->nframes);
  fflush(stdout);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  gtm->rbv = MRIallocSequence(
      gtm->rbvsegmasked->width, gtm->rbvsegmasked->height, gtm->rbvsegmasked->depth, MRI_FLOAT, gtm->nframes);
  if (gtm->rbv == NULL) {
    printf("ERROR: GTMrbv() could not alloc rbv\n");
    return (1);
  }
  MRIcopyHeader(gtm->rbvsegmasked, gtm->rbv);
  MRIcopyPulseParameters(gtm->yvol, gtm->rbv);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  // if(gtm->mask_rbv_to_brain) MRIfree(&segbrain);

  // must be in rbvseg space
  yseg = MRIallocSequence(gtm->rbvseg->width, gtm->rbvseg->height, gtm->rbvseg->depth, MRI_FLOAT, 1);
  if (yseg == NULL) {
    printf("ERROR: GTMrbv() could not alloc yseg\n");
    return (1);
  }
  MRIcopyHeader(gtm->rbvseg, yseg);
  MRIcopyPulseParameters(gtm->yvol, yseg);

  // Keep track of segmeans in RBV for QA
  gtm->rbvsegmean = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->nframes);

  printf("RBV looping over %d frames, t = %4.2f min \n", gtm->nframes, mytimer.minutes());
  fflush(stdout);
  for (f = 0; f < gtm->nframes; f++) {
    printf("   f=%d t = %4.2f\n", f, mytimer.minutes());
    fflush(stdout);
    yframe = fMRIframe(gtm->yvol, f, yframe);

    printf("   Synthesizing unsmoothed input in seg space %4.2f \n", mytimer.minutes());
    fflush(stdout);
    if (Gdiag_no > 0) PrintMemUsage(stdout);
    yhat0seg = GTMsegSynth(gtm, f, yhat0seg);
    if (yhat0seg == NULL) {
      printf("ERROR: GTMrbv() could not synthesize yhat0seg\n");
      return (1);
    }

    printf("   Smoothing synthesized in seg space %4.2f \n", mytimer.minutes());
    fflush(stdout);
    if (Gdiag_no > 0) PrintMemUsage(stdout);
    yhatseg = MRIgaussianSmoothNI(yhat0seg, gtm->cStd, gtm->rStd, gtm->sStd, yhatseg);
    if (yhatseg == NULL) {
      printf("ERROR: GTMrbv() could not smooth yhatseg\n");
      return (1);
    }

    printf("   Sampling input to seg space with trilin %4.2f \n", mytimer.minutes());
    fflush(stdout);
    if (Gdiag_no > 0) PrintMemUsage(stdout);
    lta = LTAcopy(gtm->rbvseg2pet, NULL);
    LTAchangeType(lta, LINEAR_VOX_TO_VOX);
    MRIvol2Vol(yframe, yseg, (lta->xforms[0].m_L), SAMPLE_TRILINEAR, 0.0);
    LTAfree(&lta);

    printf("   Computing RBV %4.2f \n", mytimer.minutes());
    fflush(stdout);
    if (Gdiag_no > 0) PrintMemUsage(stdout);
    nhits = MatrixAlloc(gtm->beta->rows, 1, MATRIX_REAL);
    for (c = 0; c < gtm->rbvseg->width; c++) {  // crs order not important
      for (r = 0; r < gtm->rbvseg->height; r++) {
        for (s = 0; s < gtm->rbvseg->depth; s++) {
          segid = MRIgetVoxVal(gtm->rbvseg, c, r, s, 0);
          if (segid < 0.5) continue;

          if (gtm->mask_rbv_to_brain) {
            if (c < region->x || c >= region->x + region->dx) continue;
            if (r < region->y || r >= region->y + region->dy) continue;
            if (s < region->z || s >= region->z + region->dz) continue;
          }

          for (nthseg = 0; nthseg < gtm->nsegs; nthseg++)
            if (gtm->segidlist[nthseg] == segid) break;
          if (f == 0) nhits->rptr[nthseg + 1][1]++;

          v = MRIgetVoxVal(yseg, c, r, s, 0);
          vhat0 = MRIgetVoxVal(yhat0seg, c, r, s, 0);
          vhat = MRIgetVoxVal(yhatseg, c, r, s, 0);
          val = v * vhat0 / (vhat + FLT_EPSILON);  // RBV equation
          if (gtm->mask_rbv_to_brain)
            MRIsetVoxVal(gtm->rbv, c - region->x, r - region->y, s - region->z, f, val);
          else
            MRIsetVoxVal(gtm->rbv, c, r, s, f, val);

          // track seg means for QA. Head Segs won't reflect QA if masking
          v2 = MRIgetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f);
          MRIsetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f, v2 + val);
        }
      }
    }
  }
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  printf("  t = %4.2f min\n", mytimer.minutes());
  fflush(stdout);
  MRIfree(&yseg);
  MRIfree(&yhat0seg);
  MRIfree(&yhatseg);
  MRIfree(&yframe);
  if (gtm->mask_rbv_to_brain) free(region);

  // track seg means for QA
  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    for (f = 0; f < gtm->nframes; f++) {
      val = MRIgetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f) / nhits->rptr[nthseg + 1][1];
      MRIsetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f, val);
    }
  }
  MatrixFree(&nhits);

  PrintMemUsage(stdout);
  printf("  RBV took %4.2f min\n", mytimer.minutes());

  return (0);
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
  LTA *new2seg;  // *seg2new; //

  if (gtm->rbvsegres <= 0 || gtm->rbvsegres == gtm->anatseg->xsize) {
    printf("Not changing res of seg for RBV\n");
    gtm->rbvseg = gtm->anatseg;
    gtm->rbvseg2pet = gtm->seg2pet;
    return (0);
  }

  res = gtm->rbvsegres;
  printf("Changing res of seg for RBV to %g\n", res);
  gtm->rbvseg = MRIchangeSegRes(gtm->anatseg, res, res, res, gtm->ctGTMSeg, &new2seg);
  if (gtm->rbvseg == NULL) return (1);
  ;
  gtm->rbvseg2pet = LTAconcat2(new2seg, gtm->seg2pet, 1);
  if (gtm->rbvseg2pet == NULL) exit(1);
  LTAfree(&new2seg);

  return (0);
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
  int c, r, s, f, nthseg, segid;
  double val, v, vhat0, vhat, v2;
  LTA *lta;
  MATRIX *nhits;
  MRI *yseg;      // source volume trilin resampled to seg space (used with RBV)
  MRI *yhat0seg;  // unsmoothed yhat created in seg space (used with RBV)
  MRI *yhatseg;   // smoothed yhat in seg space (used with RBV)

  Timer mytimer;

  printf("   Synthesizing unsmoothed input in seg space... ");
  fflush(stdout);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  yhat0seg = GTMsegSynth(gtm, -1, NULL);
  if (yhat0seg == NULL) {
    printf("ERROR: GTMrbv0() could not synthesize yhat0seg\n");
    return (1);
  }
  printf("  t = %4.2f min\n", mytimer.minutes());
  fflush(stdout);

  printf("   Smoothing synthesized in seg space... ");
  fflush(stdout);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  yhatseg = MRIgaussianSmoothNI(yhat0seg, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  if (yhatseg == NULL) {
    printf("ERROR: GTMrbv0() could not smooth yhatseg\n");
    return (1);
  }
  printf("  t = %4.2f min\n", mytimer.minutes());
  fflush(stdout);

  printf("   Sampling input to seg space with trilin... ");
  fflush(stdout);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  yseg =
      MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth, MRI_FLOAT, gtm->yvol->nframes);
  if (yseg == NULL) {
    printf("ERROR: GTMrbv0() could not alloc yseg\n");
    return (1);
  }
  MRIcopyHeader(gtm->anatseg, yseg);
  MRIcopyPulseParameters(gtm->yvol, yseg);
  printf("  t = %4.2f min\n", mytimer.minutes());
  fflush(stdout);

  lta = LTAcopy(gtm->rbvseg2pet, NULL);
  LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  MRIvol2Vol(gtm->yvol, yseg, (lta->xforms[0].m_L), SAMPLE_TRILINEAR, 0.0);
  LTAfree(&lta);

  printf("   Computing RBV ... ");
  fflush(stdout);
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  gtm->rbv =
      MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth, MRI_FLOAT, gtm->yvol->nframes);
  if (gtm->rbv == NULL) {
    printf("ERROR: GTMrbv0() could not alloc rbv\n");
    return (1);
  }
  MRIcopyHeader(gtm->anatseg, gtm->rbv);
  MRIcopyPulseParameters(gtm->yvol, gtm->rbv);

  if (Gdiag_no > 0) PrintMemUsage(stdout);
  gtm->rbvsegmean = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
  nhits = MatrixAlloc(gtm->beta->rows, 1, MATRIX_REAL);
  for (s = 0; s < gtm->anatseg->depth; s++) {  // crs order not important
    for (c = 0; c < gtm->anatseg->width; c++) {
      for (r = 0; r < gtm->anatseg->height; r++) {
        segid = MRIgetVoxVal(gtm->anatseg, c, r, s, 0);
        if (segid < 0.5) continue;
        for (nthseg = 0; nthseg < gtm->nsegs; nthseg++)
          if (gtm->segidlist[nthseg] == segid) break;
        nhits->rptr[nthseg + 1][1]++;
        for (f = 0; f < gtm->yvol->nframes; f++) {
          v = MRIgetVoxVal(yseg, c, r, s, f);
          vhat0 = MRIgetVoxVal(yhat0seg, c, r, s, f);
          vhat = MRIgetVoxVal(yhatseg, c, r, s, f);
          val = v * vhat0 / (vhat + FLT_EPSILON);  // RBV equation
          MRIsetVoxVal(gtm->rbv, c, r, s, f, val);
          v2 = MRIgetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f);  // track seg means for QA
          MRIsetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f, v2 + val);
        }
      }
    }
  }
  if (Gdiag_no > 0) PrintMemUsage(stdout);
  MRIfree(&yseg);
  MRIfree(&yhat0seg);
  MRIfree(&yhatseg);
  printf("  t = %4.2f min\n", mytimer.minutes());
  fflush(stdout);

  // track seg means for QA
  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    for (f = 0; f < gtm->yvol->nframes; f++) {
      val = MRIgetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f) / nhits->rptr[nthseg + 1][1];
      MRIsetVoxVal(gtm->rbvsegmean, nthseg, 0, 0, f, val);
    }
  }
  MatrixFree(&nhits);

  if (gtm->mask_rbv_to_brain) {
    // Reduce RBV to a tight mask around the brain. This can greatly reduce
    // memory requirements. The RAS space is still that of the anatseg
    // (and so also that of the conformed anat) so no new registration is necessary
    printf("   masking RBV to brain\n");
    if (Gdiag_no > 0) PrintMemUsage(stdout);
    int n, nReplace, ReplaceThis[1000], WithThat[1000];
    MRI *segtmp, *rbvtmp;
    MRI_REGION *region;
    nReplace = 0;
    for (n = 0; n < gtm->ctGTMSeg->nentries; n++) {
      if (gtm->ctGTMSeg->entries[n] == NULL) continue;
      if (gtm->ctGTMSeg->entries[n]->TissueType != 5) continue;  // should not hard-code
      ReplaceThis[nReplace] = n;
      WithThat[nReplace] = 0;
      nReplace++;
    }
    printf("  replacing head voxels with 0\n");
    segtmp = MRIreplaceList(gtm->anatseg, ReplaceThis, WithThat, nReplace, NULL, NULL);
    printf("  computing bounding box  ");
    region = REGIONgetBoundingBox(segtmp, 10);
    REGIONprint(stdout, region);
    printf("  extracting bounding box\n");
    rbvtmp = MRIextractRegion(gtm->rbv, NULL, region);
    free(region);
    region = NULL;
    MRIfree(&segtmp);
    MRIfree(&gtm->rbv);
    gtm->rbv = rbvtmp;
  }
  PrintMemUsage(stdout);
  printf("  RBV took %4.2f min\n", mytimer.minutes());

  return (0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMmgxpvc(GTM *gtm)
  \brief GLM-based modified Muller-Gartner PVC.  Residualizes input
  wrt the GLM estimate of the non-target tissue, the divides the
  residual by the fraction of the target in a voxel. 
   Target=1 cortex only
   Target=2 subcorticalgm only
   Target=3 all GM 
   Target=4 left hemi cortex
   Target=5 right hemi cortex
   Target=6 left hemi subcort gm
   Target=7 right hemi subcort gm
   Target=8 non-lateralized subcort gm
   Targets 4-8 require lateralized tissue type ctab, eg, see
    TissueTypeSchemaLat()
 */
MRI *GTMmgxpvc(GTM *gtm, int Target)
{
  int nthseg, segid, r, f, tt;
  MATRIX *betaNotTarg, *yNotTarg, *ydiff;
  double sum;
  MRI *mgx=NULL;
  COLOR_TABLE_ENTRY *cte;
  //COLOR_TABLE *ttctab = gtm->ctGTMSeg->ctabTissueType;

  // Set beta values to 0 if they are not in the target tissue type(s)
  betaNotTarg = MatrixAlloc(gtm->beta->rows, gtm->beta->cols, MATRIX_REAL);
  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    segid = gtm->segidlist[nthseg];
    tt = gtm->ctGTMSeg->entries[segid]->TissueType;
    cte = gtm->ctGTMSeg->ctabTissueType->entries[tt];
    if(Target == 1 || Target == 3){
      if(strcmp("cortex",cte->name)==0) continue;
      if(strcmp("cortex-lh",cte->name)==0) continue;
      if(strcmp("cortex-rh",cte->name)==0) continue;
    }
    if(Target == 2 || Target == 3){
      if(strcmp("subcort_gm",cte->name)==0) continue;
      if(strcmp("subcort_gm-lh",cte->name)==0) continue;
      if(strcmp("subcort_gm-rh",cte->name)==0) continue;
      if(strcmp("subcort_gm-mid",cte->name)==0) continue;
    }
    if(Target == 4 && strcmp("cortex-lh",cte->name)==0) continue;
    if(Target == 5 && strcmp("cortex-rh",cte->name)==0) continue;
    if(Target == 6 && strcmp("subcort_gm-lh",cte->name)==0) continue;
    if(Target == 7 && strcmp("subcort_gm-rh",cte->name)==0) continue;
    if(Target == 8 && strcmp("subcort_gm-mid",cte->name)==0) continue;
    // otherwise
    for (f = 0; f < gtm->nframes; f++) betaNotTarg->rptr[nthseg + 1][f + 1] = gtm->beta->rptr[nthseg + 1][f + 1];
  }

  // Compute the estimate of the image without the target
  yNotTarg = MatrixMultiplyD(gtm->X, betaNotTarg, NULL);
  // Subtract to resdiualize the PET wrt the non-target tissue
  ydiff = MatrixSubtract(gtm->y, yNotTarg, NULL);

  // Scale by the fraction of target tissue type in voxel
  for (r = 0; r < gtm->X->rows; r++) {
    sum = 0;
    for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      tt = gtm->ctGTMSeg->entries[segid]->TissueType;
      cte = gtm->ctGTMSeg->ctabTissueType->entries[tt];
      if(Target == 1){ // asking for cortex
	if(strcmp("cortex",cte->name)!=0 &&
	   strcmp("cortex-lh",cte->name)!=0 &&
	   strcmp("cortex-rh",cte->name)!=0) continue; // but this is not cortex
      }
      if(Target == 2){ // asking for subcort
	if(strcmp("subcort_gm",cte->name)!=0 && 
	   strcmp("subcort_gm-lh",cte->name)!=0 &&
	   strcmp("subcort_gm-rh",cte->name)!=0) continue; // but this is not subcort
      }
      if(Target == 3){ // asking for any GM
	if(strcmp("cortex",cte->name)!=0 &&
	   strcmp("cortex-lh",cte->name)!=0 &&
	   strcmp("cortex-rh",cte->name)!=0 &&
	   strcmp("subcort_gm",cte->name)!=0 &&
	   strcmp("subcort_gm-lh",cte->name)!=0 &&
	   strcmp("subcort_gm-rh",cte->name)!=0 &&
	   strcmp("subcort_gm-mid",cte->name)!=0) continue; // but this is not GM
      }
      if(Target == 4 && strcmp("cortex-lh",cte->name)!=0) continue;
      if(Target == 5 && strcmp("cortex-rh",cte->name)!=0) continue;
      if(Target == 6 && strcmp("subcort_gm-lh",cte->name)!=0) continue;
      if(Target == 7 && strcmp("subcort_gm-rh",cte->name)!=0) continue;
      if(Target == 8 && strcmp("subcort_gm-mid",cte->name)!=0) continue;

      // otherwise
      sum += gtm->X->rptr[r+1][nthseg+1];
    }
    if (sum < gtm->mgx_gmthresh)
      for (f = 0; f < gtm->nframes; f++) ydiff->rptr[r + 1][f + 1] = 0;
    else
      for (f = 0; f < gtm->nframes; f++) ydiff->rptr[r + 1][f + 1] /= sum;
  }

  mgx = GTMmat2vol(gtm, ydiff, NULL);


  MatrixFree(&betaNotTarg);
  MatrixFree(&yNotTarg);
  MatrixFree(&ydiff);

  return(mgx);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMmgpvc(GTM *gtm)
  \brief Performs Muller-Gartner PVC. Hardcodes tissue type IDs to
  be frame0=cortex, frame1=subcortexgm, frame2=WM.
 */
int GTMmgpvc(GTM *gtm)
{
  int c, r, s, f;
  double vgmpsf, vwmpsf, vwmtac, vtac, vmgtac;
  // MRI *ctxpvf, *subctxpvf, *gmpvf;
  MRI *wmpvf, *wmpvfpsf;

  if (gtm->mg) MRIfree(&gtm->mg);
  gtm->mg = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT, gtm->yvol->nframes);
  if (gtm->mg == NULL) return (1);
  MRIcopyHeader(gtm->yvol, gtm->mg);

  // This is now done outside
  // Compute gray matter PVF with smoothing
  // ctxpvf    = fMRIframe(gtm->ttpvf,0,NULL); // cortex PVF
  // subctxpvf = fMRIframe(gtm->ttpvf,1,NULL); // subcortex GM PVF
  // gmpvf = MRIadd(ctxpvf,subctxpvf,NULL); // All GM PVF
  // Smooth GM PVF by PSF, need to include MB
  // gtm->gmpvfpsf = MRIgaussianSmoothNI(gmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  // Need to add MB here

  // WM PVF
  wmpvf = fMRIframe(gtm->ttpvf, 2, NULL);
  // Smooth WM PVF by PSF
  wmpvfpsf = MRIgaussianSmoothNI(wmpvf, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  if (gtm->UseMBrad) {
    MB2D *mb;
    MRI *mritmp;
    mb = MB2Dcopy(gtm->mbrad, 0, NULL);
    mb->cR = 0;
    mb->rR = 0;
    mritmp = MRImotionBlur2D(wmpvfpsf, mb, NULL);
    MRIfree(&wmpvfpsf);
    wmpvfpsf = mritmp;
    MB2Dfree(&mb);
  }
  if (gtm->UseMBtan) {
    MB2D *mb;
    MRI *mritmp;
    mb = MB2Dcopy(gtm->mbtan, 0, NULL);
    mb->cR = 0;
    mb->rR = 0;
    mritmp = MRImotionBlur2D(wmpvfpsf, mb, NULL);
    MRIfree(&wmpvfpsf);
    wmpvfpsf = mritmp;
    MB2Dfree(&mb);
  }

  // Finally, do the actual MG correction
  for (c = 0; c < gtm->yvol->width; c++) {  // crs order not important
    for (r = 0; r < gtm->yvol->height; r++) {
      for (s = 0; s < gtm->yvol->depth; s++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        vgmpsf = MRIgetVoxVal(gtm->gmpvfpsf, c, r, s, 0);
        if (vgmpsf < gtm->mg_gmthresh) continue;
        vwmpsf = MRIgetVoxVal(wmpvfpsf, c, r, s, 0);
        for (f = 0; f < gtm->yvol->nframes; f++) {
          vwmtac = gtm->mg_reftac->rptr[f + 1][1];
          vtac = MRIgetVoxVal(gtm->yvol, c, r, s, f);
          vmgtac = (vtac - vwmpsf * vwmtac) / vgmpsf;
          MRIsetVoxVal(gtm->mg, c, r, s, f, vmgtac);
        }
      }
    }
  }
  // MRIfree(&ctxpvf);  MRIfree(&subctxpvf);  MRIfree(&gmpvf);
  MRIfree(&wmpvf);
  MRIfree(&wmpvfpsf);
  return (0);
}

// Compute the MG reference TAC
int GTMmgRefTAC(GTM *gtm)
{
  int f, nhits, n, found, nthseg, segid;
  double sum;

  gtm->mg_reftac = MatrixAlloc(gtm->yvol->nframes, 1, MATRIX_REAL);
  for (f = 0; f < gtm->yvol->nframes; f++) {
    nhits = 0;
    sum = 0;
    for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      found = 0;
      for (n = 0; n < gtm->n_mg_refids; n++) {
        if (segid == gtm->mg_refids[n]) {
          found = 1;
          nhits++;
          break;
        }
      }
      if (!found) continue;
      sum += gtm->beta->rptr[nthseg + 1][f + 1];
      printf("   n=%d, nthseg=%d %g\n", n, nthseg, gtm->beta->rptr[nthseg + 1][f + 1]);
    }
    gtm->mg_reftac->rptr[f + 1][1] = sum / nhits;
    printf("   wm tac %2d %2d %g\n", f, nhits, gtm->mg_reftac->rptr[f + 1][1]);
  }
  return (0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMmeltzerpvc(GTM *gtm)
  \brief Performs Meltzer PVC. Hardcodes tissue type IDs to
  be 0=cortex, 1=subcortexgm, 2=WM.
 */
int GTMmeltzerpvc(GTM *gtm)
{
  int c, r, s, f, k, segid, nhits;
  double vgmwmpsf, v, sum;
  MRI *ctxpvf, *subctxpvf, *wmpvf, *gmwmpvf, *gmwmpvfpsf, *mritmp, *nhitseg;

  if (gtm->meltzer) MRIfree(&gtm->meltzer);
  gtm->meltzer = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT, gtm->yvol->nframes);
  if (gtm->meltzer == NULL) return (1);
  MRIcopyHeader(gtm->yvol, gtm->meltzer);

  // Compute gray+white matter PVF with smoothing
  ctxpvf = fMRIframe(gtm->ttpvf, 0, NULL);     // cortex PVF
  subctxpvf = fMRIframe(gtm->ttpvf, 1, NULL);  // subcortex GM PVF
  wmpvf = fMRIframe(gtm->ttpvf, 2, NULL);      // WM PVF
  gmwmpvf = MRIadd(ctxpvf, subctxpvf, NULL);   // All GM PVF
  MRIadd(gmwmpvf, wmpvf, gmwmpvf);             // All GM+WM PVF

  // Setting BinThresh to 0 turns off binarization
  if (gtm->MeltzerBinThresh > 0.0) {
    printf("Binarizing melzter mask before smoothing %lf\n", gtm->MeltzerBinThresh);
    mritmp = MRIbinarize(gmwmpvf, NULL, gtm->MeltzerBinThresh, 0, 1);
    MRIfree(&gmwmpvf);
    gmwmpvf = mritmp;
    if (gtm->MeltzerNDil > 0) {
      printf("Dilating melzter mask %d\n", gtm->MeltzerNDil);
      for (k = 0; k < gtm->MeltzerNDil; k++) MRIdilate(gmwmpvf, gmwmpvf);
    }
  }

  // Smooth GMWM PVF by PSF
  gmwmpvfpsf = MRIgaussianSmoothNI(gmwmpvf, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  // Need to add MB here

  // Finally, do the actual Meltzer correction
  gtm->mzseg = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->nframes);
  nhitseg = MRIallocSequence(gtm->nsegs, 1, 1, MRI_INT, gtm->nframes);
  for (c = 0; c < gtm->yvol->width; c++) {  // crs order not important
    for (r = 0; r < gtm->yvol->height; r++) {
      for (s = 0; s < gtm->yvol->depth; s++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        segid = MRIgetVoxVal(gtm->gtmseg, c, r, s, 0);
        for (k = 0; k < gtm->nsegs; k++)
          if (segid == gtm->segidlist[k]) break;
        // check if not in list, can happen if PET space>AnatSpace and 0
        // not defined in anatseg
        if (k == gtm->nsegs) continue;
        nhits = MRIgetVoxVal(nhitseg, k, 0, 0, 0) + 1;
        vgmwmpsf = MRIgetVoxVal(gmwmpvfpsf, c, r, s, 0);
        if (vgmwmpsf < gtm->MeltzerMaskThresh) continue;
        // count as a hit only if seg is inside the div mask
        // should probably only count it if it is inside the bin mask
        // but if bin thresh is low, this probably won't make a difference
        MRIsetVoxVal(nhitseg, k, 0, 0, 0, nhits);
        for (f = 0; f < gtm->yvol->nframes; f++) {
          v = MRIgetVoxVal(gtm->yvol, c, r, s, f);
          MRIsetVoxVal(gtm->meltzer, c, r, s, f, v / vgmwmpsf);
          sum = MRIgetVoxVal(gtm->mzseg, k, 0, 0, f) + v / vgmwmpsf;
          MRIsetVoxVal(gtm->mzseg, k, 0, 0, f, sum);
        }
      }
    }
  }
  for (k = 0; k < gtm->nsegs; k++) {
    nhits = MRIgetVoxVal(nhitseg, k, 0, 0, 0);
    if (nhits == 0) continue;  // prob neither GM or WM
    for (f = 0; f < gtm->yvol->nframes; f++) {
      sum = MRIgetVoxVal(gtm->mzseg, k, 0, 0, f);
      MRIsetVoxVal(gtm->mzseg, k, 0, 0, f, sum / nhits);
    }
  }
  MRIfree(&ctxpvf);
  MRIfree(&subctxpvf);
  MRIfree(&wmpvf);
  MRIfree(&gmwmpvf);
  MRIfree(&gmwmpvfpsf);
  MRIfree(&nhitseg);
  return (0);
}

/*------------------------------------------------------------------*/
/*
  \fn int GTMsynth(GTM *gtm, int NoiseSeed, int nReps)
  \brief Synthesizes the unsmoothed PET image by computing
   ysynth = X0*beta and then re-packing the result into a volume.
   This is then smoothed with GTMsmoothSynth() to give a full
   synthesis of the input (same as X*beta). If NoiseSeed > 0
   then exponentially distributed noise is added. The noisy
   volume is replicated nReps times with different noise in each rep.
 */
int GTMsynth(GTM *gtm, int NoiseSeed, int nReps)
{
  MATRIX *yhat;
  MRI *mritmp;

  if (gtm->ysynth) MRIfree(&gtm->ysynth);
  gtm->ysynth =
      MRIallocSequence(gtm->gtmseg->width, gtm->gtmseg->height, gtm->gtmseg->depth, MRI_FLOAT, gtm->beta->cols);
  if (gtm->yvol) {
    MRIcopyHeader(gtm->yvol, gtm->ysynth);
    MRIcopyPulseParameters(gtm->yvol, gtm->ysynth);
  }
  yhat = MatrixMultiply(gtm->X0, gtm->beta, NULL);
  GTMmat2vol(gtm, yhat, gtm->ysynth);
  MatrixFree(&yhat);

  if (NoiseSeed > 0) {
    printf("GTMsynth(): adding noise, seed=%d, nReps=%d\n", NoiseSeed, nReps);
    mritmp = MRIrandexp(gtm->ysynth, gtm->mask, NoiseSeed, nReps, NULL);
    MRIfree(&gtm->ysynth);
    gtm->ysynth = mritmp;
    gtm->nframes = nReps;
  }

  return (0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMsmoothSynth(GTM *gtm)
  \brief Smooths the synthesized volume (ysynth) to create the ysynthsm vol.
  Should give the same result as X*beta.
 */
int GTMsmoothSynth(GTM *gtm)
{
  if (gtm->ysynth == NULL) GTMsynth(gtm, 0, 0);
  gtm->ysynthsm = MRIgaussianSmoothNI(gtm->ysynth, gtm->cStd, gtm->rStd, gtm->sStd, gtm->ysynthsm);
  if (gtm->UseMBrad) {
    MB2D *mb;
    MRI *mritmp;
    mb = MB2Dcopy(gtm->mbrad, 0, NULL);
    mb->cR = 0;
    mb->rR = 0;
    mritmp = MRImotionBlur2D(gtm->ysynthsm, mb, NULL);
    MRIfree(&gtm->ysynthsm);
    gtm->ysynthsm = mritmp;
    MB2Dfree(&mb);
  }
  if (gtm->UseMBtan) {
    MB2D *mb;
    MRI *mritmp;
    mb = MB2Dcopy(gtm->mbtan, 0, NULL);
    mb->cR = 0;
    mb->rR = 0;
    mritmp = MRImotionBlur2D(gtm->ysynthsm, mb, NULL);
    MRIfree(&gtm->ysynthsm);
    gtm->ysynthsm = mritmp;
    MB2Dfree(&mb);
  }
  return (0);
}

/*------------------------------------------------------------------*/
/*
  \fn int GTMcheckX(MATRIX *X)
  \brief Checks that all colums sum to 1
 */
int GTMcheckX(MATRIX *X)
{
  int r, c, count;
  double sum, d, dmax;

  count = 0;
  dmax = -1;
  for (r = 0; r < X->rows; r++) {
    sum = 0;
    for (c = 0; c < X->cols; c++) {
      sum += X->rptr[r + 1][c + 1];
    }
    d = abs(sum - 1);
    if (d > .00001) count++;
    if (dmax < d) dmax = d;
  }

  printf("GTMcheckX: count=%d, dmax=%g\n", count, dmax);
  return (count);
}
/*------------------------------------------------------------------------------*/
/*
  \fn int GTMbuildX(GTM *gtm)
  \brief Builds the GTM design matrix both with (X) and without (X0) PSF.  If
  gtm->DoVoxFracCor=1 then corrects for volume fraction effect.
*/
int GTMbuildX(GTM *gtm)
{
  int nthseg, err;

  if (gtm->X == NULL || gtm->X->rows != gtm->nmask || gtm->X->cols != gtm->nsegs) {
    // Alloc or realloc X
    if (gtm->X) MatrixFree(&gtm->X);
    gtm->X = MatrixAlloc(gtm->nmask, gtm->nsegs, MATRIX_REAL);
    if (gtm->X == NULL) {
      printf("ERROR: GTMbuildX(): could not alloc X %d %d\n", gtm->nmask, gtm->nsegs);
      return (1);
    }
  }
  if (gtm->X0 == NULL || gtm->X0->rows != gtm->nmask || gtm->X0->cols != gtm->nsegs) {
    if (gtm->X0) MatrixFree(&gtm->X0);
    gtm->X0 = MatrixAlloc(gtm->nmask, gtm->nsegs, MATRIX_REAL);
    if (gtm->X0 == NULL) {
      printf("ERROR: GTMbuildX(): could not alloc X0 %d %d\n", gtm->nmask, gtm->nsegs);
      return (1);
    }
  }
  gtm->dof = gtm->X->rows - gtm->X->cols;

  Timer timer;

  err = 0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) reduction(+ : err)
#endif
  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    ROMP_PFLB_begin
    
    int segid, k, c, r, s;
    MRI *nthsegpvf = NULL, *nthsegpvfbb = NULL, *nthsegpvfbbsm = NULL, *nthsegpvfbbsmmb = NULL;
    MRI_REGION *region;
    MB2D *mb;
    segid = gtm->segidlist[nthseg];
    if (gtm->DoVoxFracCor)
      nthsegpvf = fMRIframe(gtm->segpvf, nthseg, NULL);  // extract PVF for this seg
    else
      nthsegpvf = MRIbinarizeMatch(gtm->gtmseg, &segid, 1, 0, NULL);  // or get binary mask
    // Extract a region for the seg. This speeds up smoothing.
    region = REGIONgetBoundingBox(nthsegpvf, gtm->nPad);  // tight+pad bounding box
    if (region->dx < 0) {
      printf(
          "ERROR: creating region for nthseg=%d, segid=%d, %s\n", nthseg, segid, gtm->ctGTMSeg->entries[segid]->name);
      printf(
          "It may be that there are no voxels for this seg when mapped "
          "into the input space. \nCheck %s/aux/seg.nii.gz and the registration\n",
          gtm->OutDir);
      err++;
      continue;
    }
    nthsegpvfbb = MRIextractRegion(nthsegpvf, NULL, region);  // extract BB
    if (nthsegpvfbb == NULL) {
      printf("ERROR: extracting nthseg=%d, segid=%d, %s\n", nthseg, segid, gtm->ctGTMSeg->entries[segid]->name);
      printf(
          "It may be that there are no voxels for this seg when mapped "
          "into the input space. \nCheck %s/aux/seg.nii.gz and the registration\n",
          gtm->OutDir);
      err++;
      continue;
    }
    nthsegpvfbbsm = MRIgaussianSmoothNI(nthsegpvfbb, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
    if (gtm->UseMBrad) {
      // Order of operations should not matter
      mb = MB2Dcopy(gtm->mbrad, 0, NULL);
      mb->cR = region->x;
      mb->rR = region->y;
      nthsegpvfbbsmmb = MRImotionBlur2D(nthsegpvfbbsm, mb, NULL);
      MRIfree(&nthsegpvfbbsm);
      nthsegpvfbbsm = nthsegpvfbbsmmb;
      MB2Dfree(&mb);
    }
    if (gtm->UseMBtan) {
      // Order of operations should not matter
      mb = MB2Dcopy(gtm->mbtan, 0, NULL);
      mb->cR = region->x;
      mb->rR = region->y;
      nthsegpvfbbsmmb = MRImotionBlur2D(nthsegpvfbbsm, mb, NULL);
      MRIfree(&nthsegpvfbbsm);
      nthsegpvfbbsm = nthsegpvfbbsmmb;
      MB2Dfree(&mb);
    }
    // Fill X, creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way. See GTMvol2mat()
    k = 0;
    for (s = 0; s < gtm->yvol->depth; s++) {
      for (c = 0; c < gtm->yvol->width; c++) {
        for (r = 0; r < gtm->yvol->height; r++) {
          if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
          k++;  // have to incr here in case continue below
          if (c < region->x || c >= region->x + region->dx) continue;
          if (r < region->y || r >= region->y + region->dy) continue;
          if (s < region->z || s >= region->z + region->dz) continue;
          // do not use k+1 here because it has already been incr above
          if (!gtm->Optimizing)
            gtm->X0->rptr[k][nthseg + 1] = MRIgetVoxVal(nthsegpvfbb, c - region->x, r - region->y, s - region->z, 0);

          gtm->X->rptr[k][nthseg + 1] = MRIgetVoxVal(nthsegpvfbbsm, c - region->x, r - region->y, s - region->z, 0);
        }
      }
    }
    MRIfree(&nthsegpvf);
    MRIfree(&nthsegpvfbb);
    MRIfree(&nthsegpvfbbsm);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  if (!gtm->Optimizing) printf(" Build time %6.4f, err = %d\n", timer.seconds(), err);
  fflush(stdout);
  if (err) gtm->X = NULL;

  return (0);
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
  int c, r, s, f, segid, segno, nframes;

  if (frame < 0)
    nframes = gtm->nframes;
  else
    nframes = 1;

  if (synth == NULL) {
    synth = MRIallocSequence(gtm->rbvseg->width, gtm->rbvseg->height, gtm->rbvseg->depth, MRI_FLOAT, nframes);
    if (synth == NULL) {
      printf("ERROR: GTMsegSynth(): could not alloc %d frames\n", nframes);
      return (NULL);
    }
    MRIcopyHeader(gtm->rbvseg, synth);
    MRIcopyPulseParameters(gtm->yvol, synth);
  }

  for (c = 0; c < gtm->rbvseg->width; c++) {  // crs order does not matter here
    for (r = 0; r < gtm->rbvseg->height; r++) {
      for (s = 0; s < gtm->rbvseg->depth; s++) {
        segid = MRIgetVoxVal(gtm->rbvseg, c, r, s, 0);
        if (segid == 0) continue;
        for (segno = 0; segno < gtm->nsegs; segno++)
          if (segid == gtm->segidlist[segno]) break;
        if (segno == gtm->nsegs) {
          printf("ERROR: GTMsegSynth(): could not find a match for segid=%d\n", segid);
          for (segno = 0; segno < gtm->nsegs; segno++) printf("%3d %5d\n", segno, gtm->segidlist[segno]);
          return (NULL);
        }
        if (frame < 0) {
          for (f = 0; f < gtm->beta->cols; f++) MRIsetVoxVal(synth, c, r, s, f, gtm->beta->rptr[segno + 1][f + 1]);
        }
        else
          MRIsetVoxVal(synth, c, r, s, 0, gtm->beta->rptr[segno + 1][frame + 1]);
      }
    }
  }

  return (synth);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
  \brief Prints the Muller-Gartner WM reference TAC to the given file pointer
*/
int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
{
  int f;
  for (f = 0; f < gtm->yvol->nframes; f++) fprintf(fp, "%3d %10.5f\n", f, gtm->mg_reftac->rptr[f + 1][1]);

  return (0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMwriteMGRefTAC(GTM *gtm, char *filename)
  \brief Writes the Muller-Gartner WM reference TAC to the given file
*/
int GTMwriteMGRefTAC(GTM *gtm, char *filename)
{
  FILE *fp;
  fp = fopen(filename, "w");
  GTMprintMGRefTAC(gtm, fp);
  fclose(fp);
  return (0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrescale(GTM *gtm)
  \brief Computes global rescaling factor and applies it to yvol, beta, and y.
  The factor = 100/mean(beta_i) where i is the list of scale seg IDs (scale_ref_ids)
*/
int GTMrescale(GTM *gtm)
{
  int f, n, nthseg, segid, found, nhits;
  double sum;

  // global rescaling
  nhits = 0;
  sum = 0;
  for (f = 0; f < gtm->yvol->nframes; f++) {
    for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      found = 0;
      for (n = 0; n < gtm->n_scale_refids; n++) {
        if (segid == gtm->scale_refids[n]) {
          found = 1;
          nhits++;
          break;
        }
      }
      if (!found) continue;
      sum += gtm->beta->rptr[nthseg + 1][f + 1];
    }
  }
  gtm->scale = gtm->scale_refval / (sum / nhits);
  printf("gtm multiplicative scale %g\n", gtm->scale);

  MRImultiplyConst(gtm->yvol, gtm->scale, gtm->yvol);
  MatrixScalarMul(gtm->beta, gtm->scale, gtm->beta);
  MatrixScalarMul(gtm->y, gtm->scale, gtm->y);

  return (0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMsteadyState(GTM *gtm)
  \brief out = (in - ref)*dcf/(scale*bpc)
*/
int GTMsteadyState(GTM *gtm)
{
  int f, n, c, r, s;
  double v, y, ref;

  if (gtm->DoSteadyState == 0) return (1);

  v = gtm->ss_dcf / (gtm->ss_scale * gtm->ss_bpc);
  printf("SteadyState: v = %g\n", v);

  for (f = 0; f < gtm->nframes; f++) {
    ref = gtm->km_reftac->rptr[f + 1][1];
    for (n = 0; n < gtm->beta->rows; n++) gtm->beta->rptr[n + 1][f + 1] = v * (gtm->beta->rptr[n + 1][f + 1] - ref);
    for (n = 0; n < gtm->y->rows; n++) gtm->y->rptr[n + 1][f + 1] = v * (gtm->y->rptr[n + 1][f + 1] - ref);
    for (c = 0; c < gtm->yvol->width; c++) {  // crs order not important
      for (r = 0; r < gtm->yvol->height; r++) {
        for (s = 0; s < gtm->yvol->depth; s++) {
          if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
          y = MRIgetVoxVal(gtm->yvol, c, r, s, f);
          MRIsetVoxVal(gtm->yvol, c, r, s, f, v * (y - ref));
        }
      }
    }
  }

  return (0);
}

/*-------------------------------------------------------------------------------*/
/*
  \fn int GTMttest(GTM *gtm)
  \brief Computes a t-test for each contrast. This includes computing gamma,
  gammavar, t, and p.
*/
int GTMttest(GTM *gtm)
{
  MATRIX *Ct, *CiXtX, *CiXtXCt;
  int n, f, nframes;
  GTMCON *gtmc;
  double F;

  nframes = gtm->beta->cols;

  for (n = 0; n < gtm->nContrasts; n++) {
    gtmc = gtm->contrasts[n];
    gtmc->gamma = MatrixMultiply(gtmc->C, gtm->beta, NULL);
    gtmc->gammavar = MatrixAlloc(1, nframes, MATRIX_REAL);
    gtmc->t = MatrixAlloc(1, nframes, MATRIX_REAL);
    gtmc->p = MatrixAlloc(1, nframes, MATRIX_REAL);
    Ct = MatrixTranspose(gtmc->C, NULL);
    CiXtX = MatrixMultiply(gtmc->C, gtm->iXtX, NULL);
    CiXtXCt = MatrixMultiply(CiXtX, Ct, NULL);
    for (f = 0; f < nframes; f++) {
      gtmc->gammavar->rptr[1][f + 1] = CiXtXCt->rptr[1][1] * gtm->rvar->rptr[1][f + 1];
      gtmc->t->rptr[1][f + 1] = gtmc->gamma->rptr[1][f + 1] / sqrt(gtmc->gammavar->rptr[1][f + 1]);
      F = gtmc->t->rptr[1][f + 1] * gtmc->t->rptr[1][f + 1];
      gtmc->p->rptr[1][f + 1] = sc_cdf_fdist_Q(F, gtmc->C->rows, gtm->dof);
    }
    MatrixFree(&Ct);
    MatrixFree(&CiXtX);
    MatrixFree(&CiXtXCt);
  }
  return (0);
}
/*------------------------------------------------------------------------------*/
/*
  \fn int GTMwriteContrasts(GTM *GTM)
  \brief Creates an ASCII file in the output folder for each contrast with gamma,
  gammavar, t, and p
 */
int GTMwriteContrasts(GTM *gtm)
{
  int n, nframes, f;
  GTMCON *gtmc;
  char tmpstr[5000];
  FILE *fp;

  nframes = gtm->beta->cols;

  for (n = 0; n < gtm->nContrasts; n++) {
    gtmc = gtm->contrasts[n];
    sprintf(tmpstr, "%s/%s.mat", gtm->OutDir, gtmc->name);
    MatlabWrite(gtmc->C, tmpstr, "C");
    sprintf(tmpstr, "%s/%s.dat", gtm->OutDir, gtmc->name);
    fp = fopen(tmpstr, "w");
    for (f = 0; f < nframes; f++) {
      fprintf(fp,
              "%2d %20.15f %20.15f %20.15f %20.15f\n",
              f,
              gtmc->gamma->rptr[1][1],
              gtmc->gammavar->rptr[1][1],
              gtmc->t->rptr[1][1],
              gtmc->p->rptr[1][1]);
    }
    fclose(fp);
  }
  return (0);
}

/*-------------------------------------------------------------------------*/
/*
  \fn int GTMcheckRefIds(GTM *gtm)
  \brief Checks the segmentation IDs used for rescaling, MG, KM Ref, and KM HB
  to make sure that they are in the segmentation.
 */
int GTMcheckRefIds(GTM *gtm)
{
  int n, m, ok;

  if (gtm->rescale) {
    for (n = 0; n < gtm->n_scale_refids; n++) {
      ok = 0;
      for (m = 0; m < gtm->nsegs; m++) {
        if (gtm->segidlist[m] == gtm->scale_refids[n]) {
          ok = 1;
          break;
        }
      }
      if (!ok) {
        printf("ERROR: scale reference id %d cannot be found in seg id list\n", gtm->scale_refids[n]);
        fprintf(gtm->logfp, "ERROR: scale reference id %d cannot be found in seg id list\n", gtm->scale_refids[n]);
        return (1);
      }
    }
  }

  if (gtm->DoMGPVC) {
    for (n = 0; n < gtm->n_mg_refids; n++) {
      ok = 0;
      for (m = 0; m < gtm->nsegs; m++) {
        if (gtm->segidlist[m] == gtm->mg_refids[n]) {
          ok = 1;
          break;
        }
      }
      if (!ok) {
        printf("ERROR: MG reference id %d cannot be found in seg id list\n", gtm->mg_refids[n]);
        fprintf(gtm->logfp, "ERROR: scale reference id %d cannot be found in seg id list\n", gtm->mg_refids[n]);
        return (1);
      }
    }
  }

  if (gtm->DoKMRef) {
    for (n = 0; n < gtm->n_km_refids; n++) {
      ok = 0;
      for (m = 0; m < gtm->nsegs; m++) {
        if (gtm->segidlist[m] == gtm->km_refids[n]) {
          ok = 1;
          break;
        }
      }
      if (!ok) {
        printf("ERROR: KM reference id %d cannot be found in seg id list\n", gtm->km_refids[n]);
        fprintf(gtm->logfp, "ERROR: scale reference id %d cannot be found in seg id list\n", gtm->km_refids[n]);
        return (1);
      }
    }
  }

  if (gtm->DoKMHB) {
    for (n = 0; n < gtm->n_km_hbids; n++) {
      ok = 0;
      for (m = 0; m < gtm->nsegs; m++) {
        if (gtm->segidlist[m] == gtm->km_hbids[n]) {
          ok = 1;
          break;
        }
      }
      if (!ok) {
        printf("ERROR: KM high binding id %d cannot be found in seg id list\n", gtm->km_hbids[n]);
        fprintf(gtm->logfp, "ERROR: scale high binding id %d cannot be found in seg id list\n", gtm->km_hbids[n]);
        return (1);
      }
    }
  }

  return (0);
}

/*
  \fn int GTMprintRefIds(GTM *gtm, FILE *fp)
  \brief Prints the segmentation IDs used for rescaling, MG, KM Ref, and KM HB
  to the given FILE pointer.
 */
int GTMprintRefIds(GTM *gtm, FILE *fp)
{
  int n;

  if (gtm->rescale) {
    fprintf(fp, "Segmentations used for rescaling\n");
    for (n = 0; n < gtm->n_scale_refids; n++)
      fprintf(fp, "%4d %s\n", gtm->scale_refids[n], gtm->ctGTMSeg->entries[gtm->scale_refids[n]]->name);
  }

  if (gtm->DoMGPVC) {
    fprintf(fp, "Segmentations used for MG PVC WM reference\n");
    for (n = 0; n < gtm->n_mg_refids; n++)
      fprintf(fp, "%4d %s\n", gtm->mg_refids[n], gtm->ctGTMSeg->entries[gtm->mg_refids[n]]->name);
  }

  if (gtm->DoKMRef) {
    fprintf(fp, "Segmentations used for KM reference TAC\n");
    for (n = 0; n < gtm->n_km_refids; n++)
      fprintf(fp, "%4d %s\n", gtm->km_refids[n], gtm->ctGTMSeg->entries[gtm->km_refids[n]]->name);
  }

  if (gtm->DoKMHB) {
    fprintf(fp, "Segmentations used for KM High Binding TAC\n");
    for (n = 0; n < gtm->n_km_hbids; n++)
      fprintf(fp, "%4d %s\n", gtm->km_hbids[n], gtm->ctGTMSeg->entries[gtm->km_hbids[n]]->name);
  }
  fflush(fp);
  return (0);
}
/*
  \fn int GTMrefTAC(GTM *gtm)
  \brief Computes KM reference and hibinding TACs. Also writes them
  to the output folder. The HB TAC is written as both an ascii file
  and a nii.gz (the later needed for KM analysis with mri_glmfit)
 */
int GTMrefTAC(GTM *gtm)
{
  int f, n, nthseg, segid;
  double sum;
  char tmpstr[5000];
  FILE *fp;

  if (gtm->DoKMRef) {
    sprintf(tmpstr, "%s/km.ref.tac.dat", gtm->OutDir);
    fp = fopen(tmpstr, "w");
    if (fp == NULL) return (1);
    gtm->km_reftac = MatrixAlloc(gtm->yvol->nframes, 1, MATRIX_REAL);
    for (f = 0; f < gtm->yvol->nframes; f++) {
      sum = 0;
      for (n = 0; n < gtm->n_km_refids; n++) {
        for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
          segid = gtm->segidlist[nthseg];
          if (segid == gtm->km_refids[n]) break;
        }
        sum += gtm->beta->rptr[nthseg + 1][f + 1];
      }
      gtm->km_reftac->rptr[f + 1][1] = sum / gtm->n_km_refids;
      fprintf(fp, "%20.15lf\n", sum / gtm->n_km_refids);
    }
    fclose(fp);
  }

  if (gtm->DoKMHB) {
    sprintf(tmpstr, "%s/km.hb.tac.dat", gtm->OutDir);
    fp = fopen(tmpstr, "w");
    if (fp == NULL) return (1);
    gtm->km_hbtac = MatrixAlloc(gtm->yvol->nframes, 1, MATRIX_REAL);
    for (f = 0; f < gtm->yvol->nframes; f++) {
      sum = 0;
      for (n = 0; n < gtm->n_km_hbids; n++) {
        for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
          segid = gtm->segidlist[nthseg];
          if (segid == gtm->km_hbids[n]) break;
        }
        sum += gtm->beta->rptr[nthseg + 1][f + 1];
      }
      gtm->km_hbtac->rptr[f + 1][1] = sum / gtm->n_km_hbids;
      fprintf(fp, "%20.15lf\n", gtm->km_hbtac->rptr[f + 1][1]);
    }
    fclose(fp);

    MRI *mritmp = MRIallocSequence(1, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
    for (f = 0; f < gtm->yvol->nframes; f++) MRIsetVoxVal(mritmp, 0, 0, 0, f, gtm->km_hbtac->rptr[f + 1][1]);
    sprintf(tmpstr, "%s/km.hb.tac.nii.gz", gtm->OutDir);
    MRIwrite(mritmp, tmpstr);
    MRIfree(&mritmp);
  }

  return (0);
}

/*
  \fn int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat)
  \brief Prints the replacement list to the FILE pointer.
 */
int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat)
{
  int n;
  for (n = 0; n < nReplace; n++) fprintf(fp, "%5d %5d\n", ReplaceThis[n], WithThat[n]);
  return (0);
}

/*
  \fn int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat)
  \brief Checks replacement list to make sure that no item in ReplaceThis list appears in
  the WithThat list.
 */
int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat)
{
  int n, m;
  for (n = 0; n < nReplace; n++) {
    for (m = 0; m < nReplace; m++) {
      if (ReplaceThis[n] == WithThat[m]) {
        printf("ERROR: item %d appears as both source and target seg id in replacement list\n", ReplaceThis[n]);
        return (1);
      }
    }
  }
  return (0);
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
  int nlist, r, nth;
  char tmpstr[1001], *s;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    int e = errno;
    printf("ERROR: could not open %s %d\n", fname, e);
    printf("%s\n", strerror(e));
    return (1);
  }
  nlist = *nReplace;

  nth = -1;
  while (1) {
    nth++;
    s = fgets(tmpstr, 1000, fp);
    if (s == NULL) break;
    if (tmpstr[0] == '#') continue;
    r = CountItemsInString(tmpstr);
    if (r != 2) {
      printf("ERROR: line %d in %s has the wrong format (has %d elements, expecting 2)\n", nth, fname, r);
      return (1);
    }
    sscanf(tmpstr, "%d %d", &ReplaceThis[nlist], &WithThat[nlist]);
    nlist++;
  }
  if (nlist == 0) {
    printf("ERROR: could not find any replacement items in %s\n", fname);
    return (1);
  }
  printf("Read in %d replacement items from %s\n", nlist, fname);
  *nReplace += nlist;
  return (0);
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
  LTA *lta;
  double std;

  gtm->mask = MRIalloc(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT);
  MRIcopyHeader(gtm->yvol, gtm->mask);
  if (LTAmriIsSource(gtm->seg2pet, gtm->yvol))
    lta = LTAcopy(gtm->seg2pet, NULL);
  else
    lta = LTAinvert(gtm->seg2pet, NULL);
  LTAchangeType(lta, LINEAR_VOX_TO_VOX);

  // Sample the seg into the pet space
  MRIvol2Vol(gtm->anatseg, gtm->mask, lta->xforms[0].m_L, SAMPLE_NEAREST, 0);
  LTAfree(&lta);
  // Threshold it at 0.5 to give the mask
  MRIbinarize(gtm->mask, gtm->mask, 0.5, 0, 1);
  // Smooth binary mask - MB not super important here
  std = gtm->automask_fwhm / sqrt(log(256.0));  // convert fwhm to std
  MRIgaussianSmoothNI(gtm->mask, std, std, std, gtm->mask);
  // Binarize again to get final mask
  MRIbinarize(gtm->mask, gtm->mask, gtm->automask_thresh, 0, 1);

  return (0);
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
  int f, s, c, r, n, nhits, nhitsb, segid, tt;
  double sum, sumb;
  FILE *fp;
  char tmpstr[2000];

  if (gtm->rvargm == NULL) gtm->rvargm = MatrixAlloc(1, gtm->res->cols, MATRIX_REAL);
  if (gtm->rvarbrain == NULL) gtm->rvarbrain = MatrixAlloc(1, gtm->res->cols, MATRIX_REAL);
  ct = gtm->ctGTMSeg;

  for (f = 0; f < gtm->res->cols; f++) {
    sum = 0;
    sumb = 0;
    n = -1;
    nhits = 0;
    nhitsb = 0;
    // slice, col, row order is important here as is skipping the mask
    for (s = 0; s < gtm->yvol->depth; s++) {
      for (c = 0; c < gtm->yvol->width; c++) {
        for (r = 0; r < gtm->yvol->height; r++) {
          if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
          n++;
          segid = MRIgetVoxVal(gtm->gtmseg, c, r, s, 0);
          // check if not in list, can happen if PET space>AnatSpace and 0
          // not defined in anatseg
          if (ct->entries[segid] == NULL) {
            printf("WARNING: GTMrvarGM(): segid=%d (%d,%d,%d) does not exist in ctab\n", segid, c, r, s);
            // CTABwriteFileASCIItt(gtm->ctGTMSeg,"GTMrvarGM.ctab");
            // MRIwrite(gtm->gtmseg,"GTMrvarGM.gtmseg.mgz");
            // exit(1);
            continue;
          }
          tt = ct->entries[segid]->TissueType;
          if (tt == 1 || tt == 2) {
            // GM (hardcoded)
            sum += ((double)gtm->res->rptr[n + 1][f + 1] * gtm->res->rptr[n + 1][f + 1]);
            nhits++;
          }
          if (tt == 1 || tt == 2 || tt == 3) {
            // GM  or WM (hardcoded)
            sumb += ((double)gtm->res->rptr[n + 1][f + 1] * gtm->res->rptr[n + 1][f + 1]);
            nhitsb++;
          }
        }
      }
    }
    gtm->rvargm->rptr[1][f + 1] = sum / nhits;
    gtm->rvarbrain->rptr[1][f + 1] = sumb / nhitsb;
    if (f == 0 && !gtm->Optimizing)
      printf("rvargm %2d %6.4f %6.4f\n", f, gtm->rvargm->rptr[1][f + 1], gtm->rvarbrain->rptr[1][f + 1]);
  }

  if (!gtm->Optimizing) {
    sprintf(tmpstr, "%s/rvar.gm.dat", gtm->AuxDir);
    fp = fopen(tmpstr, "w");
    for (f = 0; f < gtm->res->cols; f++) fprintf(fp, "%30.20f\n", gtm->rvargm->rptr[1][f + 1]);
    fclose(fp);
    sprintf(tmpstr, "%s/rvar.brain.dat", gtm->AuxDir);
    fp = fopen(tmpstr, "w");
    for (f = 0; f < gtm->res->cols; f++) fprintf(fp, "%30.20f\n", gtm->rvarbrain->rptr[1][f + 1]);
    fclose(fp);
  }

  if (gtm->rescale) {
    sprintf(tmpstr, "%s/rvar.gm.unscaled.dat", gtm->AuxDir);
    fp = fopen(tmpstr, "w");
    for (f = 0; f < gtm->res->cols; f++)
      fprintf(fp, "%30.20f\n", gtm->rvargm->rptr[1][f + 1] / (gtm->scale * gtm->scale));
    fclose(fp);
    sprintf(tmpstr, "%s/rvar.brain.unscaled.dat", gtm->AuxDir);
    fp = fopen(tmpstr, "w");
    for (f = 0; f < gtm->res->cols; f++)
      fprintf(fp, "%30.20f\n", gtm->rvarbrain->rptr[1][f + 1] / (gtm->scale * gtm->scale));
    fclose(fp);
  }

  return (0);
}

/*
  \fn MRI **GTMlocal(MRI *src, MRI *pvf, MRI *mask, int nrad, MRI **pvc)
  \brief Performs "localized" GTM correction by estimating the contribution
  from each tissue type based on solving a GLM of neighbors around each
  target voxel. This is a map-based, instead of ROI-based, PVC.
  \param src is the source data
  \param pvf is partial volume fraction smoothed by the PSF, each
  tissue type in its own frame
  \param nrad defines neighborhood as (2*nrad+1)^3 voxels
*/
MRI **GTMlocal(GTM *gtm, MRI **pvc)
{
  int const nTT = gtm->ttpvf->nframes;

  MRI * const pvfpsf = MRIgaussianSmoothNI(gtm->ttpvf, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  // Need to add MB here

  if (pvc == NULL) {
    pvc = (MRI **)calloc(sizeof(MRI *), nTT);
    int tt;
    for (tt = 0; tt < nTT; tt++) {
      pvc[tt] = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT, gtm->yvol->nframes);
      if (pvc[tt] == NULL) {
        printf("ERROR: GTMlocal(): could not alloc %d\n", tt);
        return (NULL);
      }
      MRIcopyHeader(gtm->yvol, pvc[tt]);
      MRIcopyPulseParameters(gtm->yvol, pvc[tt]);
    }
  }
  if (MRIdimMismatch(gtm->yvol, pvc[0], 1)) {
    printf("ERROR: GTMlocal(): gtm->yvol-pvc dim mismatch\n");
    return (NULL);
  }

  if (gtm->lgtm->res == NULL) {
    gtm->lgtm->res =
        MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT, gtm->yvol->nframes);
    MRIcopyHeader(gtm->yvol, gtm->lgtm->res);
    MRIcopyPulseParameters(gtm->yvol, gtm->lgtm->res);
  }
  if (MRIdimMismatch(gtm->yvol, gtm->lgtm->res, 1)) {
    printf("ERROR: GTMlocal(): gtm->yvol-res dim mismatch\n");
    return (NULL);
  }

  if (gtm->lgtm->rvar == NULL) {
    gtm->lgtm->rvar =
        MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth, MRI_FLOAT, gtm->yvol->nframes);
    MRIcopyHeader(gtm->yvol, gtm->lgtm->rvar);
    MRIcopyPulseParameters(gtm->yvol, gtm->lgtm->rvar);
  }
  if (MRIdimMismatch(gtm->yvol, gtm->lgtm->rvar, 0)) {
    printf("ERROR: GTMlocal(): gtm->yvol-rvar dim mismatch\n");
    return (NULL);
  }

  int const nvmax = (2 * gtm->lgtm->nrad + 1) * (2 * gtm->lgtm->nrad + 1) * (2 * gtm->lgtm->nrad + 1);
  printf("GTMlocal(): nrad = %d, nvmax = %d, nTT=%d, Xthresh %f\n", gtm->lgtm->nrad, nvmax, nTT, gtm->lgtm->Xthresh);
  Timer timer;

  int c;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  printf("     nthreads = %d\n", omp_get_max_threads());
#pragma omp parallel for if_ROMP(assume_reproducible)
#endif
  for (c = 0; c < gtm->yvol->width; c++) {
    ROMP_PFLB_begin
    
    MATRIX *X, *y, *beta = NULL, *Xt = NULL, *XtX = NULL, *Xty = NULL, *iXtX = NULL, *Xsum, *ytmp, *Xtmp;
    MATRIX *yhat, *eres;
    int r, s, f, dc, dr, ds, tt, nth, nv, nkeep, indkeep[100], nthbeta;
    double v;
    for (r = 0; r < gtm->yvol->height; r++) {
      for (s = 0; s < gtm->yvol->depth; s++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        y = MatrixAlloc(nvmax, gtm->yvol->nframes, MATRIX_REAL);
        X = MatrixAlloc(nvmax, pvfpsf->nframes, MATRIX_REAL);

        nth = 0;  // row of y or X
        for (f = 0; f < gtm->yvol->nframes; f++) y->rptr[nth + 1][f + 1] = MRIgetVoxVal(gtm->yvol, c, r, s, f);
        for (tt = 0; tt < pvfpsf->nframes; tt++) X->rptr[nth + 1][tt + 1] = MRIgetVoxVal(pvfpsf, c, r, s, tt);
        nth++;
        for (dc = -gtm->lgtm->nrad; dc <= gtm->lgtm->nrad; dc++) {
          if (c + dc < 0 || c + dc >= gtm->yvol->width) continue;
          for (dr = -gtm->lgtm->nrad; dr <= gtm->lgtm->nrad; dr++) {
            if (r + dr < 0 || r + dr >= gtm->yvol->height) continue;
            for (ds = -gtm->lgtm->nrad; ds <= gtm->lgtm->nrad; ds++) {
              if (s + ds < 0 || s + ds >= gtm->yvol->depth) continue;
              if (dc == 0 && dr == 0 && ds == 0) continue;
              for (f = 0; f < gtm->yvol->nframes; f++)
                y->rptr[nth + 1][f + 1] = MRIgetVoxVal(gtm->yvol, c + dc, r + dr, s + ds, f);
              for (tt = 0; tt < pvfpsf->nframes; tt++) {
                v = MRIgetVoxVal(pvfpsf, c + dc, r + dr, s + ds, tt);
                X->rptr[nth + 1][tt + 1] = v;
              }
              nth++;
            }      // ds
          }        // dr
        }          // dc
        nv = nth;  // nvox in y/X. May not be nvmax if near edge

        // Each ttype needs a minmal representation
        // Sum the columns of X
        Xsum = MatrixSum(X, 1, NULL);
        nkeep = 0;
        for (tt = 0; tt < pvfpsf->nframes; tt++) {
          indkeep[tt] = 0;
          if (Xsum->rptr[1][tt + 1] > gtm->lgtm->Xthresh) {
            indkeep[tt] = 1;
            nkeep++;
          }
        }
        if (nkeep == 0) continue;

        if (nkeep != pvfpsf->nframes || nv != nvmax) {
          // Either not all tissue types will be kept or not all nvmax were avail
          ytmp = MatrixAlloc(nv, gtm->yvol->nframes, MATRIX_REAL);
          Xtmp = MatrixAlloc(nv, nkeep, MATRIX_REAL);
          for (nth = 0; nth < nv; nth++) {
            for (f = 0; f < gtm->yvol->nframes; f++) ytmp->rptr[nth + 1][f + 1] = y->rptr[nth + 1][f + 1];
            nkeep = 0;
            for (tt = 0; tt < pvfpsf->nframes; tt++) {
              if (Xsum->rptr[1][tt + 1] > gtm->lgtm->Xthresh) {
                Xtmp->rptr[nth + 1][nkeep + 1] = X->rptr[nth + 1][tt + 1];
                nkeep++;
              }
            }
          }
          MatrixFree(&y);
          MatrixFree(&X);
          y = ytmp;
          X = Xtmp;
        }
        MatrixFree(&Xsum);
        Xt = MatrixTranspose(X, NULL);
        XtX = MatrixMultiplyD(Xt, X, NULL);
        iXtX = MatrixInverse(XtX, NULL);
        if (iXtX == NULL) {
          MatrixFree(&y);
          MatrixFree(&X);
          MatrixFree(&Xt);
          MatrixFree(&XtX);
          continue;
        }

        Xty = MatrixMultiplyD(Xt, y, NULL);
        beta = MatrixMultiplyD(iXtX, Xty, NULL);
        nthbeta = 0;
        for (tt = 0; tt < pvfpsf->nframes; tt++) {
          if (indkeep[tt] == 0) continue;
          for (f = 0; f < gtm->yvol->nframes; f++) MRIsetVoxVal(pvc[tt], c, r, s, f, beta->rptr[nthbeta + 1][f + 1]);
          nthbeta++;
        }

        yhat = MatrixMultiply(X, beta, NULL);
        eres = MatrixSubtract(y, yhat, NULL);
        for (f = 0; f < gtm->yvol->nframes; f++) MRIsetVoxVal(gtm->lgtm->res, c, r, s, f, eres->rptr[1][f + 1]);
        for (f = 0; f < gtm->yvol->nframes; f++) {
          v = 0;
          for (nth = 0; nth < nv; nth++) v += (eres->rptr[nth + 1][f + 1] * eres->rptr[nth + 1][f + 1]);
          v = v / (X->rows - X->cols);
          MRIsetVoxVal(gtm->lgtm->rvar, c, r, s, f, v);
        }

        if (c == 56 && r == 55 && s == 33) {
          printf("c=%d; r=%d; s=%d;\n", c + 1, r + 1, s + 1);
          printf("nv=%d, nkeep=%d, nf=%d\n", nv, nkeep, pvfpsf->nframes);
          // MatrixPrint(stdout,y);
          // MatrixPrint(stdout,X);
          MatrixPrint(stdout, beta);
          fflush(stdout);
          MatlabWrite(y, "y.mat", "y");
          MatlabWrite(X, "X.mat", "X");
          MatlabWrite(beta, "b.mat", "beta");
          // exit(1);
        }
        MatrixFree(&beta);

        MatrixFree(&y);
        MatrixFree(&yhat);
        MatrixFree(&eres);
        MatrixFree(&X);
        MatrixFree(&Xt);
        MatrixFree(&XtX);
        MatrixFree(&Xty);
        MatrixFree(&iXtX);
      }  // s
    }    // r
    ROMP_PFLB_end
  }      // c
  ROMP_PF_end
  
  MRI * pvfpsf_nonconst = pvfpsf;     // don't use pvfpsf after here!
  MRIfree(&pvfpsf_nonconst);
  
  printf("\n");
  // printf("nNull=%d\n",nNull);
  printf("t=%6.4f\n", timer.seconds());
  printf("GTMlocal(); done\n");
  fflush(stdout);

  return (pvc);
}

/*
  \fn MRI int GTMttPercent(GTM *gtm)
  \brief Computes the percent of each tissue type contributing to
   each segmentation.
*/
int GTMttPercent(GTM *gtm)
{
  int nTT, k, s, c, r, segid, nthseg, mthseg, mthsegid, tt;
  double sum;

  nTT = gtm->ttpvf->nframes;
  if (gtm->ttpct != NULL) MatrixFree(&gtm->ttpct);
  gtm->ttpct = MatrixAlloc(gtm->nsegs, nTT, MATRIX_REAL);

  // Must be done in same order as GTMbuildX()
  k = 0;
  for (s = 0; s < gtm->yvol->depth; s++) {
    for (c = 0; c < gtm->yvol->width; c++) {
      for (r = 0; r < gtm->yvol->height; r++) {
        if (gtm->mask && MRIgetVoxVal(gtm->mask, c, r, s, 0) < 0.5) continue;
        segid = MRIgetVoxVal(gtm->gtmseg, c, r, s, 0);
        k++;  // have to do this here
        if (segid == 0) continue;
        for (nthseg = 0; nthseg < gtm->nsegs; nthseg++)
          if (segid == gtm->segidlist[nthseg]) break;
        for (mthseg = 0; mthseg < gtm->nsegs; mthseg++) {
          mthsegid = gtm->segidlist[mthseg];
          tt = gtm->ctGTMSeg->entries[mthsegid]->TissueType;
          // printf("k=%d, segid = %d, nthseg = %d, mthsegid = %d, mthseg = %d, tt=%d\n",
          // k,segid,nthseg,mthsegid,mthseg,tt);
          fflush(stdout);
          gtm->ttpct->rptr[nthseg + 1][tt] +=  // not tt+1
              (gtm->X->rptr[k][mthseg + 1] * gtm->beta->rptr[mthseg + 1][1]);
        }
      }
    }
  }

  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    sum = 0;
    for (tt = 0; tt < nTT; tt++) sum += gtm->ttpct->rptr[nthseg + 1][tt + 1];  // yes, tt+1
    for (tt = 0; tt < nTT; tt++)
      gtm->ttpct->rptr[nthseg + 1][tt + 1] = 100 * gtm->ttpct->rptr[nthseg + 1][tt + 1] / sum;
  }

  return (0);
}

/*
  \fn int GTMsegid2nthseg(GTM *gtm, int segid)
  \breif Returns the nthseg of the given segid
 */
int GTMsegid2nthseg(GTM *gtm, int segid)
{
  int nthseg, ok;
  ok = 0;
  for (nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
    if (segid == gtm->segidlist[nthseg]) {
      ok = 1;
      break;
    }
  }
  if (!ok) return (-1);
  return (nthseg);
}

/*!
  \fn int GTMwriteText(GTM *gtm, char *OutDir, int DeMean)
  \brief Save GTM values out to text files named after the seg. If
  DeMean==1 then the frames are demeaned. Does not create the 
  output dir. Returns 0 if no error.
 */
int GTMwriteText(GTM *gtm, char *OutDir, int DeMean)
{
  int nthseg, f, segid;
  FILE *fp;
  std::string fname;
  double mean;

  for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
    segid = gtm->segidlist[nthseg];
    mean = 0;
    if(DeMean){
      for(f=0; f < gtm->beta->cols; f++) mean += gtm->beta->rptr[nthseg+1][f+1];
      mean /= gtm->beta->cols;
    }
    fname = std::string(OutDir) + '/' +
      std::string(gtm->ctGTMSeg->entries[segid]->name) + ".dat";
    fp = fopen(fname.c_str(),"w");
    if(fp==NULL){
      printf("ERROR: GTMwriteText(): could not open %s\n",fname.c_str());
      return(1);
    }
    for(f=0; f < gtm->beta->cols; f++) 
      fprintf(fp,"%12.6f\n",gtm->beta->rptr[nthseg+1][f+1]-mean);
    fclose(fp);
  }
  return(0);
}
