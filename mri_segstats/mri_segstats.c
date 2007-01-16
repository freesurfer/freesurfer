/**
 * @file  mri_segstats.c
 * @brief Computes statistics from a segmentation.
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Dougas N Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/01/16 02:32:40 $
 *    $Revision: 1.25 $
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
#include <sys/utsname.h>

#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "version.h"
#include "cma.h"
#include "gca.h"
#include "fsenv.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);

typedef struct {
  int id;
  char name[1000];
  int nhits;
  float vol;
  float min, max, range, mean, std;
}
STATSUMENTRY;

int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg);
int *MRIsegIdList(MRI *seg, int *nlist, int frame);
int MRIsegCount(MRI *seg, int id, int frame);
int MRIsegStats(MRI *seg, int segid, MRI *mri,  int frame,
                float *min, float *max, float *range,
                float *mean, float *std);
int compare_ints(const void *v1,const void *v2);
int nunqiue_int_list(int *idlist, int nlist);
int *unqiue_int_list(int *idlist, int nlist, int *nunique);
STATSUMENTRY *LoadStatSumFile(char *fname, int *nsegid);
int DumpStatSumTable(STATSUMENTRY *StatSumTable, int nsegid);


int main(int argc, char *argv[]) ;

static char vcid[] =
  "$Id: mri_segstats.c,v 1.25 2007/01/16 02:32:40 greve Exp $";
char *Progname = NULL, *SUBJECTS_DIR = NULL, *FREESURFER_HOME=NULL;
char *SegVolFile = NULL;
char *InVolFile = NULL;
char *InIntensityName = "";
char *InIntensityUnits = "unknown";
char *MaskVolFile = NULL;
char *PVVolFile = NULL;
char *BrainMaskFile = NULL;
char *StatTableFile = NULL;
char *FrameAvgFile = NULL;
char *FrameAvgVolFile = NULL;
char *SpatFrameAvgFile = NULL;
int DoFrameAvg = 0;
int frame = 0;
int synth = 0;
int debug = 0;
int dontrun = 0;
long seed = 0;
MRI *seg, *invol, *famri, *maskvol, *pvvol, *brainvol;
int nsegid0, *segidlist0;
int nsegid, *segidlist;
int NonEmptyOnly = 0;
int UserSegIdList[1000];
int nUserSegIdList = 0;
int DoExclSegId = 0, ExclSegId = 0;
int DoExclCtxGMWM= 0;
int DoSurfCtxGMWM = 0;
double lhwhitevol, lhpialvol, lhctxvol;
double rhwhitevol, rhpialvol, rhctxvol;

char *gcafile = NULL;
GCA *gca;

float maskthresh = 0.5;
int   maskinvert = 0, maskframe = 0;
char *masksign=NULL;
int   nmaskhits;
int   nbrainsegvoxels = 0;
double brainsegvolume = 0;
int   nbrainmaskvoxels = 0;
double brainmaskvolume = 0;
int   BrainVolFromSeg = 0;
int   DoETIV = 0;

char *ctabfile = NULL;
COLOR_TABLE *ctab = NULL;
STATSUMENTRY *StatSumTable = NULL;
STATSUMENTRY *StatSumTable2 = NULL;

MRIS *mris;
char *subject = NULL;
char *hemi    = NULL;
char *annot   = NULL;


/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, n, n0, skip, nhits, f, nsegidrep, ind, nthsegid;
  int c,r,s,err;
  float voxelvolume,vol;
  float min, max, range, mean, std;
  FILE *fp;
  double  **favg, *favgmn;
  struct utsname uts;
  char *cmdline;
  char tmpstr[1000];
  double atlas_icv=0;
  int ntotalsegid=0;
  int valid;
  int usersegid=0;
  nhits = 0;
  vol = 0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  cmdline = argv2cmdline(argc,argv);
  uname(&uts);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if (subject != NULL) {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if (DoETIV) {
    // calc total intracranial volume estimation
    sprintf
    (tmpstr,
     "%s/%s/mri/transforms/talairach_with_skull.lta",
     SUBJECTS_DIR,
     subject);
    atlas_icv = MRIestimateTIV(tmpstr,NULL,NULL);
    printf("atlas_icv = %g\n",atlas_icv);
  }

  /* Make sure we can open the output summary table file*/
  fp = fopen(StatTableFile,"w");
  if (fp == NULL) {
    printf("ERROR: could not open %s for writing\n",StatTableFile);
    exit(1);
  }
  fclose(fp);

  /* Make sure we can open the output frame average file*/
  if (FrameAvgFile != NULL) {
    fp = fopen(FrameAvgFile,"w");
    if (fp == NULL) {
      printf("ERROR: could not open %s for writing\n",FrameAvgFile);
      exit(1);
    }
    fclose(fp);
  }

  /* Load the segmentation */
  if (SegVolFile) {
    printf("Loading %s\n",SegVolFile);
    seg = MRIread(SegVolFile);
    if (seg == NULL) {
      printf("ERROR: loading %s\n",SegVolFile);
      exit(1);
    }
  } else {
    printf("Constructing seg from annotation\n");
    sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
    mris = MRISread(tmpstr);
    if (mris==NULL) exit(1);
    sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot);
    err = MRISreadAnnotation(mris, tmpstr);
    if (err) exit(1);
    seg = MRISannotIndex2Seg(mris);
    // Now create a colortable in a temp location to be read out below (hokey)
    if (mris->ct) {
      sprintf(tmpstr,"/tmp/mri_segstats.tmp.%s.%s.%d.ctab",subject,hemi,
              nint(randomNumber(0, 255)));
      ctabfile = strcpyalloc(tmpstr);
      CTABwriteFileASCII(mris->ct,ctabfile);
    }
  }
  if (ctabfile != NULL) {
    /* Load the color table file */
    ctab = CTABreadASCII(ctabfile);
    if (ctab == NULL) {
      printf("ERROR: reading %s\n",ctabfile);
      exit(1);
    }
  }
  if (gcafile != NULL) {
    gca = GCAread(gcafile);
    if (gca == NULL) {
      printf("ERROR: reading %s\n",gcafile);
      exit(1);
    }
    ctab = GCAcolorTableCMA(gca);
  }

  if(DoSurfCtxGMWM){
    printf("Getting Cerebral GM and WM volumes from surfaces\n");

    sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lhwhitevol = MRISvolumeInSurf(mris);
    MRISfree(&mris);

    sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lhpialvol = MRISvolumeInSurf(mris);
    lhctxvol = lhpialvol - lhwhitevol;
    MRISfree(&mris);

    sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rhwhitevol = MRISvolumeInSurf(mris);
    MRISfree(&mris);

    sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rhpialvol = MRISvolumeInSurf(mris);
    rhctxvol = rhpialvol - rhwhitevol;
    MRISfree(&mris);
    mris = NULL;

    printf("lh surface-based volumes (mm3): w = %lf,  p = %lf c = %lf \n",
	   lhwhitevol,lhpialvol,lhctxvol);
    printf("rh surface-based volumes (mm3): w = %lf,  p = %lf c = %lf \n",
	   rhwhitevol,rhpialvol,rhctxvol);
    fflush(stdout);
  }


  /* Load the input volume */
  if (InVolFile != NULL) {
    printf("Loading %s\n",InVolFile);
    invol = MRIread(InVolFile);
    if (invol == NULL) {
      printf("ERROR: loading %s\n",InVolFile);
      exit(1);
    }
    if (frame >= invol->nframes) {
      printf("ERROR: input frame = %d, input volume only has %d frames\n",
             frame,invol->nframes);
      exit(1);
    }
    /* Should check that invol the same dim as seg, etc*/
  }

  /* Load the partial volume mri */
  if (PVVolFile != NULL) {
    printf("Loading %s\n",PVVolFile);
    pvvol = MRIread(PVVolFile);
    if (pvvol == NULL) {
      printf("ERROR: loading %s\n",PVVolFile);
      exit(1);
    }
    /* Should check that invol the same dim as seg, etc*/
  }

  /* Load the brain volume */
  if (BrainMaskFile != NULL) {
    printf("Loading %s\n",BrainMaskFile);
    brainvol = MRIread(BrainMaskFile);
    if (brainvol == NULL) {
      printf("ERROR: loading %s\n",BrainMaskFile);
      exit(1);
    }
    /* Should check that invol the same dim as seg, etc*/
    nbrainmaskvoxels = MRItotalVoxelsOn(brainvol, WM_MIN_VAL) ;
    brainmaskvolume =
      nbrainmaskvoxels * brainvol->xsize * brainvol->ysize * brainvol->zsize;
    MRIfree(&brainvol) ;
    printf("# nbrainmaskvoxels %d\n",nbrainmaskvoxels);
    printf("# brainmaskvolume %10.1lf\n",brainmaskvolume);
  }

  if (BrainVolFromSeg) {
    nbrainsegvoxels = 0;
    for (n = 0 ; n <= MAX_CMA_LABEL ; n++) {
      if (!IS_BRAIN(n)) continue ;
      nbrainsegvoxels += MRIvoxelsInLabel(seg, n) ;
    }
    brainsegvolume = nbrainsegvoxels * seg->xsize * seg->ysize * seg->zsize;
    printf("# nbrainsegvoxels %d\n",nbrainsegvoxels);
    printf("# brainsegvolume %10.1lf\n",brainsegvolume);
  }

  /* Load the mask volume */
  if (MaskVolFile != NULL) {
    printf("Loading %s\n",MaskVolFile);
    maskvol = MRIread(MaskVolFile);
    if (maskvol == NULL) {
      printf("ERROR: loading %s\n",MaskVolFile);
      exit(1);
    }
    if (maskframe >= maskvol->nframes) {
      printf("ERROR: mask frame = %d, mask volume only has %d frames\n",
             maskframe,maskvol->nframes);
      exit(1);
    }
    /* Should check that maskvol the same dim as seg, etc*/
    mri_binarize(maskvol, maskthresh, masksign, maskinvert,
                 maskvol, &nmaskhits);
    if (nmaskhits == 0) {
      printf("ERROR: no voxels in mask meet thresholding criteria\n");
      exit(1);
    }
    printf("There were %d voxels in the mask\n",nmaskhits);

    /* perform the masking */
    for (c=0; c < seg->width; c++) {
      for (r=0; r < seg->height; r++) {
        for (s=0; s < seg->depth; s++) {
          // Set voxels out of the mask to 0
          if (! (int)MRIgetVoxVal(maskvol,c,r,s,maskframe))
            MRIsetVoxVal(seg,c,r,s,0,0);
        }
      }
    }
  }

  if (!mris) {
    voxelvolume = seg->xsize * seg->ysize * seg->zsize;
    printf("Voxel Volume is %g mm^3\n",voxelvolume);
  } else {
    if (mris->group_avg_surface_area > 0)
      voxelvolume = mris->group_avg_surface_area/mris->nvertices;
    else
      voxelvolume = mris->total_area/mris->nvertices;
    printf("Vertex Area is %g mm^3\n",voxelvolume);
  }

  /* There are three ways that the list of segmentations
     can be specified:
     1. User does not specify, then get all from the seg itself
     2. User specfies with --id (can be multiple)
     3. User supplies a color table
     If the user specficies a color table and --id, then the
     segs from --id are used ang the color table is only
     used to determine the name of the segmentation.
  */

  printf("Generating list of segmentation ids\n");
  segidlist0 = MRIsegIdList(seg, &nsegid0,0);

  if (ctab == NULL && nUserSegIdList == 0) {
    /* Must get list of segmentation ids from segmentation itself*/
    segidlist = segidlist0;
    nsegid = nsegid0;
    StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
    for (n=0; n < nsegid; n++) {
      StatSumTable[n].id = segidlist[n];
      strcpy(StatSumTable[n].name, "\0");
    }
  } else { /* Get from user or color table */
    if (ctab != NULL) {
      if (nUserSegIdList == 0) {
        /* User has not spec anything, so use all the ids in the color table */
        /* We want to fill StatSumTable with all the valid entries
           from the color table. So we'll get the number of valid
           entries and create StatSumTable with that many
           elements. Then walk through the entirity of the ctab and
           skip past the invalid entries. Copy the valid entries into
           StatSumTable. We use a separate index with StatSumTable
           that only goes from 0->the number of valid entries.*/
        CTABgetNumberOfValidEntries(ctab,&nsegid);
        CTABgetNumberOfTotalEntries(ctab,&ntotalsegid);
        StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
        usersegid=0;
        for (n=0; n < ntotalsegid; n++) {
          CTABisEntryValid(ctab,n,&valid);
          if (!valid)
            continue;
          StatSumTable[usersegid].id = n;
          CTABcopyName(ctab,n,StatSumTable[usersegid].name,
                       sizeof(StatSumTable[usersegid].name));
          usersegid++;
        }
      } else {
        /* User has specified --id, use those and get names from ctab */
        nsegid = nUserSegIdList;
        StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
        for (n=0; n < nsegid; n++) {
          StatSumTable[n].id = UserSegIdList[n];
          /* Here ind should be the same as the ctab entry, but make
             sure it's valid. */
          ind = StatSumTable[n].id;
          CTABisEntryValid(ctab,ind,&valid);
          if (!valid) {
            printf("ERROR: cannot find seg id %d in %s\n",
                   StatSumTable[n].id,ctabfile);
            exit(1);
          }
          CTABcopyName(ctab,ind,StatSumTable[n].name,
                       sizeof(StatSumTable[n].name));
        }
      }
    } else { /* User specified ids, but no color table */
      nsegid = nUserSegIdList;
      StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
      for (n=0; n < nsegid; n++)
        StatSumTable[n].id = UserSegIdList[n];
    }
  }

  printf("Found %3d segmentations\n",nsegid);
  printf("Computing statistics for each segmentation\n");
  fflush(stdout);
  for (n=0; n < nsegid; n++) {
    if(DoExclSegId && StatSumTable[n].id == ExclSegId) continue;
    if(DoExclCtxGMWM && 
       (StatSumTable[n].id ==  2 || StatSumTable[n].id == 3 ||
	StatSumTable[n].id == 41 || StatSumTable[n].id == 42) )
      continue; // excludes Cortical GM and WM

    printf("%3d   %3d  %s ",n,StatSumTable[n].id,StatSumTable[n].name);
    fflush(stdout);

    // Skip ones that are not represented
    skip = 1;
    for (n0=0; n0 < nsegid0; n0++)
      if (StatSumTable[n].id == segidlist0[n0]) skip = 0;
    if (skip) {
      printf(" 0\n");
      continue;
    }

    if (!dontrun) {
      if (!mris) {
        if (pvvol == NULL)
          nhits = MRIsegCount(seg, StatSumTable[n].id, 0);
        else
          nhits =
            MRIvoxelsInLabelWithPartialVolumeEffects
            (seg, pvvol, StatSumTable[n].id);
        vol = nhits*voxelvolume;
      } else {
        // Compute area here
        nhits = 0;
        vol = 0;
        for (c=0; c < mris->nvertices; c++) {
          if (MRIgetVoxVal(seg,c,0,0,0)==StatSumTable[n].id) {
            nhits++;
            if (mris->group_avg_vtxarea_loaded)
              vol += mris->vertices[c].group_avg_area;
            else
              vol += mris->vertices[c].area;
          }
        }
      }
    } else  nhits = n;

    printf("%4d  %g\n",nhits,vol);
    fflush(stdout);
    StatSumTable[n].nhits = nhits;
    StatSumTable[n].vol = vol;
    if (InVolFile != NULL && !dontrun) {
      if (nhits > 0) {
        MRIsegStats(seg, StatSumTable[n].id, invol, frame,
                    &min, &max, &range, &mean, &std);
      } else {
        min=0;
        max=0;
        range=0;
        mean=0;
        std=0;
      }
      StatSumTable[n].min   = min;
      StatSumTable[n].max   = max;
      StatSumTable[n].range = range;
      StatSumTable[n].mean  = mean;
      StatSumTable[n].std   = std;
    }
  }
  printf("\n");


  /* Remove empty segmentations, if desired */
  if (NonEmptyOnly || DoExclSegId) {
    // Count the number of nonempty segmentations
    nsegidrep = 0;
    for (n=0; n < nsegid; n++) {
      if (NonEmptyOnly && StatSumTable[n].nhits==0) continue;
      if (DoExclSegId && StatSumTable[n].id==ExclSegId) continue;
      if(DoExclCtxGMWM && 
	 (StatSumTable[n].id ==  2 || StatSumTable[n].id == 3 ||
	  StatSumTable[n].id == 41 || StatSumTable[n].id == 42) )
	continue; // excludes Cortical GM and WM
      nsegidrep ++;
    }
    StatSumTable2 = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegidrep);
    nthsegid = 0;
    for (n=0; n < nsegid; n++) {
      if(NonEmptyOnly && StatSumTable[n].nhits==0) continue;
      if(DoExclSegId && StatSumTable[n].id==ExclSegId) continue;
      if(DoExclCtxGMWM && 
	 (StatSumTable[n].id ==  2 || StatSumTable[n].id == 3 ||
	  StatSumTable[n].id == 41 || StatSumTable[n].id == 42) )
	continue; // excludes Cortical GM and WM
      StatSumTable2[nthsegid].id    = StatSumTable[n].id;
      StatSumTable2[nthsegid].nhits = StatSumTable[n].nhits;
      StatSumTable2[nthsegid].vol   = StatSumTable[n].vol;
      StatSumTable2[nthsegid].min   = StatSumTable[n].min;
      StatSumTable2[nthsegid].max   = StatSumTable[n].max;
      StatSumTable2[nthsegid].range = StatSumTable[n].range;
      StatSumTable2[nthsegid].mean  = StatSumTable[n].mean;
      StatSumTable2[nthsegid].std   = StatSumTable[n].std;
      strcpy(StatSumTable2[nthsegid].name,StatSumTable[n].name);
      nthsegid++;
    }
    free(StatSumTable);
    StatSumTable = StatSumTable2;
    nsegid = nsegidrep;
  }
  printf("Reporting on %3d segmentations\n",nsegid);

  /* Dump the table to the screen */
  if (debug) {
    for (n=0; n < nsegid; n++) {
      printf("%3d  %8d %10.1f  ", StatSumTable[n].id,StatSumTable[n].nhits,
             StatSumTable[n].vol);
      if (ctab != NULL) printf("%-30s ",StatSumTable[n].name);
      if (InVolFile != NULL)
        printf("%10.4f %10.4f %10.4f %10.4f %10.4f ",
               StatSumTable[n].min, StatSumTable[n].max,
               StatSumTable[n].range, StatSumTable[n].mean,
               StatSumTable[n].std);
      printf("\n");
    }
  }

  /* Print the table to the output file */
  if (StatTableFile != NULL) {
    fp = fopen(StatTableFile,"w");
    fprintf(fp,"# Title Segmentation Statistics \n");
    fprintf(fp,"# \n");
    fprintf(fp,"# generating_program %s\n",Progname);
    fprintf(fp,"# cvs_version %s\n",vcid);
    fprintf(fp,"# cmdline %s\n",cmdline);
    fprintf(fp,"# sysname  %s\n",uts.sysname);
    fprintf(fp,"# hostname %s\n",uts.nodename);
    fprintf(fp,"# machine  %s\n",uts.machine);
    fprintf(fp,"# user     %s\n",VERuser());
    if (mris) fprintf(fp,"# anatomy_type surface\n");
    else     fprintf(fp,"# anatomy_type volume\n");
    fprintf(fp,"# \n");
    if (subject != NULL) {
      fprintf(fp,"# SUBJECTS_DIR %s\n",SUBJECTS_DIR);
      fprintf(fp,"# subjectname %s\n",subject);
    }
    if (BrainMaskFile) {
      fprintf(fp,"# BrainMaskFile  %s \n",BrainMaskFile);
      fprintf(fp,"# BrainMaskFileTimeStamp  %s \n",
              VERfileTimeStamp(BrainMaskFile));
      fprintf(fp,"# Measure BrainMask, BrainMaskNVox, "
              "Number of Brain Mask Voxels, %7d, unitless\n",
              nbrainmaskvoxels);
      fprintf(fp,"# Measure BrainMask, BrainMaskVol, "
              "Brain Mask Volume, %f, mm^3\n",brainmaskvolume);
    }
    if (BrainVolFromSeg) {
      fprintf(fp,"# Measure BrainSeg, BrainSegNVox, "
              "Number of Brain Segmentation Voxels, %7d, unitless\n",
              nbrainsegvoxels);
      fprintf(fp,"# Measure BrainSeg, BrainSegVol, "
              "Brain Segmentation Volume, %f, mm^3\n",
              brainsegvolume);
    }
    if (DoETIV) {
      fprintf(fp,"# Measure IntraCranialVol, ICV, "
              "Intracranial Volume, %f, mm^3\n",atlas_icv);
    }
    if (SegVolFile) {
      fprintf(fp,"# SegVolFile %s \n",SegVolFile);
      fprintf(fp,"# SegVolFileTimeStamp  %s \n",VERfileTimeStamp(SegVolFile));
    }
    if (annot) fprintf(fp,"# Annot %s %s %s\n",subject,hemi,annot);
    if (ctabfile) {
      fprintf(fp,"# ColorTable %s \n",ctabfile);
      fprintf(fp,"# ColorTableTimeStamp %s \n",VERfileTimeStamp(ctabfile));
    }
    if (gcafile) {
      fprintf(fp,"# ColorTableFromGCA %s \n",gcafile);
      fprintf(fp,"# GCATimeStamp %s \n",VERfileTimeStamp(gcafile));
    }
    if (MaskVolFile) {
      fprintf(fp,"# MaskVolFile  %s \n",MaskVolFile);
      fprintf(fp,"#   MaskVolFileTimeStamp  %s \n",
              VERfileTimeStamp(MaskVolFile));
      fprintf(fp,"#   MaskThresh %f \n",maskthresh);
      fprintf(fp,"#   MaskSign   %s \n",masksign);
      fprintf(fp,"#   MaskFrame  %d \n",maskframe);
      fprintf(fp,"#   MaskInvert %d \n",maskframe);
    }
    if (InVolFile) {
      fprintf(fp,"# InVolFile  %s \n",InVolFile);
      fprintf(fp,"# InVolFileTimeStamp  %s \n",VERfileTimeStamp(InVolFile));
      fprintf(fp,"# InVolFrame %d \n",frame);
    }
    if (PVVolFile) {
      fprintf(fp,"# PVVolFile  %s \n",PVVolFile);
      fprintf(fp,"# PVVolFileTimeStamp  %s \n",VERfileTimeStamp(PVVolFile));
    }
    if(DoSurfCtxGMWM){
      fprintf(fp,"# surface-based-volume mm3 lh-cerebral-white-matter %lf\n",lhwhitevol);
      fprintf(fp,"# surface-based-volume mm3 lh-cerebral-cortex       %lf\n",lhctxvol);
      fprintf(fp,"# surface-based-volume mm3 rh-cerebral-white-matter %lf\n",rhwhitevol);
      fprintf(fp,"# surface-based-volume mm3 rh-cerebral-cortex       %lf\n",rhctxvol);
    }
    if(DoExclSegId)  fprintf(fp,"# ExcludeSegId %d \n",ExclSegId);
    if(DoExclCtxGMWM)fprintf(fp,"# Excluding Cortical Gray and White Matter\n");
    if(NonEmptyOnly) fprintf(fp,"# Only reporting non-empty segmentations\n");
    if(!mris) fprintf(fp,"# VoxelVolume_mm3 %g \n",voxelvolume);
    else      fprintf(fp,"# VertexArea_mm2 %g \n",voxelvolume);
    fprintf(fp,"# TableCol  1 ColHeader Index \n");
    fprintf(fp,"# TableCol  1 FieldName Index \n");
    fprintf(fp,"# TableCol  1 Units     NA \n");

    fprintf(fp,"# TableCol  2 ColHeader SegId \n");
    fprintf(fp,"# TableCol  2 FieldName Segmentation Id\n");
    fprintf(fp,"# TableCol  2 Units     NA\n");
    if (!mris) {
      fprintf(fp,"# TableCol  3 ColHeader NVoxels \n");
      fprintf(fp,"# TableCol  3 FieldName Number of Voxels\n");
      fprintf(fp,"# TableCol  3 Units     unitless\n");
      fprintf(fp,"# TableCol  4 ColHeader Volume_mm3\n");
      fprintf(fp,"# TableCol  4 FieldName Volume\n");
      fprintf(fp,"# TableCol  4 Units     mm^3\n");
    } else {
      fprintf(fp,"# TableCol  3 ColHeader NVertices \n");
      fprintf(fp,"# TableCol  3 FieldName Number of Vertices\n");
      fprintf(fp,"# TableCol  3 Units     unitless\n");
      fprintf(fp,"# TableCol  4 ColHeader Area_mm2\n");
      fprintf(fp,"# TableCol  4 FieldName Area\n");
      fprintf(fp,"# TableCol  4 Units     mm^2\n");
    }
    n = 5;
    if (ctab) {
      fprintf(fp,"# TableCol %2d ColHeader StructName\n",n);
      fprintf(fp,"# TableCol %2d FieldName Structure Name\n",n);
      fprintf(fp,"# TableCol %2d Units     NA\n",n);
      n++;
    }

    if (InVolFile) {
      fprintf(fp,"# TableCol %2d ColHeader %sMean \n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Intensity %sMean\n",
              n,InIntensityName);
      fprintf(fp,"# TableCol %2d Units     %s\n",n,InIntensityUnits);
      n++;

      fprintf(fp,"# TableCol %2d ColHeader %sStdDev\n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Itensity %sStdDev\n",
              n,InIntensityName);
      fprintf(fp,"# TableCol %2d Units     %s\n",n,InIntensityUnits);
      n++;

      fprintf(fp,"# TableCol %2d ColHeader %sMin\n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Intensity %sMin\n",
              n,InIntensityName);
      fprintf(fp,"# TableCol %2d Units     %s\n",n,InIntensityUnits);
      n++;

      fprintf(fp,"# TableCol %2d ColHeader %sMax\n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Intensity %sMax\n",
              n,InIntensityName);
      fprintf(fp,"# TableCol %2d Units     %s\n",n,InIntensityUnits);
      n++;

      fprintf(fp,"# TableCol %2d ColHeader %sRange\n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Intensity %sRange\n",
              n,InIntensityName);
      fprintf(fp,"# TableCol %2d Units     %s\n",n,InIntensityUnits);
      n++;

    }
    fprintf(fp,"# NRows %d \n",nsegid);
    fprintf(fp,"# NTableCols %d \n",n-1);

    fprintf(fp,"# ColHeaders  Index SegId ");
    if (!mris) fprintf(fp,"NVoxels Volume_mm3 ");
    else      fprintf(fp,"NVertices Area_mm2 ");
    if (ctab) fprintf(fp,"StructName ");
    if (InVolFile) fprintf(fp,"%sMean %sStdDev %sMin %sMax %sRange  ",
                             InIntensityName,
                             InIntensityName,
                             InIntensityName,
                             InIntensityName,
                             InIntensityName);
    fprintf(fp,"\n");

    for (n=0; n < nsegid; n++) {
      fprintf(fp,"%3d %3d  %8d %10.1f  ", n+1, StatSumTable[n].id,
              StatSumTable[n].nhits, StatSumTable[n].vol);
      if (ctab != NULL) fprintf(fp,"%-30s ",StatSumTable[n].name);
      if (InVolFile != NULL)
        fprintf(fp,"%10.4f %10.4f %10.4f %10.4f %10.4f ",
                StatSumTable[n].mean, StatSumTable[n].std,
                StatSumTable[n].min, StatSumTable[n].max,
                StatSumTable[n].range);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  // Average input across space to create a waveform
  // for each segmentation
  if (DoFrameAvg) {
    printf("Computing spatial average of each frame\n");
    favg = (double **) calloc(sizeof(double *),nsegid);
    for (n=0; n < nsegid; n++)
      favg[n] = (double *) calloc(sizeof(double),invol->nframes);
    favgmn = (double *) calloc(sizeof(double *),nsegid);
    for (n=0; n < nsegid; n++) {
      printf("%3d",n);
      if (n%20 == 19) printf("\n");
      fflush(stdout);
      MRIsegFrameAvg(seg, StatSumTable[n].id, invol, favg[n]);
      favgmn[n] = 0.0;
      for(f=0; f < invol->nframes; f++) favgmn[n] += favg[n][f];
      favgmn[n] /= invol->nframes;
    }
    printf("\n");

    // Save mean over space and frames in simple text file
    // Each seg on a separate line
    if(SpatFrameAvgFile) {
      printf("Writing to %s\n",SpatFrameAvgFile);
      fp = fopen(SpatFrameAvgFile,"w");
      for (n=0; n < nsegid; n++) fprintf(fp,"%g\n",favgmn[n]);
      fclose(fp);
    }

    // Save as a simple text file
    if (FrameAvgFile) {
      printf("Writing to %s\n",FrameAvgFile);
      fp = fopen(FrameAvgFile,"w");
      for (f=0; f < invol->nframes; f++) {
        fprintf(fp,"%3d %7.3f ",f,f*invol->tr/1000);
        for (n=0; n < nsegid; n++) fprintf(fp,"%g ",favg[n][f]);
        fprintf(fp,"\n");
      }
      fclose(fp);
    }

    // Save as an MRI "volume"
    if (FrameAvgVolFile) {
      printf("Writing to %s\n",FrameAvgVolFile);
      famri = MRIallocSequence(nsegid,1,1,MRI_FLOAT,invol->nframes);
      for (f=0; f < invol->nframes; f++) {
        for (n=0; n < nsegid; n++)
          MRIsetVoxVal(famri,n,0,0,f,(float)favg[n][f]);
      }
      MRIwrite(famri,FrameAvgVolFile);
    }
  }// Done with Frame Average

  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--dontrun"))   dontrun = 1;
    else if (!strcasecmp(option, "--nonempty")) NonEmptyOnly = 1;
    else if ( !strcmp(option, "--brain-vol-from-seg") ) BrainVolFromSeg = 1;
    else if ( !strcmp(option, "--etiv") ) DoETIV = 1;
    else if ( !strcmp(option, "--excl-ctxgmwm") ) DoExclCtxGMWM = 1;
    else if ( !strcmp(option, "--surf-ctxgmwm") ) DoSurfCtxGMWM = 1;

    else if ( !strcmp(option, "--sd") ) {
      if(nargc < 1) argnerr(option,1);
      FSENVsetSUBJECTS_DIR(pargv[0]);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--ctab-default") ) {
      FREESURFER_HOME = getenv("FREESURFER_HOME");
      ctabfile = (char *) calloc(sizeof(char),1000);
      sprintf(ctabfile,"%s/FreeSurferColorLUT.txt",FREESURFER_HOME);
      printf("Using defalt ctab %s\n",ctabfile);
    } else if ( !strcmp(option, "--ctab-gca") ) {
      if (nargc < 1) argnerr(option,1);
      gcafile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--seg") ) {
      if (nargc < 1) argnerr(option,1);
      SegVolFile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--in") ) {
      if (nargc < 1) argnerr(option,1);
      InVolFile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--in-intensity-name") ) {
      if (nargc < 1) argnerr(option,1);
      InIntensityName = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--in-intensity-units") ) {
      if (nargc < 1) argnerr(option,1);
      InIntensityUnits = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--brainmask") ) {
      if (nargc < 1) argnerr(option,1);
      BrainMaskFile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--id") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&UserSegIdList[nUserSegIdList]);
      nUserSegIdList++;
      nargsused = 1;
    } else if ( !strcmp(option, "--excludeid") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&ExclSegId);
      DoExclSegId = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--mask") ) {
      if (nargc < 1) argnerr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--masksign")) {
      if (nargc < 1) argnerr(option,1);
      masksign = pargv[0];
      nargsused = 1;
      if (strncasecmp(masksign,"abs",3) &&
          strncasecmp(masksign,"pos",3) &&
          strncasecmp(masksign,"neg",3)) {
        fprintf(stderr,"ERROR: mask sign = %s, must be abs, pos, or neg\n",
                masksign);
        exit(1);
      }
    } else if (!strcmp(option, "--maskthresh")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&maskthresh);
      nargsused = 1;
    } else if (!strcmp(option, "--maskframe")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&maskframe);
      nargsused = 1;
    } else if (!strcasecmp(option, "--maskinvert"))  maskinvert = 1;

    else if ( !strcmp(option, "--sum") ) {
      if (nargc < 1) argnerr(option,1);
      StatTableFile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--sum-in") ) {
      if (nargc < 1) argnerr(option,1);
      StatTableFile = pargv[0];
      StatSumTable = LoadStatSumFile(StatTableFile,&nsegid);
      printf("Found %d\n",nsegid);
      DumpStatSumTable(StatSumTable,nsegid);
      exit(1);
      nargsused = 1;
    } else if ( !strcmp(option, "--avgwf") ) {
      if (nargc < 1) argnerr(option,1);
      FrameAvgFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--sfavg") ) {
      if (nargc < 1) argnerr(option,1);
      SpatFrameAvgFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--avgwfvol") ) {
      if (nargc < 1) argnerr(option,1);
      FrameAvgVolFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--ctab") ) {
      if (nargc < 1) argnerr(option,1);
      ctabfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--frame")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    } else if (!strcmp(option, "--subject")) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annot")) {
      if (nargc < 3) argnerr(option,1);
      subject = pargv[0];
      hemi    = pargv[1];
      annot   = pargv[2];
      nargsused = 3;
    } else if (!strcmp(option, "--synth")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      synth = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--pv") ) {
      if (nargc < 1) argnerr(option,1);
      PVVolFile = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --seg segvol : segmentation volume path \n");
  printf("   --annot subject hemi annot : not fully tested yet\n");
  printf("\n");
  printf("   --sum file   : stats summary table file \n");
  printf("\n");
  printf(" Other Options\n");
  printf("   --pv pvvol : use pvvol to compensate for partial voluming\n");
  printf("   --in invol : report more stats on the input volume\n");
  printf("   --frame frame : report stats on nth frame of input volume\n");
  printf("\n");
  printf("   --ctab ctabfile : color table file with seg id names\n");
  printf("   --ctab-default: use $FREESURFER_HOME/FreeSurferColorLUT.txt\n");
  printf("   --ctab-gca gcafile: get color table from GCA (CMA)\n");
  printf("   --id segid <--id segid> : manually specify seg ids\n");
  printf("   --excludeid segid : exclude seg id from report\n");
  printf("   --excl-ctxgmwm : exclude cortical gray and white matter\n");
  printf("   --surf-ctxgmwm : compute coritcal gray and white volumes from surf\n");
  printf("   --nonempty : only report non-empty segmentations\n");
  printf("\n");
  printf("Masking options\n");
  printf("   --mask maskvol : must be same size as seg \n");
  printf("   --maskthresh thresh : binarize mask with this threshold <0.5>\n");
  printf("   --masksign sign : <abs>,pos,neg\n");
  printf("   --maskframe frame : 0-based frame number <0>\n");
  printf("   --maskinvert : invert mask \n");
  printf("\n");
  printf("Brain volume options\n");
  printf("   --brain-vol-from-seg : "
         "get brain volume from brain segmentations\n");
  printf("   --brainmask brainmask: "
         "compute volume from non-zero vox in brain mask\n");
  printf("   --etiv : compute intracranial volume "
         "from subject/mri/transfomrs/talairach_with_skull.lta\n");
  printf("\n");
  printf("Average waveform options\n");
  printf("   --avgwf textfile  : save into an ascii file\n");
  printf("   --avgwfvol mrivol : save as a binary mri 'volume'\n");
  printf("   --sfavg textfile  : save mean across space and frame\n");
  printf("\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf
  (
    "\n"
    "Help Outline:\n"
    "  - SUMMARY\n"
    "  - COMMAND-LINE ARGUMENTS\n"
    "  - SPECIFYING SEGMENTATION IDS\n"
    "  - SUMMARY FILE FORMAT\n"
    "  - EXAMPLES\n"
    "  - SEE ALSO\n"
    "\n"
    "SUMMARY\n"
    "\n"
    "This program will comute statistics on segmented volumes. In its\n"
    "simplist invocation, it will report on the number of voxels and\n"
    "volume in each segmentation. However, it can also compute statistics\n"
    "on the segmentation based on the values from another volume. This\n"
    "includes computing waveforms averaged inside each segmentation.\n"
    "\n"
    "COMMAND-LINE ARGUMENTS\n"
    "\n"
    "--seg segvol\n"
    "\n"
    "Input segmentation volume. A segmentation is a volume whose voxel\n"
    "values indicate a segmentation or class. This can be as complicaated\n"
    "as a FreeSurfer automatic cortical or subcortial segmentation or as\n"
    "simple as a binary mask. The format of segvol can be anything that\n"
    "mri_convert accepts as input (eg, analyze, nifti, mgh, bhdr, bshort, \n"
    "bfloat).\n"
    "\n"
    "--sum summaryfile\n"
    "\n"
    "ASCII file in which summary statistics are saved. See SUMMARY FILE\n"
    "below for more information.\n"
    "\n"
    "--pv pvvol\n"
    "\n"
    "Use pvvol to compensate for partial voluming. This should result in\n"
    "more accurate volumes. Usually, this is only done when computing \n"
    "anatomical statistics. Usually, the mri/norm.mgz volume is used.\n"
    "\n"
    "--in invol\n"
    "\n"
    "Input volume from which to compute more statistics, including min,\n"
    "max, range, average, and standard deviation as measured spatially\n"
    "across each segmentation. The input volume must be the same size\n"
    "and dimension as the segmentation volume.\n"
    "\n"
    "--frame frame\n"
    "\n"
    "Report statistics of the input volume at the 0-based frame number.\n"
    "frame is 0 be default.\n"
    "\n"
    "--ctab ctabfile\n"
    "\n"
    "FreeSurfer color table file. This is a file used by FreeSurfer to \n"
    "specify how each segmentation index is mapped to a segmentation\n"
    "name and color. See $FREESURFER_HOME/FreeSurferColorLUT.txt "
    "for example.\n"
    "The ctab can be used to specify the segmentations to report on or\n"
    "simply to supply human-readable names to segmentations chosen with\n"
    "--id. See SPECIFYING SEGMENTATION IDS below.\n"
    "\n"
    "--ctab-default\n"
    "\n"
    "Same as --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt\n"
    "\n"
    "--ctab-gca gcafile\n"
    "\n"
    "Get color table from the given GCA file. Eg,\n"
    "   $FREESURFER_HOME/average/RB_all_2006-02-15.gca\n"
    "This can be convenient when the seg file is that produced by\n"
    "mri_ca_label (ie, aseg.mgz) as it will only report on those \n"
    "segmentations that were actually considered during mri_ca_label.\n"
    "Note that there can still be some labels do not have any voxels \n"
    "in the report.\n"
    "\n"
    "--id segid1 <--id segid2>\n"
    "\n"
    "Specify numeric segmentation ids. Multiple ids can be specified with\n"
    "multiple --id invocations. SPECIFYING SEGMENTATION IDS.\n"
    "\n"
    "--excludeid segid\n"
    "\n"
    "Exclude the given segmentation id from report. This can be convenient\n"
    "for removing id=0. Only one segid can be targeted for exclusion.\n"
    "\n"
    "--excl-ctxgmwm\n"
    "\n"
    "Exclude cortical gray and white matter. These are assumed to be IDs\n"
    "2, 3, 41, and 42. The volume structures are more accurately measured\n"
    "using surface-based methods (see mris_volume).\n"
    "\n"
    "--surf-ctxgmwm\n"
    "\n"
    "Compute cortical gray and white matter volumes based on the white and\n"
    "pial surfaces. This is more accurate than from the aseg. The aseg \n"
    "values for these are still reported in the table, but there will be\n"
    "the following lines in the table:\n"
    "\n"
    "  # surface-based-volume mm3 lh-cerebral-white-matter 266579.428518\n"
    "  # surface-based-volume mm3 lh-cerebral-cortex       230267.243704\n"
    "  # surface-based-volume mm3 rh-cerebral-white-matter 265945.120671\n"
    "  # surface-based-volume mm3 rh-cerebral-cortex       236389.131763\n"
    "\n"
    "--nonempty\n"
    "\n"
    "Only report on segmentations that have actual representations in the\n"
    "segmentation volume.\n"
    "\n"
    "--mask maskvol\n"
    "\n"
    "Exlude voxels that are not in the mask. Voxels to be excluded are\n"
    "assigned a segid of 0. The mask volume may be binary or continuous.\n"
    "The masking criteria is set by the mask threshold, sign, frame, and\n"
    "invert parameters (see below). The mask volume must be the same\n"
    "size and dimension as the segmentation volume. If no voxels meet \n"
    "the masking criteria, then mri_segstats exits with an error.\n"
    "\n"
    "--maskthresh thresh\n"
    "\n"
    "Exlude voxels that are below thresh (for pos sign), above -thresh (for\n"
    "neg sign), or between -thresh and +thresh (for abs sign). Default\n"
    "is 0.5.\n"
    "\n"
    "--masksign sign\n"
    "\n"
    "Specify sign for masking threshold. Choices are abs, pos, and neg. \n"
    "Default is abs.\n"
    "\n"
    "--maskframe frame\n"
    "\n"
    "Derive the mask volume from the 0-based frameth frame.\n"
    "\n"
    "--maskinvert\n"
    "\n"
    "After applying all the masking criteria, invert the mask.\n"
    "\n"
    "--brain-vol-from-seg\n"
    "\n"
    "Get volume of brain as the sum of the volumes "
    "of the segmentations that\n"
    "are in the brain. Based on CMA/FreeSurferColorLUT.txt. "
    "The number of voxels\n"
    "and brain volume are stored as values in the "
    "header of the summary file\n"
    "with tags nbrainsegvoxels and brainsegvolume.\n"
    "\n"
    "--brainmask brainmask\n"
    "\n"
    "Load brain mask and compute the volume of the brain as the non-zero\n"
    "voxels in this volume. The number of voxels and "
    "brain volume are stored \n"
    "as values in the header of the summary file with "
    "tags nbrainmaskvoxels \n"
    "and brainmaskvolume.\n"
    "\n"
    "--avgwf textfile\n"
    "\n"
    "For each segmentation, compute an average waveform across all the\n"
    "voxels in the segmentation (excluding voxels masked out). The results\n"
    "are saved in an ascii text file with number of rows equal to the\n"
    "number of frames and number of columns equal to the number of\n"
    "segmentations reported plus 2. The first two columns are: (1) 0-based\n"
    "frame number and (2) 0-based frame number times TR.\n"
    "\n"
    "--avgwfvol mrivol\n"
    "\n"
    "Same as --avgwf except that the resulting waveforms are stored in a\n"
    "binary mri volume format (eg, analyze, nifti, mgh, etc) with number of\n"
    "columns equal to the number segmentations, number of rows = slices =\n"
    "1, and the number of frames equal that of the input volume. This may\n"
    "be more convenient than saving as an ascii text file.\n"
    "\n"
    "--help\n"
    "\n"
    "Don't get me started ...\n"
    "\n"
    "SPECIFYING SEGMENTATION IDS\n"
    "\n"
    "There are three ways that the list of segmentations to report on\n"
    "can be specified:\n"
    "  1. User specfies with --id.\n"
    "  2. User supplies a color table but does not specify --id. All\n"
    "     the segmentations in the color table are then reported on.\n"
    "     If the user specficies a color table and --id, then the\n"
    "     segids from --id are used and the color table is only\n"
    "     used to determine the name of the segmentation for reporint\n"
    "     purposes.\n"
    "  3. If the user does not specify either --id or a color table, then \n"
    "     all the ids from the segmentation volume are used.\n"
    "This list can be further reduced by specifying masks, --nonempty,\n"
    "and --excludeid.\n"
    "\n"
    "SUMMARY FILE FORMAT\n"
    "\n"
    "The summary file is an ascii file in which the segmentation statistics\n"
    "are reported. This file will have some 'header' information. Each\n"
    "header line begins with a '#'. There will be a row for each\n"
    "segmentation reported. The number and meaning of the columns depends\n"
    "somewhat how the program was run. The indentity of each column is\n"
    "given in the header. The first col is the row number. The second col\n"
    "is the segmentation id. The third col is the number of voxels in the\n"
    "segmentation. The fourth col is the volume of the segmentation in\n"
    "mm. If a color table was specified, then the next column will be the\n"
    "segmentation name. If an input volume was specified, then the next\n"
    "five columns will be intensity min, max, range, average, and standard\n"
    "deviation measured across the voxels in the segmentation.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "1. mri_segstats --seg $SUBJECTS_DIR/bert/mri/aseg \n"
    "    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \n"
    "    --nonempty --excludeid 0 --sum bert.aseg.sum \n"
    "\n"
    "This will compute the segmentation statistics from the automatic\n"
    "FreeSurfer subcortical segmentation for non-empty segmentations and\n"
    "excluding segmentation 0 (UNKNOWN). The results are stored in\n"
    "bert.aseg.sum.\n"
    "\n"
    "2. mri_segstats --seg $SUBJECTS_DIR/bert/mri/aseg \n"
    "    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \n"
    "    --nonempty --excludeid 0 --sum bert.aseg.sum \n"
    "    --in $SUBJECTS_DIR/bert/mri/orig\n"
    "\n"
    "Same as above but intensity statistics from the orig volume\n"
    "will also be reported for each segmentation.\n"
    "\n"
    "3. mri_segstats --seg aseg-in-func.img \n"
    "    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \n"
    "    --nonempty --excludeid 0 --in func.img \n"
    "    --mask spmT.img --maskthresh 2.3 \n"
    "    --sum bert.aseg-in-func.sum \n"
    "    --avgwf bert.avgwf.dat --avgwfvol bert.avgwf.img\n"
    "\n"
    "This will compute the segmentation statistics from the automatic\n"
    "FreeSurfer subcortical segmentation resampled into the functional\n"
    "space (see below and mri_label2vol --help). It will report intensity\n"
    "statistics from the 4D analyze volume func.img (same dimension as\n"
    "aseg-in-func.img). The segmentation is masked by thresholding the\n"
    "spmT.img map at 2.3. The average functional waveform of each\n"
    "segmentation is reported in the ascii file bert.avgwf.dat and in the\n"
    "4D analyze 'volume' bert.avgwf.img. This is not a real volume but just\n"
    "another way to save the data that may be more convenient than ascii.\n"
    "\n"
    "4. mri_label2vol --seg $SUBJECTS_DIR/bert/mri/aseg \n"
    "     --temp func.img --reg register.dat \n"
    "     --fillthresh 0.5 --o aseg-in-func.img\n"
    "\n"
    "This uses mri_label2vol to resample the automatic subcortical\n"
    "segmentation to the functional space. For more information\n"
    "see mri_label2vol --help.\n"
    "\n"
    "5. mri_label2vol --annot $SUBJECTS_DIR/bert/label/lh.aparc.annot \n"
    "     --temp func.img --reg register.dat --fillthresh 0.5 \n"
    "     --hemi lh --subject bert --proj frac 0 .1 1 \n"
    "     --o lh.aparc-in-func.img\n"
    "\n"
    "This uses mri_label2vol to resample the automatic cortical\n"
    "segmentation to the functional space. For more information\n"
    "see mri_label2vol --help.\n"
    "\n"
    "SEE ALSO:\n"
    "  mri_label2vol, tkregister2, mri_vol2roi.\n"
    "\n"
    "\n"
  );

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void) {
  if (SegVolFile == NULL && annot == NULL) {
    printf("ERROR: must specify a segmentation volume\n");
    exit(1);
  }
  if (StatTableFile == NULL) {
    printf("ERROR: must specify an output table file\n");
    exit(1);
  }
  if (DoFrameAvg && InVolFile == NULL) {
    printf("ERROR: cannot do frame average without input volume\n");
    exit(1);
  }
  if (DoETIV && subject == NULL) {
    printf("ERROR: need subject with --etiv\n");
    exit(1);
  }
  if (ctabfile != NULL && gcafile != NULL) {
    printf("ERROR: cannot specify ctab and gca\n");
    exit(1);
  }
  if(DoSurfCtxGMWM && subject == NULL){
    printf("ERROR: need --subject with --surf-gmwm\n");
    exit(1);
  }
  if (masksign == NULL) masksign = "abs";
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}

/* ----------------------------------------------------------
   MRIsegCount() - returns the number of times the given
   segmentation id appears in the volume.
   --------------------------------------------------------- */
int MRIsegCount(MRI *seg, int id, int frame) {
  int nhits, v, c,r,s;
  nhits = 0;
  for (c=0; c < seg->width; c++) {
    for (r=0; r < seg->height; r++) {
      for (s=0; s < seg->depth; s++) {
        v = (int) MRIgetVoxVal(seg,c,r,s,frame);
        if (v == id) nhits ++;
      }
    }
  }
  return(nhits);
}
/* ----------------------------------------------------------
   MRIsegIdList() - returns a list of the unique segmentation ids in
   the volume. The number in the list is *nlist. The volume need not
   be an int or char, but it is probably what it will be.
   --------------------------------------------------------- */
int *MRIsegIdList(MRI *seg, int *nlist, int frame) {
  int nvoxels,r,c,s,nth;
  int *tmplist = NULL;
  int *segidlist = NULL;

  nvoxels = seg->width * seg->height * seg->depth;
  tmplist = (int *) calloc(sizeof(int),nvoxels);

  // First, load all voxels into a list
  nth = 0;
  for (c=0; c < seg->width; c++) {
    for (r=0; r < seg->height; r++) {
      for (s=0; s < seg->depth; s++) {
        tmplist[nth] = (int) MRIgetVoxVal(seg,c,r,s,frame);
        nth++;
      }
    }
  }

  segidlist = unqiue_int_list(tmplist, nvoxels, nlist);
  free(tmplist);
  //for(nth=0; nth < *nlist; nth++)
  //printf("%3d %3d\n",nth,segidlist[nth]);
  return(segidlist);
}

/* --------------------------------------------- */
int compare_ints(const void *v1,const void *v2) {
  int i1, i2;

  i1 = *((int*)v1);
  i2 = *((int*)v2);

  if (i1 < i2) return(-1);
  if (i1 > i2) return(+1);
  return(0);
}
/* ---------------------------------------------------
   nunqiue_int_list() - counts the number of unique items
   in a list of integers. The list will be sorted.
   --------------------------------------------------- */
int nunqiue_int_list(int *idlist, int nlist) {
  int idprev, nunique, n;

  qsort(idlist,nlist,sizeof(int),compare_ints);
  nunique = 1;
  idprev = idlist[0];
  for (n=1; n<nlist; n++) {
    if (idprev != idlist[n]) {
      nunique++;
      idprev = idlist[n];
    }
  }
  return(nunique);
}
/* ---------------------------------------------------
   unqiue_int_list() - the returns the unique items
   in a list of integers. The list will be sorted.
   --------------------------------------------------- */
int *unqiue_int_list(int *idlist, int nlist, int *nunique) {
  int n, *ulist, nthu;

  /* count number of unique elements in the list,
     this also sorts the list */
  *nunique = nunqiue_int_list(idlist, nlist);

  /* alloc the unqiue list */
  ulist = (int *) calloc(sizeof(int),*nunique);

  nthu = 0;
  ulist[nthu] = idlist[0];
  for (n=1; n<nlist; n++) {
    if (ulist[nthu] != idlist[n]) {
      nthu ++;
      ulist[nthu] = idlist[n];
    }
  }
  return(ulist);
}

/*---------------------------------------------------------
  MRIsegStats() - computes statistics within a given
  segmentation. Returns the number of voxels in the
  segmentation.
  ---------------------------------------------------------*/
int MRIsegStats(MRI *seg, int segid, MRI *mri,int frame,
                float *min, float *max, float *range,
                float *mean, float *std) {
  int id,nvoxels,r,c,s;
  double val, sum, sum2;

  *min = 0;
  *max = 0;
  sum  = 0;
  sum2 = 0;
  nvoxels = 0;
  for (c=0; c < seg->width; c++) {
    for (r=0; r < seg->height; r++) {
      for (s=0; s < seg->depth; s++) {
        id = (int) MRIgetVoxVal(seg,c,r,s,0);
        if (id != segid) continue;
        val =  MRIgetVoxVal(mri,c,r,s,frame);
        nvoxels++;
        if ( nvoxels == 1 ) {
          *min = val;
          *max = val;
        }
        if (*min > val) *min = val;
        if (*max < val) *max = val;
        sum  += val;
        sum2 += (val*val);
      }
    }
  }

  *range = *max - *min;

  if (nvoxels != 0) *mean = sum/nvoxels;
  else             *mean = 0.0;

  if (nvoxels > 1)
    *std = sqrt(((nvoxels)*(*mean)*(*mean) - 2*(*mean)*sum + sum2)/
                (nvoxels-1));
  else *std = 0.0;

  return(nvoxels);
}
/*---------------------------------------------------------
  MRIsegFrameAvg() - computes the average time course withing the
  given segmentation. Returns the number of voxels in the
  segmentation. favg must be preallocated to number of
  frames. favg = (double *) calloc(sizeof(double),mri->nframes);
  ---------------------------------------------------------*/
int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg) {
  int id,nvoxels,r,c,s,f;
  double val;

  /* zero it out */
  for (f=0;f<mri->nframes;f++) favg[f] = 0;

  nvoxels = 0;
  for (c=0; c < seg->width; c++) {
    for (r=0; r < seg->height; r++) {
      for (s=0; s < seg->depth; s++) {
        id = (int) MRIgetVoxVal(seg,c,r,s,0);
        if (id != segid) continue;
        for (f=0;f<mri->nframes;f++) {
          val =  MRIgetVoxVal(mri,c,r,s,f);
          favg[f] += val;
        }
        nvoxels++;
      }
    }
  }

  if (nvoxels != 0)
    for (f=0;f<mri->nframes;f++) favg[f] /= nvoxels;

  return(nvoxels);
}

/*------------------------------------------------------------*/
STATSUMENTRY *LoadStatSumFile(char *fname, int *nsegid) {
  FILE *fp;
  char tmpstr[1000];
  STATSUMENTRY *StatSumTable, *e;

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",fname);
    exit(1);
  }

  // Count the number of entries
  *nsegid = 0;
  while (fgets(tmpstr,1000,fp) != NULL) {
    if (tmpstr[0] == '#') continue;
    (*nsegid)++;
  }
  fclose(fp);

  StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),*nsegid);

  // Now actually read it in
  fp = fopen(fname,"r");
  *nsegid = 0;
  while (fgets(tmpstr,1000,fp) != NULL) {
    if (tmpstr[0] == '#') continue;
    e = &StatSumTable[*nsegid];
    sscanf(tmpstr,"%*d %d %d %f %s %f %f %f %f %f",
           &e->id,&e->nhits,&e->vol,&e->name[0],
           &e->mean,&e->std,&e->min,&e->max,&e->range);
    (*nsegid)++;
  }
  fclose(fp);

  return(StatSumTable);
}
//----------------------------------------------------------------
int DumpStatSumTable(STATSUMENTRY *StatSumTable, int nsegid) {
  int n;
  for (n=0; n < nsegid; n++) {
    printf("%3d  %8d %10.1f  ", StatSumTable[n].id,StatSumTable[n].nhits,StatSumTable[n].vol);
    printf("%-30s ",StatSumTable[n].name);
    printf("%10.4f %10.4f %10.4f %10.4f %10.4f ",
           StatSumTable[n].min, StatSumTable[n].max,
           StatSumTable[n].range, StatSumTable[n].mean,
           StatSumTable[n].std);
    printf("\n");
  }
  return(0);
}

