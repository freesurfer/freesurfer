/**
 * @file  mri_segstats.c
 * @brief Computes statistics from a segmentation.
 *
 * This program will compute statistics on segmented volumes. In its
 * simplist invocation, it will report on the number of voxels and
 * volume in each segmentation. However, it can also compute statistics
 * on the segmentation based on the values from another volume. This
 * includes computing waveforms averaged inside each segmentation.
 */
/*
 * Original Author: Dougas N Greve
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/12/21 16:41:53 $
 *    $Revision: 1.72 $
 *
 * Copyright (C) 2006-2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

/*
  BEGINHELP

This program will comute statistics on segmented volumes. In its
simplist invocation, it will report on the number of voxels and volume
in each segmentation. However, it can also compute statistics on the
segmentation based on the values from another volume. This includes
computing waveforms averaged inside each segmentation. It can opperate
on both segmentation volumes (eg, aseg.mgz) and surface parcellations
(eg, lh.aparc.annot).

Help Outline:
  - COMMAND-LINE ARGUMENTS
  - SPECIFYING SEGMENTATION IDS
  - MEASURES OF BRAIN VOLUME
  - SUMMARY FILE FORMAT
  - EXAMPLES
  - SEE ALSO

COMMAND-LINE ARGUMENTS

--seg segvol

Input segmentation volume. A segmentation is a volume whose voxel
values indicate a segmentation or class. This can be as complicaated
as a FreeSurfer automatic cortical or subcortial segmentation or as
simple as a binary mask. The format of segvol can be anything that
mri_convert accepts as input (eg, analyze, nifti, mgh, bhdr, bshort, 
bfloat).

--annot subject hemi parc

Create a segmentation from hemi.parc.annot. If parc is aparc or
aparc.a2005s, then the segmentation numbers will match those in
$FREESURFER_HOME/FreeSurferColorLUT.txt (and so aparc+aseg.mgz). The
numbering can also be altered with --segbase. If an input is used, it
must be a surface ovelay with the same dimension as the parcellation.
This functionality makes mri_segstats partially redundant with
mris_anatomical_stats.

--label subject hemi labelfile

Create a segmentation from the given surface label. The points in 
the label are given a value of 1; 0 for outside.

--sum summaryfile

ASCII file in which summary statistics are saved. See SUMMARY FILE
below for more information.

--pv pvvol

Use pvvol to compensate for partial voluming. This should result in
more accurate volumes. Usually, this is only done when computing 
anatomical statistics. Usually, the mri/norm.mgz volume is used.
Not with --annot.

--i invol

Input volume from which to compute more statistics, including min,
max, range, average, and standard deviation as measured spatially
across each segmentation. The input volume must be the same size
and dimension as the segmentation volume.

--frame frame

Report statistics of the input volume at the 0-based frame number.
frame is 0 be default.

--sqr
--sqrt

Compute square or square root of input prior to computing stats.

--ctab ctabfile

FreeSurfer color table file. This is a file used by FreeSurfer to 
specify how each segmentation index is mapped to a segmentation
name and color. See $FREESURFER_HOME/FreeSurferColorLUT.txt for example.
The ctab can be used to specify the segmentations to report on or
simply to supply human-readable names to segmentations chosen with
--id. See SPECIFYING SEGMENTATION IDS below.

--ctab-default

Same as --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt

--ctab-gca gcafile

Get color table from the given GCA file. Eg,
   $FREESURFER_HOME/average/RB_all_YYYY-MM-DD.gca
This can be convenient when the seg file is that produced by
mri_ca_label (ie, aseg.mgz) as it will only report on those 
segmentations that were actually considered during mri_ca_label.
Note that there can still be some labels do not have any voxels 
in the report.

--ctab-out

Create an output color table (like FreeSurferColor.txt) with just
the segmentations reported in the output. 

--id segid1 <<segid2> <--id segid3>>

Specify numeric segmentation ids. Multiple ids can be specified with
multiple IDs after a single --id or with multiple --id invocations. 
SPECIFYING SEGMENTATION IDS.

--excludeid segid1 <segid2 ...>

Exclude the given segmentation id(s) from report. This can be convenient
for removing id=0. 

--excl-ctxgmwm

Exclude cortical gray and white matter. These are assumed to be IDs
2, 3, 41, and 42. The volume structures are more accurately measured
using surface-based methods (see mris_volume).

--surf-wm-vol

Compute cortical matter volume based on the volume encompassed by the 
white surface. This is more accurate than from the aseg. The aseg 
values for these are still reported in the table, but there will be
the following lines in the table:

  # surface-based-volume mm3 lh-cerebral-white-matter 266579.428518
  # surface-based-volume mm3 rh-cerebral-white-matter 265945.120671

--empty

Report on segmentations listed in the color table even if they 
are not found in the segmentation volume.

--mask maskvol

Exlude voxels that are not in the mask. Voxels to be excluded are
assigned a segid of 0. The mask volume may be binary or continuous.
The masking criteria is set by the mask threshold, sign, frame, and
invert parameters (see below). The mask volume must be the same
size and dimension as the segmentation volume. If no voxels meet 
the masking criteria, then mri_segstats exits with an error.

--maskthresh thresh

Exlude voxels that are below thresh (for pos sign), above -thresh (for
neg sign), or between -thresh and +thresh (for abs sign). Default
is 0.5.

--masksign sign

Specify sign for masking threshold. Choices are abs, pos, and neg. 
Default is abs.

--maskframe frame

Derive the mask volume from the 0-based frameth frame.

--maskinvert

After applying all the masking criteria, invert the mask.

--brain-vol-from-seg

Get volume of brain as the sum of the volumes of the segmentations that
are in the brain. Based on CMA/FreeSurferColorLUT.txt. The number of voxels
and brain volume are stored as values in the header of the summary file
with tags nbrainsegvoxels and brainsegvolume.

--brainmask brainmask

Load brain mask and compute the volume of the brain as the non-zero
voxels in this volume. The number of voxels and brain volume are stored 
as values in the header of the summary file with tags nbrainmaskvoxels 
and brainmaskvolume.

--avgwf textfile

For each segmentation, compute an average waveform across all the
voxels in the segmentation (excluding voxels masked out). The results
are saved in an ascii text file with number of rows equal to the
number of frames and number of columns equal to the number of
segmentations reported plus 2. The first row is -1 -1 then all
of the segid numbers. After that, the first two columns are: 
(1) 0-based frame number and (2) 0-based frame number times TR.

--avgwfvol mrivol

Same as --avgwf except that the resulting waveforms are stored in a
binary mri volume format (eg, analyze, nifti, mgh, etc) with number of
columns equal to the number segmentations, number of rows = slices =
1, and the number of frames equal that of the input volume. This may
be more convenient than saving as an ascii text file.

--help

As if

SPECIFYING SEGMENTATION IDS

There are three ways that the list of segmentations to report on
can be specified:
  1. User specfies with --id.
  2. User supplies a color table but does not specify --id. All
     the segmentations in the color table are then reported on.
     If the user specficies a color table and --id, then the
     segids from --id are used and the color table is only
     used to determine the name of the segmentation for reporint
     purposes.
  3. If the user does not specify either --id or a color table, then 
     all the ids from the segmentation volume are used.
This list can be further reduced by specifying masks and --excludeid.

MEASURES OF BRAIN VOLUME

There will be three measures of brain volume in the output summary file:
  (1) BrainSegNotVent - sum of the volume of the structures identified in 
      the aseg.mgz volume this will include cerebellum but not ventricles,
      CSF and dura. Includes partial volume compensation with --pv.
      This is probably the number you want to report.
  (2) BrainMask - total volume of non-zero voxels in brainmask.mgz. This will
      include cerebellum, ventricles, and possibly dura. This is probably not
      what you want to report.
  (3) BrainSeg - sum of the volume of the structures identified in the aseg.mgz
      volume. This will  include cerebellum and ventricles but should exclude
      dura. This does not include partial volume compensation, so 
      this number might be different than the sum of the segmentation volumes.
  (4) IntraCranialVol (ICV) - estimate of the intracranial volume based on the
      talairach transform. See surfer.nmr.mgh.harvard.edu/fswiki/eTIV for more
      details. This is the same measure as Estimated Total Intracranial Volume
      (eTIV).

SUMMARY FILE FORMAT

The summary file is an ascii file in which the segmentation statistics
are reported. This file will have some 'header' information. Each
header line begins with a '#'. There will be a row for each
segmentation reported. The number and meaning of the columns depends
somewhat how the program was run. The indentity of each column is
given in the header. The first col is the row number. The second col
is the segmentation id. The third col is the number of voxels in the
segmentation. The fourth col is the volume of the segmentation in
mm. If a color table was specified, then the next column will be the
segmentation name. If an input volume was specified, then the next
five columns will be intensity min, max, range, average, and standard
deviation measured across the voxels in the segmentation.

EXAMPLES

1. mri_segstats --seg $SUBJECTS_DIR/bert/mri/aseg 
    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt 
    --excludeid 0 --sum bert.aseg.sum 

This will compute the segmentation statistics from the automatic
FreeSurfer subcortical segmentation for non-empty segmentations and
excluding segmentation 0 (UNKNOWN). The results are stored in
bert.aseg.sum.

2. mri_segstats --seg $SUBJECTS_DIR/bert/mri/aseg 
    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt 
    --excludeid 0 --sum bert.aseg.sum 
    --i $SUBJECTS_DIR/bert/mri/orig

Same as above but intensity statistics from the orig volume
will also be reported for each segmentation.

3. mri_segstats --seg aseg-in-func.img 
    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt 
    --excludeid 0 --i func.img 
    --mask spmT.img --maskthresh 2.3 
    --sum bert.aseg-in-func.sum 
    --avgwf bert.avgwf.dat --avgwfvol bert.avgwf.img

This will compute the segmentation statistics from the automatic
FreeSurfer subcortical segmentation resampled into the functional
space (see below and mri_label2vol --help). It will report intensity
statistics from the 4D analyze volume func.img (same dimension as
aseg-in-func.img). The segmentation is masked by thresholding the
spmT.img map at 2.3. The average functional waveform of each
segmentation is reported in the ascii file bert.avgwf.dat and in the
4D analyze 'volume' bert.avgwf.img. This is not a real volume but just
another way to save the data that may be more convenient than ascii.

4. mri_label2vol --seg $SUBJECTS_DIR/bert/mri/aseg 
     --temp func.img --reg register.dat 
     --fillthresh 0.5 --o aseg-in-func.img

This uses mri_label2vol to resample the automatic subcortical
segmentation to the functional space. For more information
see mri_label2vol --help.

5. mri_label2vol --annot $SUBJECTS_DIR/bert/label/lh.aparc.annot 
     --temp func.img --reg register.dat --fillthresh 0.5 
     --hemi lh --subject bert --proj frac 0 .1 1 
     --o lh.aparc-in-func.img

This uses mri_label2vol to resample the automatic cortical
segmentation to the functional space. For more information
see mri_label2vol --help.

6. mri_segstats --annot bert lh aparc --i lh.thickness --sum lh.thickness.sum 

Produce a summary of the thickness in each parcellation of aparc. This 
will give the same mean thicknesses as that created by mris_anatomical_stats
and found in stats/lh.aparc.stats.


SEE ALSO:
  mri_label2vol, tkregister2, mri_vol2roi.

  ENDHELP
*/




/*
   Subcort stuff that needs to be removed from the surface-based white
   matter volume:

   16  Brain-Stem                            119  159  176    0

    4  Left-Lateral-Ventricle                120   18  134    0
   10  Left-Thalamus-Proper                    0  118   14    0
   11  Left-Caudate                          122  186  220    0
   12  Left-Putamen                          236   13  176    0
   13  Left-Pallidum                          12   48  255    0
   17  Left-Hippocampus                      220  216   20    0
   18  Left-Amygdala                         103  255  255    0
   26  Left-Accumbens-area                   255  165    0    0
   28  Left-VentralDC                        165   42   42    0

   43  Right-Lateral-Ventricle               120   18  134    0
   49  Right-Thalamus-Proper                   0  118   14    0
   50  Right-Caudate                         122  186  220    0
   51  Right-Putamen                         236   13  176    0
   52  Right-Pallidum                         13   48  255    0
   53  Right-Hippocampus                     220  216   20    0
   54  Right-Amygdala                        103  255  255    0
   58  Right-Accumbens-area                  255  165    0    0
   60  Right-VentralDC                       165   42   42    0

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/utsname.h>

#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "utils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "version.h"
#include "cma.h"
#include "gca.h"
#include "fsenv.h"
#include "annotation.h"
#include "registerio.h"
#include "cmdargs.h"

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
  int red, green, blue; // 0-255
  float min, max, range, mean, std, snr;
}
STATSUMENTRY;

int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg);
int MRIsegCount(MRI *seg, int id, int frame);
int MRIsegStats(MRI *seg, int segid, MRI *mri,  int frame,
                float *min, float *max, float *range,
                float *mean, float *std);
STATSUMENTRY *LoadStatSumFile(char *fname, int *nsegid);
int DumpStatSumTable(STATSUMENTRY *StatSumTable, int nsegid);


int main(int argc, char *argv[]) ;

static char vcid[] =
"$Id: mri_segstats.c,v 1.72 2010/12/21 16:41:53 rge21 Exp $";
char *Progname = NULL, *SUBJECTS_DIR = NULL, *FREESURFER_HOME=NULL;
char *SegVolFile = NULL;
char *InVolFile = NULL;
char *InVolRegFile = NULL;
MATRIX *InVolReg = NULL;
int InVolRegHeader = 0;
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
MRI *seg, *invol, *famri, *maskvol, *pvvol, *brainvol, *mri_aseg, *mri_ribbon;
int nsegid0, *segidlist0;
int nsegid, *segidlist;
int NonEmptyOnly = 1;
int UserSegIdList[1000];
int nUserSegIdList = 0;
int DoExclSegId = 0, nExcl = 0, ExclSegIdList[1000], ExclSegId;
int DoExclCtxGMWM= 0;
int DoSurfCtxVol = 0;
int DoSurfWMVol = 0;
double lhwhitevol;
double rhwhitevol;
double lhwhitevolTot, lhpialvolTot, lhctxvol;
double rhwhitevolTot, rhpialvolTot, rhctxvol;
int DoSupraTent = 0;
double SupraTentVol, SupraTentVolCor;

char *gcafile = NULL;
GCA *gca;

float maskthresh = 0.5;
int   maskinvert = 0, maskframe = 0;
char *masksign=NULL;
int   maskerode = 0;
int   nmaskhits;
int   nbrainsegvoxels = 0;
double brainsegvolume = 0;
double brainsegvolume2 = 0;
int DoSubCortGrayVol = 0;
double SubCortGrayVol = 0;
int DoTotalGrayVol = 0;
int   nbrainmaskvoxels = 0;
double brainmaskvolume = 0;
int   BrainVolFromSeg = 0;
int   DoETIV = 0;
int   DoETIVonly = 0;
int   DoOldETIVonly = 0;
char *talxfmfile = NULL;

char *ctabfile = NULL;
COLOR_TABLE *ctab = NULL;
STATSUMENTRY *StatSumTable = NULL;
STATSUMENTRY *StatSumTable2 = NULL;
char *ctabfileOut = NULL;

MRIS *mris;
char *subject = NULL;
char *hemi    = NULL;
char *annot   = NULL;

int Vox[3], DoVox = 0;
int  segbase = -1000;

int DoSquare = 0;
int DoSquareRoot = 0;
char *LabelFile = NULL;

int DoMultiply = 0;
double MultVal = 0;

int DoSNR = 0;
struct utsname uts;
char *cmdline, cwd[2000];

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, n, nx, n0, skip, nhits, f, nsegidrep, id, ind, nthsegid;
  int c,r,s,err,DoContinue;
  float voxelvolume,vol;
  float min, max, range, mean, std, snr;
  FILE *fp;
  double  **favg, *favgmn;
  char tmpstr[1000];
  double atlas_icv=0;
  int ntotalsegid=0;
  int valid;
  int usersegid=0;
  LABEL *label;
  MRI *tmp;
  MATRIX *vox2vox = NULL;
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

  if (DoETIV || DoETIVonly || DoOldETIVonly) {
    // calc total intracranial volume estimation
    // see this page:
    // http://surfer.nmr.mgh.harvard.edu/fswiki/eTIV
    // for info on determining the scaling factor.
    // a factor of 1948 was found to be best when using talairach.xfm
    // however when using talairach_with_skull.lta, 2150 was used
    double etiv_scale_factor = 1948.106;
    if (talxfmfile) {
      // path to talairach.xfm file spec'd on the command line
      sprintf(tmpstr,"%s",talxfmfile);
    } else {
      sprintf
      (tmpstr,
       "%s/%s/mri/transforms/talairach.xfm",
       SUBJECTS_DIR,
       subject);
    }
    if (DoOldETIVonly) {
      // back-door way to get the old way of calculating etiv, for debug
      sprintf
        (tmpstr,
         "%s/%s/mri/transforms/talairach_with_skull.lta",
         SUBJECTS_DIR,
         subject);
      etiv_scale_factor = 2150;
    }
    double determinant = 0;
    atlas_icv = MRIestimateTIV(tmpstr,etiv_scale_factor,&determinant);
    printf("atlas_icv (eTIV) = %d mm^3    (det: %3f )\n",
           (int)atlas_icv,determinant);
    if (DoETIVonly || DoOldETIVonly) exit(0);
  }

  /* Make sure we can open the output summary table file*/
  if(StatTableFile){
    fp = fopen(StatTableFile,"w");
    if (fp == NULL) {
      printf("ERROR: could not open %s for writing\n",StatTableFile);
      int err = errno;
      printf("Errno: %s\n", strerror(err) );
      exit(1);
    }
    fclose(fp);
  }

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
    if(DoVox){
      printf("Replacing seg with a single voxel at %d %d %d\n",
	     Vox[0],Vox[1],Vox[2]);
      for(c=0; c < seg->width; c++){
	for(r=0; r < seg->height; r++){
	  for(s=0; s < seg->depth; s++){
	    if(c == Vox[0] && r == Vox[1] && s == Vox[2])
	      MRIsetVoxVal(seg,c,r,s,0, 1);
	    else
	      MRIsetVoxVal(seg,c,r,s,0, 0);
	  }
	}
      }
    }
  } 
  else if(annot) {
    printf("Constructing seg from annotation\n");
    sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
    mris = MRISread(tmpstr);
    if (mris==NULL) exit(1);
    sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot);
    err = MRISreadAnnotation(mris, tmpstr);
    if (err) exit(1);
    if(segbase == -1000){
      // segbase has not been set with --segbase
      if(!strcmp(annot,"aparc")){
	if(!strcmp(hemi,"lh")) segbase = 1000;
	else                   segbase = 2000;
      }
      else if(!strcmp(annot,"aparc.a2005s")){
	if(!strcmp(hemi,"lh")) segbase = 1100;
	else                   segbase = 2100;
      }
      else segbase = 0;
    }
    printf("Seg base %d\n",segbase);
    seg = MRISannot2seg(mris,segbase);
    // Now create a colortable in a temp location to be read out below (hokey)
    mris->ct->idbase = segbase;
    if (mris->ct) {
      sprintf(tmpstr,"/tmp/mri_segstats.tmp.%s.%s.%d.ctab",subject,hemi,
              nint(randomNumber(0, 255)));
      ctabfile = strcpyalloc(tmpstr);
      CTABwriteFileASCII(mris->ct,ctabfile);
    }
  }
  else {
    printf("Constructing seg from label\n");
    label = LabelRead(NULL, LabelFile);
    if(label == NULL) exit(1);
    sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
    mris = MRISread(tmpstr);
    if (mris==NULL) exit(1);
    seg = MRIalloc(mris->nvertices,1,1,MRI_INT);
    for (n = 0; n < label->n_points; n++)
      MRIsetVoxVal(seg,label->lv[n].vno,0,0,0, 1);
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

  if(DoSurfWMVol){
    printf("Getting Cerebral WM volumes from surface\n");

    sprintf(tmpstr,"%s/%s/mri/aseg.mgz",SUBJECTS_DIR,subject);
    mri_aseg = MRIread(tmpstr);
    if(mri_aseg == NULL) exit(1);

    sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lhwhitevol = MRIScomputeWhiteVolume(mris, mri_aseg, 1.0/4.0);
    MRISfree(&mris);
    printf("lh white matter volume %g\n",lhwhitevol);

    sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rhwhitevol = MRIScomputeWhiteVolume(mris, mri_aseg, 1.0/4.0);
    MRISfree(&mris);
    printf("rh white matter volume %g\n",rhwhitevol);

    MRIfree(&mri_aseg);
  }

  if(DoSurfCtxVol){
    printf("Getting Cerebral GM and WM volumes from surfaces\n");
    // Does this include the non-cortical areas of the surface?

    sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lhwhitevolTot = MRISvolumeInSurf(mris);
    MRISfree(&mris);

    sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lhpialvolTot = MRISvolumeInSurf(mris);
    lhctxvol = lhpialvolTot - lhwhitevolTot;
    MRISfree(&mris);

    sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rhwhitevolTot = MRISvolumeInSurf(mris);
    MRISfree(&mris);

    sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rhpialvolTot = MRISvolumeInSurf(mris);
    rhctxvol = rhpialvolTot - rhwhitevolTot;
    MRISfree(&mris);
    mris = NULL;

    printf("lh surface-based volumes (mm3): wTot = %lf,  pTot = %lf c = %lf \n",
           lhwhitevolTot,lhpialvolTot,lhctxvol);
    printf("rh surface-based volumes (mm3): wTot = %lf,  pTot = %lf c = %lf \n",
           rhwhitevolTot,rhpialvolTot,rhctxvol);
    fflush(stdout);

    if(DoSupraTent) {
      printf("Computing SupraTentVolCor\n");
      sprintf(tmpstr,"%s/%s/mri/aseg.mgz",SUBJECTS_DIR,subject);
      mri_aseg = MRIread(tmpstr);
      if(mri_aseg == NULL) exit(1);
      sprintf(tmpstr,"%s/%s/mri/ribbon.mgz",SUBJECTS_DIR,subject);
      mri_ribbon = MRIread(tmpstr);
      if(mri_ribbon == NULL) exit(1);
      SupraTentVolCor = SupraTentorialVolCorrection(mri_aseg, mri_ribbon);
      SupraTentVol = SupraTentVolCor + lhpialvolTot + rhpialvolTot;
      printf("SupraTentVolCor = %8.3f\n",SupraTentVolCor);
      printf("SupraTentVol = %8.3f\n",SupraTentVol);
      MRIfree(&mri_aseg);
      MRIfree(&mri_ribbon);
    }
  }


  /* Load the input volume */
  if (InVolFile != NULL) {
    printf("Loading %s\n",InVolFile);
    fflush(stdout);
    invol = MRIread(InVolFile);
    if(invol == NULL) {
      printf("ERROR: loading %s\n",InVolFile);
      exit(1);
    }
    if(frame >= invol->nframes) {
      printf("ERROR: input frame = %d, input volume only has %d frames\n",
             frame,invol->nframes);
      exit(1);
    }
    if(InVolReg || InVolRegHeader){
      if(InVolReg){
	printf("Input Volume Registration:\n");
	MatrixPrint(stdout,InVolReg);
      }
      vox2vox = MRIvoxToVoxFromTkRegMtx(invol, seg, InVolReg);
      printf("Input Volume vox2vox:\n");
      MatrixPrint(stdout,vox2vox);
      printf("Allocating %d frames\n",invol->nframes);
      tmp = MRIcloneBySpace(seg,-1,invol->nframes);
      printf("Vol2Vol\n");
      err = MRIvol2VolVSM(invol,tmp,vox2vox,SAMPLE_NEAREST,-1,NULL);
      if(err) exit(1);
      MRIfree(&invol);
      invol = tmp;
    }
    if(MRIdimMismatch(invol,seg,0)) {
      printf("ERROR: dimension mismatch between input volume and seg\n");
      printf("  input %d %d %d\n",invol->width,invol->height,invol->depth);
      printf("  seg   %d %d %d\n",seg->width,seg->height,seg->depth);
      exit(1); 
    }
    if(DoMultiply) {
      printf("Multiplying input by %lf\n",MultVal);
      MRImultiplyConst(invol,MultVal,invol);
    }
    if(DoSquare){
      printf("Computing square of input\n");
      MRIsquare(invol, NULL, invol);
    }
    if(DoSquareRoot){
      printf("Computing square root of input\n");
      MRIsquareRoot(invol, NULL, invol);
    }
  }

  /* Load the partial volume mri */
  if (PVVolFile != NULL) {
    printf("Loading %s\n",PVVolFile);
    fflush(stdout);
    pvvol = MRIread(PVVolFile);
    if (pvvol == NULL) {
      printf("ERROR: loading %s\n",PVVolFile);
      exit(1);
    }
    if(MRIdimMismatch(pvvol,seg,0)) {
      printf("ERROR: dimension mismatch between PV volume and seg\n");
      printf("  pvvol %d %d %d\n",pvvol->width,pvvol->height,pvvol->depth);
      printf("  seg   %d %d %d\n",seg->width,seg->height,seg->depth);
      exit(1); 
    }
  }

  /* Load the brain volume */
  if (BrainMaskFile != NULL) {
    printf("Loading %s\n",BrainMaskFile);
    fflush(stdout);
    brainvol = MRIread(BrainMaskFile);
    if(brainvol == NULL) {
      printf("ERROR: loading %s\n",BrainMaskFile);
      exit(1);
    }
    if(MRIdimMismatch(brainvol,seg,0)) {
      printf("ERROR: dimension mismatch between brain volume and seg\n");
      exit(1); 
    }
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
    fflush(stdout);
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
    if(MRIdimMismatch(maskvol,seg,0)) {
      printf("ERROR: dimension mismatch between brain volume and seg\n");
      exit(1); 
    }
    mri_binarize(maskvol, maskthresh, masksign, maskinvert,
                 maskvol, &nmaskhits);
    if (nmaskhits == 0) {
      printf("WARNING: no voxels in mask meet thresholding criteria.\n");
      printf("The output table will be empty.\n");
      printf("thresh = %g, sign = %s, inv = %d\n",maskthresh,masksign,maskinvert);
      //exit(1);
    }
    printf("There were %d voxels in the orginal mask\n",nmaskhits);
    if(maskerode > 0){
      printf("Eroding %d voxels in 3d\n",maskerode);
      for(n=0; n<maskerode; n++) MRIerode(maskvol,maskvol);
    }

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
  fflush(stdout);
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
  } 
  else { /* Get from user or color table */
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
          if(!valid) continue;
          StatSumTable[usersegid].id = n;
          CTABcopyName(ctab,n,StatSumTable[usersegid].name,
                       sizeof(StatSumTable[usersegid].name));
	  CTABrgbAtIndexi(ctab, n, &StatSumTable[usersegid].red, 
			  &StatSumTable[usersegid].green, 
			  &StatSumTable[usersegid].blue);
          usersegid++;
        }
      } 
      else {
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
	  CTABrgbAtIndexi(ctab, ind, &StatSumTable[n].red, 
			  &StatSumTable[n].green, 
			  &StatSumTable[n].blue);
        }
      }
    } 
    else { /* User specified ids, but no color table */
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
    if(DoExclSegId){
      DoContinue = 0;
      for(nx=0; nx < nExcl; nx++) {
	if(StatSumTable[n].id == ExclSegIdList[nx]) {
	  DoContinue = 1;
	  break;
	}
      }
      if(DoContinue) continue;
    }
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
            (seg, pvvol, StatSumTable[n].id, NULL, NULL); 
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
	snr = mean/std;
      } else {
        min=0;
        max=0;
        range=0;
        mean=0;
        std=0;
        snr=0;
      }
      StatSumTable[n].min   = min;
      StatSumTable[n].max   = max;
      StatSumTable[n].range = range;
      StatSumTable[n].mean  = mean;
      StatSumTable[n].std   = std;
      StatSumTable[n].snr   = snr;
    }
  }
  printf("\n");


  /* Remove empty segmentations, if desired */
  if (NonEmptyOnly || DoExclSegId) {
    // Count the number of nonempty segmentations
    nsegidrep = 0;
    for (n=0; n < nsegid; n++) {
      if(NonEmptyOnly && StatSumTable[n].nhits==0) continue;
      if(DoExclSegId){
	DoContinue = 0;
	for(nx=0; nx < nExcl; nx++) {
	  if(StatSumTable[n].id == ExclSegIdList[nx]) {
	    DoContinue = 1;
	    break;
	  }
	}
	if(DoContinue) continue;
      }
      nsegidrep ++;
    }
    StatSumTable2 = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegidrep);
    nthsegid = 0;
    for (n=0; n < nsegid; n++) {
      if(NonEmptyOnly && StatSumTable[n].nhits==0) continue;

      if(DoExclSegId){
	DoContinue = 0;
	for(nx=0; nx < nExcl; nx++) {
	  if(StatSumTable[n].id == ExclSegIdList[nx]) {
	    DoContinue = 1;
	    break;
	  }
	}
	if(DoContinue) continue;
      }

      StatSumTable2[nthsegid].id    = StatSumTable[n].id;
      StatSumTable2[nthsegid].nhits = StatSumTable[n].nhits;
      StatSumTable2[nthsegid].vol   = StatSumTable[n].vol;
      StatSumTable2[nthsegid].min   = StatSumTable[n].min;
      StatSumTable2[nthsegid].max   = StatSumTable[n].max;
      StatSumTable2[nthsegid].range = StatSumTable[n].range;
      StatSumTable2[nthsegid].mean  = StatSumTable[n].mean;
      StatSumTable2[nthsegid].std   = StatSumTable[n].std;
      StatSumTable2[nthsegid].snr   = StatSumTable[n].snr;
      StatSumTable2[nthsegid].red   = StatSumTable[n].red;
      StatSumTable2[nthsegid].green = StatSumTable[n].green;
      StatSumTable2[nthsegid].blue  = StatSumTable[n].blue;
      strcpy(StatSumTable2[nthsegid].name,StatSumTable[n].name);
      nthsegid++;
    }
    free(StatSumTable);
    StatSumTable = StatSumTable2;
    nsegid = nsegidrep;
  }
  printf("Reporting on %3d segmentations\n",nsegid);

  if(BrainVolFromSeg) {
    brainsegvolume2 = 0.0;
    for(n=0; n < nsegid; n++) {
      id = StatSumTable[n].id;
      if(!IS_BRAIN(id)) continue ;
      if(IS_CSF(id) || IS_CSF_CLASS(id)) continue;
      brainsegvolume2 += StatSumTable[n].vol;
    }
  }
  if(DoSubCortGrayVol) {
    SubCortGrayVol = 0.0;
    for(n=0; n < nsegid; n++) {
      id = StatSumTable[n].id;
      if(! IsSubCorticalGray(id)) continue ;
      SubCortGrayVol += StatSumTable[n].vol;
    }
    printf("SubCortGrayVol = %g\n",SubCortGrayVol);
  }

  /* Dump the table to the screen */
  if (debug) {
    for (n=0; n < nsegid; n++) {
      printf("%3d  %8d %10.1f  ", StatSumTable[n].id,StatSumTable[n].nhits,
             StatSumTable[n].vol);
      if (ctab != NULL) printf("%-30s ",StatSumTable[n].name);
      if (InVolFile != NULL){
        printf("%10.4f %10.4f %10.4f %10.4f %10.4f ",
               StatSumTable[n].min, StatSumTable[n].max,
               StatSumTable[n].range, StatSumTable[n].mean,
               StatSumTable[n].std);
	if(DoSNR) printf("%10.4f ",StatSumTable[n].snr);
      }
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
      fprintf(fp,"# Measure BrainSegNotVent, BrainSegVolNotVent, "
              "Brain Segmentation Volume Without Ventricles, %f, mm^3\n",
              brainsegvolume2);
      fprintf(fp,"# Measure BrainSeg, BrainSegNVox, "
            "Number of Brain Segmentation Voxels, %7d, unitless\n",
            nbrainsegvoxels);
      fprintf(fp,"# Measure BrainSeg, BrainSegVol, "
              "Brain Segmentation Volume, %f, mm^3\n",
              brainsegvolume);
    }
    if(DoSurfCtxVol){
      // Does this include the non-cortical areas of the surface?
      fprintf(fp,"# Measure lhCortex, lhCortexVol, "
              "Left hemisphere cortical gray matter volume, %f, mm^3\n",lhctxvol);
      fprintf(fp,"# Measure rhCortex, rhCortexVol, "
              "Right hemisphere cortical gray matter volume, %f, mm^3\n",rhctxvol);
      fprintf(fp,"# Measure Cortex, CortexVol, "
              "Total cortical gray matter volume, %f, mm^3\n",lhctxvol+rhctxvol);
    }
    if(DoSurfWMVol){
      fprintf(fp,"# Measure lhCorticalWhiteMatter, lhCorticalWhiteMatterVol, "
              "Left hemisphere cortical white matter volume, %f, mm^3\n",lhwhitevol);
      fprintf(fp,"# Measure rhCorticalWhiteMatter, rhCorticalWhiteMatterVol, "
              "Right hemisphere cortical white matter volume, %f, mm^3\n",rhwhitevol);
      fprintf(fp,"# Measure CorticalWhiteMatter, CorticalWhiteMatterVol, "
              "Total cortical white matter volume, %f, mm^3\n",lhwhitevol+rhwhitevol);
    }
    if(DoSubCortGrayVol) {
      fprintf(fp,"# Measure SubCortGray, SubCortGrayVol, "
              "Subcortical gray matter volume, %f, mm^3\n",
              SubCortGrayVol);
    }
    if(DoTotalGrayVol) {
      fprintf(fp,"# Measure TotalGray, TotalGrayVol, "
              "Total gray matter volume, %f, mm^3\n",
              SubCortGrayVol+lhctxvol+rhctxvol);
    }
    if(DoSupraTent) {
      fprintf(fp,"# Measure SupraTentorial, SupraTentorialVol, "
              "Supratentorial volume, %f, mm^3\n",SupraTentVol);
              
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
    if (LabelFile) fprintf(fp,"# Label %s %s %s\n",subject,hemi,LabelFile);
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
    if(DoExclCtxGMWM)
      fprintf(fp,"# Excluding Cortical Gray and White Matter\n");
    if(DoExclSegId){
      fprintf(fp,"# ExcludeSegId ");
      for(nx=0; nx < nExcl; nx++) fprintf(fp,"%d ",ExclSegIdList[nx]);
      fprintf(fp,"\n");
    }
    if(NonEmptyOnly)
      fprintf(fp,"# Only reporting non-empty segmentations\n");
    if(!mris) fprintf(fp,"# VoxelVolume_mm3 %g \n",voxelvolume);
    else fprintf(fp,"# VertexArea_mm2 %g \n",voxelvolume);
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
    fprintf(fp,"StructName ");
    if(InVolFile){
      fprintf(fp,"%sMean %sStdDev %sMin %sMax %sRange  ",
                           InIntensityName, InIntensityName,
                           InIntensityName,InIntensityName,
                           InIntensityName);
      if(DoSNR) fprintf(fp,"%sSNR ",InIntensityName);
    }
    fprintf(fp,"\n");

    for (n=0; n < nsegid; n++) {
      fprintf(fp,"%3d %3d  %8d %10.1f  ", n+1, StatSumTable[n].id,
              StatSumTable[n].nhits, StatSumTable[n].vol);
      if(ctab != NULL) fprintf(fp,"%-30s ",StatSumTable[n].name);
      else             fprintf(fp,"Seg%04d ",StatSumTable[n].id);
      if (InVolFile != NULL){
        fprintf(fp,"%10.4f %10.4f %10.4f %10.4f %10.4f ",
                StatSumTable[n].mean, StatSumTable[n].std,
                StatSumTable[n].min, StatSumTable[n].max,
                StatSumTable[n].range);
	if(DoSNR) fprintf(fp,"%10.4f ",StatSumTable[n].snr);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  if(ctabfileOut != NULL){
    fp = fopen(ctabfileOut,"w");
    for (n=0; n < nsegid; n++)
      fprintf(fp,"%d %-30s %3d %3d %3d 0\n",StatSumTable[n].id,
	      StatSumTable[n].name,StatSumTable[n].red,
	      StatSumTable[n].green,StatSumTable[n].blue);
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
      for (n=0; n < nsegid; n++) {
	fprintf(fp,"%g\n",favgmn[n]);
	printf("%d %g\n",n,favgmn[n]);
      }
      fclose(fp);
    }

    // Save as a simple text file
    if(FrameAvgFile) {
      printf("Writing to %s\n",FrameAvgFile);
      fp = fopen(FrameAvgFile,"w");
      //fprintf(fp,"-1 -1 ");
      //for (n=0; n < nsegid; n++) fprintf(fp,"%4d ", StatSumTable[n].id);
      //fprintf(fp,"\n");
      for (f=0; f < invol->nframes; f++) {
        //fprintf(fp,"%3d %7.3f ",f,f*invol->tr/1000);
        for (n=0; n < nsegid; n++) fprintf(fp,"%g ",favg[n][f]);
        fprintf(fp,"\n");
      }
      fclose(fp);
    }

    // Save as an MRI "volume"
    if(FrameAvgVolFile) {
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
  int  nargc , nargsused, nth;
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

    if (!strcasecmp(option, "--help")||!strcasecmp(option, "--usage")
      ||!strcasecmp(option, "-h")||!strcasecmp(option, "-u"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--dontrun"))   dontrun = 1;
    else if (!strcasecmp(option, "--nonempty")) NonEmptyOnly = 1;
    else if (!strcasecmp(option, "--empty"))    NonEmptyOnly = 0;
    else if ( !strcmp(option, "--brain-vol-from-seg") ) BrainVolFromSeg = 1;
    else if ( !strcmp(option, "--subcortgray") ) DoSubCortGrayVol = 1;
    else if ( !strcmp(option, "--supratent") ) DoSupraTent = 1;
    else if ( !strcmp(option, "--totalgray") ) DoTotalGrayVol = 1;
    else if ( !strcmp(option, "--etiv") ) DoETIV = 1;
    else if ( !strcmp(option, "--etiv-only") ) DoETIVonly = 1;
    else if ( !strcmp(option, "--old-etiv-only") ) DoOldETIVonly = 1;
    else if ( !strcmp(option, "--surf-ctx-vol") ) DoSurfCtxVol = 1;
    else if ( !strcmp(option, "--surf-wm-vol") )  DoSurfWMVol = 1;
    else if ( !strcmp(option, "--sqr") )  DoSquare = 1;
    else if ( !strcmp(option, "--sqrt") )  DoSquareRoot = 1;
    else if ( !strcmp(option, "--snr") )  DoSNR = 1;

    else if ( !strcmp(option, "--mul") ){
      if(nargc < 1) argnerr(option,1);
      DoMultiply = 1;
      sscanf(pargv[0],"%lf",&MultVal);
      nargsused = 1;
    } else if ( !strcmp(option, "--talxfm") ) {
      if (nargc < 1) argnerr(option,1);
      talxfmfile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--sd") ) {
      if(nargc < 1) argnerr(option,1);
      FSENVsetSUBJECTS_DIR(pargv[0]);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--ctab-default") ) {
      FREESURFER_HOME = getenv("FREESURFER_HOME");
      ctabfile = (char *) calloc(sizeof(char),1000);
      sprintf(ctabfile,"%s/FreeSurferColorLUT.txt",FREESURFER_HOME);
      printf("Using defalt ctab %s\n",ctabfile);
    } 
    else if ( !strcmp(option, "--ctab-gca") ) {
      if (nargc < 1) argnerr(option,1);
      gcafile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--seg") ) {
      if (nargc < 1) argnerr(option,1);
      SegVolFile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--vox") ) {
      if (nargc < 1) argnerr(option,3);
      sscanf(pargv[0],"%d",&Vox[0]);
      sscanf(pargv[1],"%d",&Vox[1]);
      sscanf(pargv[2],"%d",&Vox[2]);
      DoVox = 1;
      nargsused = 3;
    } 
    else if ( !strcmp(option, "--in") || !strcmp(option, "--i") ) {
      if (nargc < 1) argnerr(option,1);
      InVolFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcmp(option, "--reg")) {
      if(nargc < 1) argnerr(option,1);
      InVolRegFile = pargv[0];
      InVolReg = regio_read_registermat(InVolRegFile);
      if(InVolReg == NULL) exit(1);
      nargsused = 1;
    } 
    else if(!strcmp(option, "--regheader")) {
      InVolRegHeader = 1;
    } 
    else if ( !strcmp(option, "--in-intensity-name") ) {
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
    } 
    else if ( !strcmp(option, "--id") ) {
      if(nargc < 1) argnerr(option,1);
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%d",&UserSegIdList[nUserSegIdList]);
	nUserSegIdList++;
	nth++;
      }
      nargsused = nth;
    } 
    else if ( !strcmp(option, "--excl-ctxgmwm") ) {
      DoExclSegId = 1;
      DoExclCtxGMWM = 1;
      ExclSegIdList[nExcl] =  2; nExcl++;
      ExclSegIdList[nExcl] =  3; nExcl++;
      ExclSegIdList[nExcl] = 41; nExcl++;
      ExclSegIdList[nExcl] = 42; nExcl++;
    }
    else if( !strcmp(option, "--excludeid") ||
	     !strcmp(option, "--exclude") ) {
      if(nargc < 1) argnerr(option,1);
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%d",&ExclSegIdList[nExcl]);
	nExcl ++;
	nth ++;
      }
      DoExclSegId = 1;
      nargsused = nth;
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
    } else if (!strcmp(option, "--maskerode")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&maskerode);
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
    } 
    else if ( !strcmp(option, "--ctab") ) {
      if (nargc < 1) argnerr(option,1);
      ctabfile = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--ctab-out") ) {
      if (nargc < 1) argnerr(option,1);
      ctabfileOut = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--frame")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    } else if (!strcmp(option, "--subject")) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--annot")) {
      if (nargc < 3) argnerr(option,1);
      subject = pargv[0];
      hemi    = pargv[1];
      annot   = pargv[2];
      nargsused = 3;
    } 
    else if (!strcmp(option, "--slabel")) {
      if (nargc < 3) argnerr(option,1);
      subject = pargv[0];
      hemi    = pargv[1];
      LabelFile = pargv[2];
      ExclSegId = 0;
      DoExclSegId = 1;
      nargsused = 3;
    } 
    else if (!strcmp(option, "--segbase")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&segbase);
      nargsused = 1;
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
  outputHelp(Progname);

#ifdef GREGT
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --seg segvol : segmentation volume path \n");
  printf("   --annot  subject hemi parc : use surface parcellation\n");
  printf("   --slabel subject hemi label : use surface label\n");
  printf("\n");
  printf("   --sum file   : stats summary table file \n");
  printf("\n");
  printf(" Other Options\n");
  printf("   --pv pvvol : use pvvol to compensate for partial voluming\n");
  printf("   --i invol : report more stats on the input volume\n");
  printf("   --frame frame : report stats on nth frame of input volume\n");
  printf("   --sqr  : compute the square of the input\n");
  printf("   --sqrt : compute the square root of the input\n");
  printf("   --mul val : multiply input by val\n");
  printf("   --snr : save mean/std as extra column in output table\n");
  printf("\n");
  printf("   --ctab ctabfile : color table file with seg id names\n");
  printf("   --ctab-default: use $FREESURFER_HOME/FreeSurferColorLUT.txt\n");
  printf("   --ctab-gca gcafile: get color table from GCA (CMA)\n");
  printf("   --id segid <segid2 ...> : manually specify seg ids\n");
  printf("   --excludeid segid : exclude seg id from report\n");
  printf("   --excl-ctxgmwm : exclude cortical gray and white matter\n");
  printf("   --surf-wm-vol : compute cortical white volume from surf\n");
  printf("   --surf-ctx-vol : compute cortical volumes from surf\n");
  printf("   --empty : report all segmentations in ctab, even if not in the seg\n");
  printf("   --ctab-out ctaboutput : create a ctab with only your segs\n");
  printf("\n");
  printf("Masking options\n");
  printf("   --mask maskvol : must be same size as seg \n");
  printf("   --maskthresh thresh : binarize mask with this threshold <0.5>\n");
  printf("   --masksign sign : <abs>,pos,neg\n");
  printf("   --maskframe frame : 0-based frame number <0>\n");
  printf("   --maskinvert : invert mask \n");
  printf("   --maskerode nerode : erode mask \n");
  printf("\n");
  printf("Brain volume options\n");
  printf("   --brain-vol-from-seg : "
         "get brain volume from brain segmentations\n");
  printf("   --brainmask brainmask: "
         "compute volume from non-zero vox in brain mask\n");
  printf("   --subcortgray : compute volume of subcortical gray matter");
  printf("   --totalgray : compute volume of total gray matter");
  printf("   --etiv : compute intracranial volume "
         "from subject/mri/transforms/talairach.xfm\n");
  printf("   --etiv-only : compute intracranial volume "
         "from subject/mri/transforms/talairach.xfm and exit\n");
  printf("   --old-etiv-only : compute intracranial volume "
         "from subject/mri/transforms/talairach_with_skull.lta and exit\n");
  printf("   --talxfm fname : specify path to talairach.xfm file for etiv\n");
  printf("\n");
  printf("Average waveform options\n");
  printf("   --avgwf textfile  : save into an ascii file\n");
  printf("   --avgwfvol mrivol : save as a binary mri 'volume'\n");
  printf("   --sfavg textfile  : save mean across space and frame\n");
  printf("\n");
  printf("   --vox C R S : replace seg with all 0s except at CRS\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
#endif
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
#ifdef GREGT
printf("\n");
printf("This program will comute statistics on segmented volumes. In its\n");
printf("simplist invocation, it will report on the number of voxels and volume\n");
printf("in each segmentation. However, it can also compute statistics on the\n");
printf("segmentation based on the values from another volume. This includes\n");
printf("computing waveforms averaged inside each segmentation. It can opperate\n");
printf("on both segmentation volumes (eg, aseg.mgz) and surface parcellations\n");
printf("(eg, lh.aparc.annot).\n");
printf("\n");
printf("Help Outline:\n");
printf("  - COMMAND-LINE ARGUMENTS\n");
printf("  - SPECIFYING SEGMENTATION IDS\n");
printf("  - MEASURES OF BRAIN VOLUME\n");
printf("  - SUMMARY FILE FORMAT\n");
printf("  - EXAMPLES\n");
printf("  - SEE ALSO\n");
printf("\n");
printf("COMMAND-LINE ARGUMENTS\n");
printf("\n");
printf("--seg segvol\n");
printf("\n");
printf("Input segmentation volume. A segmentation is a volume whose voxel\n");
printf("values indicate a segmentation or class. This can be as complicaated\n");
printf("as a FreeSurfer automatic cortical or subcortial segmentation or as\n");
printf("simple as a binary mask. The format of segvol can be anything that\n");
printf("mri_convert accepts as input (eg, analyze, nifti, mgh, bhdr, bshort, \n");
printf("bfloat).\n");
printf("\n");
printf("--annot subject hemi parc\n");
printf("\n");
printf("Create a segmentation from hemi.parc.annot. If parc is aparc or\n");
printf("aparc.a2005s, then the segmentation numbers will match those in\n");
printf("$FREESURFER_HOME/FreeSurferColorLUT.txt (and so aparc+aseg.mgz). The\n");
printf("numbering can also be altered with --segbase. If an input is used, it\n");
printf("must be a surface ovelay with the same dimension as the parcellation.\n");
printf("This functionality makes mri_segstats partially redundant with\n");
printf("mris_anatomical_stats.\n");
printf("\n");
printf("--label subject hemi labelfile\n");
printf("\n");
printf("Create a segmentation from the given surface label. The points in \n");
printf("the label are given a value of 1; 0 for outside.\n");
printf("\n");
printf("--sum summaryfile\n");
printf("\n");
printf("ASCII file in which summary statistics are saved. See SUMMARY FILE\n");
printf("below for more information.\n");
printf("\n");
printf("--pv pvvol\n");
printf("\n");
printf("Use pvvol to compensate for partial voluming. This should result in\n");
printf("more accurate volumes. Usually, this is only done when computing \n");
printf("anatomical statistics. Usually, the mri/norm.mgz volume is used.\n");
printf("Not with --annot.\n");
printf("\n");
printf("--i invol\n");
printf("\n");
printf("Input volume from which to compute more statistics, including min,\n");
printf("max, range, average, and standard deviation as measured spatially\n");
printf("across each segmentation. The input volume must be the same size\n");
printf("and dimension as the segmentation volume.\n");
printf("\n");
printf("--frame frame\n");
printf("\n");
printf("Report statistics of the input volume at the 0-based frame number.\n");
printf("frame is 0 be default.\n");
printf("\n");
printf("--sqr\n");
printf("--sqrt\n");
printf("\n");
printf("Compute square or square root of input prior to computing stats.\n");
printf("\n");
printf("--ctab ctabfile\n");
printf("\n");
printf("FreeSurfer color table file. This is a file used by FreeSurfer to \n");
printf("specify how each segmentation index is mapped to a segmentation\n");
printf("name and color. See $FREESURFER_HOME/FreeSurferColorLUT.txt for example.\n");
printf("The ctab can be used to specify the segmentations to report on or\n");
printf("simply to supply human-readable names to segmentations chosen with\n");
printf("--id. See SPECIFYING SEGMENTATION IDS below.\n");
printf("\n");
printf("--ctab-default\n");
printf("\n");
printf("Same as --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt\n");
printf("\n");
printf("--ctab-gca gcafile\n");
printf("\n");
printf("Get color table from the given GCA file. Eg,\n");
printf("   $FREESURFER_HOME/average/RB_all_YYYY-MM-DD.gca\n");
printf("This can be convenient when the seg file is that produced by\n");
printf("mri_ca_label (ie, aseg.mgz) as it will only report on those \n");
printf("segmentations that were actually considered during mri_ca_label.\n");
printf("Note that there can still be some labels do not have any voxels \n");
printf("in the report.\n");
printf("\n");
printf("--ctab-out\n");
printf("\n");
printf("Create an output color table (like FreeSurferColor.txt) with just\n");
printf("the segmentations reported in the output. \n");
printf("\n");
printf("--id segid1 <<segid2> <--id segid3>>\n");
printf("\n");
printf("Specify numeric segmentation ids. Multiple ids can be specified with\n");
printf("multiple IDs after a single --id or with multiple --id invocations. \n");
printf("SPECIFYING SEGMENTATION IDS.\n");
printf("\n");
printf("--excludeid segid1 <segid2 ...>\n");
printf("\n");
printf("Exclude the given segmentation id(s) from report. This can be convenient\n");
printf("for removing id=0. \n");
printf("\n");
printf("--excl-ctxgmwm\n");
printf("\n");
printf("Exclude cortical gray and white matter. These are assumed to be IDs\n");
printf("2, 3, 41, and 42. The volume structures are more accurately measured\n");
printf("using surface-based methods (see mris_volume).\n");
printf("\n");
printf("--surf-wm-vol\n");
printf("\n");
printf("Compute cortical matter volume based on the volume encompassed by the \n");
printf("white surface. This is more accurate than from the aseg. The aseg \n");
printf("values for these are still reported in the table, but there will be\n");
printf("the following lines in the table:\n");
printf("\n");
printf("  # surface-based-volume mm3 lh-cerebral-white-matter 266579.428518\n");
printf("  # surface-based-volume mm3 rh-cerebral-white-matter 265945.120671\n");
printf("\n");
printf("--empty\n");
printf("\n");
printf("Report on segmentations listed in the color table even if they \n");
printf("are not found in the segmentation volume.\n");
printf("\n");
printf("--mask maskvol\n");
printf("\n");
printf("Exlude voxels that are not in the mask. Voxels to be excluded are\n");
printf("assigned a segid of 0. The mask volume may be binary or continuous.\n");
printf("The masking criteria is set by the mask threshold, sign, frame, and\n");
printf("invert parameters (see below). The mask volume must be the same\n");
printf("size and dimension as the segmentation volume. If no voxels meet \n");
printf("the masking criteria, then mri_segstats exits with an error.\n");
printf("\n");
printf("--maskthresh thresh\n");
printf("\n");
printf("Exlude voxels that are below thresh (for pos sign), above -thresh (for\n");
printf("neg sign), or between -thresh and +thresh (for abs sign). Default\n");
printf("is 0.5.\n");
printf("\n");
printf("--masksign sign\n");
printf("\n");
printf("Specify sign for masking threshold. Choices are abs, pos, and neg. \n");
printf("Default is abs.\n");
printf("\n");
printf("--maskframe frame\n");
printf("\n");
printf("Derive the mask volume from the 0-based frameth frame.\n");
printf("\n");
printf("--maskinvert\n");
printf("\n");
printf("After applying all the masking criteria, invert the mask.\n");
printf("\n");
printf("--brain-vol-from-seg\n");
printf("\n");
printf("Get volume of brain as the sum of the volumes of the segmentations that\n");
printf("are in the brain. Based on CMA/FreeSurferColorLUT.txt. The number of voxels\n");
printf("and brain volume are stored as values in the header of the summary file\n");
printf("with tags nbrainsegvoxels and brainsegvolume.\n");
printf("\n");
printf("--brainmask brainmask\n");
printf("\n");
printf("Load brain mask and compute the volume of the brain as the non-zero\n");
printf("voxels in this volume. The number of voxels and brain volume are stored \n");
printf("as values in the header of the summary file with tags nbrainmaskvoxels \n");
printf("and brainmaskvolume.\n");
printf("\n");
printf("--avgwf textfile\n");
printf("\n");
printf("For each segmentation, compute an average waveform across all the\n");
printf("voxels in the segmentation (excluding voxels masked out). The results\n");
printf("are saved in an ascii text file with number of rows equal to the\n");
printf("number of frames and number of columns equal to the number of\n");
printf("segmentations reported plus 2. The first row is -1 -1 then all\n");
printf("of the segid numbers. After that, the first two columns are: \n");
printf("(1) 0-based frame number and (2) 0-based frame number times TR.\n");
printf("\n");
printf("--avgwfvol mrivol\n");
printf("\n");
printf("Same as --avgwf except that the resulting waveforms are stored in a\n");
printf("binary mri volume format (eg, analyze, nifti, mgh, etc) with number of\n");
printf("columns equal to the number segmentations, number of rows = slices =\n");
printf("1, and the number of frames equal that of the input volume. This may\n");
printf("be more convenient than saving as an ascii text file.\n");
printf("\n");
printf("--help\n");
printf("\n");
printf("As if\n");
printf("\n");
printf("SPECIFYING SEGMENTATION IDS\n");
printf("\n");
printf("There are three ways that the list of segmentations to report on\n");
printf("can be specified:\n");
printf("  1. User specfies with --id.\n");
printf("  2. User supplies a color table but does not specify --id. All\n");
printf("     the segmentations in the color table are then reported on.\n");
printf("     If the user specficies a color table and --id, then the\n");
printf("     segids from --id are used and the color table is only\n");
printf("     used to determine the name of the segmentation for reporint\n");
printf("     purposes.\n");
printf("  3. If the user does not specify either --id or a color table, then \n");
printf("     all the ids from the segmentation volume are used.\n");
printf("This list can be further reduced by specifying masks and --excludeid.\n");
printf("\n");
printf("MEASURES OF BRAIN VOLUME\n");
printf("\n");
printf("There will be three measures of brain volume in the output summary file:\n");
printf("  (1) BrainSegNotVent - sum of the volume of the structures identified in \n");
printf("      the aseg.mgz volume this will include cerebellum but not ventricles,\n");
printf("      CSF and dura. Includes partial volume compensation with --pv.\n");
printf("      This is probably the number you want to report.\n");
printf("  (2) BrainMask - total volume of non-zero voxels in brainmask.mgz. This will\n");
printf("      include cerebellum, ventricles, and possibly dura. This is probably not\n");
printf("      what you want to report.\n");
printf("  (3) BrainSeg - sum of the volume of the structures identified in the aseg.mgz\n");
printf("      volume. This will  include cerebellum and ventricles but should exclude\n");
printf("      dura. This does not include partial volume compensation, so \n");
printf("      this number might be different than the sum of the segmentation volumes.\n");
printf("  (4) IntraCranialVol (ICV) - estimate of the intracranial volume based on the\n");
printf("      talairach transform. See surfer.nmr.mgh.harvard.edu/fswiki/eTIV for more\n");
printf("      details. This is the same measure as Estimated Total Intracranial Volume\n");
printf("      (eTIV).\n");
printf("\n");
printf("SUMMARY FILE FORMAT\n");
printf("\n");
printf("The summary file is an ascii file in which the segmentation statistics\n");
printf("are reported. This file will have some 'header' information. Each\n");
printf("header line begins with a '#'. There will be a row for each\n");
printf("segmentation reported. The number and meaning of the columns depends\n");
printf("somewhat how the program was run. The indentity of each column is\n");
printf("given in the header. The first col is the row number. The second col\n");
printf("is the segmentation id. The third col is the number of voxels in the\n");
printf("segmentation. The fourth col is the volume of the segmentation in\n");
printf("mm. If a color table was specified, then the next column will be the\n");
printf("segmentation name. If an input volume was specified, then the next\n");
printf("five columns will be intensity min, max, range, average, and standard\n");
printf("deviation measured across the voxels in the segmentation.\n");
printf("\n");
printf("EXAMPLES\n");
printf("\n");
printf("1. mri_segstats --seg $SUBJECTS_DIR/bert/mri/aseg \n");
printf("    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \n");
printf("    --excludeid 0 --sum bert.aseg.sum \n");
printf("\n");
printf("This will compute the segmentation statistics from the automatic\n");
printf("FreeSurfer subcortical segmentation for non-empty segmentations and\n");
printf("excluding segmentation 0 (UNKNOWN). The results are stored in\n");
printf("bert.aseg.sum.\n");
printf("\n");
printf("2. mri_segstats --seg $SUBJECTS_DIR/bert/mri/aseg \n");
printf("    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \n");
printf("    --excludeid 0 --sum bert.aseg.sum \n");
printf("    --i $SUBJECTS_DIR/bert/mri/orig\n");
printf("\n");
printf("Same as above but intensity statistics from the orig volume\n");
printf("will also be reported for each segmentation.\n");
printf("\n");
printf("3. mri_segstats --seg aseg-in-func.img \n");
printf("    --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \n");
printf("    --excludeid 0 --i func.img \n");
printf("    --mask spmT.img --maskthresh 2.3 \n");
printf("    --sum bert.aseg-in-func.sum \n");
printf("    --avgwf bert.avgwf.dat --avgwfvol bert.avgwf.img\n");
printf("\n");
printf("This will compute the segmentation statistics from the automatic\n");
printf("FreeSurfer subcortical segmentation resampled into the functional\n");
printf("space (see below and mri_label2vol --help). It will report intensity\n");
printf("statistics from the 4D analyze volume func.img (same dimension as\n");
printf("aseg-in-func.img). The segmentation is masked by thresholding the\n");
printf("spmT.img map at 2.3. The average functional waveform of each\n");
printf("segmentation is reported in the ascii file bert.avgwf.dat and in the\n");
printf("4D analyze 'volume' bert.avgwf.img. This is not a real volume but just\n");
printf("another way to save the data that may be more convenient than ascii.\n");
printf("\n");
printf("4. mri_label2vol --seg $SUBJECTS_DIR/bert/mri/aseg \n");
printf("     --temp func.img --reg register.dat \n");
printf("     --fillthresh 0.5 --o aseg-in-func.img\n");
printf("\n");
printf("This uses mri_label2vol to resample the automatic subcortical\n");
printf("segmentation to the functional space. For more information\n");
printf("see mri_label2vol --help.\n");
printf("\n");
printf("5. mri_label2vol --annot $SUBJECTS_DIR/bert/label/lh.aparc.annot \n");
printf("     --temp func.img --reg register.dat --fillthresh 0.5 \n");
printf("     --hemi lh --subject bert --proj frac 0 .1 1 \n");
printf("     --o lh.aparc-in-func.img\n");
printf("\n");
printf("This uses mri_label2vol to resample the automatic cortical\n");
printf("segmentation to the functional space. For more information\n");
printf("see mri_label2vol --help.\n");
printf("\n");
printf("6. mri_segstats --annot bert lh aparc --i lh.thickness --sum lh.thickness.sum \n");
printf("\n");
printf("Produce a summary of the thickness in each parcellation of aparc. This \n");
printf("will give the same mean thicknesses as that created by mris_anatomical_stats\n");
printf("and found in stats/lh.aparc.stats.\n");
printf("\n");
printf("\n");
printf("SEE ALSO:\n");
printf("  mri_label2vol, tkregister2, mri_vol2roi.\n");
printf("\n");
#endif
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
  if (SegVolFile == NULL && annot == NULL && LabelFile == NULL && 
      DoETIVonly == 0 && DoOldETIVonly == 0) {
    printf("ERROR: must specify a segmentation volume\n");
    exit(1);
  }
  if (StatTableFile == NULL && FrameAvgFile == NULL && DoETIVonly == 0 && DoOldETIVonly == 0) {
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
  if (DoSupraTent && !DoSurfCtxVol) {
    printf("ERROR: need --surf-ctx-vol  with --supratent\n");
    exit(1);
  }
  if (ctabfile != NULL && gcafile != NULL) {
    printf("ERROR: cannot specify ctab and gca\n");
    exit(1);
  }
  if(DoSurfCtxVol && subject == NULL){
    printf("ERROR: need --subject with --surf-ctx-vol\n");
    exit(1);
  }
  if(DoTotalGrayVol && !DoSurfCtxVol){
    printf("ERROR: need --surf-ctx-vol with --totalgray\n");
    exit(1);
  }
  if(ctabfileOut && !ctabfile){
    printf("ERROR: need an input ctab to create output ctab\n");
    exit(1);
  }
  if(ctabfileOut){
    if(!strcmp(ctabfileOut,ctabfile)){
      printf("ERROR: output ctab is the same as input\n");
      exit(1);
    }
  }
  if (masksign == NULL) masksign = "abs";
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
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
    printf("%3d  %8d %10.1f  ",
           StatSumTable[n].id,
           StatSumTable[n].nhits,
           StatSumTable[n].vol);
    printf("%-30s ",StatSumTable[n].name);
    printf("%10.4f %10.4f %10.4f %10.4f %10.4f ",
           StatSumTable[n].min, StatSumTable[n].max,
           StatSumTable[n].range, StatSumTable[n].mean,
           StatSumTable[n].std);
    printf("\n");
  }
  return(0);
}
