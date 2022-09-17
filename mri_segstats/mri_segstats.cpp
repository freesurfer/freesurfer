/**
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

/*
   Subcort stuff that needs to be removed from the surface-based white
   matter volume:

   16  Brain-Stem                            119  159  176    0

    4  Left-Lateral-Ventricle                120   18  134    0
   10  Left-Thalamus                           0  118   14    0
   11  Left-Caudate                          122  186  220    0
   12  Left-Putamen                          236   13  176    0
   13  Left-Pallidum                          12   48  255    0
   17  Left-Hippocampus                      220  216   20    0
   18  Left-Amygdala                         103  255  255    0
   26  Left-Accumbens-area                   255  165    0    0
   28  Left-VentralDC                        165   42   42    0

   43  Right-Lateral-Ventricle               120   18  134    0
   49  Right-Thalamus                          0  118   14    0
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
#include <unistd.h>
#include <errno.h>

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
#include "fio.h"
#include "ctrpoints.h"
#include "stats.h"
#include "gtm.h"

#include "romp_support.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);


int MRIsegCount(MRI *seg, int id, int frame);
STATSUMENTRY *LoadStatSumFile(char *fname, int *nsegid);
int DumpStatSumTable(STATSUMENTRY *StatSumTable, int nsegid);
int CountEdits(char *subject, char *outfile);
float *WMAnatStats(const char *subject, const char *volname, int nErodes, float Pct);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *SUBJECTS_DIR = NULL, *FREESURFER_HOME=NULL;
char *SegVolFile = NULL;
char *InVolFile = NULL;
char *InVolRegFile = NULL;
MATRIX *InVolReg = NULL;
int InVolRegHeader = 0;
const char *InIntensityName = "";
const char *InIntensityUnits = "unknown";
char *MaskVolFile = NULL;
char *PVVolFile = NULL;
char *BrainMaskFile = NULL;
char *StatTableFile = NULL;
char *FrameAvgFile = NULL;
char *FrameAvgVolFile = NULL;
char *SpatFrameAvgFile = NULL;
int DoFrameAvg = 0;
int DoFrameSum = 0;
int RmFrameAvgMn = 0;
int DoAccumulate = 0;
int frame = 0;
int synth = 0;
int debug = 0;
int dontrun = 0;
long seed = 0;
MRI *seg, *invol, *famri, *maskvol, *pvvol, *brainvol, *mri_aseg, *mri_ribbon,*mritmp;
int nsegid0, *segidlist0;
int nsegid, *segidlist;
int NonEmptyOnly = 1;
int UserSegIdList[1000];
int nUserSegIdList = 0;
int nErodeSeg=0;
int DoExclSegId = 0, nExcl = 0, ExclSegIdList[1000], ExclSegId;
int DoExclCtxGMWM= 0;
int DoSurfCtxVol = 0;
int DoSurfWMVol = 0;
int DoSupraTent = 0;
double SupraTentVol, SupraTentVolCor;

char *gcafile = NULL;
GCA *gca;

float maskthresh = 0.5;
int   maskinvert = 0, maskframe = 0;
const char *masksign=NULL;
int   maskerode = 0;
int   nmaskhits;
int DoSubCortGrayVol = 0;
int DoTotalGrayVol = 0;
int BrainVolFromSeg = 0;
int   DoETIV = 0;
int   DoETIVonly = 0;
int   DoOldETIVonly = 0;
char *talxfmfile = NULL;
int SegFromInput=0;

char *ctabfile = NULL;
COLOR_TABLE *ctab = NULL;
STATSUMENTRY *StatSumTable = NULL;
STATSUMENTRY *StatSumTable2 = NULL;
char *ctabfileOut = NULL;

MRIS *mris;
char *subject = NULL;
char *hemi    = NULL;
char *annot   = NULL;
const char *whitesurfname = "white";

int Vox[3], DoVox = 0;
int  segbase = -1000;

int DoSquare = 0;
int DoSquareRoot = 0;
char *LabelFile = NULL;
double LabelThresh = 0;
int UseLabelThresh = 0;

int DoMultiply = 0;
double MultVal = 0;

int DoSNR = 0;
int UseRobust = 0;
float RobustPct = 5.0;
struct utsname uts;
char *cmdline, cwd[2000];

int DoEuler = 0;
int lheno, rheno;
int DoAbs = 0;
int UsePrintSegStat = 1; // use new way to print

int nReplace , SrcReplace[1000], TrgReplace[1000]; // for replacing segs
int GetCachedBrainVolStats = 1;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, n, nx, n0, skip, nhits, f, nsegidrep, ind, nthsegid;
  int c,r,s,err,DoContinue,nvox;
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

  nargs = handleVersionOption(argc, argv, "mri_segstats");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  cmdline = argv2cmdline(argc,argv);
  uname(&uts);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0)
  {
    usage_exit();
  }

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  if (subject != NULL)
  {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL)
    {
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if (DoETIV || DoETIVonly || DoOldETIVonly)
  {
    // calc total intracranial volume estimation
    // see this page:
    // http://surfer.nmr.mgh.harvard.edu/fswiki/eTIV
    // for info on determining the scaling factor.
    // a factor of 1948 was found to be best when using talairach.xfm
    // however when using talairach_with_skull.lta, 2150 was used
    double etiv_scale_factor = 1948.106;
    if (talxfmfile)
    {
      // path to talairach.xfm file spec'd on the command line
      int req = snprintf(tmpstr,1000,"%s",talxfmfile); 
      if( req >= 1000 ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

    }
    else
    {
      int req = snprintf(tmpstr,1000,
			 "%s/%s/mri/transforms/talairach.xfm",
			 SUBJECTS_DIR,
			 subject);
      if( req >= 1000 ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    if (DoOldETIVonly)
    {
      // back-door way to get the old way of calculating etiv, for debug
      int req = snprintf(tmpstr,1000,
			 "%s/%s/mri/transforms/talairach_with_skull.lta",
			 SUBJECTS_DIR,
			 subject);
      if( req >= 1000 ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      etiv_scale_factor = 2150;
    }
    double determinant = 0;
    atlas_icv = MRIestimateTIV(tmpstr,etiv_scale_factor,&determinant);
    printf("atlas_icv (eTIV) = %d mm^3    (det: %3f )\n",
           (int)atlas_icv,determinant);
    if (DoETIVonly || DoOldETIVonly) exit(0);
  }
  fflush(stdout);

  /* Make sure we can open the output summary table file*/
  if(StatTableFile)
  {
    fp = fopen(StatTableFile,"w");
    if (fp == NULL)
    {
      printf("ERROR: could not open %s for writing\n",StatTableFile);
      int err = errno;
      printf("Errno: %s\n", strerror(err) );
      exit(1);
    }
    fclose(fp);
    unlink(StatTableFile); // delete
  }

  /* Make sure we can open the output frame average file*/
  if(FrameAvgFile != NULL) {
    fp = fopen(FrameAvgFile,"w");
    if (fp == NULL){
      printf("ERROR: could not open %s for writing\n",FrameAvgFile);
      exit(1);
    }
    fclose(fp);
    unlink(FrameAvgFile); // delete
  }

  if(DoEuler){
    int req = snprintf(tmpstr,1000,"%s/%s/surf/lh.orig.nofix",SUBJECTS_DIR,subject); 
    if( req >= 1000 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if(!fio_FileExistsReadable(tmpstr)){
      printf("Warning: cannot find %s, not computing euler number\n",tmpstr);
      DoEuler = 0;
    }
    req = snprintf(tmpstr,1000,"%s/%s/surf/rh.orig.nofix",SUBJECTS_DIR,subject); 
    if( req >= 1000 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if(!fio_FileExistsReadable(tmpstr)){
      printf("Warning: cannot find %s, not computing euler number\n",tmpstr);
      DoEuler = 0;
    }
  }
  if(DoEuler){
    int nvertices, nfaces, nedges;
    int req = snprintf(tmpstr,1000,"%s/%s/surf/lh.orig.nofix",SUBJECTS_DIR,subject);  
    if( req >= 1000 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("Computing euler number\n");
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lheno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    MRISfree(&mris);
    sprintf(tmpstr,"%s/%s/surf/rh.orig.nofix",SUBJECTS_DIR,subject);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rheno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    MRISfree(&mris);
    printf("orig.nofix lheno = %4d, rheno = %d\n",lheno,rheno);
    printf("orig.nofix lhholes = %4d, rhholes = %d\n",1-lheno/2,1-rheno/2);
  }

  /* Load the segmentation */
  if (SegVolFile){
    printf("Loading %s\n",SegVolFile);
    seg = MRIread(SegVolFile);
    if (seg == NULL){
      printf("ERROR: loading %s\n",SegVolFile);
      exit(1);
    }
    if(nReplace > 0) {
      printf("Replacing %d\n",nReplace);
      mritmp = MRIreplaceList(seg,SrcReplace, TrgReplace, nReplace, NULL, NULL);
      MRIfree(&seg);
      seg = mritmp;
    }

    if(nErodeSeg){
      printf("Eroding seg %d times\n",nErodeSeg);
      tmp = MRIerodeSegmentation(seg, NULL, nErodeSeg, 0);
      MRIfree(&seg);
      seg = tmp;
    }
    if(DoVox){
      printf("Replacing seg with a single voxel at %d %d %d\n",
             Vox[0],Vox[1],Vox[2]);
      for(c=0; c < seg->width; c++)
      {
        for(r=0; r < seg->height; r++)
        {
          for(s=0; s < seg->depth; s++)
          {
            if(c == Vox[0] && r == Vox[1] && s == Vox[2])
            {
              MRIsetVoxVal(seg,c,r,s,0, 1);
            }
            else
            {
              MRIsetVoxVal(seg,c,r,s,0, 0);
            }
          }
        }
      }
    }
  }
  else if(annot){
    printf("Constructing seg from annotation\n");
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,whitesurfname);
    mris = MRISread(tmpstr);
    if (mris==NULL) exit(1);
    if(fio_FileExistsReadable(annot))sprintf(tmpstr,"%s",annot);
    else sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot);
    printf("\nReading annotation %s\n",tmpstr);
    err = MRISreadAnnotation(mris, tmpstr);
    if(err) {
      printf(" ... trying local annot\n");
      err = MRISreadAnnotation(mris, annot); // assume annot is full path
      if(! err) printf(" Successfully read local annot\n\n");
    }
    if (err) exit(1);

    if(segbase == -1000){
      // segbase has not been set with --segbase
      if(!strcmp(annot,"aparc"))
      {
        if(!strcmp(hemi,"lh"))
        {
          segbase = 1000;
        }
        else
        {
          segbase = 2000;
        }
      }
      else if(!strcmp(annot,"aparc.a2005s"))
      {
        if(!strcmp(hemi,"lh"))
        {
          segbase = 1100;
        }
        else
        {
          segbase = 2100;
        }
      }
      else
      {
        segbase = 0;
      }
    }
    printf("Seg base %d\n",segbase);
    seg = MRISannot2seg(mris,segbase);
    // Now create a colortable in a temp location to be read out below (hokey)
    mris->ct->idbase = segbase;
    if (mris->ct)
    {
      std::string tmpfile = makeTempFile(".ctab");
      sprintf(tmpstr, "%s", tmpfile.c_str());
      ctabfile = strcpyalloc(tmpstr);
      CTABwriteFileASCII(mris->ct,ctabfile);
    }
  }
  else if(LabelFile){
    printf("Constructing seg from label\n");
    if(UseLabelThresh) printf(" Label Threshold = %g\n",LabelThresh);
    label = LabelRead(NULL, LabelFile);
    if(label == NULL) exit(1);
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,whitesurfname);
    mris = MRISread(tmpstr);
    if (mris==NULL) exit(1);
    seg = MRIalloc(mris->nvertices,1,1,MRI_INT);
    for (n = 0; n < label->n_points; n++){
      if(UseLabelThresh && label->lv[n].stat < LabelThresh) continue;
      MRIsetVoxVal(seg,label->lv[n].vno,0,0,0, 1);
    }
  }
  else {
    printf("Creating a segmentation of all 1s from %s\n",InVolFile);
    mritmp = MRIreadHeader(InVolFile,MRI_VOLUME_TYPE_UNKNOWN);
    seg = MRIconst(mritmp->width, mritmp->height, mritmp->depth, mritmp->nframes,1,NULL);
    MRIfree(&mritmp);
  }

  if (ctabfile != NULL)
  {
    /* Load the color table file */
    ctab = CTABreadASCII(ctabfile);
    if (ctab == NULL)
    {
      printf("ERROR: reading %s\n",ctabfile);
      exit(1);
    }
  }
  else {
    if(seg->ct){
      ctab = seg->ct;
      printf("Using embedded color table (and excluding seg 0)\n");
      ExclSegIdList[nExcl] = 0;
      nExcl ++;
      DoExclSegId = 1;
    }
  }

  if (gcafile != NULL)
  {
    gca = GCAread(gcafile);
    if (gca == NULL)
    {
      printf("ERROR: reading %s\n",gcafile);
      exit(1);
    }
    ctab = GCAcolorTableCMA(gca);
  }

  std::vector<double> BrainVolStats;
  if(DoSurfWMVol || DoSurfCtxVol || DoSupraTent || BrainVolFromSeg || DoSubCortGrayVol){
    sprintf(tmpstr,"%s/%s/mri/ribbon.mgz",SUBJECTS_DIR,subject);
    if(fio_FileExistsReadable(tmpstr)){
      if(GetCachedBrainVolStats){
	printf("Getting Brain Volume Statistics\n");
	BrainVolStats = ReadCachedBrainVolumeStats(subject, SUBJECTS_DIR);
      }
      else {
	printf("Computing Brain Volume Statistics\n");
	BrainVolStats = ComputeBrainVolumeStats(subject, SUBJECTS_DIR);
      }
    }
    else{
      printf("Warning: cannot find %s, not computing whole brain stats\n",tmpstr);
      DoSurfWMVol=DoSurfCtxVol=DoSupraTent=BrainVolFromSeg=DoSubCortGrayVol=0;
      DoTotalGrayVol=0;DoSupraTent=0;
    }
  }

  /* Load the input volume */
  if (InVolFile != NULL)
  {
    printf("Loading %s\n",InVolFile);
    fflush(stdout);
    invol = MRIread(InVolFile);
    if(invol == NULL)
    {
      printf("ERROR: loading %s\n",InVolFile);
      exit(1);
    }
    if(frame >= invol->nframes)
    {
      printf("ERROR: input frame = %d, input volume only has %d frames\n",
             frame,invol->nframes);
      exit(1);
    }
    if(InVolReg || InVolRegHeader)
    {
      if(InVolReg)
      {
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
      if(err)
      {
        exit(1);
      }
      MRIfree(&invol);
      invol = tmp;
    }
    if(MRIdimMismatch(invol,seg,0))
    {
      printf("ERROR: dimension mismatch between input volume and seg\n");
      printf("  input %d %d %d\n",invol->width,invol->height,invol->depth);
      printf("  seg   %d %d %d\n",seg->width,seg->height,seg->depth);
      exit(1);
    }
    if(DoMultiply)
    {
      printf("Multiplying input by %lf\n",MultVal);
      MRImultiplyConst(invol,MultVal,invol);
    }
    if(DoAbs)
    {
      printf("Computing absolute value of input\n");
      MRIabs(invol,invol);
    }
    if(DoSquare)
    {
      printf("Computing square of input\n");
      MRIsquare(invol, NULL, invol);
    }
    if(DoSquareRoot)
    {
      printf("Computing square root of input\n");
      MRIsquareRoot(invol, NULL, invol);
    }
  }

  /* Load the partial volume mri */
  if (PVVolFile != NULL)
  {
    printf("Loading %s\n",PVVolFile);
    fflush(stdout);
    pvvol = MRIread(PVVolFile);
    if (pvvol == NULL)
    {
      printf("ERROR: loading %s\n",PVVolFile);
      exit(1);
    }
    if(MRIdimMismatch(pvvol,seg,0))
    {
      printf("ERROR: dimension mismatch between PV volume and seg\n");
      printf("  pvvol %d %d %d\n",pvvol->width,pvvol->height,pvvol->depth);
      printf("  seg   %d %d %d\n",seg->width,seg->height,seg->depth);
      exit(1);
    }
  }

  /* Load the mask volume */
  if (MaskVolFile != NULL)
  {
    printf("Loading %s\n",MaskVolFile);
    fflush(stdout);
    maskvol = MRIread(MaskVolFile);
    if (maskvol == NULL)
    {
      printf("ERROR: loading %s\n",MaskVolFile);
      exit(1);
    }
    if (maskframe >= maskvol->nframes)
    {
      printf("ERROR: mask frame = %d, mask volume only has %d frames\n",
             maskframe,maskvol->nframes);
      exit(1);
    }
    if(MRIdimMismatch(maskvol,seg,0))
    {
      printf("ERROR: dimension mismatch between brain volume and seg\n");
      exit(1);
    }
    mri_binarize(maskvol, maskthresh, masksign, maskinvert,
                 maskvol, &nmaskhits);
    if (nmaskhits == 0)
    {
      printf("WARNING: no voxels in mask meet thresholding criteria.\n");
      printf("The output table will be empty.\n");
      printf("thresh = %g, sign = %s, inv = %d\n",maskthresh,masksign,maskinvert);
      //exit(1);
    }
    printf("There were %d voxels in the orginal mask\n",nmaskhits);
    if(maskerode > 0)
    {
      printf("Eroding %d voxels in 3d\n",maskerode);
      for(n=0; n<maskerode; n++)
      {
        MRIerode(maskvol,maskvol);
      }
    }

    /* perform the masking */
    for (c=0; c < seg->width; c++)
    {
      for (r=0; r < seg->height; r++)
      {
        for (s=0; s < seg->depth; s++)
        {
          // Set voxels out of the mask to 0
          if (! (int)MRIgetVoxVal(maskvol,c,r,s,maskframe))
          {
            MRIsetVoxVal(seg,c,r,s,0,0);
          }
        }
      }
    }
  }

  if (!mris)
  {
    voxelvolume = seg->xsize * seg->ysize * seg->zsize;
    printf("Voxel Volume is %g mm^3\n",voxelvolume);
  }
  else
  {
    if (mris->group_avg_surface_area > 0)
    {
      voxelvolume = mris->group_avg_surface_area/mris->nvertices;
    }
    else
    {
      voxelvolume = mris->total_area/mris->nvertices;
    }
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

  if (ctab == NULL && nUserSegIdList == 0)
  {
    /* Must get list of segmentation ids from segmentation itself*/
    segidlist = segidlist0;
    nsegid = nsegid0;
    StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
    for (n=0; n < nsegid; n++)
    {
      StatSumTable[n].id = segidlist[n];
      strcpy(StatSumTable[n].name, "\0");
    }
  }
  else   /* Get from user or color table */
  {
    if (ctab != NULL)
    {
      if (nUserSegIdList == 0)
      {
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
        for (n=0; n < ntotalsegid; n++)
        {
          CTABisEntryValid(ctab,n,&valid);
          if(!valid)
          {
            continue;
          }
          StatSumTable[usersegid].id = n;
          CTABcopyName(ctab,n,StatSumTable[usersegid].name,
                       sizeof(StatSumTable[usersegid].name));
          CTABrgbAtIndexi(ctab, n, &StatSumTable[usersegid].red,
                          &StatSumTable[usersegid].green,
                          &StatSumTable[usersegid].blue);
          usersegid++;
        }
      }
      else
      {
        /* User has specified --id, use those and get names from ctab */
        nsegid = nUserSegIdList;
        StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
        for (n=0; n < nsegid; n++)
        {
          StatSumTable[n].id = UserSegIdList[n];
          /* Here ind should be the same as the ctab entry, but make
             sure it's valid. */
          ind = StatSumTable[n].id;
          CTABisEntryValid(ctab,ind,&valid);
          if (!valid)
          {
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
    else   /* User specified ids, but no color table */
    {
      nsegid = nUserSegIdList;
      StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
      for (n=0; n < nsegid; n++)
      {
        StatSumTable[n].id = UserSegIdList[n];
      }
    }
  }

  printf("Found %3d segmentations\n",nsegid);
  if(nsegid == 0){
    printf("ERROR: no segmentations to report\n");
    exit(1);
  }

  printf("Computing statistics for each segmentation\n");
  fflush(stdout);

  DoContinue=0;nx=0;skip=0;n0=0;vol=0;nhits=0;c=0;min=0.0;max=0.0;range=0.0;mean=0.0;std=0.0;snr=0.0;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) firstprivate(DoContinue,nx,skip,n0,vol,nhits,c,min,max,range,mean,std,snr)  schedule(guided)
#endif
  for (n=0; n < nsegid; n++)
  {
    ROMP_PFLB_begin
    if(DoExclSegId)
    {
      DoContinue = 0;
      for(nx=0; nx < nExcl; nx++)
      {
        if(StatSumTable[n].id == ExclSegIdList[nx])
        {
          DoContinue = 1;
          break;
        }
      }
      if(DoContinue)
      {
        ROMP_PFLB_continue;
      }
    }

    // Skip ones that are not represented
    skip = 1;
    for (n0=0; n0 < nsegid0; n0++)
      if (StatSumTable[n].id == segidlist0[n0])
      {
        skip = 0;
      }
    if (skip)
    {
      ROMP_PFLB_continue;
    }

    if (!dontrun)
    {
      if (!mris)
      {
        if (pvvol == NULL)
        {
          nhits = MRIsegCount(seg, StatSumTable[n].id, 0);
          vol = nhits*voxelvolume;
        }
        else
        {
          vol = MRIvoxelsInLabelWithPartialVolumeEffects(seg, pvvol, StatSumTable[n].id, NULL, NULL);
          nhits = MRIsegCount(seg, StatSumTable[n].id, 0);
//          nhits = nint(vol/voxelvolume);
        }
      }
      else
      {
        // Compute area here
        nhits = 0;
        vol = 0;
        for (c=0; c < mris->nvertices; c++)
        {
          if (MRIgetVoxVal(seg,c,0,0,0)==StatSumTable[n].id)
          {
            nhits++;
            if (mris->group_avg_vtxarea_loaded)
            {
              vol += mris->vertices[c].group_avg_area;
            }
            else
            {
              vol += mris->vertices[c].area;
            }
          }
        }
      }
    }
    else
    {
      nhits = n;
    }

    StatSumTable[n].nhits = nhits;
    StatSumTable[n].vol = vol;
    if (InVolFile != NULL && !dontrun)
    {
      if (nhits > 0)
      {
        if(UseRobust == 0)
          MRIsegStats(seg, StatSumTable[n].id, invol, frame,
            &min, &max, &range, &mean, &std);
        else
          MRIsegStatsRobust(seg, StatSumTable[n].id, invol, frame,
            &min, &max, &range, &mean, &std, RobustPct);

        snr = mean/std;
      }
      else
      {
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
      if(DoAccumulate == 0) StatSumTable[n].mean  = mean;
      if(DoAccumulate == 1) StatSumTable[n].mean  = mean*nhits;
      StatSumTable[n].std   = std;
      StatSumTable[n].snr   = snr;
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  /* print results ordered */
  for (n=0; n < nsegid; n++)
  {
    if(DoExclSegId)
    {
      DoContinue = 0;
      for(nx=0; nx < nExcl; nx++)
      {
        if(StatSumTable[n].id == ExclSegIdList[nx])
        {
          DoContinue = 1;
          break;
        }
      }
      if(DoContinue)
      {
        continue;
      }
    }
    if(Gdiag_no > 1){
      printf("%3d   %3d  %33s  %6d  %10.3f\n",n,StatSumTable[n].id,StatSumTable[n].name,
	     StatSumTable[n].nhits,StatSumTable[n].vol);
      fflush(stdout);
    }
  }
  printf("\n");


  /* Remove empty segmentations, if desired */
  if (NonEmptyOnly || DoExclSegId)
  {
    // Count the number of nonempty segmentations
    nsegidrep = 0;
    for (n=0; n < nsegid; n++)
    {
      if(NonEmptyOnly && StatSumTable[n].nhits==0)
      {
        continue;
      }
      if(DoExclSegId)
      {
        DoContinue = 0;
        for(nx=0; nx < nExcl; nx++)
        {
          if(StatSumTable[n].id == ExclSegIdList[nx])
          {
            DoContinue = 1;
            break;
          }
        }
        if(DoContinue)
        {
          continue;
        }
      }
      nsegidrep ++;
    }
    StatSumTable2 = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegidrep);
    nthsegid = 0;
    for (n=0; n < nsegid; n++)
    {
      if(NonEmptyOnly && StatSumTable[n].nhits==0)
      {
        continue;
      }

      if(DoExclSegId)
      {
        DoContinue = 0;
        for(nx=0; nx < nExcl; nx++)
        {
          if(StatSumTable[n].id == ExclSegIdList[nx])
          {
            DoContinue = 1;
            break;
          }
        }
        if(DoContinue)
        {
          continue;
        }
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

  /* Dump the table to the screen */
  if (debug)
  {
    for (n=0; n < nsegid; n++)
    {
      printf("%3d  %8d %10.1f  ", StatSumTable[n].id,StatSumTable[n].nhits,
             StatSumTable[n].vol);
      if (ctab != NULL)
      {
        printf("%-30s ",StatSumTable[n].name);
      }
      if (InVolFile != NULL)
      {
        printf("%10.4f %10.4f %10.4f %10.4f %10.4f ",
               StatSumTable[n].min, StatSumTable[n].max,
               StatSumTable[n].range, StatSumTable[n].mean,
               StatSumTable[n].std);
        if(DoSNR)
        {
          printf("%10.4f ",StatSumTable[n].snr);
        }
      }
      printf("\n");
    }
  }

  /* Print the table to the output file */
  if (StatTableFile != NULL)
  {
    fp = fopen(StatTableFile,"w");
    fprintf(fp,"# Title Segmentation Statistics \n");
    fprintf(fp,"# \n");
    fprintf(fp,"# generating_program %s\n",Progname);
    fprintf(fp,"# cvs_version %s\n",getVersion().c_str());
    fprintf(fp,"# cmdline %s\n",cmdline);
    fprintf(fp,"# sysname  %s\n",uts.sysname);
    fprintf(fp,"# hostname %s\n",uts.nodename);
    fprintf(fp,"# machine  %s\n",uts.machine);
    fprintf(fp,"# user     %s\n",VERuser());
    if (mris) fprintf(fp,"# anatomy_type surface\n");
    else      fprintf(fp,"# anatomy_type volume\n");
    fprintf(fp,"# \n");
    if (subject != NULL)
    {
      fprintf(fp,"# SUBJECTS_DIR %s\n",SUBJECTS_DIR);
      fprintf(fp,"# subjectname %s\n",subject);
    }
    if (UseRobust) fprintf(fp,"# RobustPct %g\n",RobustPct);
    if (!BrainVolStats.empty()) {
      if(fabs(voxelvolume-1)>.01){
	// This indicates that the global stats has been fixed
	fprintf(fp,"# BrainVolStatsFixed see surfer.nmr.mgh.harvard.edu/fswiki/BrainVolStatsFixed\n");
      }
      else{
	fprintf(fp,"# BrainVolStatsFixed-NotNeeded because voxelvolume=1mm3\n");
      }
    }
    if (BrainVolFromSeg)
    {
      fprintf(fp,"# Measure BrainSeg, BrainSegVol, "
              "Brain Segmentation Volume, %f, mm^3\n",
              BrainVolStats[0]);
      fprintf(fp,"# Measure BrainSegNotVent, BrainSegVolNotVent, "
              "Brain Segmentation Volume Without Ventricles, %f, mm^3\n",
              BrainVolStats[1]);
      // Not computed in ComputeBrainVolumeStats2() anymore because it is very close
      // to the voxel-based version and adds needless complexity
      //fprintf(fp,"# Measure BrainSegNotVentSurf, BrainSegVolNotVentSurf, "
      //      "Brain Segmentation Volume Without Ventricles from Surf, %f, mm^3\n",
      //      BrainVolStats[14]);
    }
    if (!BrainVolStats.empty()) {
      fprintf(fp,"# Measure VentricleChoroidVol, VentricleChoroidVol, "
	      "Volume of ventricles and choroid plexus, %f, mm^3\n",
	      BrainVolStats[15]);
    }
    if(DoSurfCtxVol)
    {
      // Does this include the non-cortical areas of the surface?
      fprintf(fp,"# Measure lhCortex, lhCortexVol, "
              "Left hemisphere cortical gray matter volume, %f, mm^3\n",BrainVolStats[5]);
      fprintf(fp,"# Measure rhCortex, rhCortexVol, "
              "Right hemisphere cortical gray matter volume, %f, mm^3\n",BrainVolStats[6]);
      fprintf(fp,"# Measure Cortex, CortexVol, "
              "Total cortical gray matter volume, %f, mm^3\n",BrainVolStats[7]);
    }
    if(DoSurfWMVol)
    {
      fprintf(fp,"# Measure lhCerebralWhiteMatter, lhCerebralWhiteMatterVol, "
              "Left hemisphere cerebral white matter volume, %f, mm^3\n",BrainVolStats[9]);
      fprintf(fp,"# Measure rhCerebralWhiteMatter, rhCerebralWhiteMatterVol, "
              "Right hemisphere cerebral white matter volume, %f, mm^3\n",BrainVolStats[10]);
      fprintf(fp,"# Measure CerebralWhiteMatter, CerebralWhiteMatterVol, "
              "Total cerebral white matter volume, %f, mm^3\n",BrainVolStats[11]);
    }
    if(DoSubCortGrayVol)
    {
      fprintf(fp,"# Measure SubCortGray, SubCortGrayVol, "
              "Subcortical gray matter volume, %f, mm^3\n",
              BrainVolStats[4]);
    }
    if(DoTotalGrayVol)
    {
      fprintf(fp,"# Measure TotalGray, TotalGrayVol, Total gray matter volume, %f, mm^3\n",
              BrainVolStats[8]);
    }
    if(DoSupraTent)
    {
      fprintf(fp,"# Measure SupraTentorial, SupraTentorialVol, "
              "Supratentorial volume, %f, mm^3\n",BrainVolStats[2]);
      fprintf(fp,"# Measure SupraTentorialNotVent, SupraTentorialVolNotVent, "
              "Supratentorial volume, %f, mm^3\n",BrainVolStats[3]);
      //fprintf(fp,"# Measure SupraTentorialNotVentVox, SupraTentorialVolNotVentVox, "
      //      "Supratentorial volume voxel count, %f, mm^3\n",BrainVolStats[13]);
    }
    if (BrainMaskFile && (!BrainVolStats.empty())) {
      //fprintf(fp,"# BrainMaskFile  %s \n",BrainMaskFile);
      //fprintf(fp,"# BrainMaskFileTimeStamp  %s \n",
      //       VERfileTimeStamp(BrainMaskFile));
      //fprintf(fp,"# Measure BrainMask, BrainMaskNVox, "
      //        "Number of Brain Mask Voxels, %7d, unitless\n",
      //        nbrainmaskvoxels);
      fprintf(fp,"# Measure Mask, MaskVol, "
              "Mask Volume, %f, mm^3\n",BrainVolStats[12]);
    }
    if(DoETIV && BrainVolFromSeg){
      fprintf(fp,"# Measure BrainSegVol-to-eTIV, BrainSegVol-to-eTIV, "
              "Ratio of BrainSegVol to eTIV, %f, unitless\n",
              BrainVolStats[0]/atlas_icv);
      fprintf(fp,"# Measure MaskVol-to-eTIV, MaskVol-to-eTIV, "
              "Ratio of MaskVol to eTIV, %f, unitless\n",
              BrainVolStats[12]/atlas_icv);
    }
    if(DoEuler){
      fprintf(fp,"# Measure lhSurfaceHoles, lhSurfaceHoles, "
              "Number of defect holes in lh surfaces prior to fixing, %d, unitless\n",
              (1-lheno/2));
      fprintf(fp,"# Measure rhSurfaceHoles, rhSurfaceHoles, "
              "Number of defect holes in rh surfaces prior to fixing, %d, unitless\n",
              (1-rheno/2));
      fprintf(fp,"# Measure SurfaceHoles, SurfaceHoles, "
              "Total number of defect holes in surfaces prior to fixing, %d, unitless\n",
              (1-lheno/2) + (1-rheno/2) );
    }
    if (DoETIV)
    {
      //fprintf(fp,"# Measure IntraCranialVol, ICV, "
      //      "Intracranial Volume, %f, mm^3\n",atlas_icv);
      fprintf(fp,"# Measure EstimatedTotalIntraCranialVol, eTIV, "
              "Estimated Total Intracranial Volume, %f, mm^3\n",atlas_icv);
    }
    if (SegVolFile)
    {
      fprintf(fp,"# SegVolFile %s \n",SegVolFile);
      fprintf(fp,"# SegVolFileTimeStamp  %s \n",VERfileTimeStamp(SegVolFile));
    }
    if (annot)
    {
      fprintf(fp,"# Annot %s %s %s\n",subject,hemi,annot);
    }
    if (LabelFile)
    {
      fprintf(fp,"# Label %s %s %s\n",subject,hemi,LabelFile);
    }
    if (ctabfile)
    {
      fprintf(fp,"# ColorTable %s \n",ctabfile);
      fprintf(fp,"# ColorTableTimeStamp %s \n",VERfileTimeStamp(ctabfile));
    }
    if (gcafile)
    {
      fprintf(fp,"# ColorTableFromGCA %s \n",gcafile);
      fprintf(fp,"# GCATimeStamp %s \n",VERfileTimeStamp(gcafile));
    }
    if (MaskVolFile)
    {
      fprintf(fp,"# MaskVolFile  %s \n",MaskVolFile);
      fprintf(fp,"#   MaskVolFileTimeStamp  %s \n",
              VERfileTimeStamp(MaskVolFile));
      fprintf(fp,"#   MaskThresh %f \n",maskthresh);
      fprintf(fp,"#   MaskSign   %s \n",masksign);
      fprintf(fp,"#   MaskFrame  %d \n",maskframe);
      fprintf(fp,"#   MaskInvert %d \n",maskframe);
    }
    if (InVolFile)
    {
      fprintf(fp,"# InVolFile  %s \n",InVolFile);
      fprintf(fp,"# InVolFileTimeStamp  %s \n",VERfileTimeStamp(InVolFile));
      fprintf(fp,"# InVolFrame %d \n",frame);
    }
    if (PVVolFile)
    {
      fprintf(fp,"# PVVolFile  %s \n",PVVolFile);
      fprintf(fp,"# PVVolFileTimeStamp  %s \n",VERfileTimeStamp(PVVolFile));
    }
    if(DoExclCtxGMWM)
    {
      fprintf(fp,"# Excluding Cortical Gray and White Matter\n");
    }
    if(DoExclSegId)
    {
      fprintf(fp,"# ExcludeSegId ");
      for(nx=0; nx < nExcl; nx++)
      {
        fprintf(fp,"%d ",ExclSegIdList[nx]);
      }
      fprintf(fp,"\n");
    }
    if(NonEmptyOnly)
    {
      fprintf(fp,"# Only reporting non-empty segmentations\n");
    }
    if(!mris)
    {
      fprintf(fp,"# VoxelVolume_mm3 %g \n",voxelvolume);
    }
    else
    {
      fprintf(fp,"# VertexArea_mm2 %g \n",voxelvolume);
    }


    if(UsePrintSegStat){
      printf("Using PrintSegStat\n");
      SEGSTAT *segstat;
      segstat = (SEGSTAT *)calloc(sizeof(SEGSTAT),1);
      segstat->nentries = nsegid;
      segstat->entry = StatSumTable;
      if(!mris) segstat->IsSurf = 0;
      else      segstat->IsSurf = 1;
      if(ctab)  segstat->UseName = 1;
      else      segstat->UseName = 0;
      if(InVolFile) segstat->DoIntensity = 1;
      else          segstat->DoIntensity = 0;
      segstat->InIntensityName = InIntensityName;
      segstat->InIntensityUnits = InIntensityUnits;
      segstat->DoSNR = DoSNR;
      PrintSegStat(fp, segstat);
    }
    else {
      printf("Not using PrintSegStat\n");

    fprintf(fp,"# TableCol  1 ColHeader Index \n");
    fprintf(fp,"# TableCol  1 FieldName Index \n");
    fprintf(fp,"# TableCol  1 Units     NA \n");

    fprintf(fp,"# TableCol  2 ColHeader SegId \n");
    fprintf(fp,"# TableCol  2 FieldName Segmentation Id\n");
    fprintf(fp,"# TableCol  2 Units     NA\n");
    if (!mris)
    {
      fprintf(fp,"# TableCol  3 ColHeader NVoxels \n");
      fprintf(fp,"# TableCol  3 FieldName Number of Voxels\n");
      fprintf(fp,"# TableCol  3 Units     unitless\n");
      fprintf(fp,"# TableCol  4 ColHeader Volume_mm3\n");
      fprintf(fp,"# TableCol  4 FieldName Volume\n");
      fprintf(fp,"# TableCol  4 Units     mm^3\n");
    }
    else
    {
      fprintf(fp,"# TableCol  3 ColHeader NVertices \n");
      fprintf(fp,"# TableCol  3 FieldName Number of Vertices\n");
      fprintf(fp,"# TableCol  3 Units     unitless\n");
      fprintf(fp,"# TableCol  4 ColHeader Area_mm2\n");
      fprintf(fp,"# TableCol  4 FieldName Area\n");
      fprintf(fp,"# TableCol  4 Units     mm^2\n");
    }
    n = 5;
    if (ctab)
    {
      fprintf(fp,"# TableCol %2d ColHeader StructName\n",n);
      fprintf(fp,"# TableCol %2d FieldName Structure Name\n",n);
      fprintf(fp,"# TableCol %2d Units     NA\n",n);
      n++;
    }

    if (InVolFile)
    {
      fprintf(fp,"# TableCol %2d ColHeader %sMean \n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Intensity %sMean\n",
              n,InIntensityName);
      fprintf(fp,"# TableCol %2d Units     %s\n",n,InIntensityUnits);
      n++;

      fprintf(fp,"# TableCol %2d ColHeader %sStdDev\n",n,InIntensityName);
      fprintf(fp,"# TableCol %2d FieldName Intensity %sStdDev\n",
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
    if (!mris)
    {
      fprintf(fp,"NVoxels Volume_mm3 ");
    }
    else
    {
      fprintf(fp,"NVertices Area_mm2 ");
    }
    fprintf(fp,"StructName ");
    if(InVolFile)
    {
      fprintf(fp,"%sMean %sStdDev %sMin %sMax %sRange  ",
              InIntensityName, InIntensityName,
              InIntensityName,InIntensityName,
              InIntensityName);
      if(DoSNR)
      {
        fprintf(fp,"%sSNR ",InIntensityName);
      }
    }
    fprintf(fp,"\n");

    for (n=0; n < nsegid; n++)
    {
      fprintf(fp,"%3d %3d  %8d %10.1f  ", n+1, StatSumTable[n].id,
              StatSumTable[n].nhits, StatSumTable[n].vol);
      if(ctab != NULL)
      {
        fprintf(fp,"%-30s ",StatSumTable[n].name);
      }
      else
      {
        fprintf(fp,"Seg%04d ",StatSumTable[n].id);
      }
      if (InVolFile != NULL)
      {
        fprintf(fp,"%10.4f %10.4f %10.4f %10.4f %10.4f ",
                StatSumTable[n].mean, StatSumTable[n].std,
                StatSumTable[n].min, StatSumTable[n].max,
                StatSumTable[n].range);
        if(DoSNR)
        {
          fprintf(fp,"%10.4f ",StatSumTable[n].snr);
        }
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
}

  if(ctabfileOut != NULL)
  {
    fp = fopen(ctabfileOut,"w");
    for (n=0; n < nsegid; n++)
      fprintf(fp,"%d %-30s %3d %3d %3d 0\n",StatSumTable[n].id,
              StatSumTable[n].name,StatSumTable[n].red,
              StatSumTable[n].green,StatSumTable[n].blue);
    fclose(fp);
  }

  // Average input across space to create a waveform
  // for each segmentation
  if (DoFrameAvg){
    printf("Computing spatial average of each frame\n");
    favg = (double **) calloc(sizeof(double *),nsegid);
    for (n=0; n < nsegid; n++)
      favg[n] = (double *) calloc(sizeof(double),invol->nframes);
    favgmn = (double *) calloc(sizeof(double *),nsegid);
    for (n=0; n < nsegid; n++) {
      if(debug){
	printf("%3d",n);
	if (n%20 == 19) printf("\n");
	fflush(stdout);
      }
      nvox = MRIsegFrameAvg(seg, StatSumTable[n].id, invol, favg[n]);
      favgmn[n] = 0.0;
      for(f=0; f < invol->nframes; f++) {
	if(DoFrameSum) favg[n][f] *= nvox; // Undo spatial average
	favgmn[n] += favg[n][f];
      }
      favgmn[n] /= invol->nframes;
      if(RmFrameAvgMn) for(f=0; f < invol->nframes; f++) favg[n][f] -= favgmn[n];
    }
    printf("\n");

    // Save mean over space and frames in simple text file
    // Each seg on a separate line
    if(SpatFrameAvgFile) {
      printf("Writing to %s\n",SpatFrameAvgFile);
      fp = fopen(SpatFrameAvgFile,"w");
      for (n=0; n < nsegid; n++){
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
      for (f=0; f < invol->nframes; f++){
        //fprintf(fp,"%3d %7.3f ",f,f*invol->tr/1000);
        for (n=0; n < nsegid; n++) fprintf(fp,"%11.5f ",favg[n][f]);
        fprintf(fp,"\n");
      }
      fclose(fp);
    }

    // Save as an MRI "volume"
    if(FrameAvgVolFile) {
      if(nsegid == 0){
	printf("ERROR: no voxels found in segmentation\n");
	exit(1);
      }
      printf("Writing to %s\n",FrameAvgVolFile);
      famri = MRIallocSequence(nsegid,1,1,MRI_FLOAT,invol->nframes);
      MRIcopyHeader(invol,famri);
      for (f=0; f < invol->nframes; f++){
        for (n=0; n < nsegid; n++)
          MRIsetVoxVal(famri,n,0,0,f,(float)favg[n][f]);
      }
      MRIwrite(famri,FrameAvgVolFile);
    }
  }// Done with Frame Average

  printf("mri_segstats done\n");
  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused, nth, err;
  char **pargv, *option ;

  if (argc < 1)
  {
    usage_exit();
  }

  setRandomSeed(4321) ;
    // It was previously using a different random sequence every time
    // it was run, and did not have a --seed option to stop this
    // and it appears in many scripts!

  nargc   = argc;
  pargv = argv;
  while (nargc > 0)
  {

    option = pargv[0];
    if (debug)
    {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help")||!strcasecmp(option, "--usage")
        ||!strcasecmp(option, "-h")||!strcasecmp(option, "-u"))
    {
      print_help() ;
    }
    else if (!strcasecmp(option, "--version"))
    {
      print_version() ;
    }
    else if (!strcasecmp(option, "--debug"))
    {
      debug = 1;
    }
    else if (!strcasecmp(option, "--newprint")) UsePrintSegStat = 1;
    else if (!strcasecmp(option, "--no-newprint")) UsePrintSegStat = 0;
    else if (!strcasecmp(option, "--dontrun"))
    {
      dontrun = 1;
    }
    else if (!strcasecmp(option, "--nonempty")||!strcasecmp(option, "--non-empty"))
      NonEmptyOnly = 1;
    else if (!strcasecmp(option, "--empty"))
    {
      NonEmptyOnly = 0;
    }
    else if (!strcasecmp(option, "--no-cached")) GetCachedBrainVolStats = 0;
    else if ( !strcmp(option, "--brain-vol-from-seg") )
    {
      BrainVolFromSeg = 1;
    }
    else if ( !strcmp(option, "--subcortgray") )
    {
      DoSubCortGrayVol = 1;
    }
    else if ( !strcmp(option, "--supratent") )
    {
      DoSupraTent = 1;
    }
    else if ( !strcmp(option, "--totalgray") )
    {
      DoTotalGrayVol = 1;
    }
    else if ( !strcmp(option, "--etiv") )
    {
      DoETIV = 1;
    }
    else if ( !strcmp(option, "--etiv-only") )
    {
      DoETIVonly = 1;
    }
    else if ( !strcmp(option, "--old-etiv-only") )
    {
      DoOldETIVonly = 1;
    }
    else if ( !strcmp(option, "--surf-ctx-vol") )
    {
      DoSurfCtxVol = 1;
    }
    else if ( !strcmp(option, "--surf-wm-vol") )
    {
      DoSurfWMVol = 1;
    }
    else if ( !strcmp(option, "--no-global-stats") ){
      DoSurfWMVol = DoSurfCtxVol = DoSupraTent = BrainVolFromSeg = DoSubCortGrayVol = DoTotalGrayVol = DoEuler = 0 ;
    }
    else if ( !strcmp(option, "--euler") )
    {
      DoEuler = 1;
    }
    else if ( !strcmp(option, "--abs") )
    {
      DoAbs = 1;
    }
    else if ( !strcmp(option, "--sqr") )
    {
      DoSquare = 1;
    }
    else if ( !strcmp(option, "--sqrt") )
    {
      DoSquareRoot = 1;
    }
    else if ( !strcmp(option, "--snr") )
    {
      DoSNR = 1;
    }
    else if ( !strcmp(option, "--acc") || !strcmp(option, "--accumulate") )
    {
      DoAccumulate=1;
    }

    else if ( !strcmp(option, "--mul") )
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      DoMultiply = 1;
      sscanf(pargv[0],"%lf",&MultVal);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--div") )
    {
      if(nargc < 1) argnerr(option,1);
      DoMultiply = 1;
      sscanf(pargv[0],"%lf",&MultVal);
      MultVal = 1.0/MultVal;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--robust"))
    {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&RobustPct);
      UseRobust = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--talxfm") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      talxfmfile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--sd") )
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      FSENVsetSUBJECTS_DIR(pargv[0]);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--ctab-default") )
    {
      FREESURFER_HOME = getenv("FREESURFER_HOME");
      ctabfile = (char *) calloc(sizeof(char),1000);
      sprintf(ctabfile,"%s/FreeSurferColorLUT.txt",FREESURFER_HOME);
      printf("Using defalt ctab %s\n",ctabfile);
    }
    else if ( !strcmp(option, "--ctab-gca") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      gcafile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--seg") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      SegVolFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--seg-from-input") )
    {
      SegFromInput=1;
      UserSegIdList[0] = 1;
      nUserSegIdList = 1;
    }
    else if (!strcmp(option, "--seg-erode")){
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nErodeSeg);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--vox") )
    {
      if (nargc < 1)
      {
        argnerr(option,3);
      }
      sscanf(pargv[0],"%d",&Vox[0]);
      sscanf(pargv[1],"%d",&Vox[1]);
      sscanf(pargv[2],"%d",&Vox[2]);
      DoVox = 1;
      nargsused = 3;
    }
    else if ( !strcmp(option, "--in") || !strcmp(option, "--i") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      InVolFile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--reg"))
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      InVolRegFile = pargv[0];
      InVolReg = regio_read_registermat(InVolRegFile);
      if(InVolReg == NULL)
      {
        exit(1);
      }
      nargsused = 1;
    }
    else if(!strcmp(option, "--regheader"))
    {
      InVolRegHeader = 1;
    }
    else if(!strcmp(option, "--xfm2etiv"))
    {
      if(nargc < 2) argnerr(option,1);
      double etiv_scale_factor = 1948.106, determinant = 0, atlas_icv;
      atlas_icv = MRIestimateTIV(pargv[0],etiv_scale_factor,&determinant);
      printf("%12.4lf\n",atlas_icv);
      if(strcmp(pargv[1],"nofile")!=0){
	FILE *fp = fopen(pargv[1],"w");
	if(fp == NULL) exit(1);
	fprintf(fp,"%12.4lf\n",atlas_icv);
	fclose(fp);
      }
      exit(0);
    }
    else if ( !strcmp(option, "--in-intensity-name") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      InIntensityName = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--in-intensity-units") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      InIntensityUnits = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--brainmask") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      BrainMaskFile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--id") )
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) )
      {
        sscanf(pargv[nth],"%d",&UserSegIdList[nUserSegIdList]);
        nUserSegIdList++;
        nth++;
      }
      nargsused = nth;
    }
    else if ( !strcmp(option, "--excl-ctxgmwm") )
    {
      DoExclSegId = 1;
      DoExclCtxGMWM = 1;
      ExclSegIdList[nExcl] =  2;
      nExcl++;
      ExclSegIdList[nExcl] =  3;
      nExcl++;
      ExclSegIdList[nExcl] = 41;
      nExcl++;
      ExclSegIdList[nExcl] = 42;
      nExcl++;
    }
    else if( !strcmp(option, "--excludeid") ||
             !strcmp(option, "--exclude") )
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) )
      {
        sscanf(pargv[nth],"%d",&ExclSegIdList[nExcl]);
        nExcl ++;
        nth ++;
      }
      DoExclSegId = 1;
      nargsused = nth;
    }
    else if ( !strcmp(option, "--mask") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      MaskVolFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--masksign"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      masksign = pargv[0];
      nargsused = 1;
      if (strncasecmp(masksign,"abs",3) &&
          strncasecmp(masksign,"pos",3) &&
          strncasecmp(masksign,"neg",3))
      {
        fprintf(stderr,"ERROR: mask sign = %s, must be abs, pos, or neg\n",
                masksign);
        exit(1);
      }
    }
    else if (!strcmp(option, "--maskthresh"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%f",&maskthresh);
      nargsused = 1;
    }
    else if (!strcmp(option, "--maskframe"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&maskframe);
      nargsused = 1;
    }
    else if (!strcmp(option, "--maskerode"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&maskerode);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--maskinvert"))
    {
      maskinvert = 1;
    }
    else if ( !strcmp(option, "--sd") ) {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      FSENVsetSUBJECTS_DIR(SUBJECTS_DIR);
      nargsused = 1;
    }
    else if( !strcmp(option, "--sum") || !strcmp(option, "--o") )
    {
      if (nargc < 1) argnerr(option,1);
      StatTableFile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--sum-in") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      StatTableFile = pargv[0];
      StatSumTable = LoadStatSumFile(StatTableFile,&nsegid);
      printf("Found %d\n",nsegid);
      DumpStatSumTable(StatSumTable,nsegid);
      exit(1);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--avgwf") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      FrameAvgFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--avgwf-remove-mean") ) RmFrameAvgMn = 1;
    else if ( !strcmp(option, "--sumwf") )
    {
      if (nargc < 1) argnerr(option,1);
      FrameAvgFile = pargv[0];
      DoFrameAvg = 1;
      DoFrameSum = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--sfavg") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      SpatFrameAvgFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--avgwfvol") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      FrameAvgVolFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--ctab") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      ctabfile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--ctab-out") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      ctabfileOut = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--frame"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    }
    else if (!strcmp(option, "--subject"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--qa-stats"))
    {
      if (nargc < 2) argnerr(option,2);
      subject = pargv[0];
      err=CountEdits(subject,pargv[1]);
      exit(err);
      nargsused = 2;
    }
    else if (!strcmp(option, "--annot"))
    {
      if (nargc < 3)
      {
        argnerr(option,1);
      }
      subject = pargv[0];
      hemi    = pargv[1];
      annot   = pargv[2];
      nargsused = 3;
    }
    else if (!strcmp(option, "--surf")){
      if (nargc < 1) argnerr(option,1);
      whitesurfname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--slabel")) {
      if(nargc < 3) argnerr(option,3);
      subject = pargv[0];
      hemi    = pargv[1];
      LabelFile = pargv[2];
      ExclSegId = 0;
      DoExclSegId = 1;
      nargsused = 3;
    }
    else if (!strcmp(option, "--label-thresh")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&LabelThresh);
      UseLabelThresh = 1;
      nargsused = 1;
    }
    else if (!strcmp(option, "--segbase"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&segbase);
      nargsused = 1;
    }
    else if (!strcmp(option, "--synth"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%ld",&seed);
      synth = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--pv") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      PVVolFile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--gtm-default-seg-merge"))
      GTMdefaultSegReplacmentList(&nReplace,&(SrcReplace[0]),&(TrgReplace[0]));
    else if(!strcasecmp(option, "--gtm-default-seg-merge-choroid")){
      GTMdefaultSegReplacmentList(&nReplace,&(SrcReplace[0]),&(TrgReplace[0]));
      nReplace -= 2;       // Last two items are choroid.
    }
    else if(!strcmp(option, "--replace-file")){
      if(nargc < 1) CMDargNErr(option,1);
      int err=GTMloadReplacmentList(pargv[0],&nReplace,&(SrcReplace[0]),&(TrgReplace[0]));
      if(err) exit(1);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--replace")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplace[nReplace]);
      sscanf(pargv[1],"%d",&TrgReplace[nReplace]);
      nReplace++;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--seed")) {
      if(nargc < 1) CMDargNErr(option,1);
      setRandomSeed(atol(pargv[0])) ;
      printf("setting seed for random number genererator to %d\n",
            atoi(pargv[0])) ;
      nargsused = 1;
    }
    else
    {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
      {
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      }
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
#include "mri_segstats.help.xml.h"
static void print_usage(void)
{
  outputHelpXml(mri_segstats_help_xml,
                mri_segstats_help_xml_len);
}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
  {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  }
  else
  {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(SegVolFile == NULL && annot == NULL && LabelFile == NULL &&
     DoETIVonly == 0 && DoOldETIVonly == 0 && SegFromInput==0){
    printf("ERROR: must specify a segmentation volume\n");
    exit(1);
  }
  if(SegFromInput==1 && InVolFile==0){
    printf("ERROR: must specify an input volume with --seg-from-input\n");
    exit(1);
  }
  if (StatTableFile == NULL && FrameAvgFile == NULL && DoETIVonly == 0 && DoOldETIVonly == 0 && FrameAvgVolFile==NULL)
  {
    printf("ERROR: must specify an output table file\n");
    exit(1);
  }
  if (DoFrameAvg && InVolFile == NULL)
  {
    printf("ERROR: cannot do frame average without input volume\n");
    exit(1);
  }
  if (DoETIV && subject == NULL)
  {
    printf("ERROR: need subject with --etiv\n");
    exit(1);
  }
  if (DoSupraTent && !DoSurfCtxVol)
  {
    printf("ERROR: need --surf-ctx-vol  with --supratent\n");
    exit(1);
  }
  if (ctabfile != NULL && gcafile != NULL)
  {
    printf("ERROR: cannot specify ctab and gca\n");
    exit(1);
  }
  if(DoSurfCtxVol && subject == NULL)
  {
    printf("ERROR: need --subject with --surf-ctx-vol\n");
    exit(1);
  }
  if(DoTotalGrayVol && !DoSurfCtxVol)
  {
    printf("ERROR: need --surf-ctx-vol with --totalgray\n");
    exit(1);
  }
  if(ctabfileOut && !ctabfile)
  {
    printf("ERROR: need an input ctab to create output ctab\n");
    exit(1);
  }
  if(ctabfileOut)
  {
    if(!strcmp(ctabfileOut,ctabfile))
    {
      printf("ERROR: output ctab is the same as input\n");
      exit(1);
    }
  }
  if (masksign == NULL)
  {
    masksign = "abs";
  }
  if(DoEuler && subject == NULL){
    printf("ERROR: need subject with --euler\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"whitesurfname  %s\n",whitesurfname);
  fprintf(fp,"UseRobust  %d\n",UseRobust);
  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if (len < 2)
  {
    return(0);
  }

  if (flag[0] == '-' && flag[1] != '-')
  {
    return(1);
  }
  return(0);
}

/* ----------------------------------------------------------
   MRIsegCount() - returns the number of times the given
   segmentation id appears in the volume.
   --------------------------------------------------------- */
int MRIsegCount(MRI *seg, int id, int frame)
{
  int nhits, v, c,r,s;
  nhits = 0;
  for (c=0; c < seg->width; c++)
  {
    for (r=0; r < seg->height; r++)
    {
      for (s=0; s < seg->depth; s++)
      {
        v = (int) MRIgetVoxVal(seg,c,r,s,frame);
        if (v == id)
        {
          nhits ++;
        }
      }
    }
  }
  return(nhits);
}
/*------------------------------------------------------------*/
STATSUMENTRY *LoadStatSumFile(char *fname, int *nsegid)
{
  FILE *fp;
  char tmpstr[1000];
  STATSUMENTRY *StatSumTable, *e;

  fp = fopen(fname,"r");
  if (fp == NULL)
  {
    printf("ERROR: cannot open %s\n",fname);
    exit(1);
  }

  // Count the number of entries
  *nsegid = 0;
  while (fgets(tmpstr,1000,fp) != NULL)
  {
    if (tmpstr[0] == '#')
    {
      continue;
    }
    (*nsegid)++;
  }
  fclose(fp);

  StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),*nsegid);

  // Now actually read it in
  fp = fopen(fname,"r");
  *nsegid = 0;
  while (fgets(tmpstr,1000,fp) != NULL)
  {
    if (tmpstr[0] == '#')
    {
      continue;
    }
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
int DumpStatSumTable(STATSUMENTRY *StatSumTable, int nsegid)
{
  int n;
  for (n=0; n < nsegid; n++)
  {
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

//-----------------------------------------------------
/*!
\fn int CountEdits(char *subject, char *outfile)
\brief Prints out number of control points and number 
of wm, brainmask, and aseg edits. Right now it just
prints to a screen and/or to a file. It would be nice
to have this in the aseg.stats file at some point.
\param subject 
\param outfile saves results in outfile
*/
int CountEdits(char *subject, char *outfile)
{
  char *SUBJECTS_DIR;
  char sd[4000],tmpstr[4000];
  MPoint *pArray = 0;
  int count = 0,useRealRAS = 0;
  int c,r,s;
  MRI *mri, *mri2;
  int nWMErase, nWMFill, nBMErase, nBMClone, nASegChanges;
  double v1,v2;
  FILE *fp;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(sd,"%s/%s",SUBJECTS_DIR,subject);

  int req = snprintf(tmpstr,STRLEN,"%s/tmp/control.dat",sd);  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  count = 0;
  if(fio_FileExistsReadable(tmpstr)){
    pArray = MRIreadControlPoints(tmpstr, &count, &useRealRAS);
    free(pArray);
  }

  req = snprintf(tmpstr,STRLEN,"%s/mri/wm.mgz",sd);  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(tmpstr);
  if(mri == NULL) return(1);
  nWMErase = 0;
  nWMFill = 0;
  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	v1 = MRIgetVoxVal(mri,c,r,s,0);
	if(v1 == WM_EDITED_OFF_VAL) nWMErase++;
	if(v1 == WM_EDITED_ON_VAL)  nWMFill++;
      }
    }
  }
  MRIfree(&mri);

  req = snprintf(tmpstr,STRLEN,"%s/mri/brainmask.mgz",sd);  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(tmpstr);
  if(mri == NULL) return(1);
  req = snprintf(tmpstr,STRLEN,"%s/mri/brainmask.auto.mgz",sd);  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri2 = MRIread(tmpstr);
  if(mri2 == NULL) return(1);
  nBMErase = 0;
  nBMClone = 0;
  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	v1 = MRIgetVoxVal(mri,c,r,s,0);
	v2 = MRIgetVoxVal(mri2,c,r,s,0);
	if(v1 == v2) continue;
	if(v1 == 1) nBMErase++;
	else        nBMClone++;
      }
    }
  }
  MRIfree(&mri);
  MRIfree(&mri2);

  req = snprintf(tmpstr,STRLEN,"%s/mri/aseg.mgz",sd);  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(tmpstr);
  if(mri == NULL) return(1);
  req = snprintf(tmpstr,STRLEN,"%s/mri/aseg.auto.mgz",sd);   
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri2 = MRIread(tmpstr);
  if(mri2 == NULL) return(1);
  nASegChanges = 0;
  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	v1 = MRIgetVoxVal(mri,c,r,s,0);
	v2 = MRIgetVoxVal(mri2,c,r,s,0);
	if(v1 == v2) continue;
	nASegChanges++;
      }
    }
  }
  MRIfree(&mri);
  MRIfree(&mri2);

  // Note: ?h.orig.nofix files might not exist in longitudinal;
  // number of holes will be 0 for long anyway
  int nvertices, nfaces, nedges;
  int lheno, rheno, lhholes, rhholes, totholes;
  MRIS *mris;
  sprintf(tmpstr,"%s/%s/surf/lh.orig.nofix",SUBJECTS_DIR,subject);
  if(fio_FileExistsReadable(tmpstr)){
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    lheno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    MRISfree(&mris);
    lhholes = 1-lheno/2;
  } else lhholes = 0;
  sprintf(tmpstr,"%s/%s/surf/rh.orig.nofix",SUBJECTS_DIR,subject);
  if(fio_FileExistsReadable(tmpstr)){
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    rheno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    MRISfree(&mris);
    rhholes = 1-rheno/2;
  }
  else rhholes = 0;
  totholes = lhholes+rhholes;

  double determinant = 0;
  double etiv_scale_factor = 1948.106;
  double atlas_icv=0;
  sprintf(tmpstr,"%s/%s/mri/transforms/talairach.xfm",SUBJECTS_DIR,subject);
  atlas_icv = MRIestimateTIV(tmpstr,etiv_scale_factor,&determinant);

  std::vector<double> BrainVolStats = ReadCachedBrainVolumeStats(subject, SUBJECTS_DIR);
  double MaskVolToETIV;
  MaskVolToETIV = BrainVolStats[12] / atlas_icv;

  float *wmstats;
  // Erode=3, trim the top and bottom 2% when computing WM mean, std, etc
  wmstats = WMAnatStats(subject, "norm.mgz", 3, 2);

  // Compute gray/white contrast, its spatial stddev, cnr = mean/std
  double gwconmeansum=0, gwconvarsum=0;
  int hemi;
  for(hemi = 0; hemi < 2; hemi++){
    LABEL *clabel;
    MRI *wgcon;
    char hemistr[3];
    if(hemi==0) memcpy(hemistr,"lh",2);
    if(hemi==1) memcpy(hemistr,"rh",2);
    sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemistr);
    mris = MRISread(tmpstr);
    if(mris==NULL) exit(1);
    sprintf(tmpstr,"%s/%s/label/%s.cortex.label",SUBJECTS_DIR,subject,hemistr);
    clabel = LabelRead(NULL,tmpstr);
    if(clabel == NULL) exit(1);
    sprintf(tmpstr,"%s/%s/surf/%s.w-g.pct.mgh",SUBJECTS_DIR,subject,hemistr);
    wgcon = MRIread(tmpstr);
    seg = MRIalloc(mris->nvertices,1,1,MRI_INT);
    int n;
    for (n = 0; n < clabel->n_points; n++){
      MRIsetVoxVal(seg,clabel->lv[n].vno,0,0,0, 1);
    }
    float min, max, range, mean, std;
    MRIsegStats(seg, 1, wgcon, 0, &min, &max, &range, &mean, &std);
    gwconmeansum += mean; gwconvarsum += (std*std);
    MRISfree(&mris);
    LabelFree(&clabel);
    MRIfree(&wgcon);
    MRIfree(&seg);
    printf(" %s cnrstats: %6.3f %6.3f %6.3f\n",hemistr,mean,std,mean/std);
  }
  double gwconmean = gwconmeansum/2.0;
  double gwconstd  = sqrt(gwconvarsum/2.0);

  printf("%s nc %3d, nWMErase %3d, nWMFill %3d, nBMErase %3d, nBMClone %3d, nASegChanges %3d, "
	 "lhholes %4d, rhholes %4d, MaskVolToETIV %7.5f\n",
	 subject,count,nWMErase,nWMFill,nBMErase,nBMClone,nASegChanges,
	 lhholes,rhholes,MaskVolToETIV);
  printf("wmstats: %6.2f %6.2f %6.2f %6.2f %6.2f\n",wmstats[0],wmstats[1],
	 wmstats[2],wmstats[3],wmstats[4]);
  printf("cnrstats: %6.3f %6.3f %6.3f\n",gwconmean,gwconstd,gwconmean/gwconstd);


  if(outfile){
    fp = fopen(outfile,"w");
    fprintf(fp,"%s %3d    %4d %4d    %4d %4d   %4d  %4d %4d %4d   %7.5f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.3f %6.3f %6.3f\n",
	    subject,count,nWMErase,nWMFill,nBMErase,nBMClone,nASegChanges,
	    lhholes,rhholes,totholes,MaskVolToETIV,wmstats[0],wmstats[1],
	    wmstats[2],wmstats[3],wmstats[4],wmstats[0]/wmstats[1],
	    gwconmean,gwconstd,gwconmean/gwconstd);
    fclose(fp);
  }


  return(0);
}

float *WMAnatStats(const char *subject, const char *volname, int nErodes, float Pct)
{
  char sd[4000],tmpstr[4000];
  float *stats,val;
  MATRIX *v;
  int wmids[12] = {2, 41, 7, 46, 251, 252, 253, 254, 255, 77, 78, 79};
  int nwmids = 12;
  int c,r,s,n,Matched,nhits;
  MRI *apas, *wmvol;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(sd,"%s/%s",SUBJECTS_DIR,subject);
  sprintf(tmpstr,"%s/%s/mri/aparc+aseg.mgz",SUBJECTS_DIR,subject);
  apas = MRIread(tmpstr);
  if(apas == NULL) return(NULL);

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,volname);
  wmvol = MRIread(tmpstr);
  if(wmvol == NULL) return(NULL);

  for(c=0; c < apas->width; c++) {
    for(r=0; r < apas->height; r++) {
      for(s=0; s < apas->depth; s++) {
	val = MRIgetVoxVal(apas,c,r,s,0);
	Matched = 0;
	for(n=0; n < nwmids; n++){
	  if(fabs(val - wmids[n]) < 2*FLT_MIN){
	    MRIsetVoxVal(apas,c,r,s,0,1);
	    Matched = 1;
	    break;
	  }
	}
	if(!Matched) MRIsetVoxVal(apas,c,r,s,0,0);
      }
    }
  }

  if(nErodes > 0){
    printf("Eroding %d voxels in 3d\n",nErodes);
    for(n=0; n<nErodes; n++) MRIerode(apas,apas);
  }

  nhits = 0;
  for(c=0; c < apas->width; c++) {
    for(r=0; r < apas->height; r++) {
      for(s=0; s < apas->depth; s++) {
	val = MRIgetVoxVal(apas,c,r,s,0);
	if(val < 0.5) continue;
	nhits ++;
      }
    }
  }

  v = MatrixAlloc(nhits,1,MATRIX_REAL);
  nhits = 0;
  for(c=0; c < apas->width; c++) {
    for(r=0; r < apas->height; r++) {
      for(s=0; s < apas->depth; s++) {
	val = MRIgetVoxVal(apas,c,r,s,0);
	if(val < 0.5) continue;
	v->rptr[nhits+1][1] = MRIgetVoxVal(wmvol,c,r,s,0);
	nhits ++;
      }
    }
  }

  stats = (float *) calloc(5,sizeof(float));
  MRIsegStatsRobust(apas, 1, wmvol, 0, &stats[2],&stats[3],&stats[4],
		    &stats[0],&stats[1],Pct);
  //stats[1] = sqrt(VectorVar(v,&stats[0]));
  //stats[4] = VectorRange(v, &stats[2], &stats[3]);

  MRIfree(&apas);
  MRIfree(&wmvol);
  MatrixFree(&v);


  return(stats);
}

