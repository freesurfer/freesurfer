/**
 * @brief fixup the wm vol based on segmentations found in aseg
 *
 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "cma.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "tags.h"
#include "mrimorph.h"
#include "const.h"
#include "transform.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"

static void
DiagBreak2() {}

#define SPACKLE_MTL 1
#if SPACKLE_MTL
static int distance_to_lateral_edge(MRI *mri_seg, int label, int x, int y, int z, int left) ;
static int distance_to_superior_edge(MRI *mri_seg, int label, int x, int y, int z) ;
static int distance_to_inferior_edge(MRI *mri_seg, int label, int x, int y, int z) ;
static int anterior_edge_of_amygdala(MRI *mri_aseg, int x, int y, int z) ;

static int medial_edge_of_hippocampus(MRI *mri_seg, int x, int y, int z, int dx);
static int MRInonfilledInNbhd(MRI *mri, int x, int y, int z, int whalf) ;
static int distance_to_nonzero(MRI *mri_wm, int x,
                               int y, int z, int dx, int dy,
                               int dz, int max_dist) ;

static int distance_to_zero(MRI *mri_wm, int x,
                            int y, int z, int dx, int dy,
                            int dz, int max_dist) ;

#endif
static int remove_medial_voxels(MRI *mri_roi, MRI *mri_aseg, int left) ;
static int remove_gray_matter_voxels(MRI *mri_roi, MRI *mri_aseg) ;
static int remove_unknown_voxels(MRI *mri_roi, MRI *mri_aseg, int left) ;
static int distance_to_anterior_edge(MRI *mri_seg, int label, int x, int y, int z) ;
static int remove_lateral_and_anterior_hippocampus(MRI *mri_roi, MRI *mri_aseg, int left) ;
static int remove_anterior_and_superior_amygdala(MRI *mri_roi, MRI *mri_aseg) ;
static int remove_paths_to_cortex(MRI *mri_wm, MRI *mri_T1, MRI *mri_aseg) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label);
static int edit_segmentation(MRI *mri_im, MRI *mri_T1, MRI *mri_seg) ;
static int spackle_wm_superior_to_mtl(MRI *mri_wm, MRI *mri_T1, MRI *mri_aseg) ;
static int distance_to_label(MRI *mri_labeled, int label, int x,
                             int y, int z, int dx, int dy,
                             int dz, int max_dist) ;


const char *Progname ;

static int lh_only = 0 ;
static int rh_only = 0 ;

static int fillven = 1 ;
static int keep_edits = 0 ;
static int keep_edits_input = 0 ;
static void usage_exit(int code) ;

static int fcd = 0 ;
int FixSCMHA = 0;
int FixSCMHANdil = 0;

double FixEntoWMLhVal,FixEntoWMRhVal;
int FixEntoWMLevel;
MRI *entowm=NULL;

int main(int argc, char *argv[])
{
  MRI    *mri_wm, *mri_aseg, *mri_T1 ;
  Timer then ;
  int    msec, nargs ;
  char *output_file_name,*input_file_name, *edits_file_name ;

  std::string cmdline = getAllInfo(argc, argv, "mri_edit_wm_with_aseg");

  nargs = handleVersionOption(argc, argv, "mri_edit_wm_with_aseg");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  for(int i=0; i<argc; i++) printf("%s ",argv[i]);
  printf("\n");
  fflush(stdout);

  then.reset() ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  Progname = argv[0] ;

  for ( ; argc > 1 && (*argv[1] == '-') ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  output_file_name = argv[4] ;
  if (argc < 5)
  {
    usage_exit(1) ;  /* will exit */
  }

  input_file_name = argv[1];
  printf("reading wm segmentation from %s...\n", input_file_name) ;
  mri_wm = MRIread(input_file_name) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read wm volume from %s",
              Progname, input_file_name) ;

  if (mri_wm->type != MRI_UCHAR)
  {
    ErrorExit(ERROR_UNSUPPORTED, "%s: volume %s must be MRI_UCHAR", Progname, argv[1]) ;
  }

  MRIaddCommandLine(mri_wm, cmdline) ;
  mri_T1 = MRIread(argv[2]) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1/brain volume from %s",
              Progname, argv[2]) ;

  if (mri_T1->type != MRI_UCHAR)
  {
    MRI *mri_tmp ;

    mri_tmp = MRIchangeType(mri_T1, MRI_UCHAR, 0, 1,1);
    MRIfree(&mri_T1) ;
    mri_T1 = mri_tmp ;
  }

  mri_aseg = MRIread(argv[3]) ;
  if (!mri_aseg)
    ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s",
              Progname, argv[3]) ;


  remove_paths_to_cortex(mri_wm, mri_T1, mri_aseg) ;
  edit_segmentation(mri_wm, mri_T1, mri_aseg) ;
  spackle_wm_superior_to_mtl(mri_wm, mri_T1, mri_aseg) ;
  if(entowm)  MRIfixEntoWM(mri_wm, entowm, FixEntoWMLevel, FixEntoWMLhVal, FixEntoWMRhVal);

  if(FixSCMHA){
    FixSubCortMassHA fscmha;
    fscmha.aseg = mri_aseg;
    fscmha.subcorticalmass = mri_wm;
    fscmha.nDilate = FixSCMHANdil;
    fscmha.FixSCM();
    MRIfree(&fscmha.mask);
  }
  if (keep_edits)
  {
    MRI *mri_old ;

    if (! keep_edits_input)
    {
      edits_file_name = output_file_name;
    }
    else
    {
      edits_file_name = input_file_name;
    }
    printf("propagating editing to output volume from %s\n",edits_file_name) ;
    mri_old = MRIread(edits_file_name) ;
    if (!mri_old)
    {
      ErrorPrintf(ERROR_NOFILE, "%s: could not read file %s to preserve edits",
                  Progname, edits_file_name) ;
      exit(1);
    }
    MRIcopyLabel(mri_old, mri_wm, WM_EDITED_ON_VAL) ;
    MRIcopyLabel(mri_old, mri_wm, WM_EDITED_OFF_VAL) ;
    MRIfree(&mri_old) ;
  }
  if (lh_only || rh_only)
  {
    int x, y, z, label ;
    for (x = 0 ; x < mri_wm->width ; x++)
      for (y = 0 ; y < mri_wm->height ; y++)
	for (z = 0 ; z < mri_wm->depth ; z++)
	{
	  label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	  if ((lh_only && IS_RH_CLASS(label)) ||
	      (rh_only && IS_LH_CLASS(label)))
	    MRIsetVoxVal(mri_wm, x, y, z, 0, 0) ;
	  
	}
  }
  printf("writing edited volume to %s....\n", output_file_name) ;
  MRIwrite(mri_wm, output_file_name) ;

  
  msec = then.milliseconds() ;
  fprintf(stderr, "auto filling took %2.2f minutes\n",
          (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "dilate")) {}
  else if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage_exit(0);
  }
  else if (!strcmp(option, "fillven"))
  {
    fillven = atoi(argv[2]) ;
    printf("%sfilling ventricles\n", fillven == 0 ? "not " : "") ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "keep"))
  {
    keep_edits = 1 ;
    fprintf(stderr, "preserving editing changes in output volume...\n");
  }
  else if (!stricmp(option, "fcd"))
  {
    fcd = 1 ;
    fprintf(stderr, "preserving focal cortical dysplasias - not filling non-wm lesions\n");
  }
  else if (!stricmp(option, "lh"))
  {
    lh_only = 1 ;
    fprintf(stderr, "assuming input is only lh\n") ;
  }
  else if (!stricmp(option, "rh"))
  {
    rh_only = 1 ;
    fprintf(stderr, "assuming input is only rh\n") ;
  }
  else if (!stricmp(option, "keep-in"))
  {
    keep_edits = 1 ;
    keep_edits_input = 1 ;
    fprintf(stderr, "preserving editing changes in input volume...\n");
  }
  else if (!stricmp(option, "fix-scm-ha"))
  {
    FixSCMHA = 1 ;
    FixSCMHANdil = atoi(argv[2]) ; // usually set to 1
    printf("FixSCM HA %d\n",FixSCMHANdil);
    nargs = 1;
  }
  else if (!stricmp(option, "fix-scm-ha-only"))
  {
    // mri_edit_wm_with_aseg -fix-scm-ha-only aseg.presurf.mgz SCM.mgz ndil out.mgz
    // SCM = wm.seg.mgz, wm.mgz, filled.mgz, etc; ndil usually = 1
    FixSubCortMassHA fscmha;
    fscmha.aseg = MRIread(argv[2]);
    if(fscmha.aseg==NULL) exit(1);
    fscmha.subcorticalmass = MRIread(argv[3]);
    if(fscmha.subcorticalmass==NULL) exit(1);
    fscmha.nDilate = atoi(argv[4]) ;    
    fscmha.FixSCM();
    MRIfree(&fscmha.mask);
    int err = MRIwrite(fscmha.subcorticalmass,argv[5]);
    exit(err);
  }
  else if(!stricmp(option, "fix-ento-wm") || !stricmp(option, "sa-fix-ento-wm"))
  {
    // This is a bit of a hack to fill in WM in the gyrus ambiens that
    // is often so thin that it looks like GM. entowm.mgz is the output of mri_entowm_seg
    // mri_edit_wm_with_aseg -fix-entowm entowm.mgz level lhsetval rhsetval
    // mri_edit_wm_with_aseg -sa-fix-entowm entowm.mgz level lhsetval rhsetval invol outvol
    // for wm.seg.mgz/wm.asegedit.mgz use level=3 and setval=250 for both lh and rh
    // for brain.finalsurfs use level=2 and setval=255 for both lh and rh
    // for filled use level=3 and setval =255 for lh and =127 for rh
    // level=1 only fill entowm without gyrus ambiens (probably not useful)
    // level=2 only fill gyrus ambiens without entowm (more conservative)
    // level=3 fill both entowm and gyrus ambiens 
    entowm = MRIread(argv[2]);
    if(entowm==NULL) exit(1);
    sscanf(argv[3],"%d",&FixEntoWMLevel);
    sscanf(argv[4],"%lf",&FixEntoWMLhVal);
    sscanf(argv[5],"%lf",&FixEntoWMRhVal);
    printf("Fixing entowm %d %g %g\n",FixEntoWMLevel,FixEntoWMLhVal,FixEntoWMRhVal);
    if(!stricmp(option, "sa-fix-ento-wm")){
      MRI *invol = MRIread(argv[6]);
      if(invol==NULL) exit(1);
      MRIfixEntoWM(invol, entowm, FixEntoWMLevel, FixEntoWMLhVal, FixEntoWMRhVal);
      int err = MRIwrite(invol,argv[7]);
      exit(err);
    }
    nargs = 4;
  }
  else if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'H':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#include "mri_edit_wm_with_aseg.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_edit_wm_with_aseg_help_xml,
                mri_edit_wm_with_aseg_help_xml_len);
  exit(code) ;
}


static int
edit_segmentation(MRI *mri_wm, MRI *mri_T1, MRI *mri_seg)
{
  int   i, width, height, depth, x, y, z, label, non, noff, xi, yi, zi,  xk, yk, zk, nchanged,
        wsize, alabel, hlabel, slabel, olabel, left;
  MRI   *mri_filled ;

  mri_filled =  MRIclone(mri_wm,  NULL);

  width = mri_wm->width ;
  height = mri_wm->height ;
  depth = mri_wm->depth ;

  non = noff = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = height-2 ; y > 0 ; y--)
    {
      for (x = 1 ; x < width-1 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;

        left = 0 ;
        switch (label)
        {
        case Unknown:
          wsize=5 ;
          if (MRIlabelsInNbhd(mri_seg, x, y, z, (wsize-1)/2, Unknown) < (wsize*wsize*wsize-1))
          {
            break ;
          }

#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          /* !!! no break - erase unknown if it is surrounded by only  unknowns */

          /* erase these  labels */
        case Left_Cerebellum_White_Matter:
        case Left_Cerebellum_Exterior:
        case Left_Cerebellum_Cortex:
        case Right_Cerebellum_White_Matter:
        case Right_Cerebellum_Exterior:
        case Right_Cerebellum_Cortex:
        case Optic_Chiasm: // added by dng 10/11/2019
#if 0
	  // I don't think these make sense
        case Right_Cerebral_Cortex:
        case Left_Cerebral_Cortex:
#endif

#if 0
          /* otherwise will never be able to find pons */
        case Brain_Stem:
        case Left_VentralDC:
        case Right_VentralDC:
        case Left_Substancia_Nigra:
        case Right_Substancia_Nigra:
#endif
          if ((neighborLabel(mri_seg, x,y,z,1,Left_Cerebral_Cortex) == 0) &&
              (neighborLabel(mri_seg, x,y,z,1,Right_Cerebral_Cortex) == 0))
          {
            if (MRIvox(mri_wm, x, y, z) >= WM_MIN_VAL)
            {
              MRIvox(mri_wm, x, y, z) = 0 ;
              noff++ ;
            }
          }
          break ;

          /* fill these */
        case non_WM_hypointensities:
        case Left_non_WM_hypointensities:
        case Right_non_WM_hypointensities:
	  if (fcd)
	    break ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Left_Lesion:
        case Right_Lesion:
        case WM_hypointensities:
        case Left_WM_hypointensities:
        case Right_WM_hypointensities:
          // only fill these if they are not adjacent to cortex
          if ((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_Cortex) == 0) &&
              (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_Cortex) == 0) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
          break ;
        case Left_choroid_plexus:
        case Right_choroid_plexus:
          // don't fill in choroid next to inf lat vent
          if (((neighborLabel(mri_seg, x, y, z,1,Left_Inf_Lat_Vent) == 0) &&
               (neighborLabel(mri_seg, x, y, z,1,Right_Inf_Lat_Vent) == 0)) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              printf("filling choroid adjacent to ventricle at (%d, %d, %d)\n", x,y,z) ;
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
            break ;
          }
          break ;
        case Left_Lateral_Ventricle:
        case Right_Lateral_Ventricle:
          if (((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_White_Matter) > 0) ||
               (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_White_Matter) > 0)) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              printf("filling ventricle adjacent to wm at (%d, %d, %d)\n", x,y,z) ;
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
            break ;
          }
          if (fillven == 0)
          {
            break ;
          }
          if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Left_Inf_Lat_Vent:
        case Right_Inf_Lat_Vent:
          xi = (label ==  Left_Inf_Lat_Vent) ?  mri_wm->xi[x+1] :  mri_wm->xi[x-1] ; // lateral
          olabel = MRIgetVoxVal(mri_seg, xi, y, z, 0) ;

          /* don't allow cortex to be directly lateral to inf-lat-vent - should be some wm there
          also don't allow it to be diagonally connected */
          if (IS_CORTEX(olabel) && (MRIvox(mri_wm, xi, y, z) < WM_MIN_VAL))
          {
            if (xi == Gx && y == Gy && z == Gz)
            {
              printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", xi, y, z);
            }
            MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, xi, y, z) = AUTO_FILL ;
            non++ ;
          }

          yi = mri_wm->yi[y+1] ; // inferior
          olabel = MRIgetVoxVal(mri_seg, xi, yi, z, 0) ;
          if (IS_CORTEX(olabel) && (MRIvox(mri_wm, xi, yi, z) < WM_MIN_VAL))
          {
            if (xi == Gx && yi == Gy && z == Gz)
            {
              printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", xi, yi, z);
            }
            MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
            MRIvox(mri_filled, xi, yi, z) = AUTO_FILL ;
            non++ ;
          }

          // for spackling, don't do it if we are too close to superior wm
          if (distance_to_nonzero(mri_wm, x, y, z, 0, -1, 0, 5) <= 3)
          {
            continue ;
          }

          // check diagonally anterior/posterior and inferior
          zi = mri_wm->zi[z-1] ; // posterior
          olabel = MRIgetVoxVal(mri_seg , x, yi, zi, 0) ;
          if (IS_CORTEX(olabel) && (MRIvox(mri_wm, x, yi, zi) < WM_MIN_VAL))
          {
            if (x == Gx && yi == Gy && zi == Gz)
            {
              printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", x, yi, zi);
            }
            MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, zi) = AUTO_FILL ;
            non++ ;
          }
          zi = mri_wm->zi[z+1] ; // anterior
          olabel = MRIgetVoxVal(mri_seg, x, yi, zi, 0) ;
          if (IS_CORTEX(olabel) && (MRIvox(mri_wm, x, yi, zi) < WM_MIN_VAL))
          {
            if (x == Gx && yi == Gy && zi == Gz)
            {
              printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", x, yi, zi);
            }
            MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, zi) = AUTO_FILL ;
            non++ ;
          }

          hlabel = ((label == Left_Lateral_Ventricle) || (label == Left_Inf_Lat_Vent)) ? Left_Hippocampus : Right_Hippocampus ;
          if (distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 10) < 10)
          {
            continue ;  // don't fill ventricular voxels superior to hippo
          }
          if ((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_Cortex) > 0) &&
              (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_Cortex) > 0) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
          yi = mri_wm->yi[y+1] ;
          label = MRIgetVoxVal(mri_seg, x,yi, z, 0) ;
          if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
               (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter) ||
               (label == Unknown))
              && (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
          {
            if (x == Gx && yi == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
            non++ ;
          }
          if (label == Left_Inf_Lat_Vent || label ==  Right_Inf_Lat_Vent)  /* fill inferior wm */
          {
            int xi ;

            xi = label ==  Left_Inf_Lat_Vent ?  mri_wm->xi[x-1] :  mri_wm->xi[x+1] ;
            olabel = MRIgetVoxVal(mri_seg, xi, y, z, 0) ;
#if 0  // no longer needed with path stuff
            /* voxel lateral to this one is not hippocampus   */
            if ((olabel != label) && (MRIvox(mri_wm, xi, y, z) < WM_MIN_VAL) && !IS_AMYGDALA(olabel) &&
                ((distance_to_label(mri_seg, label ==  Left_Inf_Lat_Vent ? Left_Cerebral_White_Matter :
                                    Right_Cerebral_White_Matter, xi, y, z, 0, 1, 0, 5) < 3) ||
                 (distance_to_label(mri_seg, label ==  Left_Inf_Lat_Vent ? Left_Cerebral_Cortex :
                                    Right_Cerebral_Cortex, xi, y, z, 0, 1, 0, 5) < 3)) &&
                !IS_LAT_VENT(olabel))
            {
              if (xi == Gx && y == Gy && z == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
              MRIvox(mri_filled, xi, y, z) = AUTO_FILL ;
              non++ ;
            }
#endif
            yi = mri_wm->yi[y+1] ;
            label = MRIgetVoxVal(mri_seg, x,yi, z, 0) ;
            if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
                 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter))
                && (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
            {
              if (x == Gx && yi == Gy && z == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
              MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
              non++ ;
              yi = mri_wm->yi[y+2] ;
              if (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL)
              {
                if (x == Gx && yi == Gy && z == Gz)
                {
                  DiagBreak2() ;
                }
                MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
                MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
                non++ ;
              }
            }
          }
          break ;
        case Left_Hippocampus:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Hippocampus:
        {
          int xi ;

          // don't mess around with voxels near the medial or superior edge of hippocampus
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1,0,0, 10) <8)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, +1,0,0, 10) <8)
          {
            continue ;
          }
#else
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
#endif
          if (distance_to_nonzero(mri_wm, x, y, z, 0, -1, 0, 5) <= 3)
          {
            continue ;
          }

          // make sure we're not close to superior edge of hippo
          olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
          {
            continue ;
          }
          olabel = left ? Left_Thalamus : Right_Thalamus ;
          if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
          {
            continue ;
          }
          olabel = left ? Left_VentralDC : Right_VentralDC ;
          if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
          {
            continue ;
          }


          xi = label == Right_Hippocampus ?  mri_wm->xi[x-1] :  mri_wm->xi[x+1] ;
          yi = mri_wm->yi[y+1] ;
          olabel = MRIgetVoxVal(mri_seg, xi, y, z, 0) ;
          /* voxel lateral to this one is not hippocampus, and not
          superior to hippocampus, and not far from wm or gm */
#if 0  // not needed anymore with path stuff in place
          if (olabel != label && (MRIvox(mri_wm, xi, y, z) < MIN_WM_VAL) &&
              (distance_to_label(mri_seg, label, xi, y, z, 0, 1, 0, 10) >= 10) &&
              ((distance_to_label(mri_seg, label == Left_Hippocampus ? Left_Cerebral_Cortex : Right_Cerebral_Cortex, xi, y, z, 0, 1, 0, 5) < 3) ||
               (distance_to_label(mri_seg, label == Left_Hippocampus ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter, xi, y, z, 0, 1, 0, 5) < 3) ||
               (distance_to_label(mri_seg, Unknown, xi, y, z, 0, 1, 0, 5) < 3)))
          {
            if (xi == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, xi, y, z) = AUTO_FILL ;
            non++ ;
          }
#endif

          yi = mri_wm->yi[y+1] ;
          olabel = MRIgetVoxVal(mri_seg, xi, yi, z, 0) ;  // diagonally lateral and inferior
          if (IS_CORTEX(olabel) && (MRIvox(mri_wm, xi, yi, z) < WM_MIN_VAL))
          {
            if (xi == Gx && yi == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
            MRIvox(mri_filled, xi, yi, z) = AUTO_FILL ;
            non++ ;
          }


          label = MRIgetVoxVal(mri_seg, x,yi, z, 0) ;

          if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
               (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter))
              && (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
          {
            if (x == Gx && yi == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
            yi = mri_wm->yi[y+2] ;
            non++ ;
            if (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL)
            {
              if (x == Gx && yi == Gy && z == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
              MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;

              non++ ;
            }
          }
          break ;
        }
        case Left_Accumbens_area:
        case Right_Accumbens_area:
        case Left_Caudate:
        case Left_vessel:
        case Right_vessel:
        case Right_Caudate:
        case Left_Putamen:
        case Right_Putamen:
        case Left_Pallidum:
        case Right_Pallidum:
        case Right_Thalamus:
        case Left_Thalamus:
        case Left_VentralDC:
        case Right_VentralDC:
          if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
          break ;
        case Left_Cerebral_White_Matter:
        case Right_Cerebral_White_Matter:
          yi = mri_wm->yi[y-1] ;
          slabel = MRIgetVoxVal(mri_seg, x, yi, z, 0) ;
          if (IS_INF_LAT_VENT(slabel) && MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
              printf("changing voxel (%d, %d, %d) to WM, due to superior inf-lat-vent\n",
                     x, y, z) ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
          break ;
        default:
          break ;
        }
      }
    }
  }


  /* fill in the borders of the ventricle - 2mm thick. This shouldn't affect the folds
  but will prevent small holes from ventricle into wm
  */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Unknown:
          if (((neighborLabel(mri_seg, x, y, z, 1, Left_Lateral_Ventricle) > 0) &&
               (neighborLabel(mri_seg, x, y, z, 1, Left_Cerebral_White_Matter) > 0)) ||
              ((neighborLabel(mri_seg, x, y, z, 1, Right_Lateral_Ventricle) > 0) &&
               (neighborLabel(mri_seg, x, y, z, 1, Right_Cerebral_White_Matter) > 0)))
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
          break ;

        case Left_Lateral_Ventricle:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Lateral_Ventricle:
          olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          if (neighborLabel(mri_seg, x, y, z, 2, olabel) > 0)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
          }
          break ;
        default:
          break ;
        }
      }
    }
  }

  /*
  fill in voxels that were labeled wm by the aseg, but not by  wmfilter, and are
  neighbors  of  voxels that  have been already been filled .
  */
  do
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          if  (MRIvox(mri_filled, x,  y, z) == 0)
          {
            continue  ;
          }
          for (xk = -1 ; xk <= 1 ; xk++)
          {
            xi = mri_filled->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_filled->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = mri_filled->zi[z+zk] ;
                label = MRIgetVoxVal(mri_seg, xi, yi, zi, 0) ;
                if (IS_WM(label) &&  (MRIvox(mri_wm, xi, yi, zi) < WM_MIN_VAL))
                {
                  if (xi == Gx && yi == Gy && zi == Gz)
                  {
                    DiagBreak2() ;
                  }
                  nchanged++ ;
                  MRIvox(mri_wm, xi, yi, zi) = AUTO_FILL ;
#if 0
                  MRIvox(mri_filled, xi, yi, zi) = AUTO_FILL ;
#endif
                  non++ ;
                }
              }
            }
          }
        }
      }
    }
    printf("%d additional wm voxels added\n", nchanged)  ;
  }
  while (nchanged >  0) ;


  // now look for voxels in which there is no wm inferior to hippo and spackle them.
  for (z = 1 ; z < depth-1 ; z++)
  {
    for (y = height-2 ; y > 0 ; y--)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL)
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Cerebral_Cortex:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Cerebral_Cortex:
        case Unknown:
          // look for voxels that are lateral to amygdala, and inf to wm. Should be filled
          if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
          {
            continue ;
          }
          xi = left ? x-1 : x+1 ;  // lateral
          olabel = left ? Left_Amygdala : Right_Amygdala ;
          if (MRIgetVoxVal(mri_seg, xi, y, z, 0) != olabel)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, olabel, x, y, z, 0, 1, 0, 8) < 8)
          {
            continue ;
          }
          non++ ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
          nchanged++ ;
          break ;
        case Left_Amygdala:
        case Left_Inf_Lat_Vent:
        case Left_Hippocampus:
          left = 1 ;
        case Right_Amygdala:
        case Right_Inf_Lat_Vent:
        case Right_Hippocampus:
#if 0  // no longer needed with path stuff
          if (MRIgetVoxVal(mri_seg, x, y+1, z, 0) == label)  // not at inferior border
          {
            continue ;
          }
          if (IS_INF_LAT_VENT(label))  // check to make sure it's not all hippo inferior
          {
            olabel = MRIgetVoxVal(mri_seg, x, y+1, z, 0) ;
            if (IS_HIPPO(olabel) || IS_AMYGDALA(olabel))
            {
              continue ;
            }
          }
          if (IS_HIPPO(label))  // check to make sure it's not all hippo inferior
          {
            olabel = MRIgetVoxVal(mri_seg, x, y+1, z, 0) ;
            if (IS_INF_LAT_VENT(olabel))
            {
              continue ;
            }
            if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
            {
              continue ;
            }

          }
          if (IS_AMYGDALA(label))  // check to make sure it's not all hippo inferior
          {
            olabel = MRIgetVoxVal(mri_seg, x, y+1, z, 0) ;
            if (IS_INF_LAT_VENT(olabel) || IS_HIPPO(olabel))
            {
              continue ;
            }
          }

          // first make sure we aren't on the medial edge of hippo
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
          {
            continue ;
          }
#endif
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }

          // see if there is any wm inferior to this label
          olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          if (distance_to_nonzero(mri_wm, x, y, z, 0, 1, 0, 4) <= 4)
          {
            continue ;  // found wm inferior
          }

          // change either this voxel or the one inferior to it to wm
          if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x, y+1, z))
          {
            yi = y ;
          }
          else
          {
            yi = y+1 ;
          }
          nchanged++ ;
          if (x == Gx && yi == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
#if 0
          MRIvox(mri_filled, xi, yi, zi) = AUTO_FILL ;
#endif
          non++ ;
#endif
          break ;
        }
      }
    }
  }

  /* more spackling. Look for hippocampal or ventricular voxels that have wm
  inferior, but are diagonally connected to non-wm inferior. Fill these.
  */
  for (i = 0 ; i < 1 ; i++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = height-2 ; y > 0 ; y--)
      {
        for (x = 2 ; x < width-2 ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          if (MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL)
          {
            continue ;
          }
          label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
          left = 0 ;
          switch (label)
          {
          case Left_Inf_Lat_Vent:
          case Left_Hippocampus:
          case Left_Amygdala:
            left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Right_Inf_Lat_Vent:
          case Right_Hippocampus:
          case Right_Amygdala:
#if 0  // don't need because of next line
            if (MRIgetVoxVal(mri_seg, x, y+1, z, 0) == label)  // not at inferior border
            {
              continue ;
            }
#endif
            if (MRIvox(mri_wm, x, y+1, z) < MIN_WM_VAL)
            {
              continue ;  // no white matter inferior
            }

            // if on top of a thick piece of wm, this isn't needed
            if (distance_to_zero(mri_wm, x, y, z, 0, 1, 0, 4) > 2)
            {
              continue ;
            }

            // only if we are close to cortex or unknowns below
            olabel = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
            if ((distance_to_label(mri_seg, Unknown, x, y, z, 0, 1, 0, 10) > 5) &&
                (distance_to_label(mri_seg, olabel, x, y, z, 0, 1, 0, 10) > 5))
            {
              continue ;
            }

            // but not if there is wm close above (can cause defects!)
            if (distance_to_nonzero(mri_wm, x, y, z, 0, -1, 0, 5) <= 3)
            {
              continue ;
            }
            olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
            if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
            {
              continue ;
            }
            olabel = left ? Left_Thalamus : Right_Thalamus ;
            if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
            {
              continue ;
            }
            olabel = left ? Left_VentralDC : Right_VentralDC ;
            if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
            {
              continue ;
            }

#if 0  // don't need to do this check, because wm is inferior
            if (IS_INF_LAT_VENT(label))  // check to make sure it's not all hippo inferior
            {
              olabel = MRIgetVoxVal(mri_seg, x, y+1, z, 0) ;
              if (IS_HIPPO(olabel) || IS_AMYGDALA(olabel))
              {
                continue ;
              }
            }
            if (IS_HIPPO(label))  // check to make sure it's not all ventricle inferior
            {
              olabel = MRIgetVoxVal(mri_seg, x, y+1, z, 0) ;
              if (IS_INF_LAT_VENT(olabel))
              {
                continue ;
              }
            }
            if (IS_AMYGDALA(label))  // check to make sure it's not all hippo inferior
            {
              olabel = MRIgetVoxVal(mri_seg, x, y+1, z, 0) ;
              if (IS_INF_LAT_VENT(olabel) || IS_HIPPO(olabel))
              {
                continue ;
              }
            }
#endif

            // make sure we aren't on the medial edge of hippo
#if 0
            if (left)
            {
              if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
              {
                continue ;
              }
            }
            else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
            {
              continue ;
            }
#else
            if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
            {
              continue ;
            }
#endif


            if (MRIvox(mri_wm, x, y+1, z-1) < MIN_WM_VAL) // fill one of these
            {
              if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x, y+1, z-1))
              {
                yi = y ;
                zi = z ;
              }
              else
              {
                yi = y+1 ;
                zi = z-1 ;
              }
              if (x == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
              nchanged++ ;
              non++ ;
            }

            if (MRIvox(mri_wm, x, y+1, z+1) < MIN_WM_VAL) // fill one of these
            {
              if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x, y+1, z+1))
              {
                yi = y ;
                zi = z ;
              }
              else
              {
                yi = y+1 ;
                zi = z+1 ;
              }
              if (x == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
              nchanged++ ;
              non++ ;
            }

            if (MRIvox(mri_wm, x-1, y+1, z) < MIN_WM_VAL) // fill one of these
            {
              if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x-1, y+1, z))
              {
                xi = x ;
                yi = y ;
              }
              else
              {
                xi = x-1 ;
                yi = y+1 ;
              }
              if (xi == Gx && yi == Gy && z == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
              nchanged++ ;
              non++ ;
            }

            if (MRIvox(mri_wm, x+1, y+1, z) < MIN_WM_VAL) // fill one of these
            {
              if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x+1, y+1, z))
              {
                xi = x ;
                yi = y ;
              }
              else
              {
                xi = x-1 ;
                yi = y+1 ;
              }
              if (xi == Gx && yi == Gy && z == Gz)
              {
                DiagBreak2() ;
              }
              MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
              nchanged++ ;
              non++ ;
            }

            break ;
          }
        }
      }
    }
  }

  // spackle the amygdala
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Cerebral_Cortex:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Cerebral_Cortex:
          if (MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL)
          {
            continue ;
          }
          // make sure we aren't on the medial edge of hippo
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
          {
            continue ;
          }
#else
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
#endif

          // look for voxels that are lateral to amygdala, and inf to wm. Should be filled
          if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
          {
            continue ;
          }
          xi = left ? x-1 : x+1 ;  // lateral
          olabel = left ? Left_Amygdala : Right_Amygdala ;
          if (MRIgetVoxVal(mri_seg, xi, y, z, 0) != olabel)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, olabel, x, y, z, 0, 1, 0, 8) < 8)
          {
            continue ;
          }
          non++ ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
          nchanged++ ;
          break ;
        }
      }
    }
  }

#if 0  // no longer needed with path stuff
  // spackle diagonal connectivity topology flaws
  for (z = 1 ; z < depth-1 ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Unknown:
          if (MRIvox(mri_wm, x, y-1,z) >= MIN_WM_VAL)
          {
            continue ;
          }
          if (MRIvox(mri_wm, x, y,z) >= MIN_WM_VAL)
          {
            continue ;
          }
          if ((MRIvox(mri_wm, x, y-1, z-1) >= MIN_WM_VAL) && (MRIvox(mri_wm, x, y-1, z+1) >= MIN_WM_VAL))
          {
            continue ;
          }
          if (!IS_HIPPO(MRIgetVoxVal(mri_seg, x, y-1, z-1, 0)) && !IS_AMYGDALA(MRIgetVoxVal(mri_seg, x, y-1, z-1, 0)) &&
              !IS_INF_LAT_VENT(MRIgetVoxVal(mri_seg, x, y-1, z-1, 0)) &&
              !IS_HIPPO(MRIgetVoxVal(mri_seg, x, y-1, z+1, 0)) && !IS_AMYGDALA(MRIgetVoxVal(mri_seg, x, y-1, z+1, 0)) &&
              !IS_INF_LAT_VENT(MRIgetVoxVal(mri_seg, x, y-1, z+1, 0)))
          {
            continue ;
          }

          if (distance_to_label(mri_seg, Left_Inf_Lat_Vent, x, y, z, 0, 1, 0, 6) < 5)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, Right_Inf_Lat_Vent, x, y, z, 0, 1, 0, 6) < 5)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, Right_Hippocampus, x, y, z, 0, 1, 0, 6) < 5)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, Right_Amygdala, x, y, z, 0, 1, 0, 6) < 5)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, Left_Hippocampus, x, y, z, 0, 1, 0, 6) < 5)
          {
            continue ;
          }
          if (distance_to_label(mri_seg, Left_Amygdala, x, y, z, 0, 1, 0, 6) < 5)
          {
            continue ;
          }
          if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 10)  // only if there is enough wm in vicinity
          {
            continue ;
          }
          if (MRIvox(mri_wm, x, y-1, z-1) < MIN_WM_VAL)
          {
            if (x == Gx && y-1 == Gy && z-1 == Gz)
            {
              DiagBreak2() ;
            }
            non++ ;
            nchanged++ ;
            MRIvox(mri_wm, x, y-1, z-1) = AUTO_FILL ;
          }
          if (MRIvox(mri_wm, x, y-1, z+1) < MIN_WM_VAL)
          {
            if (x == Gx && y-1 == Gy && z+1 == Gz)
            {
              DiagBreak2() ;
            }
            non++ ;
            nchanged++ ;
            MRIvox(mri_wm, x, y-1, z+1) = AUTO_FILL ;
          }
          break ;
        case Left_Cerebral_Cortex:
          left = 1 ;
        case Right_Cerebral_Cortex:
          // make sure we aren't on the medial edge of hippo
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
          {
            continue ;
          }
#else
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
#endif
          slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          // look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
          if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
          {
            continue ;
          }
#endif
          xi = left ? x-1 : x+1 ;  // lateral
          if (MRIvox(mri_wm, xi, y-1, z) >= MIN_WM_VAL)
          {
            continue ;  // only if diagonal voxel isn't on
          }
          alabel = left ? Left_Amygdala : Right_Amygdala ;
          hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
          if ((MRIgetVoxVal(mri_seg, xi, y-1, z, 0) != alabel) && (MRIgetVoxVal(mri_seg, xi, y-1, z, 0) != hlabel))
          {
            continue ;
          }
          // check to make sure we are at inferior border of hippo or amygdala
          if ((distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8) ||
              (distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))
          {
            continue ;
          }
          if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
          {
            continue ;
          }
          non++ ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
          nchanged++ ;
          break ;
        }
      }
    }
  }
#endif
  for (z = 1 ; z < depth-1 ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Cerebral_Cortex:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Cerebral_Cortex:
          // make sure we aren't on the medial edge of hippo
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
          {
            continue ;
          }
#else
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
#endif

          slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          // look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
          if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
          {
            continue ;
          }
#endif
          if (MRIvox(mri_wm, x, y-1, z+1) >= MIN_WM_VAL)
          {
            continue ;  // only if diagonal voxel isn't on
          }
          alabel = left ? Left_Amygdala : Right_Amygdala ;
          hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
          if ((MRIgetVoxVal(mri_seg, x, y-1, z+1, 0) != alabel) && (MRIgetVoxVal(mri_seg, x, y-1, z+1, 0) != hlabel))
          {
            continue ;
          }
          // check to make sure we are at inferior border of hippo or amygdala
          if ((distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8) ||
              (distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))
          {
            continue ;
          }
          if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
          {
            continue ;
          }
          non++ ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
          nchanged++ ;
          break ;
        }
      }
    }
  }
  for (z = 1 ; z < depth-1 ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Cerebral_Cortex:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Cerebral_Cortex:
          // make sure we aren't on the medial edge of hippo
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
          {
            continue ;
          }
#else
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
#endif
          slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          // look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
          if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
          {
            continue ;
          }
#endif
          if (MRIvox(mri_wm, x, y-1, z-1) >= MIN_WM_VAL)
          {
            continue ;  // only if diagonal voxel isn't on
          }
          alabel = left ? Left_Amygdala : Right_Amygdala ;
          hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
          if ((MRIgetVoxVal(mri_seg, x, y-1, z-1, 0) != alabel) && (MRIgetVoxVal(mri_seg, x, y-1, z-1, 0) != hlabel))
          {
            continue ;
          }
          if ((MRIgetVoxVal(mri_seg, x, y-1, z-1, 0) == alabel) &&
              (anterior_edge_of_amygdala(mri_seg, x, y-1, z-1)))
          {
            continue ;
          }
          // check to make sure we are at inferior border of hippo or amygdala
          if ((distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8) ||
              (distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8))
          {
            continue ;
          }
          if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
          {
            continue ;
          }
          non++ ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
          nchanged++ ;
          break ;
        }
      }
    }
  }


  // look for unknown voxels with hippo diagonal connectivity.
  for (z = 1 ; z < depth-1 ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Unknown:
          if (MRIgetVoxVal(mri_seg, x-1, y-1, z, 0) == Left_Hippocampus)
          {
            xi = x-1 ;
            left = 1 ;
            hlabel = Left_Hippocampus ;
          }
          else if (MRIgetVoxVal(mri_seg, x+1, y-1, z, 0) != Right_Hippocampus)
          {
            continue ;
          }
          else
          {
            xi = x+1 ;
            hlabel = Right_Hippocampus ;
          }
          if (left)
          {
            i =  distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10)  ;
            if (i < 10 && MRIvox(mri_wm, x-i, y, z) < WM_MIN_VAL)
            {
              continue ;
            }
          }
          else
          {
            i = distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) ;
            if (i < 10 && MRIvox(mri_wm, x+i, y, z) < WM_MIN_VAL)
            {
              continue ;
            }
          }

          if ((distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))  // no hippo inferior
          {
            continue ;
          }
          if (distance_to_nonzero(mri_wm, x, y, z, 0, 1, 0, 8) < 7)
          {
            continue ;
          }
          non++ ;
          nchanged++ ;
          if (xi == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, xi, y-1, z) = AUTO_FILL ;
          break ;

        case Right_Cerebral_Cortex:
          // make sure we aren't on the medial edge of hippo
#if 0
          if (left)
          {
            if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
            {
              continue ;
            }
          }
          else if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
          {
            continue ;
          }
#else
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
#endif
          slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          // look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
          if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
          {
            continue ;
          }
#endif
          if (MRIvox(mri_wm, x, y-1, z+1) >= MIN_WM_VAL)
          {
            continue ;  // only if diagonal voxel isn't on
          }
          alabel = left ? Left_Amygdala : Right_Amygdala ;
          hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
          if ((MRIgetVoxVal(mri_seg, x, y-1, z+1, 0) != alabel) && (MRIgetVoxVal(mri_seg, x, y-1, z+1, 0) != hlabel))
          {
            continue ;
          }
          if ((MRIgetVoxVal(mri_seg, x, y-1, z+1, 0) == alabel) &&
              (anterior_edge_of_amygdala(mri_seg, x, y-1, z+1)))
          {
            continue ;
          }
          // check to make sure we are at inferior border of hippo or amygdala
          if ((distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8) ||
              (distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))
          {
            continue ;
          }
          if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
          {
            continue ;
          }
          non++ ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
          nchanged++ ;
          break ;
        }
      }
    }
  }

#if 0 // no longer needed due to path stuff
  // remove uneeded filled voxels by checking for ones that are hippo/amy/inf lat and have lots of wm inferior
  for (z = 1 ; z < depth ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_wm, x, y, z) != AUTO_FILL)
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Amygdala:
        case Left_Inf_Lat_Vent:
        case Left_Hippocampus:
          left = 1 ;
        case Right_Inf_Lat_Vent:
        case Right_Hippocampus:
        case Right_Amygdala:
        {
          int erase = 1 ;

          // if either of the 2 inferior planes have all 9 voxels on, erase this one
          for (yi = y+1 ; yi<=y+2 ; yi++)
          {
            erase = 1 ;
            for (xi = x-1 ; erase && xi <= x+1 ; xi++)
            {
              for (zi = z-1 ; erase && zi <= z+1 ; zi++)
              {
                if (MRIvox(mri_wm, xi, yi, zi) < MIN_WM_VAL)
                {
                  erase = 0 ;
                  break ;
                }
              }
            }
            if (erase)  // just need one of the two planes to be full
            {
              break ;
            }
          }
          if (erase)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x, y, z) = 0 ;
            noff++ ;
          }
        }
        break ;
        }
      }
    }
  }
#endif

  // add voxels that are wm in the aseg, but not in the wm vol
  for (z = 1 ; z < depth ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        int done ;

        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if ((MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL) ||
            (!IS_WHITE_CLASS(MRIgetVoxVal(mri_seg, x, y, z, 0))))
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Cerebral_White_Matter:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Cerebral_White_Matter:
          hlabel = left ?  Left_Hippocampus : Right_Hippocampus ;
          // if there is any hippo superior to this label, turn it on

          yi = y-1 ;
          done = 0 ;
          for (xi = x-1 ; !done && xi <= x+1 ; xi++)
          {
            for (zi = z-1 ; !done && zi <= z+1 ; zi++)
            {
              if (MRIgetVoxVal(mri_seg, xi, yi, zi, 0) == hlabel)
              {
                done = 1 ;
                if (x == Gx && y == Gy && z == Gz)
                {
                  DiagBreak2() ;
                }
                MRIvox(mri_wm, x,y,z) = AUTO_FILL ;
                non++ ;
                nchanged++ ;
                break ;
              }
            }
          }
          break ;
        }
      }
    }
  }
  // don't allow cortex to be directly lateral to hippo or inf lat vent - must be some wm
  for (z = 1 ; z < depth-1 ; z++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      for (x = 2 ; x < width-2 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        left = 0 ;
        switch (label)
        {
        case Left_Inf_Lat_Vent:
        case Left_Hippocampus:
        case Left_Amygdala:
          left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case Right_Hippocampus:
        case Right_Amygdala:
        case Right_Inf_Lat_Vent:
          if ((IS_AMYGDALA(label) || IS_INF_LAT_VENT(label)) &&
              (distance_to_superior_edge(mri_seg, label, x, y, z) < 3))
          {
            continue ;
          }
          if (medial_edge_of_hippocampus(mri_seg, x, y, z, left ? -1 : 1))
          {
            continue ;
          }
          if (anterior_edge_of_amygdala(mri_seg, x, y, z))
          {
            continue ;
          }
          olabel = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
          if (left)
          {
            xi = x+1 ;
          }
          else
          {
            xi = x-1 ;
          }
          if (MRIvox(mri_wm, xi, y, z) >= WM_MIN_VAL)
          {
            continue ;  // lateral voxel already on
          }
          if ((IS_CORTEX(MRIgetVoxVal(mri_seg, xi, y, z, 0)) == 0) &&
              (IS_UNKNOWN(MRIgetVoxVal(mri_seg, xi, y, z, 0)) == 0) &&
              (IS_WHITE_CLASS(MRIgetVoxVal(mri_seg, xi, y, z, 0)) == 0))
          {
            continue ;
          }

#if 0  // no longer needed with path stuff
          if (MRIvox(mri_T1, xi, y, z) > MRIvox(mri_T1, x, y,z))
          {
            if (xi == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, xi,y,z) = AUTO_FILL ;
            non++ ;
            break ;
          }
          else
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak2() ;
            }
            MRIvox(mri_wm, x,y,z) = AUTO_FILL ;
            non++ ;
            break ;
          }
#endif
          break ;
        }
      }

    }
  }
  printf("SEG EDIT: %d voxels turned on, %d voxels turned off.\n", non, noff) ;
  MRIfree(&mri_filled) ;
  return(NO_ERROR) ;
}

static int
neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label)
{
  int xi, yi, zi, xk, yk, zk ;

  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
#if 0
        if (abs(xk)+abs(yk)+abs(zk) > 1) /* only 6-connected neighbors */
        {
          continue ;
        }
#endif
        xi = mri->xi[x+xk] ;
        if (MRIgetVoxVal(mri, xi, yi, zi,0) == label)
        {
          return(1) ;
        }
      }
    }
  }
  return(0) ;
}

#if 0
static int
MRIlabelsInNbhd(MRI *mri, int x, int y, int z, int whalf, int label)
{
  int xi, yi, zi, xk, yk, zk, count ;

  for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIgetVoxVal(mri, xi, yi, zi, 0) == label)
        {
          count++;
        }
      }
    }
  }
  return(count) ;
}
static int
MRInonzeroInNbhd(MRI *mri, int x, int y, int z, int whalf)
{
  int xi, yi, zi, xk, yk, zk, count ;

  for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIgetVoxVal(mri, xi, yi, zi, 0) >= MIN_WM_VAL)
        {
          count++;
        }
      }
    }
  }
  return(count) ;
}
#endif
static int
distance_to_label(MRI *mri_labeled, int label, int x, int y, int z, int dx,
                  int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ;
    yi = y + d * dy ;
    zi = z + d * dz ;
    xi = mri_labeled->xi[xi] ;
    yi = mri_labeled->yi[yi] ;
    zi = mri_labeled->zi[zi];
    if (MRIgetVoxVal(mri_labeled, xi, yi, zi,0) == label)
    {
      break ;
    }
  }

  return(d) ;
}
#if SPACKLE_MTL
static int
distance_to_nonzero(MRI *mri_wm, int x, int y, int z, int dx,
                    int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ;
    yi = y + d * dy ;
    zi = z + d * dz ;
    xi = mri_wm->xi[xi] ;
    yi = mri_wm->yi[yi] ;
    zi = mri_wm->zi[zi];
    if (MRIgetVoxVal(mri_wm, xi, yi, zi,0) >= MIN_WM_VAL)
    {
      break ;
    }
  }

  return(d) ;
}
static int
distance_to_zero(MRI *mri_wm, int x, int y, int z, int dx,
                 int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ;
    yi = y + d * dy ;
    zi = z + d * dz ;
    xi = mri_wm->xi[xi] ;
    yi = mri_wm->yi[yi] ;
    zi = mri_wm->zi[zi];
    if (MRIgetVoxVal(mri_wm, xi, yi, zi,0) < MIN_WM_VAL)
    {
      break ;
    }
  }

  return(d) ;
}

static int
MRInonfilledInNbhd(MRI *mri, int x, int y, int z, int whalf)
{
  int xi, yi, zi, xk, yk, zk, count ;

  for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if ((MRIgetVoxVal(mri, xi, yi, zi, 0) >= MIN_WM_VAL) && (MRIgetVoxVal(mri, xi, yi, zi,0) != AUTO_FILL))
        {
          count++;
        }
      }
    }
  }
  return(count) ;
}
static int
medial_edge_of_hippocampus(MRI *mri_seg, int x, int y, int z, int dx)
{
  int xi, xk, label, vdc ;

  vdc = dx < 0 ? Left_VentralDC : Right_VentralDC ;

  /* if there are unknowns medial to us, and no hippocampus substantially
    medial then we are at medial edge.
  */
  if (x == Gx && y == Gy && z == Gz)
  {
    DiagBreak() ;
  }
  if ((distance_to_label(mri_seg, Unknown, x, y, z, dx, 0, 0, 10) < 10) ||
      (distance_to_label(mri_seg, Brain_Stem, x, y, z, dx, 0, 0, 10) < 10) ||
      (distance_to_label(mri_seg, vdc, x, y, z, dx, 0, 0, 10) < 10))
  {
    for (xk = 5 ; xk <= 10 ; xk++)
    {
      xi = mri_seg->xi[x+xk*dx] ;
      label = MRIgetVoxVal(mri_seg, xi, y, z, 0) ;
      if (IS_HIPPO(label) || IS_AMYGDALA(label))
      {
        return(0) ;
      }
    }
    return(1) ;
  }
  return(0) ;
}
#endif
static int
remove_paths_to_cortex(MRI *mri_wm, MRI *mri_T1, MRI *mri_aseg)
{
  int   x, y, z, xgm, ygm, zgm, nfilled, xmedial, path_to_cortex, xk, yk, zk, xi, yi, zi,
        niter, label, xb, yb, zb, zanterior, x1, y1, z1, x2, y2, z2, total_filled = 0 ;
  MRI   *mri_filled = NULL, *mri_hippo, *mri_roi ;
  float  intensity, brightest_intensity ;
  MRI_REGION box ;

  xgm = ygm = zgm = -1 ;  // to get rid of warning
  mri_hippo = MRIclone(mri_aseg, NULL) ;
  mri_roi = MRIclone(mri_aseg, NULL) ;
  MRIcopyLabel(mri_aseg, mri_roi, Left_Hippocampus) ;
  MRIcopyLabel(mri_aseg, mri_roi, Left_Amygdala) ;
  MRIcopyLabel(mri_aseg, mri_roi, Left_Inf_Lat_Vent) ;
  MRIbinarize(mri_roi, mri_roi, 1, 0, 128) ;
  for (niter = 0 ; niter < 5 ; niter++)
  {
    MRIdilate(mri_roi, mri_roi) ;  // build a mask and don't let things escape from that region
  }
  remove_anterior_and_superior_amygdala(mri_roi, mri_aseg) ;
  remove_lateral_and_anterior_hippocampus(mri_roi, mri_aseg, 1) ;
  remove_medial_voxels(mri_roi, mri_aseg, 1) ;
  remove_unknown_voxels(mri_roi, mri_aseg, 1) ;
  remove_gray_matter_voxels(mri_roi, mri_aseg) ;
  MRIboundingBox(mri_roi, 1, &box) ;
  x1 = MAX(box.x,1) ;
  x2 = MIN(box.x + box.dx-1, mri_aseg->width-2) ;
  y1 = MAX(box.y,1) ;
  y2 = MIN(box.y + box.dy-1, mri_aseg->height-2) ;
  z1 = MAX(box.z,1) ;
  z2 = MIN(box.z + box.dz-1, mri_aseg->depth-2) ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_roi, "lh_roi.mgz") ;
  }

  // find the most medial hippocampal voxel then erase everything within 1 cm of it

  // LH: now fill inferior and lateral until no voxels are filled or we reach cortex
  niter = 0 ;
  do
  {
    xmedial = mri_hippo->width ;
    ;
    zanterior = 0 ;
    MRIclear(mri_hippo) ;
    for (x = x1 ; x <= x2 ; x++)
      for (y = y1 ; y <= y2 ; y++)
        for (z = z1 ; z <= z2 ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_wm, x, y, z,0) > WM_MIN_VAL)  // on in wm
          {
            continue ;
          }
          label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
          switch (label)
          {
          default:
            break ;
          case Left_Hippocampus:
          case Left_Amygdala:
          case Left_Inf_Lat_Vent:
            yi = y+1 ;
            if (MRIgetVoxVal(mri_wm, x, yi, z,0) > WM_MIN_VAL)  // inf voxel is wm
            {
              MRIvox(mri_hippo, x, y, z) = 128 ;
              if (x < xmedial)
              {
                xmedial = x ;
              }
              if (z > zanterior)
              {
                zanterior = z ;
              }
            }
          }
        }
    for (x = x1 ; x <= x2 ; x++)
      for (y = y1 ; y <= y2 ; y++)
        for (z = z1 ; z <= z2 ; z++)
          if (MRIvox(mri_hippo, x, y, z) && ((x < xmedial+10) || (z > zanterior-15)))
          {
            MRIvox(mri_hippo, x, y, z) = 0 ;
          }
    if (niter == 0 && (Gdiag & DIAG_WRITE))
    {
      MRIwrite(mri_hippo, "f1.mgz") ;
    }
    mri_filled = MRIcopy(mri_hippo, mri_filled) ;
    do
    {
      path_to_cortex = nfilled = 0 ;
      for (x = xmedial ; x <= x2 && !path_to_cortex ; x++)
        for (y = y1 ; y <= y2 && !path_to_cortex ; y++)
          for (z = 0 ; z <= z2 && !path_to_cortex ; z++)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak() ;
            }
            if (MRIgetVoxVal(mri_filled, x, y, z,0) == 0)
            {
              continue ;
            }
            if (MRIgetVoxVal(mri_roi, x, y, z,0) == 0)
            {
              continue ;
            }
            for (xk = -1 ; xk <= 1 && !path_to_cortex ; xk++)
            {
              xi = mri_filled->xi[xk+x] ;
              for (yk = 0 ; yk <= 1 && !path_to_cortex ; yk++)
              {
                yi = mri_filled->yi[yk+y] ;
                for (zk = -1 ; zk <= 1 && !path_to_cortex  ; zk++)
                {
                  zi = mri_filled->zi[zk+z] ;
#if 0
                  if (zi > zanterior-15)
                  {
                    continue ;
                  }
#endif
                  if (xi == Gx && yi == Gy && zi == Gz)
                  {
                    DiagBreak() ;
                  }
                  if (MRIgetVoxVal(mri_roi, xi, yi, zi,0) == 0)
                  {
                    continue ;
                  }
                  if (MRIgetVoxVal(mri_filled, xi, yi, zi,0))
                  {
                    continue ;
                  }
                  if (MRIgetVoxVal(mri_wm, xi, yi, zi,0) > WM_MIN_VAL)  // wm barrier there
                  {
                    continue ;
                  }
                  if (IS_CORTEX(MRIgetVoxVal(mri_aseg, xi, yi, zi, 0))) // found a path to cortex
                  {
                    path_to_cortex = 1 ;
                    xgm = xi ;
                    ygm = yi ;
                    zgm = zi ;
                  }
                  else if (xi > xmedial+10)    // don't fill medially to avoid going "around bend"
                  {
                    label = MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
                    switch (label)
                    {
                    default:
                      break ;
                    case Unknown:
                    case Left_Cerebral_Cortex:
                    case Left_Hippocampus:
                    case Left_Amygdala:
                    case Left_Inf_Lat_Vent:
                    case Left_Cerebral_White_Matter:
                      if (MRIgetVoxVal(mri_roi, xi, yi, zi,0) > 0)
                      {
                        nfilled++ ;
                        MRIvox(mri_filled, xi, yi, zi) = 128 ;
                      }
                      break ;
                    }
                  }
                }
              }
            }
          }

    }
    while (nfilled > 0 && path_to_cortex == 0) ;
#define WHALF 3
    if (path_to_cortex > 0)  // fill in some of the voxels
    {
      // find brightest voxel in filling path, and turn it on in the wm vol
      if (xgm == Gx && ygm == Gy && zgm == Gz)
      {
        DiagBreak() ;
      }
      brightest_intensity = MRIgetVoxVal(mri_T1, xgm, ygm, zgm,0) ;
      xb = xgm ;
      yb = ygm ;
      zb = zgm ;
      for (xk = -WHALF ; xk <= WHALF ; xk++)
      {
        xi = mri_filled->xi[xk+xgm] ;
        for (yk = -WHALF ; yk <= 0 ; yk++)
        {
          yi = mri_filled->yi[yk+ygm] ;
          for (zk = -WHALF ; zk <= WHALF ; zk++)
          {
            zi = mri_filled->zi[zk+zgm] ;
            if (xi == Gx && yi == Gy && zi == Gz)
            {
              DiagBreak() ;
            }
            if (MRIgetVoxVal(mri_filled, xi, yi, zi,0) == 0)  // wasn't in the path
            {
              continue ;
            }
            if (MRIgetVoxVal(mri_wm, xi, yi, zi,0) > MIN_WM_VAL)  // already filled
            {
              continue ;
            }
            if (MRIneighborsOn(mri_wm, xi, yi, zi, WM_MIN_VAL) == 0)
            {
              continue ;
            }

            intensity = MRIgetVoxVal(mri_T1, xi, yi, zi,0) ;
            if (intensity > brightest_intensity)
            {
              brightest_intensity = intensity ;
              xb = xi ;
              yb = yi ;
              zb = zi ;
            }
          }
        }
      }
      if (xb == Gx && yb == Gy && zb == Gz)
      {
        DiagBreak2() ;
      }
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      {
        printf("filling voxel (%d, %d, %d)\n", xb, yb, zb) ;
      }
      MRIvox(mri_wm, xb, yb, zb) = AUTO_FILL ;
      MRIvox(mri_hippo, xb, yb, zb) = 0 ;
      total_filled++ ;
    }

  }
  while ((path_to_cortex > 0) && (niter++ < 100000)) ;




  mri_hippo = MRIclone(mri_aseg, NULL) ;
  mri_roi = MRIclone(mri_aseg, NULL) ;
  MRIcopyLabel(mri_aseg, mri_roi, Right_Hippocampus) ;
  MRIcopyLabel(mri_aseg, mri_roi, Right_Amygdala) ;
  MRIcopyLabel(mri_aseg, mri_roi, Right_Inf_Lat_Vent) ;
  MRIbinarize(mri_roi, mri_roi, 1, 0, 128) ;
  for (niter = 0 ; niter < 5 ; niter++)
  {
    MRIdilate(mri_roi, mri_roi) ;  // build a mask and don't let things escape from that region
  }
  remove_anterior_and_superior_amygdala(mri_roi, mri_aseg) ;
  remove_lateral_and_anterior_hippocampus(mri_roi, mri_aseg, 0) ;
  remove_medial_voxels(mri_roi, mri_aseg, 0) ;
  remove_unknown_voxels(mri_roi, mri_aseg, 0) ;
  remove_gray_matter_voxels(mri_roi, mri_aseg) ;
  MRIboundingBox(mri_roi, 1, &box) ;
  x1 = MAX(box.x,1) ;
  x2 = MIN(box.x + box.dx-1, mri_aseg->width-2) ;
  y1 = MAX(box.y,1) ;
  y2 = MIN(box.y + box.dy-1, mri_aseg->height-2) ;
  z1 = MAX(box.z,1) ;
  z2 = MIN(box.z + box.dz-1, mri_aseg->depth-2) ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_roi, "rh_roi.mgz") ;
  }

  // find the most medial hippocampal voxel then erase everything within 1 cm of it

  // RH: now fill inferior and lateral until no voxels are filled or we reach cortex
  niter = 0 ;
  do
  {
    xmedial = 0 ;
    zanterior = 0 ;
    MRIclear(mri_hippo) ;
    for (x = x1 ; x <= x2 ; x++)
      for (y = y1 ; y <= y2 ; y++)
        for (z = z1 ; z <= z2 ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_wm, x, y, z,0) > WM_MIN_VAL)  // on in wm
          {
            continue ;
          }
          label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
          switch (label)
          {
          default:
            break ;
          case Right_Hippocampus:
          case Right_Amygdala:
          case Right_Inf_Lat_Vent:
            yi = y+1 ;
            if (MRIgetVoxVal(mri_wm, x, yi, z,0) > WM_MIN_VAL)  // inf voxel is wm
            {
              MRIvox(mri_hippo, x, y, z) = 128 ;
              if (x > xmedial)
              {
                xmedial = x ;
              }
              if (z > zanterior)
              {
                zanterior = z ;
              }
            }
          }
        }
    for (x = x1 ; x <= x2 ; x++)
      for (y = y1 ; y <= y2 ; y++)
        for (z = z1 ; z <= z2 ; z++)
          if (MRIvox(mri_hippo, x, y, z) && ((x > xmedial-10) || (z > zanterior-15)))
          {
            MRIvox(mri_hippo, x, y, z) = 0 ;
          }
    if (niter == 0 && (Gdiag & DIAG_WRITE))
    {
      MRIwrite(mri_hippo, "f1.mgz") ;
    }
    mri_filled = MRIcopy(mri_hippo, mri_filled) ;
    do
    {
      path_to_cortex = nfilled = 0 ;
      for (x = x1 ; x <= xmedial && !path_to_cortex ; x++)
        for (y = y1 ; y <= y2 && !path_to_cortex ; y++)
          for (z = 0 ; z <= z2 && !path_to_cortex ; z++)
          {
            if (x == Gx && y == Gy && z == Gz)
            {
              DiagBreak() ;
            }
            if (MRIvox(mri_filled, x, y, z) == 0)
            {
              continue ;
            }
            if (MRIvox(mri_roi, x, y, z) == 0)
            {
              continue ;
            }
            for (xk = -1 ; xk <= 1 && !path_to_cortex ; xk++)
            {
              xi = mri_filled->xi[xk+x] ;
              for (yk = 0 ; yk <= 1 && !path_to_cortex ; yk++)
              {
                yi = mri_filled->yi[yk+y] ;
                for (zk = -1 ; zk <= 1 && !path_to_cortex  ; zk++)
                {
                  zi = mri_filled->zi[zk+z] ;
                  if (xi == Gx && yi == Gy && zi == Gz)
                  {
                    DiagBreak() ;
                  }
                  if (MRIvox(mri_filled, xi, yi, zi))
                  {
                    continue ;
                  }
                  if (MRIvox(mri_roi, xi, yi, zi) == 0)
                  {
                    continue ;
                  }
                  if (MRIgetVoxVal(mri_wm, xi, yi, zi,0) > WM_MIN_VAL)  // wm barrier there
                  {
                    continue ;
                  }
                  if (IS_CORTEX(MRIgetVoxVal(mri_aseg, xi, yi, zi, 0))) // found a path to cortex
                  {
                    path_to_cortex = 1 ;
                    xgm = xi ;
                    ygm = yi ;
                    zgm = zi ;
                  }
                  else if (xi < xmedial-10)    // don't fill medially to avoid going "around bend"
                  {
                    if (MRIgetVoxVal(mri_roi, xi, yi, zi,0) > 0)
                    {
                      nfilled++ ;
                      MRIvox(mri_filled, xi, yi, zi) = 128 ;
                    }
                  }
                }
              }
            }
          }

    }
    while (nfilled > 0 && path_to_cortex == 0) ;
#define WHALF 3
    if (path_to_cortex > 0)  // fill in some of the voxels
    {
      // find brightest voxel in filling path, and turn it on in the wm vol
      if (xgm == Gx && ygm == Gy && zgm == Gz)
      {
        DiagBreak() ;
      }
      brightest_intensity = MRIgetVoxVal(mri_T1, xgm, ygm, zgm,0) ;
      xb = xgm ;
      yb = ygm ;
      zb = zgm ;
      for (xk = -WHALF ; xk <= WHALF ; xk++)
      {
        xi = mri_filled->xi[xk+xgm] ;
        for (yk = -WHALF ; yk <= 0 ; yk++)
        {
          yi = mri_filled->yi[yk+ygm] ;
          for (zk = -WHALF ; zk <= WHALF ; zk++)
          {
            zi = mri_filled->zi[zk+zgm] ;
            if (xi == Gx && yi == Gy && zi == Gz)
            {
              DiagBreak() ;
            }
            if (MRIgetVoxVal(mri_filled, xi, yi, zi,0) == 0)  // wasn't in the path
            {
              continue ;
            }
            if (MRIgetVoxVal(mri_wm, xi, yi, zi,0) > MIN_WM_VAL)  // already filled
            {
              continue ;
            }
            if (MRIneighborsOn(mri_wm, xi, yi, zi, WM_MIN_VAL) == 0)
            {
              continue ;
            }

            intensity = MRIgetVoxVal(mri_T1, xi, yi, zi,0) ;
            if (intensity > brightest_intensity)
            {
              brightest_intensity = intensity ;
              xb = xi ;
              yb = yi ;
              zb = zi ;
            }
          }
        }
      }
      if (xb == Gx && yb == Gy && zb == Gz)
      {
        DiagBreak2() ;
      }
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      {
        printf("filling voxel (%d, %d, %d)\n", xb, yb, zb) ;
      }
      MRIvox(mri_wm, xb, yb, zb) = AUTO_FILL ;
      MRIvox(mri_hippo, xb, yb, zb) = 0 ;
      total_filled++ ;
    }

  }
  while ((path_to_cortex > 0) && (niter++ < 100000)) ;


  printf("%d voxels added to wm to prevent paths from MTL structures to cortex\n", total_filled) ;
  return(NO_ERROR) ;
}

static int
remove_anterior_and_superior_amygdala(MRI *mri_roi, MRI *mri_aseg)
{
  int  x, y, z, label, left ;

  for (x = 0 ; x < mri_roi->width ; x++)
  {
    for (y = 1 ; y < mri_roi->height ; y++)
    {
      for (z = 0 ; z < mri_roi->depth-1 ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (!IS_AMYGDALA(label))
        {
          continue ;
        }
        left = (label == Left_Amygdala) ;
#if 0
        if ((distance_to_label(mri_aseg, label, x, y, z, 0, 0, 1, 10) < 10) &&
            (distance_to_label(mri_aseg, label, x, y, z, 0, -1, 0, 10) < 10))
        {
          continue ;
        }
#else
        if ((distance_to_anterior_edge(mri_aseg, label, x, y, z) > 3) &&
            (distance_to_superior_edge(mri_aseg, label, x, y, z) > 2))
        {
          continue ;
        }
        if (distance_to_inferior_edge(mri_aseg, label, x, y, z) <= 3)
        {
          continue ;
        }
        if (distance_to_lateral_edge(mri_aseg, label, x, y, z, left) <= 3)
        {
          continue ;
        }
#endif
        MRIvox(mri_roi, x, y, z) = 0 ;
        MRIvox(mri_roi, x, y, z+1) = 0 ;
        MRIvox(mri_roi, x, y-1, z) = 0 ;
      }
    }
  }
  return(NO_ERROR) ;
}
static int
remove_lateral_and_anterior_hippocampus(MRI *mri_roi, MRI *mri_aseg, int left)
{
  int  x, y, z, label, xl ;

  for (x = 0 ; x < mri_roi->width ; x++)
  {
    // find lateral border
    if (left)
    {
      xl = 0 ;
      for (y = 1 ; y < mri_roi->height ; y++)
      {
        for (z = 0 ; z < mri_roi->depth-1 ; z++)
        {
          label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
          if ((label == Left_Hippocampus) && (x > xl))
          {
            xl = x ;
          }

        }
      }
    }
    else
    {
      xl =  mri_roi->width-1 ;
      for (y = 1 ; y < mri_roi->height ; y++)
      {
        for (z = 0 ; z < mri_roi->depth-1 ; z++)
        {
          label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
          if ((label == Right_Hippocampus) && (x < xl))
          {
            xl = x ;
          }

        }
      }
    }
    for (y = 1 ; y < mri_roi->height ; y++)
    {
      for (z = 0 ; z < mri_roi->depth-1 ; z++)
      {
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (!IS_HIPPO(label))
        {
          continue ;
        }
        if (left && x > xl-5)
        {
          continue ;
        }
        if ((left == 0) && (x < xl+5))
        {
          continue ;
        }
        if (distance_to_anterior_edge(mri_aseg, label, x, y, z) > 3)
        {
          continue ;
        }
        MRIvox(mri_roi, x, y, z) = 0 ;
        MRIvox(mri_roi, x, y, z+1) = 0 ;
        MRIvox(mri_roi, x, y-1, z) = 0 ;
      }
    }
  }
  return(NO_ERROR) ;
}
#if SPACKLE_MTL
static int
anterior_edge_of_amygdala(MRI *mri_aseg, int x, int y, int z)
{
  int label ;
  label = MRIgetVoxVal(mri_aseg,x,y,z, 0) ;
  if (IS_AMYGDALA(label) == 0)
  {
    return(0) ;
  }
  if (distance_to_label(mri_aseg, label, x, y, z, 0, 0, 1, 10) < 10)
  {
    return(0) ;
  }
  return(1);
}
static int
distance_to_superior_edge(MRI *mri_seg, int label, int x, int y, int z)
{
  int yi, dist ;

  for (dist = 0 ; dist < mri_seg->height ; dist++)
  {
    yi = mri_seg->yi[y-dist] ;
    if (distance_to_label(mri_seg, label, x, yi, z, 0, -1, 0, 10) >= 10)
    {
      return(dist) ;
    }
  }
  return(-1) ;
}
static int
distance_to_lateral_edge(MRI *mri_seg, int label, int x, int y, int z, int left)
{
  int xi, dist ;

  if (left)
  {
    for (dist = 0 ; dist < mri_seg->width ; dist++)
    {
      xi = mri_seg->xi[x+dist] ;
      if (distance_to_label(mri_seg, label, xi, y, z, 1, 0, 0, 10) >= 10)
      {
        return(dist) ;
      }
    }
  }
  else
  {
    for (dist = 0 ; dist < mri_seg->width ; dist++)
    {
      xi = mri_seg->xi[x-dist] ;
      if (distance_to_label(mri_seg, label, xi, y, z, -1, 0, 0, 10) >= 10)
      {
        return(dist) ;
      }
    }
  }
  return(dist) ;
}
static int
distance_to_inferior_edge(MRI *mri_seg, int label, int x, int y, int z)
{
  int yi, dist ;

  for (dist = 0 ; dist < mri_seg->height ; dist++)
  {
    yi = mri_seg->yi[y+dist] ;
    if (distance_to_label(mri_seg, label, x, yi, z, 0, 1, 0, 10) >= 10)
    {
      return(dist) ;
    }
  }
  return(-1) ;
}
#endif
static int
distance_to_anterior_edge(MRI *mri_seg, int label, int x, int y, int z)
{
  int zi, dist ;

  for (dist = 0 ; dist < mri_seg->height ; dist++)
  {
    zi = mri_seg->zi[z+dist] ;
    if (distance_to_label(mri_seg, label, x, y, zi, 0, 0, 1, 10) >= 10)
    {
      return(dist) ;
    }
  }
  return(-1) ;
}
static int
remove_medial_voxels(MRI *mri_roi, MRI *mri_aseg, int left)
{
  int x, y, z, xmid, xmin, xmax, label ;

  for (z = 0 ; z < mri_roi->depth ; z++)
  {
    // find medialmost hippo for this slice
    xmin = mri_roi->width-1 ;
    xmax = 0 ;
    for (x = 0 ; x < mri_roi->width ; x++)
    {
      for (y = 0 ; y < mri_roi->height ; y++)
      {
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if ((left && (label == Left_Hippocampus ||
                      label == Left_Amygdala)) ||
            ((left == 0) &&
             ((label == Right_Amygdala) || (label == Right_Hippocampus))))
        {
          if (x < xmin)
          {
            xmin = x ;
          }
          if (x > xmax)
          {
            xmax = x ;
          }

        }
      }
    }

    xmid = left ? (2*xmin+xmax)/3 : (2*xmax + xmin)/3 ;

    // now remove gray matter from ROI if it's medial
    for (x = 0 ; x < mri_roi->width ; x++)
    {
      for (y = 0 ; y < mri_roi->height ; y++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (left)
        {
          if (((label == Left_Amygdala) ||
               (label == Left_Hippocampus) ||
               (label == Left_Cerebral_Cortex)) &&
              (x < xmid))
          {
            MRIvox(mri_roi, x, y, z) = 0 ;
          }
          if ((label == Unknown) && x < xmin+3)
          {
            MRIvox(mri_roi, x, y, z) = 0 ;
          }
        }
        else   // rh
        {
          if (((label == Right_Amygdala) ||
               (label == Right_Hippocampus) ||
               (label == Right_Cerebral_Cortex)) &&
              (x > xmid))
          {
            MRIvox(mri_roi, x, y, z) = 0 ;
          }
          if ((label == Unknown) && x > xmax-3)
          {
            MRIvox(mri_roi, x, y, z) = 0 ;
          }
        }
      }
    }
  }
  return(NO_ERROR) ;
}

static int
remove_unknown_voxels(MRI *mri_roi, MRI *mri_aseg, int left)
{
  int x, y, z, label, xmin, xmax, xmid ;

  for (z = 0 ; z < mri_roi->depth ; z++)
  {
    xmin = mri_roi->width-1 ;
    xmax = 0 ;
    for (x = 0 ; x < mri_roi->width ; x++)
    {
      for (y = 0 ; y < mri_roi->height ; y++)
      {
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if ((label == (left ? Left_Hippocampus : Right_Hippocampus)) ||
            (label == (left ? Left_Amygdala : Right_Amygdala)))
        {
          if (x < xmin)
          {
            xmin = x ;
          }
          if (x > xmax)
          {
            xmax = x ;
          }
        }
      }
    }
    // only do lateral half of hippo/amy
    xmid = (xmin + xmax) / 2 ;
    if (left)
    {
      xmax = xmid ;
    }
    else
    {
      xmin = xmid ;
    }
    for (x = xmin ; x <= xmax ; x++)
    {
      for (y = 0 ; y < mri_roi->height ; y++)
      {
        if (x == Gx && y == Gy &&  z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_roi, x, y, z) == 0)
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (label != Unknown)
        {
          continue ;
        }
        if (distance_to_nonzero(mri_aseg, x, y, z, 0, 1, 0, 6) >= 6)
        {
          MRIvox(mri_roi, x, y, z) = 0 ;
        }
      }
    }
  }

  return(0) ;
}

static int
remove_gray_matter_voxels(MRI *mri_roi, MRI *mri_aseg)
{
  int x, y, z, label, alabel ;

  for (z = 0 ; z < mri_roi->depth ; z++)
  {
    for (x = 0 ; x < mri_roi->width ; x++)
    {
      for (y = 0 ; y < mri_roi->height ; y++)
      {
        if (x == Gx && y == Gy &&  z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_roi, x, y, z) == 0)
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (IS_CORTEX(label) == 0)
        {
          continue ;
        }
        alabel = (label == Left_Cerebral_Cortex ? Left_Amygdala : Right_Amygdala) ;

        // remove gray matter voxels that are at anterior border of amygdala
        if ((distance_to_label(mri_aseg, alabel, x, y, z, 0, 0, -1, 10) <= 2) &&
            (distance_to_label(mri_aseg, alabel, x, y, z, 0, 0, 1, 10) >= 10))
        {
          MRIvox(mri_roi, x, y, z) = 0 ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

static int
spackle_wm_superior_to_mtl(MRI *mri_wm, MRI *mri_T1, MRI *mri_aseg)
{
  int    x, y, z, label, slabel, yi, xi, left ;

  for (x = 0 ; x < mri_wm->width ; x++)
  {
    for (y = 1 ; y < mri_wm->height; y++)
    {
      for (z = 0 ; z < mri_wm->depth ; z++)
      {
        if (x == Gx && y  == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (!IS_HIPPO(label) && !IS_AMYGDALA(label))
        {
          continue ;
        }
        left = (label == Left_Amygdala || label == Left_Hippocampus) ;
        if (distance_to_superior_edge(mri_aseg, label, x, y, z) > 1)
        {
          continue ;
        }
        if ((distance_to_anterior_edge(mri_aseg, label, x, y, z) < 2) &&
            (distance_to_anterior_edge(mri_aseg, label, x, y+1, z) < 2))
        {
          continue ;
        }
        yi = mri_aseg->yi[y-1] ;
        slabel = MRIgetVoxVal(mri_aseg, x, yi, z, 0) ;
        if (IS_CORTEX(slabel) && (MRIgetVoxVal(mri_wm, x, yi, z,0) < MIN_WM_VAL))
        {
          if (x == Gx && yi == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
        }
        xi = left ? x+1 : x-1 ;
        if (IS_AMYGDALA(label) && IS_CORTEX(MRIgetVoxVal(mri_aseg, xi, y, z, 0)) &&
            (MRIvox(mri_wm, xi, y, z) < MIN_WM_VAL))
        {
          if (xi == Gx && y == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
        }

        if (IS_AMYGDALA(label) && IS_CORTEX(MRIgetVoxVal(mri_aseg, xi, yi, z, 0)) &&
            (MRIvox(mri_wm, xi, yi, z) < MIN_WM_VAL))
        {
          if (xi == Gx && yi == Gy && z == Gz)
          {
            DiagBreak2() ;
          }
          MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
        }
      }
    }
  }

  return(NO_ERROR) ;
}
