/**
 * @brief compute the % of gm layers 1-6, wm and CSF in each voxel in a volume
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021
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
#include <ctype.h>

#include "mri.h"
#include "mri2.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "mrisurf.h"
#include "registerio.h"
#include "cma.h"

#define WM_VAL         1
#define CSF_VAL        (nlayers+2)
#define SUBCORT_GM_VAL (nlayers+3)
#define VENTRICLE_VAL  (nlayers+4)

#define NLAYERS        6
#define NLABELS        (nlayers+3)  // wm + cortical layers + csf + subcortical gray

static int nlayers = NLAYERS ;
static const char *LAMINAR_NAME = "gwdist";
static const char *aseg_name = "aseg.mgz" ;

#define LEFT_HEMI  0
#define RIGHT_HEMI 1
#define BOTH_HEMIS 2

static int noaseg = 0 ;
static char *subject_name = NULL ;
static char *hemi = (char *)"lh" ;
static int hemi_no = LEFT_HEMI ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int cortex_only = 1 ;
static char sdir[STRLEN] = "" ;
static double resolution = .5 ;

static int FS_names = 0 ;

static int create_synth = 0 ;
static char *synth_fname = NULL ;
void MRIeraseOtherHemi(MRI *mri_dst, char *hemi) ;
MRI *add_aseg_structures_outside_ribbon(MRI *mri_src, MRI *mri_aseg, 
                                        MRI *mri_dst,
					int wm_val, int gm_val, int csf_val, int vent_val) ;
MRI *create_synth_vol(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst, int hemi_no) ;
int MRIcomputePartialVolumeFractions(MRI *mri_src, MATRIX *m_vox2vox, 
                                     MRI *mri_seg, MRI *mri_fractions) ;
int
main(int argc, char *argv[]) 
{
  char   **av, fname[STRLEN] ;
  int    ac, nargs, i, hno, hmin, hmax ;
  char   *subject, *reg_fname, *in_fname, *out_fname, *cp ;
  int    msec, minutes, seconds, nvox, float2int, layer, width, height, depth ;
  Timer start ;
  MRI_SURFACE *mris;
  MRI         *mri_aseg, *mri_layers, *mri_tmp, *mri_in, *mri_ohemi_layers = NULL,
    *mri_interior_bottom, *mri_interior_top, *mri_fractions ;
  MATRIX      *m_regdat ;
  float       intensity, betplaneres, inplaneres ;

  nargs = handleVersionOption(argc, argv, "mri_compute_layer_fractions");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;
  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  reg_fname = argv[1] ; in_fname = argv[2] ;
  out_fname = argv[3] ; Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  printf("reading registration file %s\n", reg_fname) ;
  if (stricmp(reg_fname, "identity.nofile") == 0)
  {
    printf("using identity transform\n") ;
    m_regdat = NULL ;
    inplaneres = betplaneres = intensity = 1 ;
    float2int = 0 ;
    subject = const_cast<char*>("unknown"); // Not nice....
  }
  else
  {
    regio_read_register(reg_fname, &subject, &inplaneres,
                        &betplaneres, &intensity,  &m_regdat,
                        &float2int);

    m_regdat = regio_read_registermat(reg_fname) ;
    if (m_regdat == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load registration file from %s", Progname,reg_fname) ;
  }

  if (subject_name)
    subject = subject_name ;   // specified on command line
  int req = snprintf(fname, STRLEN,
		     "%s/%s/mri/%s", sdir, subject, aseg_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  printf("reading volume %s\n", fname) ;
  mri_aseg = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load aseg volume from %s", Progname,fname) ;

  if (mri_aseg->ct)
  {
    char fname[STRLEN], *fs ;
    COLOR_TABLE *ctab_default ;
    int         *translations, src_label, trg_label, x, y, z ;

    printf("embedded color table found in aseg - translating to default to ensure compatibility\n") ;
    fs = getenv("FREESURFER_HOME") ;
    if (fs == NULL)
      ErrorExit(ERROR_BADPARM, "FREESURFER_HOME must be defined in the environment for LUT translation\n") ;
    sprintf(fname, "%s/FreeSurferColorLUT.txt", fs) ;
    ctab_default = CTABreadASCII(fname) ;
    if (ctab_default == NULL)
      ErrorExit(ERROR_BADPARM, "could not read default lut from %s\n", fname) ;
    
    translations = (int *)calloc(mri_aseg->ct->nentries, sizeof(int)) ;
    if (translations == NULL)
      ErrorExit(ERROR_BADPARM, "could not read allocated %d-len translation table\n", mri_aseg->ct->nentries) ;

    for (src_label = 0 ; src_label < mri_aseg->ct->nentries ; src_label++)
    {
      if (mri_aseg->ct->entries[src_label] == NULL)
	continue ;
      trg_label = CTABentryNameToIndex(mri_aseg->ct->entries[src_label]->name, ctab_default);
      if (trg_label < 0)
	printf("src_label %s (%d) cannot be found in default lut\n", mri_aseg->ct->entries[src_label]->name, src_label) ;
      else
	translations[src_label] = trg_label ;
    }

    for (x = 0 ; x < mri_aseg->width ; x++)
      for (y = 0 ; y < mri_aseg->height ; y++)
	for (z = 0 ; z < mri_aseg->depth ; z++)
	{
	  src_label = (int)MRIgetVoxVal(mri_aseg, x, y,z, 0) ;
	  trg_label = translations[src_label] ;
	  MRIsetVoxVal(mri_aseg, x, y, z , 0, trg_label) ;
	}
    mri_aseg->ct = ctab_default ;
    if (Gdiag & DIAG_WRITE)
    {
      printf("writing translated aseg to at.mgz") ;
      MRIwrite(mri_aseg, "at.mgz") ;
    }
    free(translations) ;
  }

  nvox = (int)ceil(256/resolution);  
  mri_in = MRIreadHeader(in_fname, MRI_VOLUME_TYPE_UNKNOWN) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load input volume from %s", Progname,in_fname) ;
  width = (int)ceil(mri_in->width * (mri_in->xsize/resolution));  
  height = (int)ceil(mri_in->height * (mri_in->ysize/resolution));  
  depth = (int)ceil(mri_in->depth * (mri_in->zsize/resolution));  
  mri_layers = MRIalloc(width, height, depth, MRI_UCHAR) ;
  MRIsetResolution(mri_layers, resolution, resolution, resolution) ;
  mri_layers->xstart = mri_in->xstart ; mri_layers->xend = mri_in->xend ;
  mri_layers->ystart = mri_in->ystart ; mri_layers->yend = mri_in->yend ;
  mri_layers->zstart = mri_in->zstart ; mri_layers->zend = mri_in->zend ;
  mri_layers->x_r = mri_in->x_r ; mri_layers->x_a = mri_in->x_a ; mri_layers->x_s = mri_in->x_s ;
  mri_layers->y_r = mri_in->y_r ; mri_layers->y_a = mri_in->y_a ; mri_layers->y_s = mri_in->y_s ;
  mri_layers->z_r = mri_in->z_r ; mri_layers->z_a = mri_in->z_a ; mri_layers->z_s = mri_in->z_s ;
  mri_layers->c_r = mri_in->c_r ; mri_layers->c_a = mri_in->c_a ; mri_layers->c_s = mri_in->c_s ;
  MRIfree(&mri_in) ;

  MRIreInitCache(mri_layers) ; 

  if (FS_names && nlayers > 2)
    ErrorExit(ERROR_UNSUPPORTED, "%s: if specifying FS_names must use -nlayers 1", Progname) ;
  printf("reading laminar surfaces from %s.?\n", LAMINAR_NAME) ;

  if (!strcmp(hemi, "lh"))
  {
    hmin = 0 ;
    hmax = 0 ;
  }
  else if (!strcmp(hemi, "rh"))
  {
    hmin = 1 ;
    hmax = 1 ;
  }
  else   // both
  {
    hmin = 0 ;
    hmax = 1 ;
  }
  for (hno = hmin ; hno <= hmax ; hno++)
  {
    hemi = (char *)(hno == 0 ? "lh" : "rh") ;
    for (i = 0 ; i <= nlayers ; i++) {
      if (FS_names && nlayers == 1) {
	if (i == 0) {
	  int req = snprintf(fname, STRLEN,
			     "%s/%s/surf/%s.white", sdir, subject,hemi) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	} else {
	  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.pial", sdir, subject,hemi) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	}
      } else {
	int req ; 
	if (FS_names && i == 0)
	  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.white", sdir, subject,hemi) ;
	else if (FS_names && i == 2)
	  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.pial", sdir, subject,hemi) ;
	else if (FS_names)
	  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, subject,hemi,LAMINAR_NAME) ;
	else
	  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s.%d", sdir, subject,hemi,LAMINAR_NAME,i) ;
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
	if (FileExists(fname) == 0) {
	  int req = snprintf(fname, STRLEN,
			     "%s/%s/surf/%s.%s%3.3d", sdir, subject,hemi,LAMINAR_NAME,i) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	}
      }
      printf("reading surface %s\n", fname) ;
      mris = MRISread(fname) ;
      if (mris == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not load %s surface %d from %s", 
		  Progname,hemi, i, fname) ;
      
      mri_interior_top = MRIclone(mri_layers, NULL) ;
      MRISfillInterior(mris, resolution, mri_interior_top) ;
      
      if (Gdiag & DIAG_WRITE) {
	sprintf(fname, "%s.top%d.mgz", hemi, i) ;
	printf("writing layer %d interior to %s\n", i, fname) ;
	MRIwrite(mri_interior_top, fname) ;
      }
      if (i == 0)  { // fill white matter
	mri_tmp = MRIclone(mri_interior_top, NULL) ;
	MRIreplaceValuesOnly(mri_interior_top, mri_tmp, 1, WM_VAL) ;
	MRIcopyLabel(mri_tmp, mri_layers, WM_VAL) ;
	MRIfree(&mri_tmp) ;
      } else  { // fill cortical layer
	mri_tmp = MRInot(mri_interior_bottom, NULL) ;
	MRIfree(&mri_interior_bottom) ;
	MRIand(mri_interior_top, mri_tmp, mri_tmp, 1) ;
	layer = nlayers-(i-1) ;
	layer = i+1 ;
	MRIreplaceValuesOnly(mri_tmp, mri_tmp, 1, layer) ; // layer # + 1
	MRIcopyLabel(mri_tmp, mri_layers, layer) ;
	MRIfree(&mri_tmp) ;
      }
      if (Gdiag & DIAG_WRITE)
      {
	sprintf(fname, "%s.layer%d.mgz", hemi, i) ;
	printf("writing layer %d to %s\n", i, fname) ;
	MRIwrite(mri_layers, fname) ;
      }
      mri_interior_bottom = mri_interior_top ;
    }
    
    MRIfree(&mri_interior_bottom) ;
    MRIreplaceValuesOnly(mri_layers, mri_layers, 0, CSF_VAL) ;
    if (hno == 0 && hmax == 1)  // save lh when processing rh
    {
      mri_ohemi_layers = MRIcloneDifferentType(mri_layers, MRI_INT) ;
      MRIcopy(mri_layers, mri_ohemi_layers) ;
      MRIsetValues(mri_layers, 0) ;
    }
  }
  if (noaseg == 0)
  {
    if (create_synth)
    {
      MRI *mri_synth, *mri_bg_lh, *mri_bg_rh ;

      if (nlayers != 2)
	ErrorExit(ERROR_UNSUPPORTED, "%s: synth only available for nlayers = 2 (%d specified)", Progname, nlayers);
      mri_synth = MRIcloneDifferentType(mri_layers, MRI_INT) ; // layers if uchar and won't fit ctx labels
      MRIcopy(mri_layers, mri_synth) ;
      if (hemi_no == LEFT_HEMI)  // lh
      {
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL, Left_Cerebral_White_Matter) ;
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL+1, ctx_lh_infragranular) ;
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL+2, ctx_lh_supragranular) ;
      }
      else if (hemi_no == RIGHT_HEMI)  // rh
      {
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL, Right_Cerebral_White_Matter) ;
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL+1, ctx_rh_infragranular) ;
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL+2, ctx_rh_supragranular) ;
      }
      else if (hemi_no == BOTH_HEMIS)  // both hemis
      {
	// replace numbers with aseg labels
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL, Right_Cerebral_White_Matter) ;
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL+1, ctx_rh_infragranular) ;
	MRIreplaceValuesOnly(mri_layers, mri_synth, WM_VAL+2, ctx_rh_supragranular) ;
	MRIreplaceValuesOnly(mri_ohemi_layers, mri_ohemi_layers, WM_VAL+1, ctx_lh_infragranular) ;
	MRIreplaceValuesOnly(mri_ohemi_layers, mri_ohemi_layers, WM_VAL+2, ctx_lh_supragranular) ;
	MRIreplaceValuesOnly(mri_ohemi_layers, mri_ohemi_layers, WM_VAL, Left_Cerebral_White_Matter) ;

	// copy the lh into the layers vol to create a single one with both hemis
	if (0)
	{
	  MRIcopyLabel(mri_ohemi_layers, mri_synth, Left_Cerebral_White_Matter) ;
	  MRIcopyLabel(mri_ohemi_layers, mri_synth, ctx_lh_infragranular) ;
	  MRIcopyLabel(mri_ohemi_layers, mri_synth, ctx_lh_supragranular) ;
	}
	
	mri_bg_lh = MRIcopy(mri_ohemi_layers, NULL) ;
	mri_bg_rh = MRIcopy(mri_synth, NULL) ;
	MRIreplaceValuesOnly(mri_synth, mri_synth, CSF_VAL, 0) ;
	MRIreplaceValuesOnly(mri_ohemi_layers, mri_ohemi_layers, CSF_VAL, 0) ;
	
	MRImax(mri_synth, mri_ohemi_layers, mri_synth) ;
//	MRIcopyLabel(mri_bg_lh, mri_synth, CSF_VAL) ;
//	MRIcopyLabel(mri_bg_rh, mri_synth, CSF_VAL) ;
	MRIreplaceValuesOnly(mri_synth, mri_synth, 0, CSF_VAL) ;
	
	MRIfree(&mri_bg_lh) ; MRIfree(&mri_bg_rh) ;
      }
      create_synth_vol(mri_synth, mri_aseg, mri_synth, hemi_no);
      printf("writing synth volume to %s\n", synth_fname) ;
      MRIwrite(mri_synth, synth_fname) ;
      MRIfree(&mri_synth) ;
    }

    add_aseg_structures_outside_ribbon(mri_layers, mri_aseg, mri_layers, 
				       WM_VAL, SUBCORT_GM_VAL, CSF_VAL, VENTRICLE_VAL) ;
  }
  printf("reading movable volume %s\n", in_fname) ;
  mri_in = MRIread(in_fname) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load input volume from %s", Progname,in_fname) ;
  
  {
    MATRIX *m_conformed_to_epi_vox2vox, *m_seg_to_conformed_vox2vox,
      *m_seg_to_epi_vox2vox ;

    if (m_regdat == NULL)    // assume identity transform
      m_seg_to_epi_vox2vox = MRIgetVoxelToVoxelXform(mri_layers, mri_in) ;
    else    // a register.dat was specified between mri_in and the aseg/surface space
    {
      m_conformed_to_epi_vox2vox = 
	MRIvoxToVoxFromTkRegMtx(mri_in, mri_aseg, m_regdat);
      m_seg_to_conformed_vox2vox = 
	MRIgetVoxelToVoxelXform(mri_layers, mri_aseg) ;
      
      m_seg_to_epi_vox2vox = 
	MatrixMultiply(m_conformed_to_epi_vox2vox,m_seg_to_conformed_vox2vox,NULL);
      MatrixFree(&m_regdat) ; MatrixFree(&m_conformed_to_epi_vox2vox) ; 
      MatrixFree(&m_seg_to_conformed_vox2vox);
    }
    
    printf("seg to EPI vox2vox matrix:\n") ;
    MatrixPrint(Gstdout, m_seg_to_epi_vox2vox) ;
    mri_fractions = 
      MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth,
                       MRI_FLOAT, NLABELS) ;
    MRIcopyHeader(mri_in, mri_fractions) ;
    MRIfree(&mri_aseg) ;
    if (Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s.layers.mgz", hemi) ;
      printf("writing layers to %s\n", fname) ;
      MRIwrite(mri_layers, fname) ;
    }
    printf("computing partial volume fractions...\n") ;
    MRIcomputePartialVolumeFractions(mri_in, m_seg_to_epi_vox2vox, mri_layers, 
                                     mri_fractions) ;
    MatrixFree(&m_seg_to_epi_vox2vox) ;
  }

  sprintf(fname, "%s.layers.mgz", out_fname) ;
  printf("writing layer labeling to %s\n", out_fname) ;
  MRIwrite(mri_fractions, out_fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("volume fraction calculation took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;

  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  } else if (!stricmp(option, "cortex")) {
    printf("limitting gm val to cortex\n") ;
    cortex_only = 1 ;
  } else if (!stricmp(option, "FS_names")) {
    printf("using standard FS names white and pial\n") ;
    FS_names = 1 ;
  } else if (!stricmp(option, "synth")) {
    create_synth = 1 ;
    synth_fname = argv[2] ;
    printf("creating synth volume instead of laminar partial volume fractions and writing to %s\n", synth_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "noaseg")) {
    noaseg = 1 ;
    printf("not computing subcortical components\n") ;
  } else if (!stricmp(option, "nlayers")) {
    nlayers = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d input layers for laminar analysis\n", nlayers) ;
  } else if (!stricmp(option, "rh") || !stricmp(option, "lh") || !stricmp(option, "both")) {
    hemi = option ;
    hemi_no = !stricmp(hemi, "lh") ? LEFT_HEMI : (!stricmp(hemi, "rh") ? RIGHT_HEMI : BOTH_HEMIS);
    printf("processing %s hemisphere (%d)\n", hemi, hemi_no) ;
  } else switch (toupper(*option)) {
  case 'W':
    Gdiag |= DIAG_WRITE ;
    printf("turning on write diagnostics\n") ;
    break ;
  case 'N':
    LAMINAR_NAME = argv[2] ;
    printf("using %s as layer surface name\n", LAMINAR_NAME) ;
    nargs = 1 ;
    break ;
  case 'S':
    subject_name = argv[2] ;
    nargs = 1 ;
    printf("overriding subject name in .dat file with %s\n", subject_name) ;
    break ;
  case 'A':
    aseg_name = argv[2] ;
    nargs = 1 ;
    printf("using aseg named %s\n", aseg_name) ;
    break ;
  case 'R':
    resolution = atof(argv[2]) ;
    printf("setting resolution = %2.3f\n", resolution) ;
    nargs = 1 ;
    break ;
  case '?':
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
static void
usage_exit(int code) {
  printf("usage: %s [options] <reg file> <input volume> <output stem>\n",
         Progname) ;
  printf("  -SDIR SUBJECTS_DIR \n");
  printf(
         "\t\n");
  exit(code) ;
}

int
MRIcomputePartialVolumeFractions(MRI *mri_src, MATRIX *m_vox2vox, 
                                 MRI *mri_seg, MRI *mri_fractions)
{
  int    x, y, z, xs, ys, zs, label ;
  VECTOR *v1, *v2 ;
  MRI    *mri_counts ;
  float  val, count ;
  MATRIX *m_inv ;

  m_inv = MatrixInverse(m_vox2vox, NULL) ;
  if (m_inv == NULL)
  {
    MatrixPrint(stdout, m_vox2vox) ;
    ErrorExit(ERROR_BADPARM, "MRIcomputePartialVolumeFractions: non-invertible vox2vox matrix");
  }
  mri_counts = MRIcloneDifferentType(mri_src, MRI_INT) ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  for (x = 0 ; x < mri_seg->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_seg->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_seg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xs = nint(V3_X(v2)) ; ys = nint(V3_Y(v2)) ; zs = nint(V3_Z(v2)) ;
        if (xs == Gx && ys == Gy && zs == Gz)
          DiagBreak() ;
        if (xs >= 0 && ys >= 0 && zs >= 0 &&
            xs < mri_src->width && ys < mri_src->height && zs < mri_src->depth)
        {
          val = MRIgetVoxVal(mri_counts, xs, ys, zs, 0) ;
          if (val > 0)
            DiagBreak() ;
          MRIsetVoxVal(mri_counts, xs, ys, zs, 0, val+1) ;

          label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
          if (label >= 0 && label <= mri_fractions->nframes)
          {
            val = MRIgetVoxVal(mri_fractions, xs, ys, zs, label-1) ;  // labels are frame+1
            MRIsetVoxVal(mri_fractions, xs, ys, zs, label-1, val+1) ;
          }
          else
            DiagBreak() ;
        }
      }
    }
  }

  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        count = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
        if (count >= 1)
        {
          for (label = 0 ; label < mri_fractions->nframes ; label++)
          {
            if (x == Gx && y == Gy && z == Gz)
              DiagBreak() ;
            val = MRIgetVoxVal(mri_fractions, x, y, z, label) ;
            MRIsetVoxVal(mri_fractions, x, y, z, label, val/count) ;
          }
        }
        else  // sample in other direction
        {
          V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
          MatrixMultiply(m_inv, v1, v2) ;
          MatrixMultiply(m_inv, v1, v2) ;
          xs = nint(V3_X(v2)) ; ys = nint(V3_Y(v2)) ; zs = nint(V3_Z(v2)) ;
          if (xs >= 0 && ys >= 0 && zs >= 0 &&
              xs < mri_seg->width && ys < mri_seg->height && zs < mri_seg->depth)
          {
            label = MRIgetVoxVal(mri_seg, xs, ys, zs, 0) ;
            if (label >= 0 && label < mri_fractions->nframes)
              MRIsetVoxVal(mri_fractions, x, y, z, label-1, 1) ;
            else
              DiagBreak() ;
          }
        }
      }
  VectorFree(&v1) ; VectorFree(&v2) ;
  MRIfree(&mri_counts) ; MatrixFree(&m_inv) ;
  
  return(NO_ERROR) ;
}
MRI *
add_aseg_structures_outside_ribbon(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst,
                                   int wm_val, int gm_val, int csf_val, int ventricle_val)
{
  VECTOR *v1, *v2 ;
  MATRIX *m_vox2vox ;
  int    x, y, z, xa, ya, za, label, seg_label ;

  if (mri_dst == NULL)
    mri_dst = MRIcopy(mri_src, NULL) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_src, mri_aseg) ;


  for (x = 0 ; x < mri_dst->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_dst->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        seg_label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0)) ;
        if (seg_label > wm_val && seg_label < csf_val)  // already labeled as cortex 
          continue ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xa = (int)(nint(V3_X(v2))) ;
        ya = (int)(nint(V3_Y(v2))) ;
        za = (int)(nint(V3_Z(v2))) ;
        if (xa < 0 || ya < 0 || za < 0 ||
            xa >= mri_aseg->width || ya >= mri_aseg->height || za >= mri_aseg->depth)
          continue ;
        if (xa == Gx && ya == Gy && za == Gz)
          DiagBreak() ;
        label = nint(MRIgetVoxVal(mri_aseg, xa, ya, za, 0)) ;
        switch (label)
        {
        case Left_Cerebellum_White_Matter:
        case Right_Cerebellum_White_Matter:
        case Brain_Stem:
          MRIsetVoxVal(mri_dst, x, y, z, 0, wm_val) ;
          break ;
        case Left_Cerebellum_Cortex:
        case Right_Cerebellum_Cortex:
	  if (seg_label != wm_val)  // otherwise it is interior to the ribbon
	    MRIsetVoxVal(mri_dst, x, y, z, 0, gm_val) ;
          break ;
        case Left_Pallidum:
        case Right_Pallidum:
        case Left_Thalamus:
        case Right_Thalamus:
        case Right_Putamen:
        case Left_Putamen:
        case Right_Caudate:
        case Left_Caudate:
        case Left_Accumbens_area:
        case Right_Accumbens_area:  // remove them from cortex
          MRIsetVoxVal(mri_dst, x, y, z, 0, gm_val) ;
          break ;
        case Left_Inf_Lat_Vent:
        case Right_Inf_Lat_Vent:
        case Left_Lateral_Ventricle:
        case Right_Lateral_Ventricle:
	case Third_Ventricle:
	case Fourth_Ventricle:
	case CSF:
          MRIsetVoxVal(mri_dst, x, y, z, 0, ventricle_val) ;
          break ;
        default:
          break ;
        }
      }
    }
  }

  VectorFree(&v1) ; VectorFree(&v2) ; MatrixFree(&m_vox2vox) ;
  return(mri_dst) ;
}

MRI *
create_synth_vol(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst, int hemi_no)
{
  VECTOR *v1, *v2 ;
  MATRIX *m_vox2vox ;
  int    x, y, z, xa, ya, za, aseg_label, seg_label, in_mtl, dilations, dno ;
  char   *hemi = (char *)(hemi_no == LEFT_HEMI ? "lh" : (hemi_no == RIGHT_HEMI ? "rh" : "both"));
  MRI    *mri_mtl ;

  mri_mtl = MRIclone(mri_aseg, NULL) ;
  if (hemi_no == LEFT_HEMI || hemi_no == BOTH_HEMIS)
  {
    MRIcopyLabel(mri_aseg, mri_mtl, Left_Hippocampus) ;
    MRIcopyLabel(mri_aseg, mri_mtl, Left_Amygdala) ;
    MRIcopyLabel(mri_aseg, mri_mtl, Left_Inf_Lat_Vent) ;
  }
  if (hemi_no == RIGHT_HEMI || hemi_no == BOTH_HEMIS)
  {
    MRIcopyLabel(mri_aseg, mri_mtl, Right_Hippocampus) ;
    MRIcopyLabel(mri_aseg, mri_mtl, Right_Amygdala) ;
    MRIcopyLabel(mri_aseg, mri_mtl, Right_Inf_Lat_Vent) ;
  }
  MRIbinarize(mri_mtl, mri_mtl, 1, 0, 1) ;
  dilations = nint(2.0 / (mri_aseg->xsize)) ;
  for (dno = 0 ; dno < dilations ; dno++)
    MRIdilate(mri_mtl, mri_mtl) ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.l0.mgz", hemi) ;
    printf("writing input volume to synth as %s (and input aseg as a.mgz)\n", fname) ;
    MRIwrite(mri_src, fname);
    MRIwrite(mri_aseg, "a.mgz");
  }
  if (mri_dst == NULL)
  {
    mri_dst = MRIcloneDifferentType(mri_src, MRI_INT) ;
    MRIcopy(mri_src, mri_dst) ;
  }
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_src, mri_aseg) ;

  for (x = 0 ; x < mri_dst->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_dst->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xa = (int)(nint(V3_X(v2))) ;
        ya = (int)(nint(V3_Y(v2))) ;
        za = (int)(nint(V3_Z(v2))) ;
        if (xa < 0 || ya < 0 || za < 0 ||
            xa >= mri_aseg->width || ya >= mri_aseg->height || za >= mri_aseg->depth)
          continue ;
        if (xa == Gx && ya == Gy && za == Gz)
          DiagBreak() ;
        aseg_label = nint(MRIgetVoxVal(mri_aseg, xa, ya, za, 0)) ;
        seg_label = nint(MRIgetVoxVal(mri_src, x, y, z, 0)) ;
	// already labeled properly
	// check to see if we are in MTL where surfaces can't be trusted
	in_mtl = (int)MRIgetVoxVal(mri_mtl, xa, ya, za, 0) ;

        if (IS_LAYER(seg_label) && !IS_HIPPO(aseg_label) && !IS_AMYGDALA(aseg_label) &&
	    !IS_CHOROID(aseg_label) && !IS_VENTRICLE(aseg_label)) 
	{
	  if (!in_mtl)
	    continue ;
	}

        switch (aseg_label)
        {
        case Right_Cerebral_Cortex:
	case Left_Cerebral_Cortex:
	  if (IS_CEREBRAL_WM(seg_label))   // if it is inside the white surface use the aseg value
	    MRIsetVoxVal(mri_dst, x, y, z, 0, seg_label) ;
	  else if (seg_label == CSF_VAL)  // was outside the ribbon
	    MRIsetVoxVal(mri_dst, x, y, z, 0, Unknown) ;
	  else
	    MRIsetVoxVal(mri_dst, x, y, z, 0, aseg_label) ;
	  break ;
        case Left_Cerebral_White_Matter:
        case Right_Cerebral_White_Matter:
//	case Left_choroid_plexus:
//	case Right_choroid_plexus:
	  if (IS_CEREBRAL_WM(seg_label))   // if it is inside the white surface use the aseg value
	    MRIsetVoxVal(mri_dst, x, y, z, 0, aseg_label) ;
	  else if (in_mtl)  // don't trust ribbon
	      MRIsetVoxVal(mri_dst, x, y, z, 0, aseg_label) ;
	  else if (seg_label == CSF_VAL)
	      MRIsetVoxVal(mri_dst, x, y, z, 0, aseg_label) ;
          break ;
        case Left_Cerebellum_Cortex:
        case Right_Cerebellum_Cortex:
	  if (!IS_CEREBRAL_WM(seg_label))  // otherwise it is not inside the white surface use cerebellar label
	    MRIsetVoxVal(mri_dst, x, y, z, 0, aseg_label) ;
          break ;
        default:
	case Unknown:
          MRIsetVoxVal(mri_dst, x, y, z, 0, aseg_label) ;
          break ;
        }
      }
    }
  }

  if (hemi_no != BOTH_HEMIS)
    MRIeraseOtherHemi(mri_dst, hemi) ;
  VectorFree(&v1) ; VectorFree(&v2) ; MatrixFree(&m_vox2vox) ; MRIfree(&mri_mtl) ;
  return(mri_dst) ;
}


void
MRIeraseOtherHemi(MRI *mri_dst, char *hemi)
{
  int label, x, y, z ;

  if (!strcmp(hemi, "rh"))
  {
    for (x = 0 ; x < mri_dst->width ; x++)
      for (y = 0 ; y < mri_dst->height ; y++)
	for (z = 0 ; z < mri_dst->depth ; z++)
	{
	  label = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
	  switch (label)
	  {
	  case Left_Cerebral_White_Matter:
	  case Left_Cerebral_Cortex:
	  case Left_choroid_plexus:
	  case Left_Cerebellum_White_Matter:
	  case Left_Cerebellum_Cortex:
	  case Left_Pallidum:
	  case Left_Thalamus:
	  case Left_Putamen:
	  case Left_Caudate:
	  case Left_Accumbens_area:  // remove them from cortex
	  case Left_Inf_Lat_Vent:
	  case Left_Amygdala:
	  case Left_Hippocampus:
	  case Left_Lateral_Ventricle:
	  case Left_VentralDC:
	    MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
	    break ;
	  default:
	    break ;
	  }
	}
  }
  else
  {
    for (x = 0 ; x < mri_dst->width ; x++)
      for (y = 0 ; y < mri_dst->height ; y++)
	for (z = 0 ; z < mri_dst->depth ; z++)
	{
	  label = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
	  switch (label)
	  {
	  case Right_VentralDC:
	  case Right_choroid_plexus:
	  case Right_Cerebral_White_Matter:
	  case Right_Cerebellum_White_Matter:
	  case Right_Cerebellum_Cortex:
	  case Right_Cerebral_Cortex:
	  case Right_Pallidum:
	  case Right_Thalamus:
	  case Right_Putamen:
	  case Right_Caudate:
	  case Right_Accumbens_area:  // remove them from cortex
	  case Right_Inf_Lat_Vent:
	  case Right_Amygdala:
	  case Right_Hippocampus:
	  case Right_Lateral_Ventricle:
	    MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
	    break ;
	  default:
	    break ;
	  }
	}
  }
}
