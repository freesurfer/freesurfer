/**
 * @brief program for computing/optimizing registration of a surface to a volume
 *        
 *
 * Program to compute a rigid alignment between a surface and a label by maximizing the gradient
 * magnitude across the gray/white boundary, divided by its variance
 * Now supports multiple similarity functions.
 */
/*
 * Original Author: Greg Grev
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
BEGINUSAGE --------------------------------------------------------------

 mri_register_to_label
  --reg regfile
  --mov fvol
  --surf surface   : surface to read in
  --res resolution : specify the resolution used for the distance transform
  --pial pial surface name   : pial surface to read in
  --pial_only pial surface name   : pial surface to read in (don't use white in similarity)
  --targ vol
  --angle_init ax0 ay0 az0 : search around these rotation angles instead of 0
  --trans_init tx0 ty0 tz0 : search around these translations instead of 0

  --median        : use median cost function
  --L1            : use L1 norm cost function
  --patch patch   :  patch  to read in

  --cost costfile

  --interp interptype : interpolation trilinear or nearest (def is trilin)
  --profile : print out info about exec time

  --border border : size of border region to ignore

  --out-reg outreg : reg at lowest cost (updated continuously)

ENDUSAGE ---------------------------------------------------------------
*/

/*
BEGINHELP --------------------------------------------------------------

FORMATS

Data file format can be specified implicitly (through the path name)
or explicitly. All formats accepted by mri_convert can be used.

BUGS

sinc interpolation is broken except for maybe COR to COR.


BUG REPORTING

Report bugs to analysis-bugs@nmr.mgh.harvard.edu. Include the following
formatted as a list as follows: (1) command-line, (2) directory where
the program was run (for those in the MGH-NMR Center), (3) version,
(4) text output, (5) description of the problem.

SEE ALSO

mri_vol2vol mri_convert, tkregister2 mri_segreg


ENDHELP --------------------------------------------------------------

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"

#include "matrix.h"
#include "mri.h"
#include "version.h"
#include "mri2.h"
#include "mri_identify.h"
#include "MRIio_old.h"
#include "registerio.h"
#include "resample.h"
#include "gca.h"
#include "gcamorph.h"
#include "fio.h"
#include "cmdargs.h"
#include "pdf.h"
#include "timer.h"
#include "numerics.h"
#include "mri_circulars.h"
#include "romp_support.h"

#ifdef X
#undef X
#endif

static float resolution = 0.5 ;
static void check_options(void);
static int write_lta(MATRIX *m, char *fname_from_caller, MRI *mri_src, MRI *mri_dst);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
#include "tags.h"
static int istringnmatch(const char *str1, const char *str2, int n);

int main(int argc, char *argv[]) ;
static int  parse_commandline(int argc, char **argv);
static double
find_optimal_linear_xform(MRI *mri_dist, MRI *mri_src, LABEL *area,
			  MATRIX *R0,
			  float min_angle, float max_angle,
			  float min_scale, float max_scale,
			  float min_trans, float max_trans,
			  float angle_steps, float scale_steps, float trans_steps,
			  int nreductions, int cost_function,
			  double tx0, double ty0, double tz0,
			  double ax0, double ay0, double az0,
                          MRI *mri_targ_vol);

static double tx0 = 0 ;//-10 ;
static double ty0 = 0 ;
static double tz0 = 0 ;
static double ax0 = 0 ; //RADIANS(10) ;
static double ay0 = 0 ;
static double az0 = 0 ;

#define COST_RMS     0
#define COST_MEDIAN  1
#define COST_MAX     2
#define COST_L1      3

static int write_diags = 1; 
static int debug = 0  ;
static int downsample = 0 ;

static MRI *mri_targ_vol ;
static int cost_function = COST_RMS ;
static char *regfile, *surf_fname, *outregfile, *vol_fname, *targ_vol_fname, *subject, *patch_fname, *label_fname ;
static int write_iter = -1 ;
#define MAX_ROT RADIANS(10) 
#define MAX_TRANS 20
#define ANGLE_STEP_SIZE (.33333)
#define TRANS_STEP_SIZE 1
#define SCALE_STEP_SIZE .025

static double max_rot = MAX_ROT ;
static double max_trans = MAX_TRANS ;
static double min_scale = 1.0 ;
static double max_scale = 1.0 ;
static MRI_SURFACE   *mris ; 
static double angle_step_size = ANGLE_STEP_SIZE ;
static double trans_step_size = TRANS_STEP_SIZE ;
//static double scale_step_size = SCALE_STEP_SIZE ;
static double trans_steps = (2/TRANS_STEP_SIZE)*MAX_TRANS+1 ;
static double angle_steps = (2/ANGLE_STEP_SIZE)*DEGREES(MAX_ROT)+1 ;
static double scale_steps = 1 ;


const char *Progname ;

static char base_name[STRLEN] ;
static char *SUBJECTS_DIR ;
static MATRIX *R0 ;
int
main(int argc, char **argv) 
{
  char          fname[STRLEN] ;
  MRI           *mri_dist, *mri_src ;
  LABEL         *area, *ltmp ;
  int           msec, nargs ;
  Timer then ;

  nargs = handleVersionOption(argc, argv, "mris_register_to_label");
  if(nargs && argc - nargs == 1) exit (0);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= DIAG_SHOW ;

  then.reset() ;
  if(argc == 0) 
    usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  printf("Loading surface %s\n", surf_fname);
  mris = MRISread(surf_fname) ;
  if (mris == NULL)
    ErrorExit(Gerror, "") ;
  area = LabelRead(NULL, label_fname) ;
  if (area == NULL)
    ErrorExit(Gerror, "") ;
  printf("label loaded with %d points\n", area->n_points) ;
  ltmp = LabelRemoveAlmostDuplicates(area, resolution, NULL) ;
  LabelFree(&area) ; area = ltmp ;
  mri_src = MRIread(vol_fname) ;
  if (mri_src == NULL)
    ErrorExit(Gerror, "") ;

  if (targ_vol_fname)
  {
    mri_targ_vol = MRIread(targ_vol_fname) ;
    if (mri_targ_vol == NULL)
      ErrorExit(Gerror, "") ;
  }
  else
    mri_targ_vol = NULL ;

  if (regfile != NULL) 
  {
    LTA *lta = LTAread(regfile) ;
    if (lta == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load initial registration from %s", regfile) ;
    if (lta->type != LINEAR_RAS_TO_RAS)
      LTAchangeType(lta, LINEAR_RAS_TO_RAS) ;
    R0 = MatrixCopy(lta->xforms[0].m_L, NULL) ;
    LTAfree(&lta) ;
  }
  FileNameRemoveExtension(outregfile, base_name) ;
  mri_dist = MRIScomputeDistanceToSurface(mris, NULL, resolution) ;
  MRIabs(mri_dist, mri_dist) ;
  if (write_diags)
  {
    int req = snprintf(fname, STRLEN, "%s.dist.mgz", base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing distance transform to %s\n", fname) ;
    MRIwrite(mri_dist, fname) ;
  }

  if (strstr(area->space, "TkReg") != NULL)
  {
#if 0
    printf("!!!!!!!!! converting label to scanner coords from TkReg !!!!!!!!\n") ;
    LabelToScannerRAS(area, mri_src, area) ;   // not implemented due to freeview bug
#endif
    sprintf(area->space, "coords=scanner") ;
  }
  area->coords = LABEL_COORDS_SCANNER_RAS ;
  find_optimal_linear_xform(mri_dist, mri_src, area,
			    R0,
			    -max_rot, max_rot,
			    min_scale, max_scale,
			    -max_trans, max_trans,
			    angle_steps, scale_steps, trans_steps, 0, cost_function,
			    tx0, ty0, tz0, ax0, ay0, az0, mri_targ_vol) ;
  printf("writing final matrix to %s:\n", outregfile) ;
  MatrixPrint(stdout, R0) ;
  if (TransformFileNameType(outregfile) == REGISTER_DAT)
    regio_write_surfacexform_to_register_dat(R0, outregfile, mris, mri_src, subject, FLT2INT_ROUND) ;
  else
    write_lta(R0,  outregfile, mri_dist, mri_src) ;
  msec = then.milliseconds() ;
  fprintf(stderr,
          "registration took %2.1f minutes\n", (float)msec/(60*1000.0f));
  exit(0) ;
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
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

    if (!strcasecmp(option,      "--help"))     print_help() ;
    else if (!strcasecmp(option, "--version"))  print_version() ;
    else if (!strcasecmp(option, "--debug"))    debug = 1;
    else if (istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--out-reg",0)) {
      if (nargc < 1) argnerr(option,1);
      outregfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--res",0)) {
      if (nargc < 1) argnerr(option,1);
      resolution = atof(pargv[0]);
      nargsused = 1;
    } else if (istringnmatch(option, "--downsample",0)) {
      if (nargc < 1) argnerr(option,1);
      downsample = atoi(pargv[0]);
      nargsused = 1;
    } else if (istringnmatch(option, "--mov",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      vol_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--targ",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      targ_vol_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--angle_init",0)) {
      if (nargc < 3) 
        argnerr(option,1);
      ax0 = RADIANS(atof(pargv[0])) ;
      ay0 = RADIANS(atof(pargv[1])) ;
      az0 = RADIANS(atof(pargv[2])) ;
      printf("performing angular search around (%2.1f, %2.1f, %2.1f) degrees\n",
	     DEGREES(ax0), DEGREES(ay0), DEGREES(az0)) ;
      nargsused = 3;
    } else if (istringnmatch(option, "--trans_init",0)) {
      if (nargc < 3) 
        argnerr(option,1);
      tx0 = (atof(pargv[0])) ;
      ty0 = (atof(pargv[1])) ;
      tz0 = (atof(pargv[2])) ;
      printf("performing translation search around (%2.1f, %2.1f, %2.1f) \n", tx0, ty0, tz0) ;
      nargsused = 3;
    } else if (istringnmatch(option, "--surf",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      surf_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--s",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--patch",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      patch_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--label",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      label_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--sdir",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--w")) {
      write_iter = atoi(pargv[0]) ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--max")) {
      cost_function = COST_MAX ;
      nargsused = 0;
      printf("using MAX cost functional\n") ;
    } else if (!strcasecmp(option, "--rms")) {
      cost_function = COST_RMS ;
      nargsused = 0;
      printf("using RMS cost functional\n") ;
    } else if (!strcasecmp(option, "--median")) {
      cost_function = COST_MEDIAN ;
      nargsused = 0;
      printf("using median cost functional\n") ;
    } else if (!strcasecmp(option, "--L1")) {
      cost_function = COST_L1 ;
      nargsused = 0;
      printf("using L1 cost functional\n") ;
    } else if (!strcasecmp(option, "--max_rot")) {
      if (nargc < 1) CMDargNErr(option,1);
      max_rot = RADIANS(atof(pargv[0])) ;
      angle_steps = (2/angle_step_size)*DEGREES(max_rot)+1 ;
      printf("setting max_rot = %2.1f, angle_steps = %d\n", DEGREES(max_rot), nint(angle_steps)) ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--max_trans")) {
      if (nargc < 1) CMDargNErr(option,1);
      max_trans = atof(pargv[0]) ;
      trans_steps = (2/trans_step_size)*max_trans+1 ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--scale")) {
      if (nargc < 2) CMDargNErr(option,1);
      sscanf(pargv[0], "%lf:%lf:%lf", &min_scale, &scale_steps, &max_scale) ;
      printf("searching scale range %2.4f --> %2.4f in %2.0f steps\n", min_scale, max_scale, scale_steps) ;
      nargsused = 1;
    } else if ( !strcmp(option, "--trans_step_size") ) {
      if (nargc < 1) argnerr(option,1);
      trans_step_size = atof(pargv[0]) ;
      printf("using trans step size %2.1f\n", trans_step_size) ;
      nargsused = 1;
    } else if ( !strcmp(option, "--angle_step_size") ) {
      if (nargc < 1) argnerr(option,1);
      angle_step_size = atof(pargv[0]) ;
      printf("using angle step size %2.1f\n", angle_step_size) ;
      nargsused = 1;
    } else if ( !strcmp(option, "--gdiagno") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
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
printf("\n");
printf("mris_register_to_label\n");
printf("  --surf surface\n");
printf("  --reg regfile\n");
printf("  --mri_reg fvol\n");
printf("\n");
printf("\n");
printf("  --cost costfile\n");
printf("\n");
printf("  --mov volfname : volume on which label points are specified\n");
printf("  --res resolution\n");
printf("  --max_rot angle : max angle (degrees) to search over\n");
printf("  --max_trans dist :max translation (mm) to search over\n");
printf("  --s subject     : specify name of subject (for register.dat file)\n");
printf("  --label <label>   : load label <label> and limit calculations to it\n") ;
printf("\n");
printf("  --out-reg outreg : reg at lowest cost (updated continuously)\n");
printf("  --downsample N : downsample input volume by a factor of N\n");
printf("\n");

}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n%s\n\n",getVersion().c_str());
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (SUBJECTS_DIR == NULL)
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    printf("ERROR: SUBJECTS_DIR undefined.\n");
    exit(1);
  }
  if (vol_fname == NULL) {
    printf("ERROR: No mov volume supplied.\n");
    exit(1);
  }
  if (surf_fname == NULL) {
    printf("ERROR: No surface supplied.\n");
    exit(1);
  }

  if (label_fname == NULL) {
    printf("ERROR: No label supplied.\n");
    exit(1);
  }
  if (regfile == NULL) {
    printf("using identity as initial registration\n") ;
    R0 = MatrixIdentity(4, NULL) ;
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) 
{
  //  int n;
  fprintf(fp,"movvol %s\n",vol_fname);
  fprintf(fp,"surface %s\n",surf_fname);
  fprintf(fp,"regfile %s\n",regfile);
  if(outregfile) fprintf(fp,"outregfile %s\n",outregfile);
  if (Gdiag_no >= 0)
    fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);
  if (downsample > 0)
    printf("downsampling input volume by a factor of %d\n", downsample) ;
  return;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
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
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);
  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*------------------------------------------------------------
  istringnmatch() - compare the first n characters of two strings,
  return a 1 if they match (ignoring case), a zero otherwise. If
  n=0, then do a full comparison.
  ------------------------------------------------------------*/
static int istringnmatch(const char *str1, const char *str2, int n) {
  if (n > 0  && ! strncasecmp(str1,str2,n)) return(1);
  if (n <= 0 && ! strcasecmp(str1,str2)) return(1);
  return(0);
}

LABEL *
LabelTransformToVoxelCoords(LABEL *lsrc, MRI *mri, MATRIX *R, MATRIX *Mras2vox, LABEL *ldst)
{
//  double xv, yv, zv ;
  int    i ;
  static VECTOR *v1 = NULL,  *v2 ;
  MATRIX  *M ;

  if (v1 == NULL)  // first time
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1,4) = 1.0 ;
    VECTOR_ELT(v2,4) = 1.0 ;
  }
  
  if (ldst == NULL)
  {
    ldst = LabelClone(lsrc) ;
    ldst->n_points = lsrc->n_points ;
  }
  M = MatrixMultiply(Mras2vox, R, NULL) ;
  for (i = 0 ; i < lsrc->n_points ; i++)
  {
    V3_X(v1) = lsrc->lv[i].x ; V3_Y(v1) = lsrc->lv[i].y ; V3_Z(v1) = lsrc->lv[i].z ;
    MatrixMultiply(M, v1, v2) ;
    ldst->lv[i].x = V3_X(v2) ; ldst->lv[i].y = V3_Y(v2) ; ldst->lv[i].z = V3_Z(v2) ;
    ldst->lv[i].stat = lsrc->lv[i].stat ;
  }

  strcpy(ldst->space, "voxel") ;
  MatrixFree(&M) ;
  return(ldst) ;
}

LABEL *
LabelFromVoxelCoords(LABEL *lsrc, MRI *mri, LABEL *ldst, MATRIX *Mvox2ras)
{
  int    i ;
//  double r, a, s ;
  static VECTOR *v1 = NULL,  *v2 ;

  if (v1 == NULL)  // first time
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1,4) = 1.0 ;
    VECTOR_ELT(v2,4) = 1.0 ;
  }
  
  if (ldst == NULL)
  {
    ldst = LabelClone(lsrc) ;
    ldst->n_points = lsrc->n_points ;
  }
  for (i = 0 ; i < lsrc->n_points ; i++)
  {
    V3_X(v1) = lsrc->lv[i].x ; V3_Y(v1) = lsrc->lv[i].y ; V3_Z(v1) = lsrc->lv[i].z ;
    MatrixMultiply(Mvox2ras, v1, v2) ;
    ldst->lv[i].x = V3_X(v2) ; ldst->lv[i].y = V3_Y(v2) ;  ldst->lv[i].z = V3_Z(v2) ;
    ldst->lv[i].stat = lsrc->lv[i].stat ;
  }

  strcpy(ldst->space, "scanner") ;
  return(ldst) ;
}

static double
compute_error(MRI *mri_dist, MATRIX *R0, MATRIX *Mras2vox, LABEL *lras, LABEL *lvol, int cost, double current_error) 
{
  double sse, val, max_val = 0 ;
  int    i ;

  lvol = LabelTransformToVoxelCoords(lras, mri_dist, R0, Mras2vox, lvol) ;
  sse = 0 ;
  switch (cost)
  {
  case COST_MAX:
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) firstprivate(val, mri_dist) shared(lvol) schedule(static,1) 
#endif
    for (i = 0 ; i < lvol->n_points ; i++)
    {
      ROMP_PFLB_begin
      MRIsampleVolume(mri_dist, lvol->lv[i].x, lvol->lv[i].y, lvol->lv[i].z, &val) ;
      lvol->lv[i].stat = val ;
      if (val > current_error)
	max_val = val ; // no point in looking any furter
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    if (!FZERO(max_val))
      return(max_val) ;

    for (i = 0 ; i < lvol->n_points ; i++)
    {
      MRIsampleVolume(mri_dist, lvol->lv[i].x, lvol->lv[i].y, lvol->lv[i].z, &val) ;
      if (val > max_val)
      {
	max_val = val ;
	if (val > current_error)
	  return(val) ; // no point in looking any furter
      }
    }
    return(max_val) ;
    break ;

  default:
  case COST_RMS:
     ROMP_PF_begin
#ifdef HAVE_OPENMP
     #pragma omp parallel for if_ROMP(experimental) firstprivate(val, mri_dist) shared(lvol) schedule(static,1) reduction(+: sse)
#endif
    for (i = 0 ; i < lvol->n_points ; i++)
    {
      ROMP_PFLB_begin
      MRIsampleVolume(mri_dist, lvol->lv[i].x, lvol->lv[i].y, lvol->lv[i].z, &val) ;
      lvol->lv[i].stat = val ;
      sse += val*val ;
      ROMP_PFLB_end
    }
    ROMP_PF_end
    return(sqrt(sse/lras->n_points)) ;
    break ;
  case COST_L1:
    ROMP_PF_begin
#ifdef HAVE_OPENMP
#pragma omp parallel for if_ROMP(experimental) firstprivate(val, mri_dist) shared(lvol) schedule(static,1) reduction(+: sse)
#endif
    for (i = 0 ; i < lvol->n_points ; i++)
    {
      ROMP_PFLB_begin
      MRIsampleVolume(mri_dist, lvol->lv[i].x, lvol->lv[i].y, lvol->lv[i].z, &val) ;
      lvol->lv[i].stat = val ;
      sse += fabs(val) ;
      ROMP_PFLB_end
    }
    ROMP_PF_end
    return(sse/lras->n_points) ;
    break ;
  }

  return(0) ;  // never get here
}

static double
find_optimal_linear_xform(MRI *mri_dist, MRI *mri_src, LABEL *area,
			  MATRIX *R0,
			  float min_angle, float max_angle,
			  float min_scale, float max_scale,
			  float min_trans, float max_trans,
			  float angle_steps, float scale_steps, float trans_steps,
			  int nreductions, int cost_function,
			  double tx0, double ty0, double tz0,
			  double ax0, double ay0, double az0, MRI *mri_targ_vol)
{
  double  min_error, error ;
  LABEL   *lvol ;
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp, *m_origin, *m_origin_inv, *m_inv, 
    *m_tmp2, *m_scale, *m_trans, *m_tmp3 = NULL, *m_best, *Mras2vox, *Mvox2ras, *m_xy_rot, *m_best_inv ;
  double  x_max_rot;
  double  y_max_rot, z_max_rot, delta_rot;
  double  x_max_scale, y_max_scale, z_max_scale;
  double  delta_scale = 1, delta_trans;
  double   mean_angle, x_angle, y_angle, z_angle, x_scale, y_scale, z_scale, min_z_scale, max_z_scale, delta_z_scale;
  double  mean_scale, x_max_trans, y_max_trans, z_max_trans, mean_trans, x_trans, y_trans, z_trans ;
  int     i, nfound = 0  ;

  
  m_origin = MatrixIdentity(4, NULL) ;
  // set the translation  (center = 1)
  if (mri_targ_vol && regfile == NULL)  // if registration already specified, don't mess with it
  {
    *MATRIX_RELT(m_origin, 1, 4) = mri_targ_vol->c_r ;
    *MATRIX_RELT(m_origin, 2, 4) = mri_targ_vol->c_a ;
    *MATRIX_RELT(m_origin, 3, 4) = mri_targ_vol->c_s ;
  }
#if 0
  *MATRIX_RELT(m_origin, 1, 4) = mri_src->c_r ;
  *MATRIX_RELT(m_origin, 2, 4) = mri_src->c_a ;
  *MATRIX_RELT(m_origin, 3, 4) = mri_src->c_s ;
  *MATRIX_RELT(m_origin, 1, 4) = 0 ;
  *MATRIX_RELT(m_origin, 2, 4) = 0 ;
  *MATRIX_RELT(m_origin, 3, 4) = 0 ;
#endif
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;
  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  Mvox2ras = MRIgetVoxelToRasXform(mri_dist) ;
  Mras2vox = MatrixInverse(Mvox2ras, NULL) ;
    
  m_trans = MatrixIdentity(4, NULL) ;
  // set the translation  (center = 1)

  m_best = MatrixCopy(R0, NULL) ;
  m_best_inv = MatrixInverse(m_best, NULL) ;
  lvol = LabelTransformToVoxelCoords(area, mri_dist, R0, Mras2vox, NULL) ;
  min_error = compute_error(mri_dist, R0, Mras2vox, area, lvol, cost_function, 1e10) ;
  printf("starting error = %2.4f\n", min_error) ;

  m_L_tmp = m_x_rot = m_xy_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = m_inv = NULL ;
  x_max_trans = y_max_trans = z_max_trans = x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;

  if (mri_src->depth == 1)
  {
    min_z_scale = max_z_scale = 1.0 ; delta_z_scale = 1.0 ;
  }
  else
  {
    min_z_scale = min_scale ; max_z_scale = max_scale ; delta_z_scale = delta_scale ;
  }
  for (i = 0 ; i < nreductions+1 ; i++)
  {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
    if (scale_steps <= 1)
      delta_scale = 1 ;
    else
      delta_scale = (max_scale-min_scale) / (scale_steps-1) ;
    if (FZERO(delta_scale))
    {
      delta_scale = max_scale ;
    }

    if (angle_steps == 1)
    {
      min_angle = max_angle = 0.0 ;
      delta_rot = 1 ;
    }
    else
    {
      delta_rot = (max_angle-min_angle) / (angle_steps-1) ;
    }
    if (Gdiag & DIAG_SHOW)
    {
      printf("  scanning %2.2f degree nbhd (%2.3f)\n"
             "  scale %2.3f->%2.3f (step %2.3f), "
             "trans %2.2f->%2.2f (step %2.2f)\n",
             (float)DEGREES(max_angle), (float)DEGREES(delta_rot),
             min_scale,max_scale, delta_scale,
             min_trans, max_trans, delta_trans);
      fflush(stdout) ;
    }

  if (write_diags)
  {
    char fname[STRLEN] ;
    MATRIX *OCT_voxel_to_RAS, *OCT_voxel_to_xformed_RAS  ;
    LABEL *l = LabelFromVoxelCoords(lvol, mri_dist, NULL, Mvox2ras) ;
    MRI   *mri_tmp ;
    
    // compute scanner vox2ras
    OCT_voxel_to_RAS = MRIgetVoxelToRasXform(mri_src);
    OCT_voxel_to_xformed_RAS = MatrixMultiply(R0, OCT_voxel_to_RAS, NULL) ;
    mri_tmp = MRIcopy(mri_src, NULL) ;
    MRIsetVoxelToRasXform(mri_tmp, OCT_voxel_to_xformed_RAS)  ;

    if (debug) {
      int req = snprintf(fname, STRLEN, "%s.opt.%4.4d.label", "debug",nfound) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s.opt.%4.4d.label", base_name,nfound) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("writing new optimal label to %s, v(0): RAS=%2.1f %2.1f %2.1f, VOX=(%d %d %d)\n", 
	   fname, l->lv[0].x, l->lv[0].y, l->lv[0].z,
	   nint(lvol->lv[0].x), nint(lvol->lv[0].y), nint(lvol->lv[0].z)) ;
    LabelWrite(l, fname) ;
    if (debug) {
      int req = snprintf(fname, STRLEN, "%s.opt.%4.4d.mgz", "debug",nfound) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s.opt.%4.4d.mgz", base_name,nfound) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("saving %s\n", fname) ;
    MRIwrite(mri_tmp, fname) ;
    MRIfree(&mri_tmp) ; LabelFree(&l) ;
    MatrixFree(&OCT_voxel_to_RAS) ;  MatrixFree(&OCT_voxel_to_xformed_RAS) ;
  
    // compute scanner vox2ras for target vol if given
    if (mri_targ_vol)
    {
      MATRIX *MRI_voxel_to_RAS, *MRI_voxel_to_xformed_RAS  ;
      char   fname[STRLEN], ext[STRLEN], fname_only[STRLEN] ;
      
      MRI_voxel_to_RAS = MRIgetVoxelToRasXform(mri_targ_vol);
      MRI_voxel_to_xformed_RAS = MatrixMultiply(R0, MRI_voxel_to_RAS, NULL) ;
      mri_tmp = MRIcopy(mri_targ_vol, NULL) ;
      MRIsetVoxelToRasXform(mri_tmp, MRI_voxel_to_xformed_RAS)  ;
      int req = snprintf(fname, STRLEN, "%s.targ.%4.4d.mgz", base_name,nfound) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing reoriented target volume to %s\n", fname) ;
      MRIwrite(mri_tmp, fname) ;
      
      MRIfree(&mri_tmp) ; 
      MatrixFree(&MRI_voxel_to_RAS) ; MatrixFree(&MRI_voxel_to_xformed_RAS) ;
      
      FileNameExtension(outregfile, ext) ;
      FileNameRemoveExtension(outregfile, fname_only) ;
      req = snprintf(fname, STRLEN, "%s.%4.4d.%s", fname_only, nfound, ext) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing intermediate registration to %s\n", fname) ;
      
      if (TransformFileNameType(outregfile) == REGISTER_DAT)
	regio_write_surfacexform_to_register_dat(m_best, fname, mris, mri_src, subject, FLT2INT_ROUND) ;
      else
	write_lta(m_best,  fname, mri_dist, mri_src) ;
    }
  }

// scale /////////////////////////////////////////////////////////////
  for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta_scale)
  {
    printf("x_scale = %2.3f\n", x_scale) ;
    *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
    for (y_scale = min_scale ;
	 y_scale <= max_scale ;
	 y_scale += delta_scale)
    {
      *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
      printf("y_scale = %2.3f\n", y_scale) ;
      for (z_scale= min_z_scale ;
	   z_scale <= max_z_scale;
	   z_scale += delta_z_scale)
      {
	*MATRIX_RELT(m_scale, 3, 3) = z_scale ;
	
	/* reset scaling values */
	*MATRIX_RELT(m_scale, 1, 4) =
	  *MATRIX_RELT(m_scale, 2, 4) =
	  *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
	m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
	MatrixMultiply(m_origin, m_tmp, m_scale) ;
	
	// angle //////////////////////////////
	for (x_angle = min_angle+ax0 ;
	     x_angle <= max_angle+ax0 ;
	     x_angle += delta_rot)
	{
	  if (debug)
	    x_angle = RADIANS(-9) ;
	  printf("------------ x_angle = %2.3f ----------------\n", DEGREES(x_angle)) ;
	  m_x_rot = MatrixReallocRotation (4, x_angle, X_ROTATION, m_x_rot) ;
	  for (y_angle = min_angle+ay0 ;
	       y_angle <= max_angle+ay0 ;
	       y_angle += delta_rot)
	  {
	    if (debug)
	      y_angle = RADIANS(9.5) ;
	    m_y_rot = MatrixReallocRotation(4, y_angle, Y_ROTATION, m_y_rot);
	    m_xy_rot = MatrixMultiply(m_y_rot, m_x_rot, m_xy_rot) ;
	    for (z_angle= min_angle+az0;
		 z_angle <= max_angle+az0;
		 z_angle += delta_rot)
	    {
	      if (debug)
		z_angle = RADIANS(-9.5) ;
	      m_z_rot = MatrixReallocRotation(4, z_angle,Z_ROTATION,m_z_rot);
	      m_rot = MatrixMultiply(m_z_rot, m_xy_rot, m_rot) ;
	      MatrixMultiply(m_rot, m_origin_inv, m_tmp) ;
	      MatrixMultiply(m_origin, m_tmp, m_rot) ;
	      m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
	      m_tmp3 = MatrixMultiply(m_tmp2, R0, m_tmp3) ;
	      
	      // translation //////////
	      for (x_trans = min_trans+tx0 ;
		   x_trans <= max_trans+tx0 ;
		   x_trans += delta_trans)
	      {
		if (debug)
		  x_trans = -9 ;
		*MATRIX_RELT(m_trans, 1, 4) = x_trans ;
		for (y_trans = min_trans+ty0 ;
		     y_trans <= max_trans+ty0 ;
		     y_trans += delta_trans)
		{
		  if (debug)
		    y_trans = -5 ;
		  *MATRIX_RELT(m_trans, 2, 4) =  y_trans ;
		  for (z_trans= min_trans+tz0 ;
		       z_trans <= max_trans+tz0 ;
		       z_trans += delta_trans)
		  {
		    if (debug)
		      z_trans = -5 ;
		    *MATRIX_RELT(m_trans, 3, 4) = z_trans ;
		    
		    m_L_tmp = MatrixMultiply(m_trans, m_tmp3, m_L_tmp) ;
		    m_inv = MatrixInverse(m_L_tmp, m_inv) ;
		    error = compute_error(mri_dist, m_inv, Mras2vox, area, lvol, cost_function, min_error) ;
		    if (error < min_error || debug)
		    {
		      MatrixCopy(m_L_tmp, m_best) ;
		      MatrixCopy (m_inv, m_best_inv) ;
		      min_error = error ;
		      x_max_scale = x_scale ; y_max_scale = y_scale ; z_max_scale = z_scale ;
		      x_max_rot = x_angle ;   y_max_rot = y_angle ;   z_max_rot = z_angle ;
		      x_max_trans = x_trans ; y_max_trans = y_trans ; z_max_trans = z_trans ;
#if 1
		      printf("\noptimum #%d\n", nfound+1) ;
		      printf( "%s: min_error = %2.4f\n", __FUNCTION__, min_error );
		      printf( "%s: Translation (%4.2f, %4.2f, %4.2f)\n",
			      __FUNCTION__, x_trans, y_trans, z_trans );
		      printf( "%s: Rotation (%4.3f, %4.3f, %4.3f)\n",
			      __FUNCTION__, DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle) );
		      printf( "%s: Scale (%4.4f, %4.4f, %4.4f)\n",
			      __FUNCTION__, x_scale, y_scale, z_scale );
#endif
		      if (write_diags)
		      {
			char fname[STRLEN] ;
			MATRIX *OCT_voxel_to_RAS, *OCT_voxel_to_xformed_RAS  ;
			LABEL *l = LabelFromVoxelCoords(lvol, mri_dist, NULL, Mvox2ras) ;
			MRI   *mri_tmp ;
			
			// compute scanner vox2ras
			OCT_voxel_to_RAS = MRIgetVoxelToRasXform(mri_src);
			OCT_voxel_to_xformed_RAS = MatrixMultiply(m_best_inv, OCT_voxel_to_RAS, NULL) ;
			mri_tmp = MRIcopy(mri_src, NULL) ;
			MRIsetVoxelToRasXform(mri_tmp, OCT_voxel_to_xformed_RAS)  ;
			
			if (debug) {
			  sprintf(fname, "%s.opt.%4.4d.label", "debug",++nfound) ;
			} else {
			  int req = snprintf(fname, STRLEN, "%s.opt.%4.4d.label", base_name,++nfound) ;
			  if( req >= STRLEN ) {
			    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
			  }
			}
			printf("writing new optimal label to %s, v(0): RAS=%2.1f %2.1f %2.1f, VOX=(%d %d %d)\n", 
			       fname, l->lv[0].x, l->lv[0].y, l->lv[0].z,
			       nint(lvol->lv[0].x), nint(lvol->lv[0].y), nint(lvol->lv[0].z)) ;
			LabelWrite(l, fname) ;
			if (debug) {
			  sprintf(fname, "%s.opt.%4.4d.mgz", "debug",nfound) ;
			} else {
			  int req = snprintf(fname, STRLEN, "%s.opt.%4.4d.mgz", base_name,nfound) ;
			  if( req >= STRLEN ) {
			    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
			  }
			}
			printf("saving %s\n", fname) ;
			MRIwrite(mri_tmp, fname) ;
			if (debug || nfound == Gdiag_no)
			  MatrixPrint(stdout, m_best) ;
			MRIfree(&mri_tmp) ; LabelFree(&l) ;
			MatrixFree(&OCT_voxel_to_RAS) ;  MatrixFree(&OCT_voxel_to_xformed_RAS) ;
			
			// compute scanner vox2ras for target vol if given
			if (mri_targ_vol)
			{
			  MATRIX *MRI_voxel_to_RAS, *MRI_voxel_to_xformed_RAS  ;
			  char   fname[STRLEN], ext[STRLEN], fname_only[STRLEN] ;
			  
			  MRI_voxel_to_RAS = MRIgetVoxelToRasXform(mri_targ_vol);
			  MRI_voxel_to_xformed_RAS = MatrixMultiply(m_L_tmp, MRI_voxel_to_RAS, NULL) ;
			  mri_tmp = MRIcopy(mri_targ_vol, NULL) ;
			  MRIsetVoxelToRasXform(mri_tmp, MRI_voxel_to_xformed_RAS)  ;
			  int req = snprintf(fname, STRLEN, "%s.targ.%4.4d.mgz", base_name,nfound) ;
			  if( req >= STRLEN ) {
			    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
			  }
			  printf("writing reoriented target volume to %s\n", fname) ;
			  MRIwrite(mri_tmp, fname) ;
			  
			  MRIfree(&mri_tmp) ; 
			  MatrixFree(&MRI_voxel_to_RAS) ; MatrixFree(&MRI_voxel_to_xformed_RAS) ;
			  
			  FileNameExtension(outregfile, ext) ;
			  FileNameRemoveExtension(outregfile, fname_only) ;
			  req = snprintf(fname, STRLEN, "%s.%4.4d.%s", fname_only, nfound, ext) ;
			  if( req >= STRLEN ) {
			    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
			  }
			  printf("writing intermediate registration to %s\n", fname) ;
			  
			  if (TransformFileNameType(outregfile) == REGISTER_DAT)
			    regio_write_surfacexform_to_register_dat(m_best, fname, mris, mri_src, subject, FLT2INT_ROUND) ;
			  else
			    write_lta(m_best,  fname, mri_dist, mri_src) ;
			}
			
			if (debug)
			  exit(0) ;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  if (Gdiag & DIAG_SHOW)
  {
    printf("  min error = %2.1f @ R=(%2.3f,%2.3f,%2.3f),"
	   "S=(%2.3f,%2.3f,%2.3f), T=(%2.1f,%2.1f,%2.1f)\n",
	   min_error, DEGREES(x_max_rot), DEGREES(y_max_rot),
	   DEGREES(z_max_rot),x_max_scale, y_max_scale, z_max_scale,
	   x_max_trans, y_max_trans,z_max_trans) ;
  }
  

#if 0
    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) =
      *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
    *MATRIX_RELT(m_scale,1,1) = x_max_scale ;
    *MATRIX_RELT(m_scale,2,2) = y_max_scale ;
    *MATRIX_RELT(m_scale,3,3) = z_max_scale ;

    x_max_scale = y_max_scale = z_max_scale = 1.0 ;

    /* update L to reflect new maximum and search around it */
    m_x_rot = MatrixReallocRotation(4, x_max_rot, X_ROTATION, m_x_rot) ;
    m_y_rot = MatrixReallocRotation(4, y_max_rot, Y_ROTATION, m_y_rot) ;
    m_z_rot = MatrixReallocRotation(4, z_max_rot, Z_ROTATION, m_z_rot) ;
    m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
    m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;

    m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
    m_tmp3 = MatrixMultiply(m_tmp2, R0, m_tmp3) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max_trans ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max_trans ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max_trans ;
    m_L_tmp = MatrixMultiply(m_trans, m_tmp3, m_L_tmp) ;



    MatrixCopy(m_L_tmp, R0) ;

    x_max_trans = y_max_trans = z_max_trans = 0.0 ;
    x_max_rot = y_max_rot = z_max_rot = 0.0 ;
#endif

    mean_scale = (max_scale + min_scale) / 2 ;
    delta_scale = (max_scale-min_scale)/4 ;
    min_scale = mean_scale - delta_scale ;
    max_scale = mean_scale + delta_scale ;
    mean_trans = (max_trans + min_trans) / 2 ;
    delta_trans = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta_trans ;
    max_trans = mean_trans + delta_trans ;
    
    mean_angle = (max_angle + min_angle) / 2 ;
    delta_rot = (max_angle-min_angle)/4 ;
    min_angle = mean_angle - delta_rot ;
    max_angle = mean_angle + delta_rot ;
  }


  MatrixCopy(m_best, R0) ;

  MatrixFree(&m_xy_rot) ; MatrixFree(&m_x_rot) ; MatrixFree(&m_y_rot) ; MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;  MatrixFree(&m_scale) ; MatrixFree(&m_trans) ; 
  MatrixFree(&m_tmp) ; MatrixFree(&m_tmp2) ; MatrixFree(&m_tmp3) ;
  MatrixFree(&m_L_tmp) ; MatrixFree(&m_best) ; MatrixFree(&m_best_inv) ; MatrixFree(&m_inv) ;

  LabelFree(&lvol) ;

  return(min_error) ;
}

static int
write_lta(MATRIX *m, char *fname_from_caller, MRI *mri_src, MRI *mri_dst)
{
  char *dot, fname[STRLEN] ;
  LTA  *lta ;

  strcpy(fname, fname_from_caller) ;
  dot = strrchr(fname, '.') ;
  strcpy(dot+1, "lta") ;

  lta = LTAalloc(1, NULL) ;
  lta->type = LINEAR_RAS_TO_RAS ;

  LTAsetVolGeom(lta, mri_src, mri_dst) ;
  MatrixCopy(m, lta->xforms[0].m_L) ;

  LTAwriteEx(lta, fname) ;
  return(NO_ERROR) ;
}
