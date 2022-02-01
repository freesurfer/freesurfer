#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>

#include "version.h"
#include "diag.h"
#include "log.h"
#include "cmdargs.h"
#include "fio.h"
#include "mri.h"
#include "mriBSpline.h"
#include "GradUnwarp.h"

/* examples:
 *
 * unwarp at given crs
 * 1. mris_gradient_unwarp/mris_gradient_unwarp 
 *    --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad 
 *    --invol $FS_TEST/gradunwarp/example/orig.mgz 
 *    --crs 0,0,0
 *
 * unwarp for given ras
 * 2. mris_gradient_unwarp/mris_gradient_unwarp 
 *    --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad 
 *    --invol $FS_TEST/gradunwarp/example/orig.mgz 
 *    --ras 0.1,0.2,0.3
 *
 * save gradient unwarp transformation table
 * 3. mris_gradient_unwarp/mris_gradient_unwarp 
 *    --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad 
 *    --invol $FS_TEST/gradunwarp/example/orig.mgz 
 *    --save_transtbl gradunwarp.m3z
 *    --nthreads 10
 *
 * transform given volume, save gradient unwarp transformation table
 * 4. mris_gradient_unwarp/mris_gradient_unwarp 
 *    --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad 
 *    --invol $FS_TEST/gradunwarp/example/orig.mgz 
 *    --outvol orig.unwarped.cubic.mgz --interp cubic 
 *    --save_transtbl gradunwarp.m3z
 *    --nthreads 10
 *
 * transform given volume using input gradient unwarp transformation table
 * 5. mris_gradient_unwarp/mris_gradient_unwarp 
 *    --load_transtbl $FS_TEST/gradunwarp/transtbl/savetbl/gradunwarp.m3z
 *    --invol $FS_TEST/gradunwarp/example/orig.mgz 
 *    --outvol orig.unwarped.cubic.mgz --interp cubic 
 *    --nthreads 10
 *
 */
/****** BIG PICTURES ******/
/*
 * MATRIX *vox2ras = MRIxfmCRS2XYZ(const MRI *mri, 0)
 *
 * vol = MRIread('vol.mgz');
 * [Ax Ay Az Bx By Bz] = ReadYourCoef('a.coef');
 * for col = 0; col < vol->width; col++
 *   for row
 *     for slice
 *       RAS = vox2ras*[C R S 1]'
 *       DeltaRAS = YourFunction(Ax,Ay,Az,Bx,By,Bz,RAS,R0);
 *       DistortedRAS = RAS + DeltaRAS;
 *       DistortedCRS = inv(vox2ras)*DistortedRAS; // floating point
 *       NewVol(col,row,slice) =
 *           MRIsampleVol(vol,distcol,distrow,distslice, interpolationcode);
 *
 * MRIS *surf = MRISread('lh.white');
 * for n = 0:surf->nvertices-1
 *  VERTEX *v = surf->vertices[n]
 *     v->x, v->y, v->z // by default these are in the warped space
 *     tkRAS->rptr[1][1] = v->x;
 *     tkRAS->rptr[2][1] = v->y;
 *     tkRAS->rptr[3][1] = v->z;
 *     MATRIX *M = TkrRAS2VoxfromVolGeom(&surf->vg); // convert from tkreg space to voxel
 *     MATRIX *V = vg_getVoxelToRasXform(&surf->vg); // converts from voxel to RAS
 *     MATRIX *Q = MatrixMultiply(V,M,NULL); // convert from tkreg space to RAS
 *     DistortedRAS = MatrixMultiply(Q,tkRAS,DistortedCRS);
       spharm_evaluate(Distx, Disty, Distz, &Dx, &Dy, &Dz);
       v->x +/- Dx;
       v->y +/- Dy;
       v->z +/- Dz;
 */

static void dump_exit_codes(FILE *fp);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
static int isflag(char *flag);
static int nth_is_arg(int nargc, char **argv, int nth);

static void unwarpCRS(GradUnwarp *gradUnwarp, MATRIX *vox2ras, MATRIX *inv_vox2ras);

int debug = 0, checkoptsonly = 0;
const char *Progname = NULL;
std::string gradfile_str, origmgz_str, unwarpedmgz_str, loadtrans_str, savetrans_str;
const char *gradfile = NULL, *origmgz = NULL, *unwarpedmgz = NULL, *loadtrans = NULL, *savetrans = NULL;
int inputras = 0, inputcrs = 0;
double ras_x, ras_y, ras_z;
int crs_c = 0, crs_r = 0, crs_s = 0;
const char *interpmethod = "trilinear";
int   interpcode = -1;
int   sinchw = 0;
int   nthreads = 1;

int main(int argc, char *argv[])
{
  int nargs;
  struct utsname uts;
  char *cmdline, cwd[2000];

  interpcode = MRIinterpCode(interpmethod);

  nargs = handleVersionOption(argc, argv, "mris_gradient_unwarp");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  printf("\n");
  printf("%s\n",getVersion().c_str());
  printf("cwd %s\n",cwd);
  printf("cmdline %s\n",cmdline);
  printf("sysname  %s\n",uts.sysname);
  printf("hostname %s\n",uts.nodename);
  printf("machine  %s\n",uts.machine);

  dump_options(stdout);

  if (gradfile == NULL && loadtrans == NULL)
  {
    printf("Use --gradcoeff to specify gradient coeff file, \n");
    printf(" or --load_transtbl to specify unwarp transformation table \n");
    return 0; 
  }

  if (gradfile != NULL && loadtrans != NULL)
  {
    printf("--gradcoeff and --load_transtbl are mutually excluded. \n");
    printf("Use --gradcoeff to specify gradient coeff file, \n");
    printf(" or --load_transtbl to specify unwarp transformation table \n");
    return 0; 
  }

  if (origmgz == NULL)
  {
    printf("Use --invol to specify input volumn\n");
    return 0;
  }

  MRI *origvol = MRIread(origmgz);

  printf("\nVolume information for %s\n", origmgz);
  if (origvol->nframes > 1)
    printf("Dimensions: %d x %d x %d x %d\n", origvol->width, origvol->height, origvol->depth, origvol->nframes) ;
  else
    printf("Dimensions: %d x %d x %d\n", origvol->width, origvol->height, origvol->depth) ;
  printf("Data Type: %s (%d)\n",
         origvol->type == MRI_UCHAR   ? "UCHAR" :
         origvol->type == MRI_SHORT   ? "SHORT" :
	 origvol->type == MRI_USHRT   ? "USHRT" :
         origvol->type == MRI_INT     ? "INT" :
         origvol->type == MRI_LONG    ? "LONG" :
         origvol->type == MRI_BITMAP  ? "BITMAP" :
         origvol->type == MRI_TENSOR  ? "TENSOR" :
         origvol->type == MRI_FLOAT   ? "FLOAT" : "UNKNOWN", origvol->type) ;

  MATRIX *vox2ras_orig = extract_i_to_r(origvol);  //MRIxfmCRS2XYZ(origvol, 0);
  printf("\nvoxel to ras transform:\n"); 
  MatrixPrint(stdout, vox2ras_orig);

  MATRIX *inv_vox2ras_orig = MatrixInverse(vox2ras_orig, NULL);  //extract_r_to_i(origvol);
  printf("\nras to voxel transform:\n");
  MatrixPrint(stdout, inv_vox2ras_orig);

  MRI *unwarpedvol = NULL;
  if (unwarpedmgz != NULL)
  {
    unwarpedvol = MRIallocSequence(origvol->width, origvol->height, origvol->depth, MRI_FLOAT, origvol->nframes);
    MRIcopyHeader(origvol, unwarpedvol);
    MRIcopyPulseParameters(origvol, unwarpedvol);
  }

  GradUnwarp *gradUnwarp = new GradUnwarp(nthreads);
  if (gradfile != NULL)
  {
    gradUnwarp->read_siemens_coeff(gradfile);
    gradUnwarp->initSiemensLegendreNormfact();
    if (getenv("PRN_GRADCOEFF_ONLY"))
    {
      gradUnwarp->printCoeff();
      delete gradUnwarp;

      return 0;
    }
    else if (inputras) // for debugging only
    {
      float Dx, Dy, Dz;
      gradUnwarp->spharm_evaluate(ras_x, ras_y, ras_z, &Dx, &Dy, &Dz);

      printf(" x = %.6f,  Dx = %.6lf\n", ras_x, Dx);
      printf(" y = %.6f,  Dy = %.6lf\n", ras_y, Dy);
      printf(" z = %.6f,  Dz = %.6lf\n", ras_z, Dz);

      delete gradUnwarp;

      return 0;
    }
    else if (inputcrs) // for debugging
    {
      unwarpCRS(gradUnwarp, vox2ras_orig, inv_vox2ras_orig);
      return 0;
    }

    gradUnwarp->create_transtable(origvol, unwarpedvol, vox2ras_orig, inv_vox2ras_orig);
  }
  else
  {
    gradUnwarp->load_transtable(loadtrans);
  }

  if (unwarpedvol != NULL)
  {
    gradUnwarp->unwarp_volume(origvol, unwarpedvol, interpcode, sinchw);

    printf("Writing to %s\n", unwarpedmgz);
    int err = MRIwrite(unwarpedvol, unwarpedmgz);
    if (err) 
      printf("something went wrong\n");

    MRIfree(&unwarpedvol);
  }

  if (savetrans != NULL)
    gradUnwarp->save_transtable(savetrans);

  delete gradUnwarp;

  MatrixFree(&vox2ras_orig);
  MatrixFree(&inv_vox2ras_orig);

  MRIfree(&origvol);

  return 0;
}

/* --------------------------------------------- */
static void unwarpCRS(GradUnwarp *gradUnwarp, MATRIX *vox2ras, MATRIX *inv_vox2ras)
{
  int (*nintfunc)( double );
  nintfunc = &nint;

  MATRIX *CRS = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *RAS = MatrixAlloc(4, 1, MATRIX_REAL);;
  MATRIX *DeltaRAS = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *DistortedRAS = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *DistortedCRS = MatrixAlloc(4, 1, MATRIX_REAL);

  CRS->rptr[1][1] = crs_c;
  CRS->rptr[2][1] = crs_r;
  CRS->rptr[3][1] = crs_s;
  CRS->rptr[4][1] = 1;          

  // Convert the CRS to RAS
  RAS->rptr[4][1] = 1;
  RAS = MatrixMultiply(vox2ras, CRS, RAS);

  double x = RAS->rptr[1][1];
  double y = RAS->rptr[2][1];
  double z = RAS->rptr[3][1];
  float Dx = 0, Dy = 0, Dz = 0;
  gradUnwarp->spharm_evaluate(x, y, z, &Dx, &Dy, &Dz);

  DeltaRAS->rptr[1][1] = Dx;
  DeltaRAS->rptr[2][1] = Dy;
  DeltaRAS->rptr[3][1] = Dz;
  DeltaRAS->rptr[4][1] = 1; 
        
  DistortedRAS = MatrixAdd(RAS, DeltaRAS, DistortedRAS);
  DistortedRAS->rptr[4][1] = 1;
  DistortedCRS = MatrixMultiply(inv_vox2ras, DistortedRAS, DistortedCRS);
	  
  printf("\n");
  printf("CRS:\n");
  MatrixPrint(stdout, CRS);

  printf("RAS:\n");
  MatrixPrint(stdout, RAS);

  printf("DeltaRAS:\n");
  MatrixPrint(stdout, DeltaRAS);

  printf("DistortedRAS:\n");
  MatrixPrint(stdout, DistortedRAS);

  printf("DistortedCRS:\n");
  MatrixPrint(stdout, DistortedCRS);

  printf("\n c=%d r=%d s=%d\n", crs_c, crs_r, crs_s);
  printf(" x=%.6f  y=%.6f  z=%.6f\n", x, y, z);
  printf("Dx=%.6lf Dy=%.6lf Dz=%.6lf\n", Dx, Dy, Dz);
 
  float fcs = DistortedCRS->rptr[1][1];
  float frs = DistortedCRS->rptr[2][1];
  float fss = DistortedCRS->rptr[3][1];

  int ics =  nintfunc(fcs);
  int irs =  nintfunc(frs);
  int iss =  nintfunc(fss);

  printf("fcs = %lf (%d), frs = %lf (%d), fss = %lf (%d)\n", fcs, ics, frs, irs, fss, iss);

  MatrixFree(&CRS);
  MatrixFree(&RAS);
  MatrixFree(&DeltaRAS);
  MatrixFree(&DistortedRAS);
  MatrixFree(&DistortedCRS);
}

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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcmp(option, "--ras")) {
      inputras = 1;
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%lf,%lf,%lf", &ras_x, &ras_y, &ras_z);
      nargsused = 1;
    } else if (!strcmp(option, "--crs")) {
      inputcrs = 1;
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d,%d,%d", &crs_c, &crs_r, &crs_s);
      nargsused = 1;
    } else if (!strcmp(option, "--gradcoeff")) {
      if (nargc < 1) CMDargNErr(option,1);
      gradfile_str = fio_fullpath(pargv[0]);
      gradfile = gradfile_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--invol")) {
      if (nargc < 1) CMDargNErr(option,1);
      origmgz_str = fio_fullpath(pargv[0]);
      origmgz = origmgz_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--outvol")) {
      if (nargc < 1) CMDargNErr(option,1);
      unwarpedmgz_str = fio_fullpath(pargv[0]);
      unwarpedmgz = unwarpedmgz_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--load_transtbl")) {
      if (nargc < 1) CMDargNErr(option,1);
      loadtrans_str = fio_fullpath(pargv[0]);
      loadtrans = loadtrans_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--save_transtbl")) {
      if (nargc < 1) CMDargNErr(option,1);
      savetrans_str = fio_fullpath(pargv[0]);
      savetrans = savetrans_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--interp")) {
      if (nargc < 1) CMDargNErr(option, 1);
      interpmethod = pargv[0];
      interpcode = MRIinterpCode(interpmethod);
      printf("!!! <interpmethod=%s, interpcode=%d> !!!\n", interpmethod, interpcode);
      nargsused = 1;
      if (!strcmp(interpmethod, "sinc") && nth_is_arg(nargc, pargv, 1)) {
        sscanf(pargv[1], "%d", &sinchw);
        printf("!!! <sinchw = %d> !!!\n", sinchw);
        nargsused ++;
      }
    } else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d", &nthreads);
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(1);
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
  printf("   --gradcoeff     gradient coeff input \n");
  printf("   --load_transtbl load unwarp transform table in m3z format \n");
  printf("   --invol         input volume \n");
  printf("   --outvol        unwarped output volume \n");
  printf("   --save_transtbl save unwarp transform table in m3z format \n");
  printf("\n");
  printf("   --ras       x,y,z\n");
  printf("   --crs       c,r,s\n");
  printf("\n"); 
  printf("   --interp    interptype nearest | trilinear | sinc | cubic\n");
  printf("\n"); 
  printf("   --nthreads  nthreads : Set OPEN MP threads\n");
  printf("\n"); 
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "program description\n"
    "...\n"
    "...\n");

  printf("\n");
  dump_exit_codes(stdout);
  printf("\n");

  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void) {
  if (interpcode < 0) {
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilinear, sinc, and cubic\n");
    exit(1);
  }
  /*if (vol1File.size() == 0) {
    printf("ERROR: must specify a vol1 file\n");
    exit(1);
  }
  if (vol2File.size() == 0) {
    printf("ERROR: must specify a vol2 file\n");
    exit(1);
  }
  return;
  */
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  /*fprintf(fp,"vol1    %s\n",vol1File.c_str());
  fprintf(fp,"vol2    %s\n",vol2File.c_str());
  fprintf(fp,"pix thresh  %g\n",pixdiff_thresh);
  fprintf(fp,"vox2ras thresh %g\n",vox2ras_thresh);
  return;
  */
}

/* --------------------------------------------- */
static void dump_exit_codes(FILE *fp) {
  /*  fprintf(fp,"dimensions inconsistent   %d\n",DIMENSION_EC);
  fprintf(fp,"precision  inconsistent   %d\n",PRECISION_EC);
  fprintf(fp,"resolution inconsistent   %d\n",RESOLUTION_EC);
  fprintf(fp,"vox2ras    inconsistent   %d\n",VOX2RAS_EC);
  fprintf(fp,"pixel      inconsistent   %d\n",PIXEL_EC);
  */
}

/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}

/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if (nargc <= nth) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
}
