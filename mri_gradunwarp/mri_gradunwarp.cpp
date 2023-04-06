#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "version.h"
#include "diag.h"
#include "log.h"
#include "cmdargs.h"
#include "fio.h"
#include "mri.h"
#include "mri_identify.h"
#include "mriBSpline.h"
#include "GradUnwarp.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static int isflag(char *flag);
static int nth_is_arg(int nargc, char **argv, int nth);

static void printVolInfo(MRI *vol, const VOL_GEOM *vg, MATRIX *vox2ras_orig, MATRIX *inv_vox2ras_orig);
static void printVolGeom(const VOL_GEOM *vg, MATRIX *vox2ras_orig, MATRIX *inv_vox2ras_orig);
static void unwarpCRS(GradUnwarp *gradUnwarp, MATRIX *vox2ras, MATRIX *inv_vox2ras);

int checkoptsonly = 0;
const char *Progname = NULL;
std::string gradfile_str, inf_str, outf_str, loadtrans_str, outtrans_str, invgcamfile_str;
const char *gradfile = NULL, *inf = NULL, *outf = NULL, *loadtrans = NULL, *outtrans = NULL, *invgcamfile=NULL;
int inputras = 0, inputcrs = 0, unwarp = 0, m3zonly = 0;
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

  nargs = handleVersionOption(argc, argv, "mri_gradunwarp");
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

  MRI *origvol = NULL;
  MRIS *origsurf = NULL;
  VOL_GEOM vg;
  MATRIX *vox2ras_orig = NULL;
  MATRIX *inv_vox2ras_orig = NULL;

  int unwarpvol = 0, unwarpsurf = 0;
  if (inf != NULL)
  {
    // determine if input is a volume or surface
    unwarpvol = 1; unwarpsurf = 0;

    int voltype = mri_identify(inf);
    if (voltype == MRI_VOLUME_TYPE_UNKNOWN)
    {
      unwarpvol = 0; unwarpsurf = 1;
    }
  }

  if (unwarpvol)
  {
    origvol = MRIread(inf);
    if (origvol == NULL)
    {
      printf("ERROR: could not open volume %s\n", inf);
      exit(1);
    }

    if (m3zonly)
      unwarpvol = 0;

    getVolGeom(origvol, &vg);
  }
  else if (unwarpsurf)
  {
    origsurf = MRISread(inf);
    if (origsurf == NULL)
    {
      printf("ERROR: could not open surface %s\n", inf);
      exit(1);
    }

    if (m3zonly)
      unwarpsurf = 0;

    copyVolGeom(&origsurf->vg, &vg);
  }  

  vox2ras_orig     = vg_i_to_r(&vg);
  inv_vox2ras_orig = vg_r_to_i(&vg);

  printVolGeom(&vg, vox2ras_orig, inv_vox2ras_orig); 


  GradUnwarp *gradUnwarp = new GradUnwarp(nthreads);
  if (gradfile != NULL)
  {
    gradUnwarp->read_siemens_coeff(gradfile);
    gradUnwarp->initSiemensLegendreNormfact();
    if (getenv("GRADUNWARP_PRN_GRADCOEFF_ONLY"))
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

    if (getenv("GRADUNWARP_USE_GRADFILE") == NULL)
      gradUnwarp->create_transtable(&vg, vox2ras_orig, inv_vox2ras_orig);
  }
  else
  {
    gradUnwarp->load_transtable(loadtrans);
  }

  if (unwarpvol)
  {
    MRI *unwarpedvol= NULL;

    if (getenv("GRADUNWARP_USE_GRADFILE"))
    {
      printf("****** env GRADUNWARP_USE_GRADFILE set. ******\n");

      unwarpedvol = MRIallocSequence(origvol->width, origvol->height, origvol->depth, origvol->type, origvol->nframes);
      MRIcopyHeader(origvol, unwarpedvol);
      MRIcopyPulseParameters(origvol, unwarpedvol);

      gradUnwarp->unwarp_volume_gradfile(origvol, unwarpedvol, vox2ras_orig, inv_vox2ras_orig, interpcode, sinchw);
    }
    else
    {
      unwarpedvol = gradUnwarp->unwarp_volume(origvol, NULL, interpcode, sinchw);
    }

    printf("Writing to %s\n", outf);
    int err = MRIwrite(unwarpedvol, outf);
    if(err) exit(1);

    MRIfree(&unwarpedvol);

    MatrixFree(&vox2ras_orig);
    MatrixFree(&inv_vox2ras_orig);

    MRIfree(&origvol);
  }
  else if (unwarpsurf)
  {
    // origsurf will be updated in unwarp_surface() and unwarp_surface_gradfile()
    if (getenv("GRADUNWARP_USE_GRADFILE"))
    {
      printf("****** env GRADUNWARP_USE_GRADFILE set. ******\n");
      gradUnwarp->unwarp_surface_gradfile(origsurf, origsurf);
    }
    else
    {
      gradUnwarp->unwarp_surface(origsurf, origsurf);
    }

    //origsurf->vg.copy(&unwarpedsurf->vg);

    printf("Writing to %s\n", outf);
    int err = MRISwrite(origsurf, outf);
    if (err) 
      printf("Erros outputing surface %s\n", outf);

    MRISfree(&origsurf);
  }
  
  if(outtrans != NULL) {
    // this also covers save_transtbl_only 
    // if we used gradient file to unwarp, we need to create the m3z transform table now
    if(getenv("GRADUNWARP_USE_GRADFILE"))
      gradUnwarp->create_transtable(&vg, vox2ras_orig, inv_vox2ras_orig);
    gradUnwarp->save_transtable(outtrans);
  }

  if(0 && invgcamfile != NULL) {
    // This does not actually invert the gcam, which is what I wanted it to do
    if(outtrans == NULL) {
      if(getenv("GRADUNWARP_USE_GRADFILE"))
	gradUnwarp->create_transtable(&vg, vox2ras_orig, inv_vox2ras_orig);
    }
    gradUnwarp->invert_gcam(unwarpedvol);
    gradUnwarp->save_transtable(invgcamfile);
    // WARNING: unwarpedvol will now be warped!!
  }

  delete gradUnwarp;

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
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    /* these two options --ras, --crs are for debugging purposes only */
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
    }
    else if (!strcmp(option, "--gradcoeff")) {
      if (nargc < 1) CMDargNErr(option,1);

      struct stat stat_buf;
      if (stat(pargv[0], &stat_buf) < 0)
      {
        printf("ERROR: could not find --gradcoeff %s\n", pargv[0]);
        exit(1);
      }

      gradfile_str = fio_fullpath(pargv[0]);
      gradfile = gradfile_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);

      struct stat stat_buf;
      if (stat(pargv[0], &stat_buf) < 0)
      {
        printf("ERROR: could not find --i %s\n", pargv[0]);
        exit(1);
      }

      inf_str = fio_fullpath(pargv[0]);
      inf = inf_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outf_str = fio_fullpath(pargv[0]);
      outf = outf_str.c_str();
      nargsused = 1;
    } else if (!strcmp(option, "--load_transtbl")) {
      if (nargc < 1) CMDargNErr(option,1);

      struct stat stat_buf;
      if (stat(pargv[0], &stat_buf) < 0)
      {
        printf("ERROR: could not find --load_transtbl %s\n", pargv[0]);
        exit(1);
      }

      loadtrans_str = fio_fullpath(pargv[0]);
      loadtrans = loadtrans_str.c_str();
      nargsused = 1;
    } 
    else if(!strcmp(option, "--out_transtbl") || !strcmp(option, "--gcam")) {
      if (nargc < 1) CMDargNErr(option,1);
      outtrans_str = fio_fullpath(pargv[0]);
      outtrans = outtrans_str.c_str();
      nargsused = 1;
    } 
    else if (!strcmp(option, "--inv-gcam")) {
      if(nargc < 1) CMDargNErr(option,1);
      invgcamfile_str = fio_fullpath(pargv[0]);
      invgcamfile= invgcamfile_str.c_str();
      nargsused = 1;
    } 
    else if(!strcmp(option, "--save_transtbl_only") || !strcmp(option, "--gcam-only")) {
      m3zonly = 1; unwarp = 0;
    } 
    else if (!strcmp(option, "--interp")) {
      if (nargc < 1) CMDargNErr(option, 1);
      interpmethod = pargv[0];
      interpcode = MRIinterpCode(interpmethod);
      nargsused = 1;
      /*printf("!!! <interpmethod=%s, interpcode=%d> !!!\n", interpmethod, interpcode);
      if (!strcmp(interpmethod, "sinc") && nth_is_arg(nargc, pargv, 1)) {
        sscanf(pargv[1], "%d", &sinchw);
        printf("!!! <sinchw = %d> !!!\n", sinchw);
        nargsused ++;
	}*/
    } else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d", &nthreads);
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      print_help();
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
  printf("\n");
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --gradcoeff <gradient-file>  gradient coeff input (not with --load_transtbl)\n");
  printf("   --load_transtbl <m3z-table>  load unwarp transform table in m3z format (not with --gradcoeff)\n");
  printf("   --i <input-warped-file>       input volume, or surface \n");
  printf("   --o <output-unwarped-file>    unwarped output volume, or surface \n");
  printf("   --out_transtbl <output-m3z-table>  save unwarp transform table in m3z format (or --gcam)\n");
  printf("   --save_transtbl_only   just save unwarp transform table in m3z format, need --gradcoeff <> (or --gcam-only)\n");
  //printf("   --inv-gcam invgcam.m3z : save inverse of m3z\n"); // does not work the way I wanted it to
  printf("\n");
  //printf("   --ras       x,y,z\n");
  //printf("   --crs       c,r,s\n");
  //printf("\n"); 
  printf("   --interp <interptype>                                       interpolation method: nearest | trilinear | cubic (default is trilinear)\n");
  printf("\n"); 
  printf("   --nthreads <nthreads>                                       number of threads to run\n");
  printf("\n"); 
  printf("   --checkopts                                                 don't run anything, just check options and exit\n");
  printf("   --help                                                      print out information on how to use this program\n");
  printf("   --version                                                   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage();

  printf("\n%s\n\n",getVersion().c_str());
  printf("\n");
  printf("This program provides a tool to correct gradient non-linearity distortions in MRI images.\n");
  printf("\n");
  printf("Examples:\n");
  printf("\n");
  printf("1. create gradient unwarp m3z transformation table only from a volume using gradient file:\n");
  printf("  mri_gradunwarp \n"); 
  printf("    --gradcoeff coeff_Sonata.grad \n");
  printf("    --i invol.mgz \n");
  printf("    --out_transtbl gradunwarp.m3z \n");
  printf("    --save_transtbl_only \n");
  printf("    --nthreads 10 \n");
  printf("\n");
  printf("2. unwarp given volume using gradient file, save gradient unwarp m3z transformation table:\n");
  printf("  mri_gradunwarp \n");
  printf("    --gradcoeff coeff_Sonata.grad \n");
  printf("    --i invol.mgz \n");
  printf("    --o invol.unwarped.nearest.mgz --interp nearest\n");
  printf("    --out_transtbl gradunwarp.m3z \n");
  printf("    --nthreads 10 \n");
  printf("\n");
  printf("3. unwarp given surface using gradient unwarp m3z transformation table:\n");
  printf("  mri_gradunwarp \n");
  printf("    --load_transtbl gradunwarp.m3z \n");
  printf("    --i lh.white \n");
  printf("    --o lh.unwarped.nearest.white --interp nearest\n");
  printf("    --nthreads 10 \n");
  printf("\n");

#if 0
  printf("* unwarp at given crs (debug only): \n");
  printf("  mri_gradunwarp \n");
  printf("    --gradcoeff coeff_Sonata.grad \n");
  printf("    --i orig.mgz \n");
  printf("    --crs 0,0,0 \n");
  printf("\n");
  printf("* unwarp for given ras (debug only): \n");
  printf("  mri_gradunwarp \n");
  printf("    --gradcoeff coeff_Sonata.grad \n");
  printf("    --i orig.mgz \n");
  printf("    --ras 0.1,0.2,0.3 \n");
#endif

  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void)
{
  // debug options --crs and --ras need --gradfile <> 
  if ((m3zonly || inputcrs || inputras) && gradfile == NULL)
  {
    printf("Use --gradcoeff to specify gradient coeff file, \n");
    print_usage();
    exit(1);
  }

  if (getenv("GRADUNWARP_USE_GRADFILE") && gradfile == NULL)
  {
    printf("Environment variable GRADUNWARP_USE_GRADFILE is set.\n");
    printf("Use --gradcoeff to specify gradient coeff file\n");
    print_usage();
    exit(1);
  }

  if (gradfile == NULL && loadtrans == NULL)
  {
    printf("Use --gradcoeff to specify gradient coeff file, \n");
    printf(" or --load_transtbl to specify unwarp transformation table \n");
    print_usage();
    exit(1); 
  }
  
  if (gradfile != NULL && loadtrans != NULL)
  {
    if (m3zonly)
      printf("--gradcoeff and --load_transtbl are mutually excluded. \n");
    else
    {
      printf("--gradcoeff and --load_transtbl are mutually excluded. \n");
      printf("Use --gradcoeff to specify gradient coeff file, \n");
      printf(" or --load_transtbl to specify unwarp transformation table \n");
    }

    print_usage();
    exit(1); 
  }

  if (inf == NULL)
  {
    printf("Use --i to specify input volume or surface\n");
    print_usage();
    exit(1);
  }

  if (m3zonly && outtrans == NULL)
  {
    printf("Use --out_transtbl to specify output m3z table\n");
    print_usage();
    exit(1);
  }

  if (!m3zonly && !inputcrs && !inputras)
    unwarp = 1;

  // need unwarped output file
  if (unwarp && outf == NULL)
  {
    printf("Use --o to specify output unwarped output\n");
    print_usage();
    exit(1);
  }   

  if (interpcode != SAMPLE_NEAREST && interpcode != SAMPLE_TRILINEAR && interpcode != SAMPLE_CUBIC_BSPLINE) 
  {
    printf("ERROR: interpolation method %s unrecognized\n", interpmethod);
    printf("       legal values are nearest, trilinear, and cubic.\n");
    print_usage();
    exit(1);
  }
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

static void printVolInfo(MRI *vol, const VOL_GEOM *vg, MATRIX *vox2ras_orig, MATRIX *inv_vox2ras_orig)
{
  printf("\nVolume information for %s\n", inf);
  if (vol->nframes > 1)
    printf("Dimensions: %d x %d x %d x %d\n", vol->width, vol->height, vol->depth, vol->nframes) ;
  else
    printf("Dimensions: %d x %d x %d\n", vol->width, vol->height, vol->depth) ;

  printf("Data Type: %s (%d)\n",
         vol->type == MRI_UCHAR   ? "UCHAR"  :
         vol->type == MRI_SHORT   ? "SHORT"  :
	 vol->type == MRI_USHRT   ? "USHRT"  :
         vol->type == MRI_INT     ? "INT"    :
         vol->type == MRI_LONG    ? "LONG"   :
         vol->type == MRI_BITMAP  ? "BITMAP" :
         vol->type == MRI_TENSOR  ? "TENSOR" :
         vol->type == MRI_FLOAT   ? "FLOAT"  : "UNKNOWN", vol->type) ;

  printVolGeom(vg, vox2ras_orig, inv_vox2ras_orig);
}

static void printVolGeom(const VOL_GEOM *vg, MATRIX *vox2ras_orig, MATRIX *inv_vox2ras_orig)
{
  printf("volume geometry:\n");
  printf("extent  : (%d, %d, %d)\n", vg->width, vg->height, vg->depth);
  printf("voxel   : (%7.9f, %7.9f, %7.9f)\n", vg->xsize, vg->ysize, vg->zsize);
  printf("x_(ras) : (%7.9f, %7.9f, %7.9f)\n", vg->x_r, vg->x_a, vg->x_s);
  printf("y_(ras) : (%7.9f, %7.9f, %7.9f)\n", vg->y_r, vg->y_a, vg->y_s);
  printf("z_(ras) : (%7.9f, %7.9f, %7.9f)\n", vg->z_r, vg->z_a, vg->z_s);
  printf("c_(ras) : (%7.9f, %7.9f, %7.9f)\n", vg->c_r, vg->c_a, vg->c_s);
  printf("file    : %s\n", vg->fname);

  printf("\nvoxel to ras transform:\n"); 
  MatrixPrint(stdout, vox2ras_orig);

  printf("\nras to voxel transform:\n");
  MatrixPrint(stdout, inv_vox2ras_orig);
}
