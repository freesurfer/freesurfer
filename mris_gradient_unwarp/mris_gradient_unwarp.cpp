#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>

#include "romp_support.h"

#include "version.h"
#include "diag.h"
#include "log.h"
#include "cmdargs.h"
#include "fio.h"
#include "mri.h"
#include "mriBSpline.h"
#include "GradUnwarp.h"

/* examples:
 * 1. mris_gradient_unwarp/mris_gradient_unwarp --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad --invol $FS_TEST/gradunwarp/example/orig.mgz --outvol unwarped.mgz
 * 2. mris_gradient_unwarp/mris_gradient_unwarp --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad --invol $FS_TEST/gradunwarp/example/orig.mgz --outvol orig.dup.mgz --dupinvol
 * 3. mris_gradient_unwarp/mris_gradient_unwarp --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad --invol $FS_TEST/gradunwarp/example/orig.mgz --ras 0.1,0.2,0.3
 * 4. mris_gradient_unwarp/mris_gradient_unwarp --gradcoeff $FS_TEST/gradunwarp/example/coeff_Sonata.grad --invol $FS_TEST/gradunwarp/example/orig.mgz  --crs 0,0,0
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
 * surf = MRISread('lh.white');
 * for n = 0:surf->nvertices-1
 *  VERTEX *v = surf->vertices[n]
 *  v->x, v->y, v->z
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

static void unwarpVol(const char* origmgz, const char* unwarpedmgz, GradUnwarp* gradUnwarp);

int debug = 0, checkoptsonly = 0;
const char *Progname = NULL;
std::string gradfile0, origmgz0, unwarpedmgz0;
int dupinvol = 0, inputras = 0, inputcrs = 0;
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

  GradUnwarp *gradUnwarp = new GradUnwarp();
  gradUnwarp->read_siemens_coeff(gradfile0.c_str());
  if (getenv("PRN_GRADCOEFF_ONLY"))
  {
    gradUnwarp->printCoeff();
    return 0;
  }

  gradUnwarp->initSiemensLegendreNormfact();

  if (inputras) // for debugging only
  {
    float Dx, Dy, Dz;
    gradUnwarp->spharm_evaluate(ras_x, ras_y, ras_z, &Dx, &Dy, &Dz);

    printf(" x = %.6f,  Dx = %.6lf\n", ras_x, Dx);
    printf(" y = %.6f,  Dy = %.6lf\n", ras_y, Dy);
    printf(" z = %.6f,  Dz = %.6lf\n", ras_z, Dz);

    delete gradUnwarp;

    return 0;
  }

  unwarpVol(origmgz0.c_str(), unwarpedmgz0.c_str(), gradUnwarp); 

  delete gradUnwarp;

  return 0;
}

/* --------------------------------------------- */
static void unwarpVol(const char* origmgz, const char* unwarpedmgz, GradUnwarp* gradUnwarp)
{
  int (*nintfunc)( double );
  nintfunc = &nint;

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

  MRI *unwarpedvol = MRIallocSequence(origvol->width, origvol->height, origvol->depth, MRI_FLOAT, origvol->nframes);
  MRIcopyHeader(origvol, unwarpedvol);
  MRIcopyPulseParameters(origvol, unwarpedvol);

  MRI_BSPLINE *bspline = NULL;
  if (interpcode == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(origvol, NULL, 3);

#ifdef HAVE_OPENMP
  printf("\nSet OPEN MP NUM threads to %d\n", nthreads);
  omp_set_num_threads(nthreads);
#endif

  int c; 
  int outofrange_total = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+ : outofrange_total) 
#endif
  for (c = 0; c < origvol->width; c++)
  {
    // You could make a vector of CRS nthreads long
    MATRIX *CRS = MatrixAlloc(4, 1, MATRIX_REAL);
    MATRIX *RAS = MatrixAlloc(4, 1, MATRIX_REAL);;
    MATRIX *DeltaRAS = MatrixAlloc(4, 1, MATRIX_REAL);
    MATRIX *DistortedRAS = MatrixAlloc(4, 1, MATRIX_REAL);
    MATRIX *DistortedCRS = MatrixAlloc(4, 1, MATRIX_REAL);

    int r = 0, s = 0;
    //int outofrange_local = 0;

#if 0
#ifdef HAVE_OPENMP
    int tid = omp_get_thread_num();
    printf("hello from thread #%d (%d, %d, %d) ...\n", tid, c, r, s);
#else
    printf("hello from process %d ...\n", getpid());
#endif
#endif

    for (r = 0; r < origvol->height; r++)
    {
      for (s = 0; s < origvol->depth; s++)
      {
        if (!dupinvol)
        {
          // clear CRS, RAS, DeltaRAS, DistortedRAS, DistortedCRS
          MatrixClear(CRS);
          MatrixClear(RAS);
          MatrixClear(DeltaRAS);
          MatrixClear(DistortedRAS);
          MatrixClear(DistortedCRS);

          if (inputcrs)
          {
            CRS->rptr[1][1] = crs_c;
            CRS->rptr[2][1] = crs_r;
            CRS->rptr[3][1] = crs_s;
            CRS->rptr[4][1] = 1;          
          }
          else
	  {
            CRS->rptr[1][1] = c;
            CRS->rptr[2][1] = r;
            CRS->rptr[3][1] = s;
            CRS->rptr[4][1] = 1;
	  }

          // Convert the CRS to RAS
          RAS->rptr[4][1] = 1;
          RAS = MatrixMultiply(vox2ras_orig, CRS, RAS);

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
          DistortedCRS = MatrixMultiply(inv_vox2ras_orig, DistortedRAS, DistortedCRS);
	  
          if (inputcrs)
	  {
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
	  }
        } 
 
        float fcs = 0, frs = 0, fss = 0;
        int ics = 0, irs = 0, iss = 0;

        if (!dupinvol)
	{
          fcs = DistortedCRS->rptr[1][1];
          frs = DistortedCRS->rptr[2][1];
          fss = DistortedCRS->rptr[3][1];

          ics =  nintfunc(fcs);
          irs =  nintfunc(frs);
          iss =  nintfunc(fss);

          if (inputcrs)
	  {
            printf("fcs = %lf (%d), frs = %lf (%d), fss = %lf (%d)\n", fcs, ics, frs, irs, fss, iss);
            exit(0);
	  }

          if (ics < 0 || ics >= origvol->width  ||
              irs < 0 || irs >= origvol->height || 
              iss < 0 || iss >= origvol->depth)
          {
            outofrange_total++;
#if 0
#ifdef HAVE_OPENMP
            outofrange_local++;
#else
            outofrange_total++;
#endif
#endif
            continue;
          }
	}
        else
	{
          // this should output the same vol as orig
          fcs =  c; frs =  r; fss =  s;
          ics =  c; irs =  r; iss =  s;
	}

        /* Assign output volume values */
        if (interpcode == SAMPLE_TRILINEAR) {
          float *valvect = new float[origvol->nframes]; 
          MRIsampleSeqVolume(origvol, fcs, frs, fss, valvect, 0, origvol->nframes - 1);
          
          int f;
          for (f = 0; f < origvol->nframes; f++) 
            MRIsetVoxVal2(unwarpedvol, c, r, s, f, valvect[f]);

          free(valvect);
	} else {
          double rval = 0;

          int f;
          for (f = 0; f < origvol->nframes; f++) {
            switch (interpcode) {
              case SAMPLE_NEAREST:
                rval = MRIgetVoxVal2(origvol, ics, irs, iss, f);
                break;
              case SAMPLE_CUBIC_BSPLINE:
                MRIsampleBSpline(bspline, fcs, frs, fss, f, &rval);
                break;
              case SAMPLE_SINC: /* no multi-frame */
                MRIsincSampleVolume(origvol, fcs, frs, fss, sinchw, &rval);
                break;
              default:
                printf("ERROR: MR: interpolation method '%i' unknown\n", interpcode);
                exit(1);
            } // switch

            MRIsetVoxVal2(unwarpedvol, c, r, s, f, rval);
          } // f
        }
      }   // s
    }     // r

#if 0
#ifdef HAVE_OPENMP
#pragma omp critical 
    outofrange_total += outofrange_local; 
    //printf("update out of range voxel count: + %d = %d\n", outofrange_local, outofrange_total);
#endif
#endif

    MatrixFree(&CRS);
    MatrixFree(&RAS);
    MatrixFree(&DeltaRAS);
    MatrixFree(&DistortedRAS);
    MatrixFree(&DistortedCRS);
  }       // c

  printf("Writing to %s\n",unwarpedmgz);
  int err = MRIwrite(unwarpedvol, unwarpedmgz);
  if(err) printf("something went wrong\n");

  printf("Total %d voxels are out of range\n", outofrange_total);

  MatrixFree(&vox2ras_orig);
  MatrixFree(&inv_vox2ras_orig);

  MRIfree(&origvol);
  MRIfree(&unwarpedvol);

  if (bspline)
    MRIfreeBSpline(&bspline);
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

    else if (!strcasecmp(option, "--dupinvol"))   dupinvol = 1;
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
      gradfile0 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--invol")) {
      if (nargc < 1) CMDargNErr(option,1);
      origmgz0 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--outvol")) {
      if (nargc < 1) CMDargNErr(option,1);
      unwarpedmgz0 = fio_fullpath(pargv[0]);
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
  printf("   --gradcoeff gradient coeff input \n");
  printf("   --invol     input volume \n");
  printf("   --outvol    unwarped outout volume \n");
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
