////////////////////////////////////////////////////////////////////
// mri_info.c
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: greve $
// Revision Date  : $Date: 2004/10/04 22:18:50 $
// Revision       : $Revision: 1.32 $
//
////////////////////////////////////////////////////////////////////
char *MRI_INFO_VERSION = "$Revision: 1.32 $";
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "const.h"
#include "machine.h"
#include "fio.h"
#include "utils.h"
#include "mri.h"
#include "volume_io.h"
#include "analyze.h"
#include "mri_identify.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "mghendian.h"
#include "fio.h"

static void do_file(char *fname);
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;

static char vcid[] = "$Id: mri_info.c,v 1.32 2004/10/04 22:18:50 greve Exp $";

char *Progname ;

char *inputlist[100];
int nthinput=0;
int PrintTR=0;
int PrintTE=0;
int PrintTI=0;
int PrintFlipAngle=0;
int PrintCRes = 0;
int PrintRRes = 0;
int PrintSRes = 0;
int PrintNCols = 0;
int PrintNRows = 0;
int PrintNSlices = 0;
int PrintNFrames = 0;
int PrintFormat = 0;

int debug = 0;

/***-------------------------------------------------------****/
int main(int argc, char *argv[])
{
  int nargs, n;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  for(n=0;n<nthinput;n++) {
    if(debug) printf("%d %s ----- \n",n,inputlist[n]);
    do_file(inputlist[n]);
  }

  exit(0);

} /* end main() */

/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--tr"))   PrintTR = 1;
    else if (!strcasecmp(option, "--te"))   PrintTE = 1;
    else if (!strcasecmp(option, "--ti"))   PrintTI = 1;
    else if (!strcasecmp(option, "--fa"))   PrintFlipAngle = 1;
    else if (!strcasecmp(option, "--flip_angle"))   PrintFlipAngle = 1;

    else if (!strcasecmp(option, "--cres"))    PrintCRes = 1;
    else if (!strcasecmp(option, "--xsize"))   PrintCRes = 1;
    else if (!strcasecmp(option, "--rres"))    PrintRRes = 1;
    else if (!strcasecmp(option, "--ysize"))   PrintCRes = 1;
    else if (!strcasecmp(option, "--sres"))    PrintSRes = 1;
    else if (!strcasecmp(option, "--zsize"))   PrintCRes = 1;

    else if (!strcasecmp(option, "--ncols"))     PrintNCols = 1;
    else if (!strcasecmp(option, "--width"))     PrintNCols = 1;
    else if (!strcasecmp(option, "--nrows"))     PrintNRows = 1;
    else if (!strcasecmp(option, "--height"))    PrintNRows = 1;
    else if (!strcasecmp(option, "--nslices"))   PrintNSlices = 1;
    else if (!strcasecmp(option, "--depth"))     PrintNSlices = 1;

    else if (!strcasecmp(option, "--nframes"))   PrintNFrames = 1;
    else if (!strcasecmp(option, "--format")) PrintFormat = 1;

    else{
      inputlist[nthinput] = option;
      nthinput++;
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s fname1 <fname2> <options> \n",Progname) ;
  printf("\n");
  printf("   --tr : print TR to stdout\n");
  printf("   --te : print TE to stdout\n");
  printf("   --ti : print TI to stdout\n");
  printf("   --fa : print flip angle to stdout\n");
  printf("   --cres : print column voxel size (xsize) to stdout\n");
  printf("   --rres : print row    voxel size (ysize) to stdout\n");
  printf("   --sres : print slice  voxel size (zsize) to stdout\n");
  printf("   --ncols : print number of columns (width) to stdout\n");
  printf("   --nrows : print number of rows (height) to stdout\n");
  printf("   --nslices : print number of columns (depth) to stdout\n");
  printf("   --nframes : print number of frames to stdout\n");
  printf("   --format : file format\n");
  printf("\n");
  //printf("   --svol svol.img (structural volume)\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
"\n"
"Dumps information about the volume to stdout. Specific pieces \n"
"of information can be printed out as well by specifying the proper\n"
"flag (eg, --tr for TR). Time is in msec. Distance is in MM. Angles\n"
"are in radians.\n"
);


  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(nthinput == 0){
    printf("ERROR: no input volume supplied\n");
    exit(1);
  }
  return;
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

/***-------------------------------------------------------****/
int PrettyMatrixPrint(MATRIX *mat)
{
  int row;

  if (mat == NULL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat = NULL!")) ;

  if (mat->type != MATRIX_REAL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not Real type")) ;
 
  if (mat->rows != 4 || mat->cols != 4)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not of 4 x 4")) ;
    
  for (row=1; row < 5; ++row)
    printf("              %8.4f %8.4f %8.4f %10.4f\n",
	   mat->rptr[row][1], mat->rptr[row][2], mat->rptr[row][3], mat->rptr[row][4]);
  return (NO_ERROR);
}

/***-------------------------------------------------------****/
static void do_file(char *fname)
{
  MRI *mri ;
  MATRIX *m ;

  if(PrintFormat){
    printf("%s\n", type_to_string(mri_identify(fname)));
    return;
  }
  mri = MRIreadHeader(fname, MRI_VOLUME_TYPE_UNKNOWN) ;
  if(!mri) return;

  if(PrintTR){
    printf("%g\n",mri->tr);
    return;
  }
  if(PrintTE){
    printf("%g\n",mri->te);
    return;
  }
  if(PrintTI){
    printf("%g\n",mri->ti);
    return;
  }
  if(PrintFlipAngle){
    printf("%g\n",mri->flip_angle);
    return;
  }
  if(PrintCRes){
    printf("%g\n",mri->xsize);
    return;
  }
  if(PrintRRes){
    printf("%g\n",mri->ysize);
    return;
  }
  if(PrintSRes){
    printf("%g\n",mri->zsize);
    return;
  }
  if(PrintNCols){
    printf("%d\n",mri->width);
    return;
  }
  if(PrintNRows){
    printf("%d\n",mri->height);
    return;
  }
  if(PrintNSlices){
    printf("%d\n",mri->depth);
    return;
  }
  if(PrintNFrames){
    printf("%d\n",mri->nframes);
    return;
  }


  printf("Volume information for %s\n", fname);
  // mri_identify has been called but the result is not stored and thus I have to call it again
  printf("          type: %s\n", type_to_string(mri_identify(fname)));
  if (mri->nframes > 1)
    printf("    dimensions: %d x %d x %d x %d\n", mri->width, mri->height, mri->depth, mri->nframes) ;
  else
    printf("    dimensions: %d x %d x %d\n", mri->width, mri->height, mri->depth) ;
  printf("   voxel sizes: %6.4f, %6.4f, %6.4f\n", mri->xsize, mri->ysize, mri->zsize) ;
  printf("          type: %s (%d)\n",
	 mri->type == MRI_UCHAR   ? "UCHAR" :
	 mri->type == MRI_SHORT   ? "SHORT" :
	 mri->type == MRI_INT     ? "INT" :
	 mri->type == MRI_LONG    ? "LONG" :
	 mri->type == MRI_BITMAP  ? "BITMAP" :
	 mri->type == MRI_TENSOR  ? "TENSOR" :
	 mri->type == MRI_FLOAT   ? "FLOAT" : "UNKNOWN", mri->type) ;
  printf("           fov: %2.3f\n", mri->fov) ;
  printf("        xstart: %2.1f, xend: %2.1f\n", mri->xstart*mri->xsize, mri->xend*mri->xsize) ;
  printf("        ystart: %2.1f, yend: %2.1f\n", mri->ystart*mri->ysize, mri->yend*mri->ysize) ;
  printf("        zstart: %2.1f, zend: %2.1f\n", mri->zstart*mri->zsize, mri->zend*mri->zsize) ;
  printf("            TR: %2.2f msec, TE: %2.2f msec, TI: %2.2f msec, flip angle: %2.2f degrees\n",
	 mri->tr, mri->te, mri->ti, DEGREES(mri->flip_angle)) ;
  printf("       nframes: %d\n", mri->nframes) ;
  printf("ras xform %spresent\n", mri->ras_good_flag ? "" : "not ") ;
  printf("    xform info: x_r = %8.4f, y_r = %8.4f, z_r = %8.4f, c_r = %10.4f\n",
	 mri->x_r, mri->y_r, mri->z_r, mri->c_r);
  printf("              : x_a = %8.4f, y_a = %8.4f, z_a = %8.4f, c_a = %10.4f\n",
	 mri->x_a, mri->y_a, mri->z_a, mri->c_a);
  printf("              : x_s = %8.4f, y_s = %8.4f, z_s = %8.4f, c_s = %10.4f\n",
	 mri->x_s, mri->y_s, mri->z_s, mri->c_s);

  if (fio_IsDirectory(fname))
    printf("\ntalairach xfm : %s\n", mri->transform_fname);
  else
  {
    char *ext = 0;
    ext = fio_extension(fname);
    if (ext)
    {
      if (strcmp(ext, "mgz") == 0 || strcmp(ext, "mgh")==0)
	printf("\ntalairach xfm : %s\n", mri->transform_fname);
      free(ext);
    }
  }
  m = MRIgetVoxelToRasXform(mri) ; // extract_i_to_r(mri) (just macto)
  printf("\nvoxel to ras transform:\n") ; PrettyMatrixPrint(m) ;
  MatrixFree(&m) ;
  m = extract_r_to_i(mri);
  printf("\nras to voxel transform:\n"); PrettyMatrixPrint(m);
  MatrixFree(&m);
  MRIfree(&mri);
  
  return;

} /* end do_file */

