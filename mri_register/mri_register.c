

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "diag.h"

char         *Progname ;

static unsigned char threshold = 40 ;

/* options */
static int    dx = 0, dy = 0, dz = 0 ;
static float  xa = 0.0f, ya = 0.0f, za = 0.0f ;
static int    verbose = 0 ;

typedef struct
{
  MATRIX *evectors ;
  MATRIX *evalues ;
  MATRIX *means ;
} REG_PARMS ;

/*static MRI  *register_mri(MRI *mri_in, MRI *mri_reg) ;*/
static MRI  *register_mri(MRI *mri_in, MRI *mri_ref, MRI *mri_reg,
                          REG_PARMS *ref_parms, REG_PARMS *in_parms);
static int get_option(int argc, char *argv[]) ;
void main(int argc, char *argv[]) ;
static int load_parms(char *dir_name, REG_PARMS *parms) ;
static int write_parms(char *dir_name, REG_PARMS *parms) ;
static int free_parms(REG_PARMS *parms) ;

/* 
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/
void
main(int argc, char *argv[])
{
  char       ref_dir_name[STR_LEN], in_dir_name[STR_LEN],out_dir_name[STR_LEN];
  MRI        *mri_ref, *mri_in, *mri_reg ;
  char       **av ;
  int        ac, nargs ;
  REG_PARMS  ref_parms, in_parms ;

  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM, 
              "usage: %s <ref brain dir> <in brain dir> <out brain dir>\n",
              Progname) ;

  strcpy(ref_dir_name, argv[1]) ;
  strcpy(in_dir_name, argv[2]) ;

  if (argc < 4)  /* no output directory specified - construct name */
    sprintf(out_dir_name, "%s_reg", in_dir_name) ;
  else
    strcpy(out_dir_name, argv[3]) ;

  /*
     check to make sure output directory exists, and, if not, create it.
     */
  if (verbose)
    fprintf(stderr, "reading '%s'...", ref_dir_name) ;
  mri_ref = MRIread(ref_dir_name) ;
  if (verbose)
  {
    fprintf(stderr, "done.\n") ;
    fflush(stderr) ;
  }

  if (verbose)
    fprintf(stderr, "reading '%s'...", in_dir_name) ;
  mri_in = MRIread(in_dir_name) ;
  if (verbose)
  {
    fprintf(stderr, "done.\n") ;
    fflush(stderr) ;
  }

  load_parms(ref_dir_name, &ref_parms) ;
  mri_reg = register_mri(mri_in, mri_ref, NULL, &ref_parms, &in_parms) ;
  write_parms(out_dir_name, &in_parms) ;
  write_parms(ref_dir_name, &ref_parms) ;
  MRIwrite(mri_reg, out_dir_name) ;


#if 0
  if (!FZERO(za))
  {
    mri_tmp = MRIrotateZ(mri_in, NULL, za) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }
  if (!FZERO(ya))
  {
    mri_tmp = MRIrotateY(mri_in, NULL, ya) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }
  if (!FZERO(xa))
  {
    mri_tmp = MRIrotateX(mri_in, NULL, xa) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }

  if (dx || dy || dz)
  {
    fprintf(stderr, "translating by (%d, %d, %d)\n", dx, dy, dz) ;
    mri_tmp = MRItranslate(mri_in, NULL, -25, 5, 0) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }
#endif

  free_parms(&ref_parms) ;
  free_parms(&in_parms) ;
  if (mri_reg)
    MRIfree(&mri_reg) ;
  if (mri_ref)
    MRIfree(&mri_ref) ;
  if (mri_in)
    MRIfree(&mri_in) ;
  exit(0) ;
}

static MRI *
register_mri(MRI *mri_in, MRI *mri_ref, MRI *mri_reg, REG_PARMS *ref_parms,
             REG_PARMS *in_parms)
{
  int    dx, dy, dz, row, col, in_means[3], ref_means[3] ;
  MATRIX *mRot, *m_in_T, *mOrigin ;
  MRI    *mri_tmp ;
  float  x_angle, y_angle, z_angle, r11, r21, r31, r32, r33, cosy, dot ;
  float  ref_evalues[3], in_evalues[3] ;

  in_parms->evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
  in_parms->evalues = MatrixAlloc(3,1, MATRIX_REAL) ;
  in_parms->means = MatrixAlloc(3,1,MATRIX_REAL) ;

  if (!ref_parms->evectors)  /* not loaded from file */
  {
    ref_parms->evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
    ref_parms->evalues = MatrixAlloc(3,1, MATRIX_REAL) ;
    ref_parms->means = MatrixAlloc(3,1,MATRIX_REAL) ;
    MRIprincipleComponents(mri_ref, ref_parms->evectors, ref_evalues, 
                           ref_means, threshold);
    
    for (row = 1 ; row <= 3 ; row++)
    {
      ref_parms->evalues->rptr[row][1] = ref_evalues[row-1] ;
      ref_parms->means->rptr[row][1] = ref_means[row-1] ;
    }
  }
  else
  {
    /* copy eigenvalues and means from MATRIX structures */
    for (row = 1 ; row <= 3 ; row++)
    {
      ref_evalues[row-1] = ref_parms->evalues->rptr[row][1] ;
      ref_means[row-1]  = ref_parms->means->rptr[row][1] ;
    }
  }

  MRIprincipleComponents(mri_in, in_parms->evectors, in_evalues,in_means,
                         threshold);
  if (verbose)
    MatrixPrint(stdout, in_parms->evectors) ;

  /* copy eigenvalues and means into MATRIX structures */
  for (row = 1 ; row <= 3 ; row++)
  {
    in_parms->evalues->rptr[row][1] = in_evalues[row-1] ;
    in_parms->means->rptr[row][1] = in_means[row-1] ;
  }

  /* check to make sure eigenvectors aren't reversed */
  for (col = 1 ; col <= 3 ; col++)
  {
    for (dot = 0.0f, row = 1 ; row <= 3 ; row++)
      dot += 
        in_parms->evectors->rptr[row][col] * 
          ref_parms->evectors->rptr[row][col] ;

    if (dot < 0.0f)
    {
      dot *= -1.0f ;
      for (row = 1 ; row <= 3 ; row++)
        in_parms->evectors->rptr[row][col] *= -1.0f ;
    }
  }

  m_in_T = MatrixTranspose(in_parms->evectors, NULL) ;
  mRot = MatrixMultiply(ref_parms->evectors, m_in_T, NULL) ;

  r11 = mRot->rptr[1][1] ;
  r21 = mRot->rptr[2][1] ;
  r31 = mRot->rptr[3][1] ;
  r32 = mRot->rptr[3][2] ;
  r33 = mRot->rptr[3][3] ;
  y_angle = atan2(-r31, sqrt(r11*r11+r21*r21)) ;
  cosy = cos(y_angle) ;
  z_angle = atan2(r21 / cosy, r11 / cosy) ;
  x_angle = atan2(r32 / cosy, r33 / cosy) ;
  if (verbose)
    fprintf(stderr, "rotation: (%2.2f, %2.2f, %2.2f)\n",
            DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;
  mOrigin = VectorAlloc(3, MATRIX_REAL) ;
  mOrigin->rptr[1][1] = in_means[0] ;
  mOrigin->rptr[2][1] = in_means[1] ;
  mOrigin->rptr[3][1] = in_means[2] ;

  dx = nint(ref_means[0] - in_means[0]) ;
  dy = nint(ref_means[1] - in_means[1]) ;
  dz = nint(ref_means[2] - in_means[2]) ;

  mri_tmp = MRIrotate(mri_in, NULL, mRot, mOrigin) ;
  mri_reg = MRItranslate(mri_tmp, mri_reg, dx, dy, dz) ;

  if (mri_tmp)
    MRIfree(&mri_tmp) ;
  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(mri_reg) ;
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
  StrUpper(option) ;
  if (!strcmp(option, "DX"))
  {
    sscanf(argv[2], "%d", &dx) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DY"))
  {
    sscanf(argv[2], "%d", &dy) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DZ"))
  {
    sscanf(argv[2], "%d", &dz) ;
    nargs = 1 ;
  }
  else switch (*option)
  {
  case '?':
  case 'U':
    printf("usage: %s [input directory] [output directory]\n", argv[0]) ;
    exit(1) ;
    break ;
  case 'X':
    sscanf(argv[2], "%f", &xa) ;
    nargs = 1 ;
    fprintf(stderr, "rotating about x axis by %2.1f degrees\n", xa) ;
    xa = RADIANS(xa) ;
    break ;
  case 'Y':
    sscanf(argv[2], "%f", &ya) ;
    nargs = 1 ;
    fprintf(stderr, "rotating about y axis by %2.1f degrees\n", ya) ;
    ya = RADIANS(ya) ;
    break ;
  case 'Z':
    sscanf(argv[2], "%f", &za) ;
    nargs = 1 ;
    fprintf(stderr, "rotating about z axis by %2.1f degrees\n", za) ;
    za = RADIANS(za) ;
    break ;
  case 'V':
    verbose = 1 ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static int
load_parms(char *dir_name, REG_PARMS *parms)
{
  char fname[STR_LEN] ;

  sprintf(fname, "%s/evectors.mat", dir_name) ;
  parms->evectors = MatrixRead(fname) ;
  if (!parms->evectors)
    return(ERROR_NO_FILE) ;

  if (verbose)
    fprintf(stderr, "using previously calculated reference parameters\n") ;
  sprintf(fname, "%s/evalues.mat", dir_name) ;
  parms->evalues = MatrixRead(fname) ;

  sprintf(fname, "%s/means.mat", dir_name) ;
  parms->means = MatrixRead(fname) ;
  return(NO_ERROR) ;
}

static int
write_parms(char *dir_name, REG_PARMS *parms)
{
  char fname[STR_LEN] ;

  sprintf(fname, "%s/evectors.mat", dir_name) ;
  MatrixWrite(parms->evectors, fname, "evectors") ;
  sprintf(fname, "%s/evalues.mat", dir_name) ;
  MatrixWrite(parms->evalues, fname, "evalues") ;
  sprintf(fname, "%s/means.mat", dir_name) ;
  MatrixWrite(parms->means, fname, "means") ;
  return(NO_ERROR) ;
}

static int
free_parms(REG_PARMS *parms)
{
  if (parms->evectors)
    MatrixFree(&parms->evectors) ;
  if (parms->evalues)
    MatrixFree(&parms->evalues) ;
  if (parms->means)
    MatrixFree(&parms->means) ;
  return(NO_ERROR) ;
}

