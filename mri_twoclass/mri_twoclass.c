#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


#include "timer.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "const.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "fio.h"
#include "sig.h"
#include "label.h"
#include "cma.h"
#include "vlabels.h"

static char vcid[] = "$Id: mri_twoclass.c,v 1.5 2002/07/19 16:32:16 fischl Exp $";


/*-------------------------------- STRUCTURES ----------------------------*/


/*-------------------------------- CONSTANTS -----------------------------*/

#define STAT_T           0
#define STAT_F           1
#define STAT_MEAN        2

/*-------------------------------- PROTOTYPES ----------------------------*/

static char *read_fname1 = NULL ;
static char *read_fname2 = NULL ;
static int Gxv = -1 ;
static int Gyv = -1 ;
static int Gzv = -1 ;
static int Gxn = -1 ;
static int Gyn = -1 ;
static int Gzn = -1 ;

static float resolution = 2 ;
static float fthresh = -1 ;

int main(int argc, char *argv[]) ;

static void write_bfloats(MRI *mri, char *out_name, char *output_subject);
static int find_label_index(VL ***voxel_labels, int x, int y, int z,int label);
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  add_volume_labels_to_average(MRI *mri, VL ***voxel_labels,
                                         float resolution, TRANSFORM *transform);

static int voxel_to_node(MRI *mri, float resolution, int xv, int yv, int zv, 
                         int *pxn, int *pyn, int *pzn, TRANSFORM *transform) ;
static char *xform_fname = "talairach.lta" ;
static MRI *compute_voxel_statistics(VL ***voxel_labels_class1, 
                                     VL ***voxel_labels_class2,
                                     int width, int height,
                                     int depth, float resolution,
                                     int nclass1, int nclass2, MRI *mri_stats);

/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static char *read_dir = NULL ;

static char *test_subject = NULL ;
static char *label_name = NULL ;
static char *prefix = "" ;
static char *vl1_name = NULL, *vl2_name = NULL ;
static int output_bfloats = 1 ;
static int bonferroni = 0 ;
static char subjects_dir[STRLEN] ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  MRI          *mri, *mri_stats ;
  char         **av, fname[STRLEN], *aseg_name,
               *cp, *subject_name, *out_fname,
               **c1_subjects, **c2_subjects, *output_subject ;
  int          ac, nargs, n, num_class1, num_class2, i, width, height, depth, 
               awidth, aheight, adepth ;
  LABEL        *area ;
  struct timeb start ;
  int          msec, minutes, seconds ;
  VL           ***voxel_labels_class1 = NULL, ***voxel_labels_class2=NULL, 
               ***vls ;
  TRANSFORM    *transform ;
  VLI          *vli1 = NULL, *vli2 = NULL ;
  

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  TimerStart(&start) ;

  /* subject_name hemi surface curvature */
  if (argc < 7)
    usage_exit() ;
  
  if (!strlen(subjects_dir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname) ;

    strcpy(subjects_dir, cp) ;
  }

  aseg_name = argv[1] ;
  output_subject = argv[2] ;
  out_fname = argv[3] ;

  width = awidth = aheight = adepth = height = 0 ;

#define ARGV_OFFSET 4

  /* first determine the number of subjects in each class */
  num_class1 = 0 ; n = ARGV_OFFSET ;
  do
  {
    num_class1++ ;
    n++ ;
    if (argv[n] == NULL || n >= argc)
      ErrorExit(ERROR_BADPARM, "%s: must spectify ':' between class lists",
                Progname) ;
  } while(argv[n][0] != ':') ;



  num_class2 = 0 ; n++ ; /* skip ':' */
  if (n >= argc)
    ErrorExit(ERROR_BADPARM, "%s: class2 list empty", Progname) ;
  do
  {
    num_class2++ ;
    n++ ;
    if (n >= argc)
      break ;
  } while(argv[n] != NULL) ;

  c1_subjects = (char **)calloc(num_class1, sizeof(char *)) ;
  c2_subjects = (char **)calloc(num_class2, sizeof(char *)) ;
  fprintf(stderr, "%d subjects in class 1, %d subjects in class 2\n",
          num_class1, num_class2) ;


  /* read in subject names for the two classes */
  for (n = 0 ; n < num_class1 ; n++)
  {
    c1_subjects[n] = argv[ARGV_OFFSET+n] ;
  }
  i = n+1+ARGV_OFFSET ;  /* starting index */
  for (n = 0 ; n < num_class2 ; n++)
  {
    c2_subjects[n] = argv[i+n] ;
  }

  /* read in label if limitting calculation to an ROI */
  if (label_name)
  {
    area = LabelRead(output_subject, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
  }
  else
    area = NULL ;

  if (read_dir)  /* read in precomputed data */
  {
  }
  else
  {
    /* real all the data for both groups  */
    for (n = 0 ; n < num_class1+num_class2 ; n++)
    {
      /* transform each subject's curvature into the output subject's space */
      subject_name = n < num_class1 ? c1_subjects[n]:c2_subjects[n-num_class1];
      fprintf(stderr, "reading subject %d of %d: %s\n",
              n+1, num_class1+num_class2, subject_name) ;
      sprintf(fname, "%s/%s/mri/%s", subjects_dir,subject_name,aseg_name);
      mri = MRIread(fname) ;
      if (!mri)
        ErrorExit(ERROR_NOFILE, "%s: could not read segmentation %s",
                  Progname, fname) ;
      sprintf(fname, "%s/%s/mri/transforms/%s", subjects_dir,subject_name,
              xform_fname);
      transform = TransformRead(fname) ;
      if (!transform)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform %s",
                  Progname, fname) ;
      TransformInvert(transform, mri) ;


      if (!width)
      {
        width = mri->width ; height = mri->height ; depth = mri->depth ;
        awidth = width/resolution ;
        aheight = height/resolution ;
        adepth = depth/resolution ;
        vli1 = VLalloc(awidth, aheight, adepth, resolution) ;
        vli2 = VLalloc(awidth, aheight, adepth, resolution) ;
        if (!vli1 || !vli2)
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate VL ***", Progname);
        
        voxel_labels_class1 = vli1->vl ;
        voxel_labels_class2 = vli2->vl ;
      }

      vls = n < num_class1 ? voxel_labels_class1 : voxel_labels_class2 ;
      add_volume_labels_to_average(mri, vls, resolution, transform) ;
      MRIfree(&mri) ; TransformFree(&transform) ;
    }
  }


  mri_stats = 
    compute_voxel_statistics(voxel_labels_class1, voxel_labels_class2,
                             awidth, aheight, adepth, resolution,
                             num_class1, num_class2, NULL) ;
  sprintf(fname, "%s/%s/mri/%s", subjects_dir,output_subject,out_fname);
  printf("writing stats to %s...\n", fname) ;

  if (!output_bfloats)
    MRIwrite(mri_stats, fname) ;
  else
    write_bfloats(mri_stats, fname, output_subject) ;

  if (vl1_name && stricmp(vl1_name, "none"))
  {
    printf("writing voxel labels for group 1 to %s...\n", vl1_name) ;
    VLwrite(vli1, vl1_name) ;
  }
  if (vl2_name && stricmp(vl2_name, "none"))
  {
    printf("writing voxel labels for group 2 to %s...\n", vl2_name) ;
    VLwrite(vli2, vl2_name) ;
  }

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "cross-subject statistics took %d minutes and %d seconds.\n",
          minutes, seconds) ;
  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "test"))
  {
    test_subject = argv[2] ;
    fprintf(stderr, "writing test.dat for subject %s\n", test_subject) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gxv = atoi(argv[2]) ;
    Gyv = atoi(argv[3]) ;
    Gzv = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gxv,Gyv,Gzv) ;
  }
  else if (!stricmp(option, "DEBUG_NODE"))
  {
    Gxn = atoi(argv[2]) ;
    Gyn = atoi(argv[3]) ;
    Gzn = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gxn,Gyn,Gzn) ;
  }
  else if (!stricmp(option, "WRITE_LABELS"))
  {
    vl1_name = argv[2] ;
    vl2_name = argv[3] ;
    nargs = 2 ;
    printf("writing voxel labels to %s and %s\n", vl1_name, vl2_name) ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(subjects_dir, argv[2]) ;
    fprintf(stderr, "using subjects_dir %s\n", subjects_dir) ; 
    nargs = 1 ;
  }
  else if (!stricmp(option, "read"))
  {
    read_fname1 = argv[2] ;
    read_fname2 = argv[3] ;
    printf("reading precomputed voxel labels from %s and %s\n",
           read_fname1, read_fname2) ;
    nargs = 2 ;
  }
  else switch (toupper(*option))
  {
  case 'X':
    xform_fname = argv[2] ;
    nargs = 1 ;
    printf("using xform fname %s\n", xform_fname) ;
    break ;
  case 'R':
    resolution = atof(argv[2]) ;
    printf("using resolution %2.1f\n", resolution) ;
    nargs = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'L':
    label_name = argv[2] ;
    fprintf(stderr, "masking label %s\n", label_name) ;
    nargs = 1 ;
    break ;
  case 'P':
    prefix = argv[2] ;
    fprintf(stderr, "using label prefix %s\n", prefix) ;
    nargs = 1 ;
    break ;
  case 'T':
    fthresh = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using F snr threshold of %2.2f...\n", fthresh) ;
    break ;
  case 'B':
    bonferroni = 1 ;
    fprintf(stderr, "doing bonferroni correction of SNR values...\n") ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  printf("usage: %s [options] \n"
         "\t<segmentation volume> <output subject> <output volume> \n\t<c1_subject1> "
         "<c1_subject2>... : \n"
         "\t<c2_subject1> <c2_subject2>...\n",
         Progname) ;
  printf("where surf must be a spherical surface suitable for "
         "computing geodesics\n") ;
  printf("The <c1_subject> ... is a list of subjects from one class\n"
         "and the <c2_subject>... is a list of subjects from another "
         "class.\n");
  printf("valid options are:\n") ;
  printf("\t-t <fthresh>       -  specify F threshold\n"
         "\t-o <subject>       -  specify output subject name\n"
         "\t-b                 -  perform bonferroni correction\n") ;
}

static void
print_help(void)
{
  print_usage() ;
  printf( 
          "\nThis program will compute the cross-subject statistics of two "
          "sets of labels.\n") ;
  printf( "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  printf( "%s\n", vcid) ;
  exit(1) ;
}

static int
add_volume_labels_to_average(MRI *mri, VL ***voxel_labels,float resolution, 
                             TRANSFORM *transform)
{
  int          x, y, z, width, height, depth, index, label, xv, yv, zv ;
  VOXEL_LABELS *vl ;

  width = mri->width; height = mri->height; depth = mri->height;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri, x, y, z) ;
        voxel_to_node(mri, resolution, x, y, z, &xv, &yv, &zv, transform) ;
        index = find_label_index(voxel_labels, xv, yv, zv,label) ;
        if (xv == Gxn && yv == Gyn && zv == Gzn)
          printf("(%d, %d, %d) --> (%d,%d,%d), label = %d\n",
                 x, y, z, xv, yv, zv, label) ;

        if (xv == 96 && yv == 66 && (zv == 56 || zv == 57))
          DiagBreak() ;
        vl = &voxel_labels[xv][yv][zv] ;
        if (index < 0)  /* allocate space for it */
        {
          if (vl->nlabels > 0)
          {
            unsigned char *old_labels ;
            unsigned short *old_counts ;

            if (xv == Gxn && yv == Gyn && zv == Gzn)
              DiagBreak() ;
            old_labels = vl->labels ; old_counts = vl->counts ;
            vl->labels = (unsigned char *)calloc(vl->nlabels+1, 
                                                 sizeof(unsigned char));
            vl->counts = (unsigned short *)calloc(vl->nlabels+1,
                                                  sizeof(unsigned short));
            if (!vl->labels)
              ErrorExit(ERROR_NOMEMORY, "%s: could not reallocate %d labels "
                        "at %d,%d,%d", Progname, vl->nlabels+1, xv, yv, zv) ;
            memmove(vl->labels, old_labels,vl->nlabels*sizeof(unsigned char));
            memmove(vl->counts, old_counts,vl->nlabels*sizeof(unsigned short));
            index = vl->nlabels++ ;
            free(old_labels) ; free(old_counts) ;
          }
          else
          {
            vl->labels = 
              (unsigned char *)calloc(vl->nlabels+1, sizeof(unsigned char));
            vl->counts = 
              (unsigned short *)calloc(vl->nlabels+1,sizeof(unsigned short));
            if (!vl->labels || !vl->counts)
              ErrorExit(ERROR_NOMEMORY, "%s: could not reallocate %d labels "
                        "at %d,%d,%d", Progname, vl->nlabels+1, xv, yv, zv) ;
            if (!vl->labels)
              ErrorExit(ERROR_NOMEMORY, "%s: could not reallocate %d labels "
                        "at %d,%d,%d", Progname, vl->nlabels+1, xv, yv, zv) ;
            index = vl->nlabels++ ;
          }
        }
        if (xv == Gxn && yv == Gyn && zv == Gzn)
          DiagBreak() ;
        if (xv == Gxn && yv == Gyn && zv == Gzn && vl->nlabels > 1)
          DiagBreak() ;
        vl->labels[index] = label ;
        vl->counts[index]++ ;
        if (Gxn > 0 && voxel_labels[Gxn][Gyn][Gzn].nlabels > 1)
          DiagBreak() ;
      }
    }
  }

  return(NO_ERROR) ;
}

static int
voxel_to_node(MRI *mri, float resolution, int xv, int yv, int zv, 
              int *pxn, int *pyn, int *pzn, TRANSFORM *transform)
{
#if 1
  float  xt, yt, zt, xscale, yscale, zscale ;
  int    width, height, depth ;

  TransformSample(transform, xv*mri->xsize, yv*mri->ysize, zv*mri->zsize, &xt, &yt, &zt) ;

  xt = nint(xt/mri->xsize) ; yt = nint(yt/mri->ysize) ; zt = nint(zt/mri->zsize); 
  width = mri->width/resolution ;
  height = mri->height/resolution ;
  depth = mri->depth/resolution ;
  xscale = mri->xsize / resolution ; 
  yscale = mri->ysize / resolution ;
  zscale = mri->zsize / resolution ;
  *pxn = nint((xt) * xscale) ;
  if (*pxn < 0)
    *pxn = 0 ;
  else if (*pxn >= width)
    *pxn = width-1 ;

  *pyn = nint((yt) * yscale) ;
  if (*pyn < 0)
    *pyn = 0 ;
  else if (*pyn >= height)
    *pyn = height-1 ;

  *pzn = nint((zt) * zscale) ;
  if (*pzn < 0)
    *pzn = 0 ;
  else if (*pzn >= depth)
    *pzn = depth-1 ;
#else
  static VECTOR *v_input, *v_canon = NULL ;
  int    width, height, depth ;


  width = mri->width * xscale ;
  height = mri->height * yscale ;
  depth = mri->depth * zscale ;

  if (!v_canon)
  {
    v_input = VectorAlloc(4, MATRIX_REAL) ;
    v_canon = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v_input, 4, 1) = 1.0 ;
    *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
  }

  /* scale voxel coords  up to mm */
  V3_X(v_input) = (float)xv*mri->xsize ; 
  V3_Y(v_input) = (float)yv*mri->ysize ; 
  V3_Z(v_input) = (float)zv*mri->zsize ; 
  MatrixMultiply(lta->xforms[0].m_L, v_input, v_canon) ;

  /* scale mm back to voxel coords */
  xt = (V3_X(v_canon)/resolution) ; 
  yt = (V3_Y(v_canon)/resolution) ; 
  zt = (V3_Z(v_canon)/resolution) ; 

  *pxn = nint(xt) ;
  if (*pxn < 0)
    *pxn = 0 ;
  else if (*pxn >= width)
    *pxn = width-1 ;

  *pyn = nint(yt) ;
  if (*pyn < 0)
    *pyn = 0 ;
  else if (*pyn >= height)
    *pyn = height-1 ;

  *pzn = nint(zt) ;
  if (*pzn < 0)
    *pzn = 0 ;
  else if (*pzn >= depth)
    *pzn = depth-1 ;
#endif
  return(NO_ERROR) ;
}

static int
find_label_index(VL ***voxel_labels, int x, int y, int z, int label)
{
  int  index ;
  VL   *vl ;

  vl = &voxel_labels[x][y][z] ;
  for (index = 0 ; index < vl->nlabels ; index++)
    if (vl->labels[index] == label)
      return(index) ;
  return(-1) ;
}

static MRI *
compute_voxel_statistics(VL ***voxel_labels_class1, VL ***voxel_labels_class2,
                         int width, int height, int depth,
                         float resolution, int nclass1, int nclass2, 
                         MRI *mri_stats)
{
  int          x, y, z, label_counts_c1[256], label_counts_c2[256], index, dof;
  int          xc, yc, zc, xp, yp, zp ;
  VOXEL_LABELS *vl ;
  double       chisq, pi, xi1, xi2, numer, denom, n1, n2, ntotal, p, 
               max_chisq, min_p ;

  n1 = (double)nclass1 ; n2 = (double)nclass2 ;
  ntotal = n1 + n2 ;

  if (!mri_stats)
  {
    mri_stats = MRIalloc(width, height, depth, MRI_FLOAT) ;
    if (!mri_stats)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %dx%dx%d stats volume",
                Progname, width, height, depth) ;
    mri_stats->xsize = mri_stats->ysize = mri_stats->zsize = resolution ;
  }
  max_chisq = 0.0 ; min_p = 1.0 ;
  xp = yp = zp = xc = yc = zc = 0 ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gxn && y == Gyn && z == Gzn)
          DiagBreak() ;

        if (x == 96 && y == 66 && (z == 56 || z == 57))
          DiagBreak() ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        memset(label_counts_c1, 0, sizeof(label_counts_c1)) ;
        vl = &voxel_labels_class1[x][y][z] ;
        for (index = 0 ; index < vl->nlabels ; index++)
          label_counts_c1[vl->labels[index]] += vl->counts[index] ;

        memset(label_counts_c2, 0, sizeof(label_counts_c2)) ;
        vl = &voxel_labels_class2[x][y][z] ;
        for (index = 0 ; index < vl->nlabels ; index++)
          label_counts_c2[vl->labels[index]] += vl->counts[index] ;
        
        chisq = 0.0 ; dof = 0 ;
        for (n1 = n2 = 0.0, index = 0 ; index <= MAX_CMA_LABEL ; index++)
        {
          n1 += (double)label_counts_c1[index] ;
          n2 += (double)label_counts_c2[index] ;
        }
        ntotal = n1 + n2 ;

        for (index = 0 ; index <= MAX_CMA_LABEL ; index++)
        {
          if (label_counts_c1[index] == 0 && label_counts_c2[index] == 0)
            continue ;
          dof++ ;
          xi1 = (double)label_counts_c1[index] ;
          xi2 = (double)label_counts_c2[index] ;
          pi = (xi1 + xi2) / ntotal ;

          numer = (xi1 - (n1 * pi)) ; numer *= numer ; denom = n1 * pi ;
          if (!FZERO(denom))
            chisq += numer/denom ;
          numer = (xi2 - (n2 * pi)) ; numer *= numer ; denom = n2 * pi ;
          if (!FZERO(denom))
            chisq += numer/denom ;
        }
        if (dof > 1)
          DiagBreak() ;
        p = sigchisq(chisq, dof-1) ;
        if (DZERO(p) || p < 0)
        {
          p = 1e-20 ;
          MRIFvox(mri_stats, x, y, z) = 20 ;
        }
        else
        {
          if (!finite(-log(p)))
            DiagBreak() ;
          MRIFvox(mri_stats, x, y, z) = -log(p) ;
        }
        if (chisq > max_chisq)
        {
          xc = x  ; yc = y ; zc = z ;
          max_chisq = chisq ;
        }
        if (p < min_p)
        {
          xp = x  ; yp = y ; zp = z ;
          min_p = p ;
        }
        if (x == Gxn && y == Gyn && z == Gzn)
        {
          printf("class 1 @ (%d, %d, %d):\n", x, y, z) ;
          vl = &voxel_labels_class1[x][y][z] ;
          for (index = 0 ; index < MAX_CMA_LABEL ; index++)
            if (label_counts_c1[index] > 0)
              printf("\tlabel %03d: %d\n", 
                     index, label_counts_c1[index]) ;
          printf("class 2:\n") ;
          vl = &voxel_labels_class2[x][y][z] ;
          for (index = 0 ; index < MAX_CMA_LABEL ; index++)
            if (label_counts_c2[index] > 0)
              printf("\tlabel %03d: %d\n", 
                     index, label_counts_c2[index]) ;
          printf("\tchisq = %2.2f, p = %2.2e, dof = %d\n",
                 chisq, p, dof) ;
        }
      }
    }
  }

  printf("max chisq = %2.2f at (%d, %d, %d)\n", max_chisq, xc, yc, zc) ;
  printf("min p     = %2.2e at (%d, %d, %d)\n", min_p, xp, yp, zp) ;
  return(mri_stats) ;
}

#include "stats.h"
#include "matrix.h"
static void
write_bfloats(MRI *mri, char *out_name, char *output_subject)
{
  STAT_VOLUME  *sv ;
  fMRI_REG     *reg ;
  
  sv = StatAllocVolume(NULL, 1, mri->width, mri->height, mri->depth,1,0);
  strcpy(sv->reg->name, output_subject) ;
  MRIfree(&sv->mri_avgs[0]) ;
  sv->mri_avgs[0] = mri ;
  sv->nslices = mri->depth ;
  sv->slice_width = mri->width ;
  sv->slice_height = mri->height ;
  sv->voltype = 0 ;  /* Raw */
  sv->nevents = 1 ; sv->time_per_event = 1 ; 
  reg = sv->reg ;
  
  reg->in_plane_res = mri->xsize ;
  reg->slice_thickness = mri->zsize ;
  reg->brightness_scale = 1.0 ;
  
  StatWriteVolume(sv, out_name) ;
  StatFree(&sv) ;
}

