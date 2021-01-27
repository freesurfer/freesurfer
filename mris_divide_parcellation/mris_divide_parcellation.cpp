/*
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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



/*!
\author Bruce Fischl

*/

/*
BEGINUSAGE --------------------------------------------------------------

mris_divide_parcellation [options] subject hemi sourceannot [splitfile|areathresh] outannot

ENDUSAGE --------------------------------------------------------------

BEGINHELP --------------------------------------------------------------

This program divides one or more parcellations into divisions
perpendicular to the long axis of the label.  The number of divisions
can be specified in one of two ways, depending upon the nature of the
fourth argument.

First, a splitfile can be specified as the fourth argument. The
splitfile is a text file with two columns. The first column is the
name of the label in the source annotation, and the second column is
the number of units that label should be divided into. The names of
the labels depends upon the source parcellation.  For aparc.annot and
aparc.a2005.annot, the names can be found in
$FREESURFER_HOME/FreeSurferColorLUT.txt.  For aparc.annot, the labels
are between the ranges of 1000-1034.  For aparc.a2005s.annot, the
labels are between the ranges of 1100-1181.  The name for the label is
the name of the segmentation without the 'ctx-lh'. Note that the name
included in the splitfile does not indicate the hemisphere. For
example, 1023 is 'ctx-lh-posteriorcingulate'.  You should put
'posteriorcingulate' in the splitfile. Eg, to divide it into three
segments, the following line should appear in the splitfile:

posteriorcingulate 3

Only labels that should be split need be specified in the splitfile.

The second method is to specify an area threshold (in mm^2) as the
fourth argument, in which case each label is divided until each
subdivision is smaller than the area threshold.

The output label name will be the original name with _divN appended,
where N is the division number. N will go from 2 to the number of
divisions. The first division has the same name as the original label.


EXAMPLES:

Method 1:
cd $SUBJECTS_DIR/mysubj001/label

Create a split file to divide the superior frontal gyrus into 4
segements and the precentral gyrus into 3 segments:

echo superiorfrontal 4 >  splittable.txt
echo precentral      3 >> splittable.txt

Run program:

mris_divide_parcellation mysubj001 rh aparc.annot splittable.txt rh.aparc.split-sfg+pc

This reads in rh.aparc.annot, splits the SFG and PC, and and saves the
result to rh.aparc.split-sfg+pc

View

tksurfer mysubj001 rh inflated -annot aparc.split-sfg+pc

The SFG divisions will have the following names: superiorfrontal,
superiorfrontal_div2, superiorfrontal_div3, superiorfrontal_div4. The
PC divisions will be precentral, precentral_div2, precentral_div3.

Method 2:
Run program (from the label directory):

mris_divide_parcellation subj001 rh aparc.annot 100 rh.aparc.split.100mmsquared, and
saves the result to rh.aparc.split.100mmsquared 

This reads in rh.aparac.annot, splits every parcellation until each subdivision is less than 100 mm^2.
  

ENDHELP ----------------------------------------------------------------

*/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "version.h"
#include "utils.h"
#include "mrisurf.h"
#include "error.h"
#include "diag.h"
#include "annotation.h"
#include "matrix.h"

static int get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int main(int argc, char *argv[]) ;

const char *Progname = NULL;

static char sdir[STRLEN] = "" ;

static int rgb_scale = 30 ;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int         ac, nargs, *nunits, index ;
  MRI_SURFACE *mris ;
  char        **av, fname[STRLEN], *annot_name, *out_fname, *cp,
  *subject, *hemi ;
  float       area_thresh ;

  nargs = handleVersionOption(argc, argv, "mris_divide_parcellation");
  if (nargs && argc - nargs == 1)  exit (0);
  argc -= nargs;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)  usage_exit() ;
  subject = argv[1] ;
  hemi = argv[2] ;
  annot_name = argv[3] ;
  area_thresh = atof(argv[4]) ;
  out_fname = argv[5] ;

  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.sphere", sdir, subject, hemi) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s",
              Progname, fname) ;
  if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read annotation %s",
              Progname, annot_name) ;
  if (mris->ct == NULL)
    ErrorExit(ERROR_NOFILE, "%s: annotation %s must have embedded color table",
              Progname, annot_name) ;

  nunits = (int *)calloc(mris->ct->nentries, sizeof(int)) ;
  if (!nunits)
    ErrorExit(ERROR_BADPARM, "%s: could not allocate %d nunits table", Progname,mris->ct->nentries) ;

  MRIScomputeMetricProperties(mris) ;
  if (isdigit(*argv[4]))  // area threshold specified
  {
    int    vno ;
    VERTEX *v ;
    float  *area ;
    area = (float *)calloc(mris->ct->nentries, sizeof(float)) ;
    if (!area)
      ErrorExit(ERROR_BADPARM, "%s: could not allocate %d area table", Progname,mris->ct->nentries) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      CTABfindAnnotation(mris->ct, v->annotation, &index) ;
      if (index < 0 || v->ripflag)
        continue ;
      area[index] += v->area ;
    }
    for (index = 0 ; index < mris->ct->nentries ; index++)
      nunits[index] =  (int)(area[index] / area_thresh)+1 ;
    free(area) ;
  } else                 // interpret it as a file with parcellation names and # of units
  {
    char  line[STRLEN], *cp, name[STRLEN] ;
    FILE  *fp ;
    int   num ;

    printf("interpreting 4th command line arg as split file name\n") ;
    fp = fopen(argv[4], "r") ;
    if (fp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open parcellation division file %s", Progname,argv[4]) ;
    while ((cp = fgetl(line, STRLEN-1, fp)) != NULL) {
      if (sscanf(line, "%s %d", name, &num) != 2)
        ErrorExit(ERROR_BADFILE, "%s: could not parse name/num pair from '%s'", Progname, line) ;
      CTABfindName(mris->ct, name, &index) ;
      if (index < 0)
        ErrorExit(ERROR_BADFILE, "%s: could not find name '%s' in color table", Progname, name) ;
      nunits[index] = num ;
    }
  }

  for (index = 0 ; index < mris->ct->nentries ; index++)
    if (nunits[index] > 1)
      printf("dividing %s (%d) into %d parts\n", mris->ct->entries[index]->name, index, nunits[index]) ;

  MRISdivideAnnotation(mris, nunits) ;

  free(nunits) ;
  printf("saving annotation to %s\n", out_fname) ;
  MRISwriteAnnotation(mris, out_fname) ;
  return 0;
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
  if (!stricmp(option, "-help") || !stricmp(option, "help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "scale")) {
    printf("ERROR: cannot scale rgb values\n");
    rgb_scale = atoi(argv[2]) ;
    printf("scaling rgb values by %d\n", rgb_scale) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      usage_exit() ;
      break ;
    }
  return(nargs) ;
}
static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  printf("\n");
  printf("mris_divide_parcellation [options] subject hemi sourceannot [splitfile|areathresh] outannot\n");
  printf("\n");
  printf("\n");
  printf("options\n");
  printf( "  -help\n");
  printf( "  -scale <scale>   specify offset scaling for rgb values (default=20)\n");
  printf( "  -l <label name>  only process the label <label name> (not implemented yet)\n");
}

static void
print_help(void) {
  print_usage() ;
printf("\n");
printf("This program divides one or more parcellations into divisions\n");
printf("perpendicular to the long axis of the label.  The number of divisions\n");
printf("can be specified in one of two ways, depending upon the nature of the\n");
printf("fourth argument.\n");
printf("\n");
printf("First, a splitfile can be specified as the fourth argument. The\n");
printf("splitfile is a text file with two columns. The first column is the\n");
printf("name of the label in the source annotation, and the second column is\n");
printf("the number of units that label should be divided into. The names of\n");
printf("the labels depends upon the source parcellation.  For aparc.annot and\n");
printf("aparc.a2005.annot, the names can be found in\n");
printf("$FREESURFER_HOME/FreeSurferColorLUT.txt.  For aparc.annot, the labels\n");
printf("are between the ranges of 1000-1034.  For aparc.a2005s.annot, the\n");
printf("labels are between the ranges of 1100-1181.  The name for the label is\n");
printf("the name of the segmentation without the 'ctx-lh'. Note that the name\n");
printf("included in the splitfile does not indicate the hemisphere. For\n");
printf("example, 1023 is 'ctx-lh-posteriorcingulate'.  You should put\n");
printf("'posteriorcingulate' in the splitfile. Eg, to divide it into three\n");
printf("segments, the following line should appear in the splitfile:\n");
printf("\n");
printf("posteriorcingulate 3\n");
printf("\n");
printf("Only labels that should be split need be specified in the splitfile.\n");
printf("\n");
printf("The second method is to specify an area threshold (in mm^2) as the\n");
printf("fourth argument, in which case each label is divided until each\n");
printf("subdivision is smaller than the area threshold.\n");
printf("\n");
printf("The output label name will be the original name with _divN appended,\n");
printf("where N is the division number. N will go from 2 to the number of\n");
printf("divisions. The first division has the same name as the original label.\n");
printf("\n");
printf("\n");
printf("EXAMPLE:\n");
printf("\n");
printf("cd $SUBJECTS_DIR/mysubj001/label\n");
printf("\n");
printf("Create a split file to divide the superior frontal gyrus into 4\n");
printf("segements and the precentral gyrus into 3 segments:\n");
printf("\n");
printf("echo superiorfrontal 4 >  splittable.txt\n");
printf("echo precentral      3 >> splittable.txt\n");
printf("\n");
printf("Run program:\n");
printf("\n");
printf("mris_divide_parcellation mysubj001 rh aparc.annot splittable.txt rh.aparc.split-sfg+pc\n");
printf("\n");
printf("This reads in rh.aparc.annot, splits the SFG and PC, and and saves the\n");
printf("result to rh.aparc.split-sfg+pc\n");
printf("\n");
printf("View\n");
printf("\n");
printf("tksurfer mysubj001 rh inflated -annot aparc.split-sfg+pc\n");
printf("\n");
printf("The SFG divisions will have the following names: superiorfrontal,\n");
printf("superiorfrontal_div2, superiorfrontal_div3, superiorfrontal_div4. The\n");
printf("PC divisions will be precentral, precentral_div2, precentral_div3.\n");
printf("\n");
printf("\n");
printf("\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

