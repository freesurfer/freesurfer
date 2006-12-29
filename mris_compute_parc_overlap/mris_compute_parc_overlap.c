/**
 * @file  mris_compute_parc_overlap.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:10 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
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


/*
 * This is a simple program to compare two parcellated (annotated) surfaces
 * to compute the label Dice coefficients.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "mrisurf.h"
#include "annotation.h"
#include "version.h"
#include "hipsu.h"
#include "error.h"

static void usage(int exit_val);
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void dump_options(FILE *fp);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);

char *Progname;
static char vcid[] =
  "$Id: mris_compute_parc_overlap.c,v 1.2 2006/12/29 02:09:10 nicks Exp $";
static char *SUBJECTS_DIR = NULL;
static char *subject = NULL;
static char *hemi = NULL;
static char *annot1 = NULL;
static char *annot2 = NULL;
static int debug = 0;

int main(int argc, char *argv[]) {
  MRIS *white1, *white2;
  int nargs,err,n;
  char tmpstr[2000];
  char annotfile1[2000];
  char annotfile2[2000];

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc--;
  argv++;
  if (argc == 0) usage(1);
  if (argc < 8) usage(1);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if (SUBJECTS_DIR==NULL) {
    printf
    ("ERROR: SUBJECTS_DIR not defined in environment or on command-line\n");
    exit(1);
  }

  /* ------ Load subject's white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  printf("\nReading %s white surface \n %s\n",hemi,tmpstr);
  white1 = MRISread(tmpstr);
  if (white1 == NULL) {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  white2 = MRISread(tmpstr);

  /* ------ Load subject's native annotation ------ */
  sprintf(annotfile1,"%s/%s/label/%s.%s.annot",
          SUBJECTS_DIR,subject,hemi,annot1);
  printf("\nLoading %s annotations from %s\n",hemi,annotfile1);
  err = MRISreadAnnotation(white1, annotfile1);
  if (err) {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile1);
    exit(1);
  }

  /* ------ Load the resampled annotation (from mri_surf2surf) ------ */
  sprintf(annotfile2,"%s/%s/label/%s.%s.annot",
          SUBJECTS_DIR,subject,hemi,annot2);
  printf("\nLoading rh annotations from %s\n",annotfile2);
  err = MRISreadAnnotation(white2, annotfile2);
  if (err) {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile2);
    exit(1);
  }

  /* make sure files contain embedded color look-up table. */
  if (white1->ct == NULL) {
    printf("ERROR: %s does not contain color table!\n",annot1);
    exit(1);
  }
  if (white2->ct == NULL) {
    printf("ERROR: %s does not contain color table!\n",annot2);
    exit(1);
  }
  char *ctabname1 = white1->ct->fname;
  char *ctabname2 = white2->ct->fname;
  if (strcmp(ctabname1,ctabname2) != 0) {
    printf("ERROR: annotation files based on different colortables:\n");
    printf("\tannot1: %s\n",ctabname1);
    printf("\tannot2: %s\n",ctabname2);
    exit(1);
  }
  if (white1->nvertices != white2->nvertices) {
    printf("ERROR: (white1->nvertices=%d != white2->nvertices=%d)",
           white1->nvertices,white2->nvertices);
    exit(1);
  }

  /*
     get annot info on the label 'unknown', so that we can use it to
     create an .annot file containing vertices were perfect overlap occurs
     which are labeled as unknown, making it easy to see where
     the failure in overlap occurs.
  */
  int unknown_index=0;
  int unknown_annot=0;
  CTABfindName(white1->ct, "unknown", &unknown_index);
  if (unknown_index == -1) {
    printf
    ("ERROR: could not retrieve index for label 'unknown'\n");
    exit(1);
  }//else printf("label 'unknown' has index %d\n",unknown_index);
  err=CTABannotationAtIndex(white1->ct,unknown_index,&unknown_annot);
  if (err != NO_ERROR) {
    printf
    ("ERROR: could not retrieve annotation for label 'unknown'\n");
    exit(1);
  }//else printf("label 'unknown' has annotation 0x%8.8X\n",unknown_annot);

  /* ------- Do vertex-by-vertex comparison ------- */
  if (debug) printf("checking %d vertices...\n",white1->nvertices);
  fflush(stdout);
  err=0;
  int surfannot_overlap = 0;
  int surfannot1 = 0;
  int surfannot2 = 0;
  for (n = 0; n < white1->nvertices; n++) {
    /*
     * surface 1
     */
    int structure1 = 0;
    VERTEX *v1 = &white1->vertices[n];
    /* Get the color and then the index. */
    if ( 0 != v1->annotation ) {
      //printf("v1->annotation=0x%8.8X\n",v1->annotation);
      CTABfindAnnotation (white1->ct, v1->annotation, &structure1);
    }
    /*
     * surface 2
     */
    int structure2 = 0;
    VERTEX *v2 = &white2->vertices[n];
    /* Get the color and then the index. */
    if ( 0 != v2->annotation ) {
      CTABfindAnnotation (white2->ct, v2->annotation, &structure2);
    }
    /*
     * compare
     */
    //printf("structure1=%d, structure2=%d\n", structure1, structure2);
    if (structure1 != structure2) {
      if (debug)
        printf("structure1=%d != structure2=%d\n", structure1, structure2);
      err++;
    }
    if (structure1 != 0) {
      if (structure1 == structure2) surfannot_overlap++;
      surfannot1++;
      surfannot2++;
    }

    /*
     * create an overlap file, were vertices where structures are equal
     * are relabelled as unknown, thus leaving unequal structures labelled
     * as-is, allowing visualization of overlap failure.
     */
    if (structure1 == structure2) {
      v2->annotation = unknown_annot;
    }
  }

  printf("Found %d mismatches out of %d vertices\n",err,white1->nvertices);

  sprintf(tmpstr,"%s/%s/label/%s.overlap.annot",
          SUBJECTS_DIR,subject,hemi);
  printf("Writing %s\n",tmpstr);
  MRISwriteAnnotation(white2,tmpstr);

  printf("Overall Dice = %1.4f \n",
         surfannot_overlap*2.0/(float)(surfannot1 + surfannot2));

  exit(0);

}  /*  end main()  */





/* --------------------------------------------- */
static void usage(int exit_val) {
  FILE *fout;
  char progname[] = "mris_compute_parc_overlap";

  fout = (exit_val ? stderr : stdout);

  fprintf
  (fout,
   "Compares two parcellated (annotated) surfaces\n"
   "and computes and overall Dice coefficient.\n\n") ;
  fprintf
  (fout,
   "Usage:\n"
   "  %s --s subject --hemi hemi \\ \n"
   "    --annot1 annotfile --annot2 annotfile\n\n",
   progname);
  fprintf(fout, "Required:\n");
  fprintf(fout, "  --s subject          subject to check\n");
  fprintf(fout, "  --hemi hemi          hemisphere: rh or lh\n");
  fprintf(fout, "  --annot1 annotfile   first .annot file\n");
  fprintf(fout, "  --annot2 annotfile   second .annot file\n");
  fprintf(fout, "\nOptional:\n");
  fprintf(fout, "  --sd subj_dir        set SUBJECTS_DIR\n");
  fprintf(fout, "  --version            version info\n");
  fprintf(fout, "  --help               this usage info\n");
  fprintf(fout, "\nExample:\n");
  fprintf
  (fout,
   "  %s --s bert --hemi lh \\"
   "\n    --annot1 aparc --annot2 aparc.ernie\n",
   progname);
  fprintf
  (fout,
   "\nIn this example, the annotation file named lh.aparc.annot, which is\n"
   "created by the utility mris_ca_label (executed during the -autorecon3\n"
   "stage of recon-all), is compared against the annotation file named \n"
   "lh.aparc.ernie.annot.  This second annotation file is created by the\n"
   "utility mri_surf2surf, which resamples one surface onto another. This\n"
   "resampling is necessary so that a vertex-by-vertex comparison is\n"
   "meaningful.  An example command-line is: \n"
   "  mri_surf2surf --srcsubject ernie --trgsubject bert --hemi lh \\ \n"
   "    --sval-annot $SUBJECTS_DIR/ernie/label/lh.aparc.annot \\ \n"
   "    --tval       $SUBJECTS_DIR/bert/label/lh.ernie.aparc.annot\n\n"
   "Note that the resampling output file, lh.ernie.annot, is deposited\n"
   "in the label directory of the subject (bert) supplied as the input\n"
   "to the %s utility.  Supply --help to\n"
   "mri_surf2surf for its usage information.\n\n"
   "The output of %s is an overall Dice coefficient.\n"
   "A value of 1 indicates perfect overlap.\n\n"
   "A file called ?h.overlap.annot is created in the subjects label\n"
   "directory, and is a copy of the annot2 input, except wherever the\n"
   "labels are identical with annot1, that label is replaced with the\n"
   "label 'unknown', thus leaving any mismatches labeled as they were\n"
   "in the annot2 file.  If the Dice coefficient is less than 1, then this\n"
   "file is a handy way to visualize the mismatches, which are typically\n"
   "around the label borders."
   "\n",progname,progname);

  exit(exit_val);
}  /*  end usage()  */


/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage(1);

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
    else if (!strcmp(option, "--sd")) {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
      nargsused = 1;
    } else if (!strcmp(option, "--s")) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annot1")) {
      if (nargc < 1) argnerr(option,1);
      annot1 = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annot2")) {
      if (nargc < 1) argnerr(option,1);
      annot2 = pargv[0];
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


/* --------------------------------------------- */
static void print_help(void) {
  usage(1);
}


/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}


/* --------------------------------------------- */
static void check_options(void) {
  if (subject == NULL) {
    printf("ERROR: must specify a subject\n");
    exit(1);
  }
  if (hemi == NULL) {
    printf("ERROR: must specify a hemisphere (rh or lh)\n");
    exit(1);
  }
  if (annot1 == NULL) {
    printf("ERROR: must specify an annotation file\n");
    exit(1);
  }
  if (annot2 == NULL) {
    printf("ERROR: must specify the second annotation file\n");
    exit(1);
  }
  return;
}


/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\nSUBJECTS_DIR: %s\n",SUBJECTS_DIR);
  fprintf(fp,"subject:      %s\n",subject);
  fprintf(fp,"hemi:         %s\n",hemi);
  fprintf(fp,"annot1:       %s\n",annot1);
  fprintf(fp,"annot2:       %s\n",annot2);
  return;
}


/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}

/*  EOF  */
