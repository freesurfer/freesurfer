/*----------------------------------------------------------
  Name: mri_annotation2label.c
  $Id: mri_annotation2label.c,v 1.9.8.1 2006/09/19 17:15:08 greve Exp $
  Author: Douglas Greve
  Purpose: Converts an annotation to a labels.

  -----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string.h>

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "label.h"
#include "annotation.h"
#include "version.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_annotation2label.c,v 1.9.8.1 2006/09/19 17:15:08 greve Exp $";
char *Progname = NULL;

char  *subject   = NULL;
char  *annotation = "aparc";
char  *hemi       = NULL;
char  *surfacename = "white";

char  *outdir = NULL;
char  *labelbase = NULL;

MRI_SURFACE *Surf;
LABEL *label     = NULL;

int debug = 0;

char *SUBJECTS_DIR = NULL;
FILE *fp;

char tmpstr[2000];
char annotfile[1000];
char labelfile[1000];
int  nperannot[1000];


/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
int main(int argc, char **argv)
{
  VERTEX *vtx;
  int nthpoint,err,vtxno,ano,ani,vtxani,animax;
  int annotnum;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_annotation2label.c,v 1.9.8.1 2006/09/19 17:15:08 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  annotation_table_file = NULL; /* defined in annoation.h */

  printf("\n");

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  /*--- Get environment variables ------*/
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  /* ------ Load subject's surface ------ */
  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfacename); 
  printf("Reading surface \n %s\n",tmpstr);
  Surf = MRISread(tmpstr);
  if(Surf == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }

  /* ------ Load annotation ------ */
  sprintf(annotfile,"%s/%s/label/%s.%s.annot",
    SUBJECTS_DIR,subject,hemi,annotation);
  printf("Loading annotations from %s\n",annotfile);
  err = MRISreadAnnotation(Surf, annotfile);
  if(err){
    printf("INFO: could not load from %s, trying ",annotfile);
    sprintf(annotfile,"%s/%s/label/%s_%s.annot",
	    SUBJECTS_DIR,subject,hemi,annotation);
    printf("%s\n",annotfile);
    err = MRISreadAnnotation(Surf, annotfile);
    if(err){
      printf("ERROR: MRISreadAnnotation() failed\n");
      exit(1);
    }
    printf("OK, that worked.\n");
  }

  /* Loop through each vertex */
  animax = -1;
  for(vtxno = 0; vtxno < Surf->nvertices; vtxno++){
    vtx = &(Surf->vertices[vtxno]);
    ano = Surf->vertices[vtxno].annotation;
    vtxani = annotation_to_index(ano);
    nperannot[vtxani] ++;    
    if(animax < vtxani) animax = vtxani;
  }

  printf("animax = %d\n",animax);
  for(ani=0; ani <= animax; ani++){

    if(nperannot[ani] == 0){
      printf("%3d  %5d   --- skipping \n",ani,nperannot[ani]);
      continue;
    }

    if(labelbase != NULL) 
      sprintf(labelfile,"%s-%03d.label",labelbase,ani);
    if(outdir != NULL){
      annotnum = index_to_annotation(ani);
      sprintf(labelfile,"%s/%s.%s.label",outdir,hemi,
	      annotation_to_name(annotnum,NULL));
    }

    printf("%3d  %5d %s\n",ani,nperannot[ani],labelfile);

    label = LabelAlloc(nperannot[ani],subject,labelfile);
    label->n_points = nperannot[ani];

    nthpoint = 0;
    for(vtxno = 0; vtxno < Surf->nvertices; vtxno++){
      vtx = &(Surf->vertices[vtxno]);
      ano = Surf->vertices[vtxno].annotation;
      vtxani = annotation_to_index(ano);
      if(vtxani == ani){
	label->lv[nthpoint].vno = vtxno;
	label->lv[nthpoint].x = vtx->x;
	label->lv[nthpoint].y = vtx->y;
	label->lv[nthpoint].z = vtx->z;
	nthpoint ++;
      }
    }

    LabelWrite(label,labelfile);
    LabelFree(&label);
  }


  return(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcmp(option, "--help"))  print_help() ;

    else if (!strcmp(option, "--version")) print_version() ;

    else if (!strcmp(option, "--debug"))   debug = 1;

    /* -------- source inputs ------ */
    else if (!strcmp(option, "--subject")){
      if(nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--surface") || !strcmp(option, "--surf")){
      if(nargc < 1) argnerr(option,1);
      surfacename = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelbase")){
      if(nargc < 1) argnerr(option,1);
      labelbase = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--outdir")){
      if(nargc < 1) argnerr(option,1);
      outdir = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annotation")){
      if(nargc < 1) argnerr(option,1);
      annotation = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--table")){
      if(nargc < 1) argnerr(option,1);
      annotation_table_file = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--hemi")){
      if(nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
  fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --subject    source subject\n");
  printf("   --hemi       hemisphere (lh or rh) (with surface)\n");
  printf("   --table      annotation table file\n");  
  printf("   --labelbase  output will be base-XXX.label \n");
  printf("   --outdir dir :  output will be dir/hemi.name.label \n");
  printf("\n");
  printf("   --annotation as found in SUBJDIR/labels <aparc>\n");
  printf("   --surface    name of surface <white>\n");  

  printf("   --help       display help\n");  
  printf("   --version    display version\n");  
  printf("\n");
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"subject = %s\n",subject);
  fprintf(fp,"annotation = %s\n",annotation);
  fprintf(fp,"hemi = %s\n",hemi);
  if(labelbase) fprintf(fp,"labelbase = %s\n",labelbase);
  if(outdir) fprintf(fp,"outdir = %s\n",labelbase);
  fprintf(fp,"surface   = %s\n",surfacename);
  fprintf(fp,"\n");
  return;
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("This program will convert an annotation into multiple label files.\n");

  printf(
"User specifies the subject, hemisphere, label base, and (optionally)\n"
"the annotation base and surface. By default, the annotation base is\n"
"aparc. The program will retrieves the annotations from\n"
"SUBJECTS_DIR/subject/label/hemi_annotbase.annot. A separate label file\n"
"is created for each annotation index. The output file names can take\n"
"one of two forms: (1) If --outdir dir is used, then the output will\n"
"be dir/hemi.name.lable, where name is the corresponding name in \n"
"the table. (2) If --labelbase is used, name of the file conforms to\n"
"labelbase-XXX.label, where XXX is the zero-padded 3 digit annotation\n"
"index. If labelbase includes a directory path, that directory must\n"
"already exist. If there are no points in the annotation file for a\n"
"particular index, no label file is created. The xyz coordinates in the\n"
"label file are obtained from the values at that vertex of the\n"
"specified surface. The default surface is 'white'. Other options\n"
"include 'pial' and 'orig'.\n"
"\n"
"The human-readable names that correspond\n"
"to the annotation indices for aparc depend upon how the annotations\n"
"were created. They are created with the program mris_ca_label, and\n"
"the human readable names are in the file given as the argument to \n"
"the -t flag. Unfortunately, this information is not maintained inside\n"
"the annotation file, and it must be supplied to mri_annotation2label \n"
"through the --table flag.\n"
" \n"
"Bugs:\n"
"\n"
"  If the name of the label base does not include a forward slash (ie, '/') \n"
"  then the program will attempt to put the label files in \n"
"  $SUBJECTS_DIR/subject/label.  So, if you want the labels to go into the \n"
"  current directory, make sure to put a './' in front of the label base.\n"
"\n"
"Example:\n"
"\n"
"  mri_annotation2label --subject LW --hemi rh \n"
"        --labelbase ./labels/aparc-rh \n"
"        --table $FREESURFER_HOME/Simple_surface_labels2002.txt \n"
"\n"
"  This will get annotations from $SUBJECTS_DIR/LW/label/rh_aparc.annot\n"
"  and then create about 94 label files: aparc-rh-001.label, \n"
"  aparc-rh-002.label, ... Note that the directory 'labels' must already\n"
"  exist. \n"
"\n"
"Example:\n"
"\n"
"  mri_annotation2label --subject LW --hemi rh \n"
"        --outdir ./labels \n"
"        --table $FREESURFER_HOME/Simple_surface_labels2002.txt \n"
"\n"
"  This will do the same thing as above except that the output files\n"
"  will have names of the form lh.S_occipital_anterior.label\n"
"\n"
"Testing:\n"
"\n"
"  1. Start tksurfer:  \n"
"       tksurfer -LW lh inflated\n"
"       read_annotations lh_aparc\n"
"     When a point is clicked on, it prints out a lot of info, including\n"
"     something like:\n"
"       annot = S_temporalis_sup (93, 3988701) (221, 220, 60)\n"
"     This indicates that annotion number 93 was hit. Save this point.\n"
"   \n"
"  2. Start another tksurfer and load the label:\n"
"       tksurfer -LW lh inflated\n"
"       [edit label field and hit the 'Read' button]\n"
"     Verify that label pattern looks like the annotation as seen in\n"
"     the tksurfer window from step 1.\n"
"\n"
"  3. Load label into tkmedit\n"
"       tkmedit LW T1\n"
"       [Load the label]\n"
"      [Goto the point saved from step 1] \n\n\n");

  //printf("Annotation Key\n");
  //print_annotation_table(stdout);

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(subject == NULL){
    fprintf(stderr,"ERROR: no source subject specified\n");
    exit(1);
  }
  if(hemi == NULL){
    fprintf(stderr,"ERROR: No hemisphere specified\n");
    exit(1);
  }
  if(annotation_table_file == NULL){
    fprintf(stderr,"ERROR: No annotation table specified\n");
    exit(1);
  }
  if(outdir == NULL && labelbase == NULL){
    fprintf(stderr,"ERROR: no output specified\n");
    exit(1);
  }
  if(outdir != NULL && labelbase != NULL){
    fprintf(stderr,"ERROR: cannot specify both --outdir and --labelbase\n");
    exit(1);
  }

  return;
}

/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
