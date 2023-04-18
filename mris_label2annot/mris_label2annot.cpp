/**
 * @brief rogram to convert one or more labels into an annotation
 *
 * Converts a set of surface labels to an annotation file.
 */
/*
 * Original Author: Doug Greve
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
  BEGINHELP

Converts a set of surface labels to an annotation file.

--s subject

Name of FreeSurfer subject.

--h hemi

Hemisphere (lh or rh).

--ctab colortablefile

File that defines the structure names, their indices, and their
color. This must have the same format as that found in
$FREESUFER_HOME/FreeSurferColorLUT.txt. This can be used to generate
the names of the label files (see --l below).

--l labelfile1 <--l labelfile2 ...>

List of label files. If no label files are specified, then the label
file name is constructed from the list in the ctab as
hemi.parcname.label.  The labels should be defined on the surface (eg,
with tksurfer). The first label will be mapped to index 1 in the color
table file. The next label will be mapped to index 2, etc. Verticies
that are not mapped to a label are assigned index 0. If --no-unknown
is specified, then the first label is mapped to index 0, etc, and
unhit vertices are not mapped.

--dilate_into_unknown label

dilates <label> into bordering vertices labeled unknown 

--ldir labeldir

When getting label file names from the ctab, find the actual label files
in ldir. If a label file based on the name in the ctab does not exist,
it is skipped.

--a annotname

Name of the annotation to create. The actual file will be called
hemi.annotname.annot, and it will be created in subject/label.
If this file exists, then mris_label2annot exits immediately
with an error message. It is then up to the user to manually
delete this file (this is so annotations are not accidentally
deleted, which could be a huge inconvenience).

--nhits nhitsfile

Overlay file with the number of labels for each vertex. Ideally, each
vertex would have only one label, but there is not constraint that
forces this. If there is more than one label for a vertex, the vertex
will be assigned to the last label as specified on the cmd line. This
can then be loaded as an overlay in tksurfer (set fthresh to 1.5). This
is mainly good for debugging.

--no-unknown

Start label numbering at index 0 instead of index 1. Do not map unhit
vertices (ie, vertices without a label) to 0.

--thresh threshold

Require that the stat field of the vertex in the label be greather 
than threshold.

EXAMPLE:

mris_label2annot --s bert --h lh --ctab aparc.annot.ctab \\
  --a myaparc --l lh.unknown.label --l lh.bankssts.label \\
  --l lh.caudalanteriorcingulate.label --nhits nhits.mgh

This will create lh.myaparc.annot in bert/labels using the three
labels specified. Any vertices that have multiple labels will then
be stored in nhits.mgh (as a volume-encoded surface file).

To view, run:

tksurfer bert lh inflated -overlay nhits.mgh -fthresh 1.5

Then File->Label->ImportAnnotation and select lh.myaparc.annot.

EXAMPLE: 

To create an annotation with a few labels from the aparc, run

cd $SUBJECTS_DIR/yoursubject/label
mri_annotation2label --hemi lh --subject yoursubject --outdir deleteme
rm deleteme/lh.superiortemporal.label   # remove superior temporal, etc
mris_label2annot --hemi lh --subject yoursubject --ctab aparc.annot.ctab
   --ldir deletme --no-unknown --a myaparc
tksurferfv yoursubject lh inflated -annot myaparc
rm -r deletme

  ENDHELP
*/

/*
  BEGINUSAGE
  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static void label2annotation();
static void label2annotationV2();

static int dilate_label_into_unknown(MRI_SURFACE *mris, int annot) ;
static char *dilate_label_name = NULL ;
static int dilate_label_index = -1 ;
static int dilate_label_annot = 0 ;
const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
int verbose=1; // print overlap warnings and maxstat overrides
struct utsname uts;

char tmpstr[1000];
char *subject, *hemi, *SUBJECTS_DIR;
char *LabelFiles[1000];
int nlabels = 0;
char *CTabFile = NULL;
char *newCTabFile = NULL;
char *AnnotName=NULL, *AnnotPath=NULL;
MRIS *mris;
LABEL *label;
COLOR_TABLE *ctab = NULL, *ctab2 = NULL;
MRI *nhits;
char *NHitsFile=NULL;
MRI *maxstat;
int maxstatwinner=0;
int MapUnhitToUnknown=1;
const char *labeldir=NULL;
int labeldirdefault=0;
int DoLabelThresh = 0;
double LabelThresh = 0;
const char *surfname = "orig";
int IndexOffset=0;
/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs = handleVersionOption(argc, argv, "mris_label2annot");
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
  dump_options(stdout);

  // Make sure subject exists
  sprintf(tmpstr,"%s/%s",SUBJECTS_DIR,subject);
  if (!fio_IsDirectory(tmpstr)) {
    printf("ERROR: cannot find %s\n",tmpstr);
    exit(1);
  }

  if(AnnotPath == NULL){
    // Get path to annot, make sure it does not exist
    sprintf(tmpstr,"%s/%s/label/%s.%s.annot",
	    SUBJECTS_DIR,subject,hemi,AnnotName);
    if (fio_FileExistsReadable(tmpstr)) {
      printf("ERROR: %s already exists\n",tmpstr);
      exit(1);
    }
    AnnotPath = strcpyalloc(tmpstr);
  }

  // Read the surf
  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  printf("Loading %s\n",tmpstr);
  mris = MRISread(tmpstr);
  if (mris == NULL) exit(1);

  // Set up color table
  set_atable_from_ctable(ctab);
  mris->ct = ctab;
  //CTABwriteFileASCII(mris->ct, "new2.ctab");

  // Set up something to keep track of nhits
  nhits = MRIalloc(mris->nvertices,1,1,MRI_INT);

  // Set up something to keep track of max stat for that vertex
  if (maxstatwinner) maxstat = MRIalloc(mris->nvertices,1,1,MRI_FLOAT);

  label2annotation();

  int nunhit = 0;
  if (MapUnhitToUnknown) {
    printf("Mapping unhit to unknown\n");
    for (int vtxno = 0; vtxno < mris->nvertices; vtxno++) {
      if (MRIgetVoxVal(nhits,vtxno,0,0,0) == 0) {
        int ano = index_to_annotation(0);
        mris->vertices[vtxno].annotation = ano;
        nunhit ++;
      }
    }
    printf("Found %d unhit vertices\n",nunhit);
  }

  if (dilate_label_name)
  {
    dilate_label_into_unknown(mris, dilate_label_annot) ;
  }
  printf("Writing annot to %s\n",AnnotPath);
  MRISwriteAnnotation(mris, AnnotPath);

  if (NHitsFile != NULL) MRIwrite(nhits,NHitsFile);

  return 0;
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

    if (!strcmp(option, "--help"))  print_help() ;
    else if (!strcmp(option, "--version")) print_version() ;
    else if (!strcmp(option, "--debug"))   debug = 1;
    else if (!strcmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcmp(option, "--no-unknown")) MapUnhitToUnknown = 0;
    else if (!strcmp(option, "--ldir-default")) labeldirdefault = 1;
    else if (!strcmp(option, "--noverbose")) verbose = 0;
    else if (!strcmp(option, "--maxstatwinner")) maxstatwinner = 1;
    else if (!strcmp(option, "--dilate-into-unknown")) 
    {
      dilate_label_name = pargv[0] ;
      printf("dilating label %s into bordering unknown vertices\n", 
             dilate_label_name) ;
      nargsused = 1;
    }
    else if (!strcmp(option, "--s") || !strcmp(option, "--subject")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--h") || !strcmp(option, "--hemi")) {
      if (nargc < 1) CMDargNErr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--ctab")) {
      if (nargc < 1) CMDargNErr(option,1);
      if (!fio_FileExistsReadable(pargv[0])) {
        printf("ERROR: cannot find or read %s\n",pargv[0]);
        exit(1);
      }
      CTabFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--prn-ctab")) {
      if (nargc < 1) CMDargNErr(option,1);
      //if (fio_FileExistsReadable(pargv[0])) {
      //  printf("ERROR: %s already exists\n", pargv[0]);
      //  exit(1);
      //}
      newCTabFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--l")) {
      if (nargc < 1) CMDargNErr(option,1);
      if (!fio_FileExistsReadable(pargv[0])) {
        printf("ERROR: cannot find or read %s\n",pargv[0]);
        exit(1);
      }
      LabelFiles[nlabels] = pargv[0];
      nlabels++;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&LabelThresh);
      DoLabelThresh = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--offset")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&IndexOffset);
      DoLabelThresh = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sd")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--ldir")) {
      if (nargc < 1) CMDargNErr(option,1);
      labeldir = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--a") || !strcmp(option, "--annot")) {
      if (nargc < 1) CMDargNErr(option,1);
      AnnotName = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annot-path")) {
      if (nargc < 1) CMDargNErr(option,1);
      AnnotPath = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--nhits")) {
      if (nargc < 1) CMDargNErr(option,1);
      NHitsFile = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --s subject : FreeSurfer subject \n");
  printf("   --h hemi : hemisphere (lh or rh)\n");
  printf("   --ctab ctabfile : colortable (like FreeSurferColorLUT.txt)\n");
  printf("   --offset IndexOffset : add to label number to get CTAB index\n");
  printf("   --l label1 <--l label 2 ...> : label file(s)\n");
  printf("   --a annotname : output annotation file (hemi.annotname.annot)\n");
  printf("   --annot-path annotpath : full name/path of annotation file\n");
  printf("   --ldir labeldir : when not using --l\n");
  printf("   --ldir-default : use subject/labels as labeldir\n");
  printf("   --no-unknown : do not map unhit labels to index 0\n");
  printf("   --thresh thresh : threshold label by stats field\n");
  printf("   --maxstatwinner : keep label with highest 'stat' value\n");
  printf("   --surf surfname : default is orig\n");
  printf("   --sd SUBJECTS_DIR\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --noverbose turn off overlap and stat override messages\n");
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
printf("\n");
printf("Converts a set of surface labels to an annotation file.\n");
printf("\n");
printf("--s subject\n");
printf("\n");
printf("Name of FreeSurfer subject.\n");
printf("\n");
printf("--h hemi\n");
printf("\n");
printf("Hemisphere (lh or rh).\n");
printf("\n");
printf("--ctab colortablefile\n");
printf("\n");
printf("File that defines the structure names, their indices, and their\n");
printf("color. This must have the same format as that found in\n");
printf("$FREESUFER_HOME/FreeSurferColorLUT.txt. This can be used to generate\n");
printf("the names of the label files (see --l below).\n");
printf("\n");
printf("--l labelfile1 <--l labelfile2 ...>\n");
printf("\n");
printf("List of label files. If no label files are specified, then the label\n");
printf("file name is constructed from the list in the ctab as\n");
printf("hemi.parcname.label.  The labels should be defined on the surface (eg,\n");
printf("with tksurfer). The first label will be mapped to index 1 in the color\n");
printf("table file. The next label will be mapped to index 2, etc. Verticies\n");
printf("that are not mapped to a label are assigned index 0. If --no-unknown\n");
printf("is specified, then the first label is mapped to index 0, etc, and\n");
printf("unhit vertices are not mapped.\n");
printf("\n");
printf("--dilate_into_unknown label\n");
printf("\n");
printf("dilates <label> into bordering vertices labeled unknown \n");
printf("\n");
printf("--ldir labeldir\n");
printf("\n");
printf("When getting label file names from the ctab, find the actual label files\n");
printf("in ldir. If a label file based on the name in the ctab does not exist,\n");
printf("it is skipped.\n");
printf("\n");
printf("--a annotname\n");
printf("\n");
printf("Name of the annotation to create. The actual file will be called\n");
printf("hemi.annotname.annot, and it will be created in subject/label.\n");
printf("If this file exists, then mris_label2annot exits immediately\n");
printf("with an error message. It is then up to the user to manually\n");
printf("delete this file (this is so annotations are not accidentally\n");
printf("deleted, which could be a huge inconvenience).\n");
printf("\n");
printf("--nhits nhitsfile\n");
printf("\n");
printf("Overlay file with the number of labels for each vertex. Ideally, each\n");
printf("vertex would have only one label, but there is not constraint that\n");
printf("forces this. If there is more than one label for a vertex, the vertex\n");
printf("will be assigned to the last label as specified on the cmd line. This\n");
printf("can then be loaded as an overlay in tksurfer (set fthresh to 1.5). This\n");
printf("is mainly good for debugging.\n");
printf("\n");
printf("--no-unknown\n");
printf("\n");
printf("Start label numbering at index 0 instead of index 1. Do not map unhit\n");
printf("vertices (ie, vertices without a label) to 0.\n");
printf("\n");
printf("--thresh threshold\n");
printf("\n");
printf("Require that the stat field of the vertex in the label be greather \n");
printf("than threshold.\n");
printf("\n");
printf("EXAMPLE:\n");
printf("\n");
printf("mris_label2annot --s bert --h lh --ctab aparc.annot.ctab \\\n");
printf("  --a myaparc --l lh.unknown.label --l lh.bankssts.label \\\n");
printf("  --l lh.caudalanteriorcingulate.label --nhits nhits.mgh\n");
printf("\n");
printf("This will create lh.myaparc.annot in bert/labels using the three\n");
printf("labels specified. Any vertices that have multiple labels will then\n");
printf("be stored in nhits.mgh (as a volume-encoded surface file).\n");
printf("\n");
printf("To view, run:\n");
printf("\n");
printf("tksurfer bert lh inflated -overlay nhits.mgh -fthresh 1.5\n");
printf("\n");
printf("Then File->Label->ImportAnnotation and select lh.myaparc.annot.\n");
printf("\n");
printf("EXAMPLE: \n");
printf("\n");
printf("To create an annotation with a few labels from the aparc, run\n");
printf("\n");
printf("cd $SUBJECTS_DIR/yoursubject/label\n");
printf("mri_annotation2label --hemi lh --subject yoursubject --outdir deleteme\n");
printf("rm deleteme/lh.superiortemporal.label   # remove superior temporal, etc\n");
printf("mris_label2annot --hemi lh --subject yoursubject --ctab aparc.annot.ctab \n");
printf("   --ldir deletme --no-unknown --a myaparc\n");
printf("tksurferfv yoursubject lh inflated -annot myaparc\n");
printf("rm -r deletme\n");
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
  int n;

  if (subject == NULL) {
    printf("ERROR: need to specify subject\n");
    exit(1);
  }
  if (hemi == NULL) {
    printf("ERROR: need to specify hemi\n");
    exit(1);
  }
  if (CTabFile == NULL) {
    printf("ERROR: need to specify color table file\n");
    exit(1);
  }
  if(AnnotName && AnnotPath) {
    printf("ERROR: cannot spec both --annot and --annot-path\n");
    exit(1);
  }
  if(AnnotName == NULL && AnnotPath == NULL) {
    printf("ERROR: need to specify annotation name\n");
    exit(1);
  }

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  // Read the color table
  printf("Reading ctab %s\n",CTabFile);
  ctab = CTABreadASCII(CTabFile);
  if (ctab == NULL) {
    printf("ERROR: reading %s\n",CTabFile);
    exit(1);
  }
  if (dilate_label_name)
  {
    CTABfindName(ctab, dilate_label_name, &dilate_label_index) ;
    CTABannotationAtIndex(ctab, dilate_label_index, &dilate_label_annot);
    printf("label %s maps to index %d, annot %x\n",
           dilate_label_name, dilate_label_index, dilate_label_annot) ;
  }
  printf("Number of ctab entries %d\n",ctab->nentries);
  if (nlabels == 0) {
    printf("INFO: no labels specified, generating from ctab\n");
    if (labeldir == NULL) labeldir = ".";
    if (labeldirdefault) {
      sprintf(tmpstr,"%s/%s/label",SUBJECTS_DIR,subject);
      labeldir = strcpyalloc(tmpstr);
    }
    // ctab2 is a shorter version of ctab
    // each label in ctab2 has its corresponding label file in labeldir
    ctab2 = CTABalloc(ctab->nentries);
    nlabels = 0;  // ??? should we start from 1 ??? 0 is reserved for unknown/annotation=0
    for (n=0; n<ctab->nentries; n++) {
      if(ctab->entries[n] == NULL) continue;
      if (strlen(ctab->entries[n]->name) == 0) continue;
      int req = snprintf(tmpstr,1000,"%s/%s.%s.label",labeldir,hemi,ctab->entries[n]->name); 
      if( req >= 1000 ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if(!fio_FileExistsReadable(tmpstr)) continue;
      printf("%2d %s\n",n,tmpstr);
      LabelFiles[nlabels] = strcpyalloc(tmpstr);
      strcpy(ctab2->entries[nlabels]->name,ctab->entries[n]->name);
      ctab2->entries[nlabels]->ri = ctab->entries[n]->ri;
      ctab2->entries[nlabels]->gi = ctab->entries[n]->gi;
      ctab2->entries[nlabels]->bi = ctab->entries[n]->bi;
      ctab2->entries[nlabels]->ai = ctab->entries[n]->ai;
      ctab2->entries[nlabels]->rf = ctab->entries[n]->rf;
      ctab2->entries[nlabels]->gf = ctab->entries[n]->gf;
      ctab2->entries[nlabels]->bf = ctab->entries[n]->bf;
      ctab2->entries[nlabels]->af = ctab->entries[n]->af;
      nlabels ++;
    }
    CTABfree(&ctab);
    ctab = ctab2;
    //CTABwriteFileASCII(ctab, "new.ctab");
  }
  return;
}


/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"\n");
  fprintf(fp,"subject %s\n",subject);
  fprintf(fp,"hemi    %s\n",hemi);
  fprintf(fp,"SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"ColorTable %s\n",CTabFile);
  if(AnnotName) fprintf(fp,"AnnotName  %s\n",AnnotName);
  if(AnnotPath) fprintf(fp,"AnnotPath  %s\n",AnnotPath);
  if (NHitsFile) fprintf(fp,"NHitsFile %s\n",NHitsFile);
  fprintf(fp,"nlables %d\n",nlabels);
  fprintf(fp,"LabelThresh %d %lf\n",DoLabelThresh,LabelThresh);
  return;
}
static int
dilate_label_into_unknown(MRI_SURFACE *mris, int annot) 
{
  int    vno, n, nchanged, iter ;

  iter = 0 ;
  do
  {
    nchanged = 0 ;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->annotation != annot)
        continue ;
      for (n = 0 ; n < vt->vnum ; n++)
      {
        VERTEX * const vn = &mris->vertices[vt->v[n]] ;
        if (vn->ripflag || vn->annotation != 0)
          continue ;

        nchanged++ ;
        vn->annotation = annot ;
      }
    }
    printf("iter %d: %d annots changed\n", iter+1, nchanged) ;
  } while (nchanged > 0 && iter++ < 1000) ;
  return(NO_ERROR) ;
}


/* The implementation assumes the label files are specified with '-l' in the same order in LUT.
 * There will be problems if label indexes in LUT are not continous.
 */
static void label2annotation()
{
  // Go thru each label
  printf("Index Offset %d\n",IndexOffset);
  for (int nthlabel = 0; nthlabel < nlabels; nthlabel ++) {
    printf("%d reading %s\n",nthlabel,LabelFiles[nthlabel]);    
    label = LabelRead(subject,LabelFiles[nthlabel]);
    if (label == NULL) {
      printf("ERROR: reading %s\n",LabelFiles[nthlabel]);
      exit(1);
    }
    int index = nthlabel+IndexOffset;
    if (MapUnhitToUnknown) index ++;
    int ano = index_to_annotation(index);
    printf("%2d %2d %s\n",index,ano,index_to_name(index));

    for (int n = 0; n < label->n_points; n++) {
      int vtxno = label->lv[n].vno;
      if (vtxno < 0 || vtxno > mris->nvertices) {
        printf("ERROR: %s, n=%d, vertex %d out of range\n",
               LabelFiles[nthlabel],n,vtxno);
        exit(1);
      }
      if(DoLabelThresh && label->lv[n].stat < LabelThresh) continue;

      if (maxstatwinner) {
        float stat = MRIgetVoxVal(maxstat,vtxno,0,0,0);
        if (label->lv[n].stat < stat) {
          if (verbose) {
            printf("Keeping prior label for vtxno %d "
                   "(old_stat=%f > this_stat=%f)\n",
                   vtxno,stat,label->lv[n].stat);
          }
          continue;
        }
        MRIsetVoxVal(maxstat,vtxno,0,0,0,label->lv[n].stat);
      }

      if (verbose) {
        if (MRIgetVoxVal(nhits,vtxno,0,0,0) > 0) {
          printf
            ("WARNING: vertex %d maps to multiple labels! (overwriting)\n",
             vtxno);
        }
      }
      MRIsetVoxVal(nhits,vtxno,0,0,0,MRIgetVoxVal(nhits,vtxno,0,0,0)+1);

      mris->vertices[vtxno].annotation = ano;
      //printf("%5d %2d %2d %s\n",vtxno,segid,ano,index_to_name(segid));
    } // label ponts
    LabelFree(&label);
  }// Label
}


/*
 * The implementation assumes the given label files [?h.]label-name.label indicate the labels intended.
 * When retrieving the label name, [?h.] is discarded if it is present.
 * If a label doesn't exist in CTAB, add it to CTAB with randomly generated RGB.
 */
static void label2annotationV2()
{
  printf("convert labels to annotation with label2annotationV2() ...\n");

  // create MRI (nvertices x (nlabels+4) x 1 x 1) to track labels assigned to each vertex
  // the MRI is used a 2D int array with these columns:
  //      col 0      col 1            col 2               col 3       col 4    ...   col (4+nlabels) 
  //    annotation  labelno.   nlabel-mapped-to-vertex   maxstat     labelno.  ...     labelno.
  //   
  //   col 0 and col 1 is the final annotation and labelno. assigned to the vertex;
  //   col 4 to (4+nlabels) record all labelno. assigned to the vertex
  //   * col 3 should be float
  //
  std::vector<int> shape{mris->nvertices, (nlabels+4), 1, 1};
  MRI *labelStat =  new MRI(shape, MRI_INT);
  if (labelStat == NULL)
  {
    printf("ERROR: failed to create labelStat MRI struct\n");
    exit(1);
  }

  // initialize col 4 to (4+nlabels) to -1 (not assigned to any labels)
  for (int r = 4; r < labelStat->height; r++)
    for (int c = 0; c < labelStat->width; c++)
      MRIsetVoxVal(labelStat, c, r, 0, 0, -1);

  // Go thru each label
  for (int nthlabel = 0; nthlabel < nlabels; nthlabel++) {
    printf("%d reading %s\n", nthlabel, LabelFiles[nthlabel]);    
    label = LabelRead(subject, LabelFiles[nthlabel]);
    if (label == NULL) {
      printf("ERROR: reading %s\n", LabelFiles[nthlabel]);
      exit(1);
    }

    // it assumes the given label files indicate the labels intended
    char *labelname = basename(LabelFiles[nthlabel]);
    char *ext = strrchr(labelname, '.');
    if (ext != NULL)
      *ext = '\0';

    // discard [?h.] if it is present
    char hemistr[16] = {'\0'};
    sprintf(hemistr, "%s.", hemi);
    if (strncmp(labelname, hemistr, strlen(hemistr)) == 0)
      labelname += strlen(hemistr);

    // look up the CTAB by label name instead of index that matching LabelFiles[] index
    int label_annot = 0, label_index = -1;
    CTABfindName(ctab, labelname, &label_index);
    if (label_index >= 0)
      CTABannotationAtIndex(ctab, label_index, &label_annot);
    else
    {
      // the label doesn't exist in CTAB
      // add it to CTAB with randomly generated RGB

      // set environment variable FREESURFER_SEED to have reproducible outcome
      // setenv("FREESURFER_SEED", "12", 1);
      setRandomSeed(12);  // gifti.cpp is using seed 12
      CTABaddUniqueEntryAtEnd(ctab, labelname, &label_index);
      if (label_index < -1)
      {
        printf("WARN: failed to add %s to colortab. It is not included in generated .annot!\n", labelname);
        continue;
      }

      printf("INFO: Added label %s to colortab index %d\n", labelname, label_index);
      CTABannotationAtIndex(ctab, label_index, &label_annot);
    }

    printf("%2d %2d %s (%d %d %d %d)\n",
           label_index, label_annot, ctab->entries[label_index]->name,
           ctab->entries[label_index]->ri, ctab->entries[label_index]->gi, ctab->entries[label_index]->bi, ctab->entries[label_index]->ai);

    for (int n = 0; n < label->n_points; n++) {
      int vtxno = label->lv[n].vno;
      if (vtxno < 0 || vtxno > mris->nvertices) {
        printf("ERROR: %s, n=%d, vertex %d out of range\n", LabelFiles[nthlabel], n, vtxno);
        exit(1);
      }

      if(DoLabelThresh && label->lv[n].stat < LabelThresh) continue;

      if (maxstatwinner) {
        float stat = MRIgetVoxVal(maxstat, vtxno, 0, 0, 0);
        if (label->lv[n].stat < stat) {
          if (verbose)
            printf("Keeping prior label for vtxno %d (old_stat=%f > this_stat=%f)\n", vtxno,stat,label->lv[n].stat);

          continue;
        }

        MRIsetVoxVal(maxstat, vtxno, 0, 0, 0, label->lv[n].stat);
      }

      //if (verbose) {
      //  if (MRIgetVoxVal(nhits,vtxno,0,0,0) > 0)
      //    printf("WARNING: vertex %d maps to multiple labels! (overwriting)\n", vtxno);
      //}

      MRIsetVoxVal(nhits, vtxno, 0, 0, 0, MRIgetVoxVal(nhits, vtxno, 0, 0, 0)+1);

      if (labelStat != NULL)
      {
        MRIsetVoxVal(labelStat, vtxno, 0, 0, 0, label_annot);
        MRIsetVoxVal(labelStat, vtxno, 1, 0, 0, nthlabel);
        MRIsetVoxVal(labelStat, vtxno, 2, 0, 0, MRIgetVoxVal(labelStat, vtxno, 2, 0, 0)+1);

        // this should be a float
        int stat =  MRIgetVoxVal(labelStat, vtxno, 3, 0, 0);
        MRIsetVoxVal(labelStat, vtxno, 3, 0, 0, (label->lv[n].stat > stat) ?  label->lv[n].stat : stat);

        MRIsetVoxVal(labelStat, vtxno, 4+nthlabel, 0, 0, nthlabel);
      }

      mris->vertices[vtxno].annotation = label_annot;
    } // label ponts
    LabelFree(&label);
  }// Label

  // debug
  if (newCTabFile != NULL)
    CTABwriteFileASCII(ctab, newCTabFile);

  if (labelStat != NULL && verbose)
  {
    int nunmapped = 0;

    // generate report from MRI* labelStat for vertex mapped to multiple labels
    printf("\n\n****** vertex mapped to multiple labels ******\n");
    for (int vtxno = 0; vtxno < mris->nvertices; vtxno++)
    {
      int nmapped = MRIgetVoxVal(labelStat, vtxno, 2, 0, 0);
      if (nmapped == 0)
      {
        nunmapped++;
      }
      else if (nmapped > 1)
      {
        int mapped_count = 0;
        char report[1024] = {'\0'};
        sprintf(report, "vertex #%8d mapped to %d labels: ", vtxno, nmapped);
        for (int nthlabel = 0; nthlabel < nlabels; nthlabel++)
	{
          int labelmapped = MRIgetVoxVal(labelStat, vtxno, 4+nthlabel, 0, 0);
          if (labelmapped >= 0)
	  {
            sprintf(report, "%s%s%d", report, (mapped_count) ? ", " : " ", labelmapped);
            mapped_count++;
          }
        }

        printf("%s\n", report);
      }
    }
    printf("Found %d unmapped vertices out of %d\n", nunmapped, mris->nvertices);
    printf("**********************************************\n\n");

    MRIfree(&labelStat);
  }
}
