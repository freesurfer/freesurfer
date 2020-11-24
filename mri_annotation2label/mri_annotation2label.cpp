/**
 * @brief Converts an annotation to labels files
 *
 */
/*
 * Original Author: Douglas Greve
 *
 * Copyright © 2011-2016 The General Hospital Corporation (Boston, MA) "MGH"
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "label.h"
#include "annotation.h"
#include "version.h"
#include "mri.h"
#include "fio.h"
#include "cma.h"

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

const char *Progname = NULL;

char  *subject   = NULL;
const char  *annotation = "aparc";
char  *hemi       = NULL;
const char  *surfacename = "white";

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

char *segfile=NULL;
MRI  *seg;
int  segbase = -1000;
char *ctabfile = NULL;

char *borderfile=NULL;
char *BorderAnnotFile=NULL;
char *LobesFile=NULL;
char *StatFile=NULL;
MRI  *Stat=NULL;
static int label_index = -1 ;  // if >= 0 only extract this label

// The global 'gi_lobarDivision'
typedef enum _lobarDivision {
  e_default, e_strict, e_strict_phcg
} e_LOBARDIVISION;
e_LOBARDIVISION Ge_lobarDivision = e_default;

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
int main(int argc, char **argv)
{
  int err,vtxno,ano,ani,vtxani,animax;
  int nargs,IsBorder,k;
  MRI *border;

  nargs = handleVersionOption(argc, argv, "mri_annotation2label");
  if (nargs && argc - nargs == 1) {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) {
    usage_exit();
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if (outdir != NULL) {
    // Create the output directory
    err = mkdir(outdir,0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",outdir);
      perror(NULL);
      return(1);
    }
  }

  /*--- Get environment variables ------*/
  if (SUBJECTS_DIR == NULL) { // otherwise specified on cmdline
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  /* ------ Load subject's surface ------ */
  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfacename);
  printf("Reading surface \n %s\n",tmpstr);
  Surf = MRISread(tmpstr);
  if (Surf == NULL) {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }

  /* ------ Load annotation ------ */
  if(fio_FileExistsReadable(annotation)) {
    strcpy(annotfile,annotation);
  } else  sprintf(annotfile,"%s/%s/label/%s.%s.annot",
                    SUBJECTS_DIR,subject,hemi,annotation);
  printf("Loading annotations from %s\n",annotfile);
  err = MRISreadAnnotation(Surf, annotfile);
  if (err) {
    printf("INFO: could not load from %s, trying ",annotfile);
    sprintf(annotfile,"%s/%s/label/%s_%s.annot",
            SUBJECTS_DIR,subject,hemi,annotation);
    printf("%s\n",annotfile);
    err = MRISreadAnnotation(Surf, annotfile);
    if (err) {
      printf("ERROR: MRISreadAnnotation() failed\n");
      exit(1);
    }
    printf("OK, that worked.\n");
  }
  if (Surf->ct == NULL) {
    printf("ERROR: cannot find embedded color table\n");
    exit(1);
  }

  if(segbase == -1000) {
    // segbase has not been set with --segbase
    if(!strcmp(annotation,"aparc")) {
      if(!strcmp(hemi,"lh")) {
        segbase = 1000;
      } else {
        segbase = 2000;
      }
    } else if(!strcmp(annotation,"aparc.a2005s")) {
      if(!strcmp(hemi,"lh")) {
        segbase = 1100;
      } else {
        segbase = 2100;
      }
    } else if(!strcmp(annotation,"aparc.a2009s")) {
      if(!strcmp(hemi,"lh")) {
        segbase = 11100;
      } else {
        segbase = 12100;
      }
    } else {
      segbase = 0;
    }
  }
  printf("Seg base %d\n",segbase);

  if(LobesFile) {
    MRISaparc2lobes(Surf, (int) Ge_lobarDivision);
    MRISwriteAnnotation(Surf,LobesFile);
    if(ctabfile != NULL) {
      Surf->ct->idbase = segbase;
      CTABwriteFileASCII(Surf->ct,ctabfile);
    }
    exit(0);
  }
  if(ctabfile != NULL) {
    Surf->ct->idbase = segbase;
    CTABwriteFileASCII(Surf->ct,ctabfile);
  }


  if(borderfile || BorderAnnotFile) {
    printf("Computing annot border\n");
    border = MRISannot2border(Surf);
    if(BorderAnnotFile) {
      for (vtxno = 0; vtxno < Surf->nvertices; vtxno++) {
        IsBorder = MRIgetVoxVal(border,vtxno,0,0,0);
        if(IsBorder) {
          continue;
        }
        Surf->vertices[vtxno].annotation = 0;
      }
      err = MRISwriteAnnotation(Surf,BorderAnnotFile);
      if(err) {
        printf("ERROR: cannot write %s\n",BorderAnnotFile);
        exit(1);
      }
    }
    if(borderfile) {
      MRIwrite(border,borderfile);
      if(err) {
        printf("ERROR: cannot write %s\n",BorderAnnotFile);
        exit(1);
      }
    }
    exit(0);
  }


  if(segfile != NULL) {
    printf("Converting to a segmentation\n");
    seg = MRISannot2seg(Surf,segbase);
    MRIwrite(seg,segfile);
    exit(0);
  }

  /* Determine which indices are present in the file by
     examining the index at each vertex */
  animax = -1;
  for (vtxno = 0; vtxno < Surf->nvertices; vtxno++) {
    // get the annotation (rgb packed into an int)
    ano = Surf->vertices[vtxno].annotation;
    // From annotation, get index into color table
    CTABfindAnnotation(Surf->ct, ano, &vtxani);
    nperannot[vtxani] ++;
    if (animax < vtxani) {
      animax = vtxani;
    }
  }
  printf("max index = %d\n",animax);

  // Load statistic file
  if(StatFile) {
    printf("Loading stat file %s\n",StatFile);
    Stat = MRIread(StatFile);
    if(Stat == NULL) {
      exit(1);
    }
    if(Stat->width != Surf->nvertices) {
      printf("ERROR: dimension mismatch between surface (%d) and stat (%d)\n",
             Surf->nvertices,Stat->width);
      exit(1);
    }
  }

  // Go thru each index present and save a label
  for (ani=0; ani <= animax; ani++) {

    if (label_index >= 0 && ani != label_index) {
      continue ;
    }
    if (nperannot[ani] == 0) {
      printf("%3d  %5d  empty --- skipping \n",ani,nperannot[ani]);
      continue;
    }

    if (labelbase != NULL) {
      sprintf(labelfile,"%s-%03d.label",labelbase,ani);
    }
    if (outdir != NULL)  {
      int req = snprintf(labelfile,STRLEN,"%s/%s.%s.label",outdir,hemi,
			 Surf->ct->entries[ani]->name);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }

    printf("%3d  %5d %s\n",ani,nperannot[ani],labelfile);
    label = annotation2label(ani, Surf);
    if(label == NULL) {
      ErrorPrintf(ERROR_BADPARM, "index %d not found, cannot write %s - skipping",
                  ani, labelfile) ;
      continue ;
    }
    strcpy(label->subject_name,subject);

    if(Stat) {
      for(k=0; k < label->n_points; k++) {
        vtxno = label->lv[k].vno;
        label->lv[k].stat = MRIgetVoxVal(Stat,vtxno,0,0,0);
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

  if (argc < 1) {
    usage_exit();
  }

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcmp(option, "--help")) {
      print_help() ;
    }

    else if (!strcmp(option, "--version")) {
      print_version() ;
    }

    else if (!strcmp(option, "--debug")) {
      debug = 1;
    } else if (!strcmp(option, "--a2005s")) {
      annotation = "aparc.a2005s";
    }

    /* -------- source inputs ------ */
    else if(!strcmp(option, "--subject") || !strcmp(option, "--s")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sd")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SUBJECTS_DIR = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--surface") || !strcmp(option, "--surf")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      surfacename = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--stat")) {
      if(nargc < 1) {
        argnerr(option,1);
      }
      StatFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--labelbase")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      labelbase = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--label")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      label_index = atoi(pargv[0]);
      printf("only extracting label %s (%d)\n",
             cma_label_to_name(label_index), label_index) ;
      nargsused = 1;
    } else if (!strcmp(option, "--outdir")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      outdir = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annotation")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      annotation = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--table")) {
      printf("ERROR: this program no long accepts --table as\n");
      printf("an argument. The colortable is now read directly\n");
      printf("from the annotation.\n");
      exit(1);
    } else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--seg")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      segfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--segbase")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&segbase);
      nargsused = 1;
    } else if (!strcmp(option, "--ctab")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      ctabfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--border")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      borderfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--border-annot")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      BorderAnnotFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--lobes")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      LobesFile = pargv[0];
      nargsused = 1;
      Ge_lobarDivision = e_default;
    } else if (!strcmp(option, "--lobesStrict")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      LobesFile = pargv[0];
      nargsused = 1;
      Ge_lobarDivision = e_strict;
    } else if (!strcmp(option, "--lobesStrictPHCG")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      LobesFile = pargv[0];
      nargsused = 1;
      Ge_lobarDivision = e_strict_phcg;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option)) {
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      }
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
  printf("   --lobes <LobesFile>\n");
  printf("   	Create an annotation based on cortical lobes.\n");
  printf("   	Note that the precentral and postcentral labels are not\n");
  printf("   	included as part of the 'frontal' and 'parietal' lobes.\n");
  printf("   	The lobar annotation is saved to <LobesFile>.\n");
  printf("   --lobesStrict <LobesFile>\n");
  printf("   	Use a slightly stricter lobe definition that adds the\n");
  printf("   	precentral to the 'frontal' and includes the postcentral\n");
  printf("   	with the 'parietal' lobe.\n");
  printf("   	The lobar annotation is saved to <LobesFile>.\n");
  printf("   --lobesStrictPHCG <LobesFile>\n");
  printf("   	Use stricter lobe definition, and adds an additional lobe called\n");
  printf("   	parahippocampalgyrus, which includes parahippocampal, entorhinal,\n");
  printf("   	temporalpole and fusiform.\n");
  printf("   	The lobar annotation is saved to <LobesFile>.\n");
  printf("\n");
  printf("Output options:\n");
  printf("   --label <int> : extract only single label \n");
  printf("   --labelbase  output will be base-XXX.label \n");
  printf("   --outdir dir :  output will be dir/hemi.name.label \n");
  printf("   --seg segfile : output will be a segmentation 'volume'\n");
  printf("   --segbase base : add base to the annotation number to get seg value\n");
  printf("   --ctab colortable : colortable like FreeSurferColorLUT.txt\n");
  printf("   --border borderfile : output will be a binary overlay of the parc borders \n");
  printf("   --border-annot borderannot : default goes in subject/label\n");
  printf("\n");
  printf("   --annotation as found in SUBJDIR/labels <aparc>\n");
  printf("   --sd <SUBJECTS_DIR>  specify SUBJECTS_DIR on cmdline\n") ;
  printf("   --surface    name of surface <white>. Only affect xyz in label.\n");
  printf("   --stat statfile : surface overlay file (curv or volume format).\n");
  printf("\n");
  printf("   --table : obsolete. Now gets from annotation file\n");
  printf("\n");
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
  if(labelbase) {
    fprintf(fp,"labelbase = %s\n",labelbase);
  }
  if(outdir) {
    fprintf(fp,"outdir = %s\n",outdir);
  }
  if(segfile) {
    fprintf(fp,"segfile = %s\n",segfile);
  }
  fprintf(fp,"surface   = %s\n",surfacename);
  fprintf(fp,"\n");
  return;
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("This program will convert an annotation into multiple label files\n");
  printf("or into a segmentaion 'volume'. It can also create a border overlay.\n\n");

  printf(
    "User specifies the subject, hemisphere, label base, and (optionally)\n"
    "the annotation base and surface. By default, the annotation base is\n"
    "aparc. The program will retrieves the annotations from\n"
    "SUBJECTS_DIR/subject/label/hemi_annotbase.annot. A separate label file\n"
    "is created for each annotation index. The output file names can take\n"
    "one of two forms: (1) If --outdir dir is used, then the output will\n"
    "be dir/hemi.name.lable, where name is the corresponding name in \n"
    "the embedded color table. (2) If --labelbase is used, name of the file conforms to\n"
    "labelbase-XXX.label, where XXX is the zero-padded 3 digit annotation\n"
    "index. If labelbase includes a directory path, that directory must\n"
    "already exist. If there are no points in the annotation file for a\n"
    "particular index, no label file is created. The xyz coordinates in the\n"
    "label file are obtained from the values at that vertex of the\n"
    "specified surface. The default surface is 'white'. Other options\n"
    "include 'pial' and 'orig'.\n"
    "\n"
    "The human-readable names that correspond to the annotation indices for \n"
    "aparc are embedded in the annotation file itself. It is no longer\n"
    "necessary (or possible) to specify the table explicitly with the\n"
    "--table option.\n"
    " \n"
    "--seg segfile \n"
    "--segbase segbase \n"
    " \n"
    "Convert annotation into a volume-encoded surface segmentation file. This \n"
    "is a volume format where the value at each 'voxel' is an integer index.\n"
    "The value of the index depends on several factors. By default, it should \n"
    "match the index for the label as found in $FREESURFER_HOME/FreeSurferColorLUT.txt;\n"
    "this requires that the annotation be either aparc or aparc.a2005s. If \n"
    "aparc and hemi=lh, then segbase=1000, etc. This makes the index match\n"
    "that found in aparc+aseg.mgz. If the annotation is neither \n"
    "aparc nor aparc.a2005s, then segbase=0. This behavior can be overridden \n"
    "by manually specifying a segbase with --segbase.\n"
    "\n"
    " \n"
    "--border borderfile \n"
    "\n"
    "Creates an overlay file in which the boundaries of the parcellations are \n"
    "set to 1 and everything else is 0. This can then be loaded as an overlay\n"
    "in tksurfer.\n"
    " \n"
    " --stat StatFile\n"
    " \n"
    " Put the value from the StatFile into the Stat field in the label. StatFile\n"
    " must be a curv or surface overlay (eg, mgh) format (eg, lh.thickness).\n"
    " \n"
    " \n"
    "Convert annotation into a volume-encoded surface segmentation file. This \n"
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
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1) {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  } else {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if (subject == NULL) {
    fprintf(stderr,"ERROR: no source subject specified\n");
    exit(1);
  }
  if (hemi == NULL) {
    fprintf(stderr,"ERROR: No hemisphere specified\n");
    exit(1);
  }
  if(outdir == NULL && labelbase == NULL && segfile == NULL &&
      borderfile == NULL && LobesFile == NULL && ctabfile == NULL) {
    fprintf(stderr,"ERROR: no output specified\n");
    exit(1);
  }
  if(outdir != NULL && labelbase != NULL) {
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
  if (len < 2) {
    return(0);
  }

  if (flag[0] == '-' && flag[1] != '-') {
    return(1);
  }
  return(0);
}

