/*
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


// mri_aparc2aseg.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
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

#undef X

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
char *SUBJECTS_DIR = NULL;
char *subject = NULL;
char *WMSegFile = NULL;
MRI *ASeg, *mritmp;
MRI *WMSeg;
MRIS *lhwhite, *rhwhite;
MHT *lhwhite_hash, *rhwhite_hash;
MHT *lhpial_hash, *rhpial_hash;
static struct { float x,y,z; } vtx;
int  lhwvtx, rhwvtx;
MATRIX *Vox2RAS, *CRS, *RAS;
float dlhw, drhw;
char annotfile[1000];


int debug = 0;
char tmpstr[2000];


/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, err, asegid, c, r, s, annot, hemioffset;
  int annotid;
  struct utsname uts;
  char *cmdline, cwd[2000];

  nargs = handleVersionOption(argc, argv, "mri_aparc2wmseg");
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
  dump_options(stdout);

  printf("\n");
  printf("%s\n",getVersion().c_str());
  printf("cwd %s\n",cwd);
  printf("cmdline %s\n",cmdline);
  printf("sysname  %s\n",uts.sysname);
  printf("hostname %s\n",uts.nodename);
  printf("machine  %s\n",uts.machine);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  printf("SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  printf("subject %s\n",subject);
  printf("\n");
  fflush(stdout);

  /* ------ Load subject's lh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  printf("Reading lh white surface \n %s\n",tmpstr);
  lhwhite = MRISread(tmpstr);
  if (lhwhite == NULL) {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }

  printf("Building hash of lh white\n");
  lhwhite_hash = MHTcreateVertexTable_Resolution(lhwhite, CURRENT_VERTICES,16);

  /* ------ Load lh annotation ------ */
  sprintf(annotfile,"%s/%s/label/lh.aparc.annot",SUBJECTS_DIR,subject);
  printf("Loading lh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(lhwhite, annotfile);
  if (err) {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  /* ------ Load subject's rh surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
  printf("Reading rh white surface \n %s\n",tmpstr);
  rhwhite = MRISread(tmpstr);
  if (rhwhite == NULL) {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  if (debug) printf("Building hash of rh white\n");
  rhwhite_hash = MHTcreateVertexTable_Resolution(rhwhite, CURRENT_VERTICES,16);

  /* ------ Load rh annotation ------ */
  sprintf(annotfile,"%s/%s/label/rh.aparc.annot",SUBJECTS_DIR,subject);
  printf("Loading rh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(rhwhite, annotfile);
  if (err) {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  if (debug && lhwhite->ct) printf("Have color table for annotation\n");
  if (debug) print_annotation_table(stdout);

  /* ------ Load ASeg ------ */
  sprintf(tmpstr,"%s/%s/mri/aseg.mgz",SUBJECTS_DIR,subject);
  if (!fio_FileExistsReadable(tmpstr)) {
    sprintf(tmpstr,"%s/%s/mri/aseg.mgh",SUBJECTS_DIR,subject);
    if (!fio_FileExistsReadable(tmpstr)) {
      sprintf(tmpstr,"%s/%s/mri/aseg/COR-.info",SUBJECTS_DIR,subject);
      if (!fio_FileExistsReadable(tmpstr)) {
        printf("ERROR: cannot find aseg\n");
        exit(1);
      } else
        sprintf(tmpstr,"%s/%s/mri/aseg/",SUBJECTS_DIR,subject);
    }
  }

  printf("Loading aseg from %s\n",tmpstr);
  ASeg = MRIread(tmpstr);
  if (ASeg == NULL) {
    printf("ERROR: loading aseg %s\n",tmpstr);
    exit(1);
  }
  mritmp = MRIchangeType(ASeg,MRI_INT,0,0,0);
  MRIfree(&ASeg);
  ASeg = mritmp;
  WMSeg = MRIclone(ASeg,NULL);

  Vox2RAS = MRIxfmCRS2XYZtkreg(ASeg);
  if (debug) {
    printf("ASeg Vox2RAS: -----------\n");
    MatrixPrint(stdout,Vox2RAS);
    printf("-------------------------\n");
  }
  CRS = MatrixAlloc(4,1,MATRIX_REAL);
  CRS->rptr[4][1] = 1;
  RAS = MatrixAlloc(4,1,MATRIX_REAL);
  RAS->rptr[4][1] = 1;

  // Go through each voxel in the aseg
  printf("\n");
  printf("Labeling WM\n");
  for (c=0; c < ASeg->width; c++) {
    if (debug) printf("%3d ",c);
    if (debug && c%20 ==19) printf("\n");
    for (r=0; r < ASeg->height; r++) {
      for (s=0; s < ASeg->depth; s++) {

        // If it's not labeled as white matter in the aseg, set
        // seg value to that from the aseg and skip the rest
        asegid = MRIgetVoxVal(ASeg,c,r,s,0);
        if (asegid != 2 && asegid != 41) {
          MRIsetVoxVal(WMSeg,c,r,s,0,asegid);
          continue;
        }

        // Convert the CRS to RAS
        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(Vox2RAS,CRS,RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        // Get the index of the closest vertex in the
        // lh.white, rh.white
        lhwvtx = MHTfindClosestVertexNoXYZ(lhwhite_hash,lhwhite,vtx.x,vtx.y,vtx.z,&dlhw);
        rhwvtx = MHTfindClosestVertexNoXYZ(rhwhite_hash,rhwhite,vtx.x,vtx.y,vtx.z,&drhw);

        if ( (lhwvtx < 0) && (rhwvtx < 0) ) {
          printf("ERROR: could not map to any surface.\n");
          printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
                 c,r,s,vtx.x,vtx.y,vtx.z);
          exit(1);
        }

        if (lhwvtx < 0) dlhw = 1000000000000000.0;
        if (rhwvtx < 0) drhw = 1000000000000000.0;

        if (dlhw < drhw) {
          // Left hemi is closer than the right
          annot = lhwhite->vertices[lhwvtx].annotation;
          hemioffset = 1000;
          if (lhwhite->ct)
            CTABfindAnnotation(lhwhite->ct, annot, &annotid);
          else
            annotid = annotation_to_index(annot);
        } else {
          // Right hemi is closer than the left
          annot = rhwhite->vertices[rhwvtx].annotation;
          hemioffset = 2000;
          if (rhwhite->ct)
            CTABfindAnnotation(lhwhite->ct, annot, &annotid);
          else
            annotid = annotation_to_index(annot);
        }

        MRIsetVoxVal(WMSeg,c,r,s,0,annotid+hemioffset);
      }
    }
  }

  printf("\nWriting output wmseg to %s\n",WMSegFile);
  MRIwrite(WMSeg,WMSegFile);

  printf("mri_aparc2wmseg done\n");
  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

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
    else if (!strcmp(option, "--s")) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--wmseg")) {
      if (nargc < 1) argnerr(option,1);
      WMSegFile = pargv[0];
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
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --s subject \n");
  printf("   --wmseg wmsegfile\n");
  printf("\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
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
  if (WMSegFile == NULL) {
    printf("ERROR: must specify a wm seg file\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
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


