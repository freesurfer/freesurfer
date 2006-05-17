// mri_aparc2aseg.c 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
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

static char vcid[] = "$Id: mri_aparc2aseg.c,v 1.7 2006/05/17 22:11:35 greve Exp $";
char *Progname = NULL;
char *SUBJECTS_DIR = NULL;
char *subject = NULL;
char *OutASegFile = NULL;
char *OutAParcFile = NULL;
char *OutDistFile = NULL;
int debug = 0;
int UseRibbon = 0;
MRI *ASeg, *mritmp;
MRI *AParc;
MRI *Dist;
MRI *lhRibbon,*rhRibbon;
MRIS *lhwhite, *rhwhite;
MRIS *lhpial, *rhpial;
MHT *lhwhite_hash, *rhwhite_hash;
MHT *lhpial_hash, *rhpial_hash;
VERTEX vtx;
int  lhwvtx, lhpvtx, rhwvtx, rhpvtx;
MATRIX *Vox2RAS, *CRS, *RAS;
float dlhw, dlhp, drhw, drhp;
float dminctx = 5.0;
int LabelWM=0;
int LabelHypoAsWM=0;

char tmpstr[2000];
char annotfile[1000];
char *annotname = "aparc";
int baseoffset = 0;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, err, asegid, c, r, s, nctx, annot;
  int annotid, IsCortex=0, IsWM=0, IsHypo=0, hemi=0, segval=0;
  float dmin=0.0, lhRibbonVal=0, rhRibbonVal=0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  /* ------ Load subject's lh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  lhwhite = MRISread(tmpstr);
  if(lhwhite == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }

  printf("\n");
  printf("Building hash of lh white\n");
  lhwhite_hash = MHTfillVertexTableRes(lhwhite, NULL,CURRENT_VERTICES,16);

  /* ------ Load subject's lh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
  printf("\nReading lh pial surface \n %s\n",tmpstr);
  lhpial = MRISread(tmpstr);
  if(lhpial == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  printf("Building hash of lh pial\n");
  lhpial_hash = MHTfillVertexTableRes(lhpial, NULL,CURRENT_VERTICES,16);

  if(lhwhite->nvertices != lhpial->nvertices){
    printf("ERROR: lh white and pial have a different number of vertices (%d,%d)\n",
	   lhwhite->nvertices,lhpial->nvertices);
    exit(1);
  }

  /* ------ Load lh annotation ------ */
  sprintf(annotfile,"%s/%s/label/lh.%s.annot",SUBJECTS_DIR,subject,annotname);
  printf("\nLoading lh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(lhwhite, annotfile);
  if(err){
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  /* ------ Load subject's rh surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
  printf("\nReading rh white surface \n %s\n",tmpstr);
  rhwhite = MRISread(tmpstr);
  if(rhwhite == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  printf("Building hash of rh white\n");
  rhwhite_hash = MHTfillVertexTableRes(rhwhite, NULL,CURRENT_VERTICES,16);

  /* ------ Load subject's rh surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
  printf("\nReading rh pial surface \n %s\n",tmpstr);
  rhpial = MRISread(tmpstr);
  if(rhpial == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  printf("Building hash of rh pial\n");
  rhpial_hash = MHTfillVertexTableRes(rhpial, NULL,CURRENT_VERTICES,16);

  if(rhwhite->nvertices != rhpial->nvertices){
    printf("ERROR: rh white and pial have a different number of vertices (%d,%d)\n",
	   rhwhite->nvertices,rhpial->nvertices);
    exit(1);
  }

  /* ------ Load rh annotation ------ */
  sprintf(annotfile,"%s/%s/label/rh.%s.annot",SUBJECTS_DIR,subject,annotname);
  printf("\nLoading rh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(rhwhite, annotfile);
  if(err){
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  if(lhwhite->ct) printf("Have color table for lh white annotation\n");
  if(rhwhite->ct) printf("Have color table for rh white annotation\n");
  //print_annotation_table(stdout);

  if(UseRibbon){
    sprintf(tmpstr,"%s/%s/mri/lh.ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading lh ribbon mask from %s\n",tmpstr);
    lhRibbon = MRIread(tmpstr);
    if(lhRibbon == NULL){
      printf("ERROR: loading aseg %s\n",tmpstr);
      exit(1);
    }
    sprintf(tmpstr,"%s/%s/mri/rh.ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading rh ribbon mask from %s\n",tmpstr);
    rhRibbon = MRIread(tmpstr);
    if(rhRibbon == NULL){
      printf("ERROR: loading aseg %s\n",tmpstr);
      exit(1);
    }
  }

  /* ------ Load ASeg ------ */
  sprintf(tmpstr,"%s/%s/mri/aseg.mgz",SUBJECTS_DIR,subject);
  if(!fio_FileExistsReadable(tmpstr)){
    sprintf(tmpstr,"%s/%s/mri/aseg.mgh",SUBJECTS_DIR,subject);
    if(!fio_FileExistsReadable(tmpstr)){
      sprintf(tmpstr,"%s/%s/mri/aseg/COR-.info",SUBJECTS_DIR,subject);
      if(!fio_FileExistsReadable(tmpstr)){
	printf("ERROR: cannot find aseg\n");
	exit(1);
      }
      else
	sprintf(tmpstr,"%s/%s/mri/aseg/",SUBJECTS_DIR,subject);
    }
  }

  printf("\nLoading aseg from %s\n",tmpstr);
  ASeg = MRIread(tmpstr);
  if(ASeg == NULL){
    printf("ERROR: loading aseg %s\n",tmpstr);
    exit(1);
  }
  mritmp = MRIchangeType(ASeg,MRI_INT,0,0,0);
  MRIfree(&ASeg);
  ASeg = mritmp;

  AParc = MRIclone(ASeg,NULL);
  if(OutDistFile != NULL){
    Dist = MRIclone(ASeg,NULL);
    mritmp = MRIchangeType(Dist,MRI_FLOAT,0,0,0);
    if(mritmp == NULL){
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&Dist);
    Dist = mritmp;
  }

  Vox2RAS = MRIxfmCRS2XYZtkreg(ASeg);
  printf("ASeg Vox2RAS: -----------\n");  
  MatrixPrint(stdout,Vox2RAS);
  printf("-------------------------\n");  
  CRS = MatrixAlloc(4,1,MATRIX_REAL);
  CRS->rptr[4][1] = 1;
  RAS = MatrixAlloc(4,1,MATRIX_REAL);
  RAS->rptr[4][1] = 1;

  printf("\nLabeling Slice\n");
  nctx = 0;
  annot = 0;
  annotid = 0;

  // Go through each voxel in the aseg
  for(c=0; c < ASeg->width; c++){
    printf("%3d ",c);
    if(c%20 ==19) printf("\n");
    for(r=0; r < ASeg->height; r++){
      for(s=0; s < ASeg->depth; s++){

	asegid = MRIgetVoxVal(ASeg,c,r,s,0);
	if(asegid == 3 || asegid == 42)  IsCortex = 1;
	else                             IsCortex = 0;
	if(asegid >= 77 && asegid <= 82) IsHypo = 1;
	else                             IsHypo = 0;
	if(asegid == 2 || asegid == 41)  IsWM = 1;
	else                             IsWM = 0;
	if(LabelHypoAsWM && IsHypo )     IsWM = 1;

	// If it's not labeled as cortex or wm in the aseg, skip
	if(!IsCortex && !IsWM) continue;

	// If it's wm but not labeling wm, skip
	if(IsWM && !LabelWM) continue;

	// Check whether this point is in the ribbon
	if(UseRibbon){
	  lhRibbonVal = MRIgetVoxVal(lhRibbon,c,r,s,0);
	  rhRibbonVal = MRIgetVoxVal(rhRibbon,c,r,s,0);
	  if(IsCortex){
	    // ASeg says it's in cortex
	    if(lhRibbonVal < 0.5 && rhRibbonVal < 0.5){
	      // but it is not part of the ribbon,
	      // so set it to unknown (0) and go to the next voxel.
	      MRIsetVoxVal(ASeg,c,r,s,0,0);
	      continue;
	    }
	  }
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
	// lh.white, lh.pial, rh.white, rh.pial
	lhwvtx = MHTfindClosestVertexNo(lhwhite_hash,lhwhite,&vtx,&dlhw);
	lhpvtx = MHTfindClosestVertexNo(lhpial_hash, lhpial, &vtx,&dlhp);
	rhwvtx = MHTfindClosestVertexNo(rhwhite_hash,rhwhite,&vtx,&drhw);
	rhpvtx = MHTfindClosestVertexNo(rhpial_hash, rhpial, &vtx,&drhp);

	if(lhwvtx < 0 && lhpvtx < 0 && rhwvtx < 0 && rhpvtx < 0){
	  printf("ERROR: could not map to any surface.\n");
	  printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
		 c,r,s,vtx.x,vtx.y,vtx.z);
	  exit(1);
	}

	if(lhwvtx < 0) dlhw = 1000000000000000.0;
	if(lhpvtx < 0) dlhp = 1000000000000000.0;
	if(rhwvtx < 0) drhw = 1000000000000000.0;
	if(rhpvtx < 0) drhp = 1000000000000000.0;

	if(dlhw < dlhp && dlhw < drhw && dlhw < drhp){
	  annot = lhwhite->vertices[lhwvtx].annotation;
	  hemi = 1;
	  if(lhwhite->ct)
	    annotid = CTABannotationToIndex(lhwhite->ct, annot);
	  else
	    annotid = annotation_to_index(annot);
	  dmin = dlhw;
	}
	if(dlhp < dlhw && dlhp < drhw && dlhp < drhp){
	  annot = lhwhite->vertices[lhpvtx].annotation;
	  hemi = 1;
	  if(lhwhite->ct)
	    annotid = CTABannotationToIndex(lhwhite->ct, annot);
	  else
	    annotid = annotation_to_index(annot);
	  dmin = dlhp;
	}

	if(drhw < dlhp && drhw < dlhw && drhw < drhp){
	  annot = rhwhite->vertices[rhwvtx].annotation;
	  hemi = 2;
	  if(rhwhite->ct)
	    annotid = CTABannotationToIndex(rhwhite->ct, annot);
	  else
	    annotid = annotation_to_index(annot);
	  dmin = drhw;
	}
	if(drhp < dlhp && drhp < drhw && drhp < dlhw){
	  annot = rhwhite->vertices[rhpvtx].annotation;
	  hemi = 2;
	  if(rhwhite->ct)
	    annotid = CTABannotationToIndex(rhwhite->ct, annot);
	  else
	    annotid = annotation_to_index(annot);
	  dmin = drhp;
	}

	if( IsCortex && hemi == 1) segval = annotid+1000 + baseoffset;
	if( IsCortex && hemi == 2) segval = annotid+2000 + baseoffset;
	if(!IsCortex && hemi == 1) segval = annotid+3000 + baseoffset;
	if(!IsCortex && hemi == 2) segval = annotid+4000 + baseoffset;

	if(!IsCortex && dmin > dminctx && hemi == 1) segval = 5001;
	if(!IsCortex && dmin > dminctx && hemi == 2) segval = 5002;

	MRIsetVoxVal(ASeg,c,r,s,0,segval);
	MRIsetVoxVal(AParc,c,r,s,0,annot);
	if(OutDistFile != NULL) MRIsetVoxVal(Dist,c,r,s,0,dmin);

	if(debug || annotid == -1){
	  printf("\n");
	  printf("aseg id = %d\n",asegid);
	  printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
		 c,r,s,vtx.x,vtx.y,vtx.z);
	  if(lhwvtx > 0) printf("lhw  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
		 lhwvtx, dlhw,
		 lhwhite->vertices[lhwvtx].x,
		 lhwhite->vertices[lhwvtx].y,
		 lhwhite->vertices[lhwvtx].z);
	  if(lhpvtx > 0) printf("lhp  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
		 lhpvtx, dlhp,
		 lhpial->vertices[lhpvtx].x,
		 lhpial->vertices[lhpvtx].y,
		 lhpial->vertices[lhpvtx].z);
	  if(rhwvtx > 0) printf("rhw  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
		 rhwvtx, drhw,
		 rhwhite->vertices[rhwvtx].x,
		 rhwhite->vertices[rhwvtx].y,
		 rhwhite->vertices[rhwvtx].z);
	  if(rhpvtx > 0) printf("rhp  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
		 rhpvtx, drhp,
		 rhpial->vertices[rhpvtx].x,
		 rhpial->vertices[rhpvtx].y,
		 rhpial->vertices[rhpvtx].z);
	  printf("annot = %d, annotid = %d\n",annot,annotid);
	  CTABprint(stdout,lhwhite->ct);
	  exit(1);
	}

	nctx++;
      }
    }
  }
  printf("nctx = %d\n",nctx);

  printf("Writing output aseg to %s\n",OutASegFile);
  MRIwrite(ASeg,OutASegFile);

  if(OutAParcFile != NULL){
    printf("Writing output aparc to %s\n",OutAParcFile);
    MRIwrite(AParc,OutAParcFile);
  }
  if(OutDistFile != NULL){
    printf("Writing output dist file to %s\n",OutDistFile);
    MRIwrite(Dist,OutDistFile);
  }


  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

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

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--ribbon"))  UseRibbon = 1;
    else if (!strcasecmp(option, "--noribbon"))  UseRibbon = 0;
    else if (!strcasecmp(option, "--labelwm"))  LabelWM = 1;
    else if (!strcasecmp(option, "--hypo-as-wm"))  LabelHypoAsWM = 1;
    else if (!strcmp(option, "--s")){
      if(nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaseg") || !strcmp(option, "--o")){
      if(nargc < 1) argnerr(option,1);
      OutASegFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--a2005s")){
      annotname = "aparc.a2005s";
      baseoffset = 100;
    }
    else if (!strcmp(option, "--annot")){
      if(nargc < 1) argnerr(option,1);
      annotname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaparc")){
      if(nargc < 1) argnerr(option,1);
      OutAParcFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--dist")){
      if(nargc < 1) argnerr(option,1);
      OutDistFile = pargv[0];
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
  printf("   --s subject \n");
  printf("   --o volfile : output aparc+aseg volume file\n");
  //printf("   --oaparc file : output aparc-only volume file\n");
  printf("   --ribbon : use mri/hemi.ribbon.mgz as a mask for ctx.\n");
  printf("\n");
  printf("   --a2005 : use aparc.a2005 instead of aparc\n");
  printf("   --annot annotname : use annotname instead of aparc\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
"Maps the cortical labels from the automatic cortical parcellation (aparc)\n"
"to the automatic segmentation volume (aseg). The result can be used as\n"
"the aseg would. The algorithm is to find each aseg voxel labeled as\n"
"cortex (3 and 42) and assign it the label of the closest cortical vertex.\n"
"If the voxel is not in the ribbon (as defined by mri/lh.ribbon and \n"
"rh.ribbon), then the voxel is marked as unknown (0). This can be turned\n"
"off with --noribbon. The cortical parcellation is obtained from \n"
"subject/label/hemi.aparc.annot which should be based on the \n"
"curvature.buckner40.filled.desikan_killiany.gcs atlas \n"
"The aseg is obtained from subject/mri/aseg.mgz and should be based on\n"
"the RB40_talairach_2005-07-20.gca atlas. If these atlases are used, then\n"
"the segmentations can be viewed with tkmedit and the FreeSurferColorLUT.txt\n"
"color table found in $FREESURFER_HOME. These are the default atlases\n"
"used by recon-all.\n"
"\n"
"--s subject\n"
"\n"
"Name of the subject as found in the SUBJECTS_DIR.\n"
"\n"
"--o volfile\n"
"\n"
"Full path of file to save the output segmentation in. Default\n"
"is mri/aparc+aseg.mgz\n"
"\n"
"--ribbon\n"
"\n"
"Mask cortical voxels with mri/hemi.ribbon.mgz. \n"
"\n"
"--a2005s\n"
"\n"
"Use ?h.aparc.a2005s.annot. Output will be aparc.a2005s+aseg.mgz.   \n"
"Creates index numbers that match a2005s entries in FreeSurferColorsLUT.txt\n"
"\n"
"--annot annotname\n"
"\n"
"Use annotname surface annotation. By default, uses ?h.aparc.annot. \n"
"With this option, it will load ?h.annotname.annot. \n"
"The output file will be set to annotname+aseg.mgz, but this can be \n"
"changed with --o. Note: running --annot aparc.a2005s is NOT the\n"
"same as running --a2005s. The index numbers will be different.\n"
"\n"
"EXAMPLE:\n"
"\n"
"  mri_aparc2aseg --s bert\n"
"  tkmedit bert orig.mgz -segmentation mri/aparc+aseg.mgz \n"
"       $FREESURFER_HOME/FreeSurferColorLUT.txt\n"
"\n"
"NOTES AND BUGS:\n"
"\n"
"The volumes of the cortical labels will be different than reported by\n"
"mris_anatomical_stats because partial volume information is lost when\n"
"mapping the surface to the volume. The values reported by  \n"
"mris_anatomical_stats will be more accurate than the volumes from\n"
"the aparc+aseg volume.\n"
"\n"
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
    printf("ERROR: must specify a subject\n");
    exit(1);
  }
  if(OutASegFile == NULL){
    sprintf(tmpstr,"%s/%s/mri/%s+aseg.mgz",SUBJECTS_DIR,subject,annotname);
    OutASegFile = strcpyalloc(tmpstr);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"subject %s\n",subject);
  fprintf(fp,"outvol %s\n",OutASegFile);
  fprintf(fp,"useribbon %d\n",UseRibbon);
  fprintf(fp,"baseoffset %d\n",baseoffset);
  if(LabelWM){
    printf("labeling wm\n");
    if(LabelHypoAsWM) printf("labeling hypo-intensities as wm\n");
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
