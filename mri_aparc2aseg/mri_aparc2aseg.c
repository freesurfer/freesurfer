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

static char vcid[] = "$Id: mri_aparc2aseg.c,v 1.2 2005/09/12 22:35:32 greve Exp $";
char *Progname = NULL;
char *SUBJECTS_DIR = NULL;
char *subject = NULL;
char *OutASegFile = NULL;
char *OutAParcFile = NULL;
char *OutDistFile = NULL;
int debug = 0;
MRI *ASeg, *mritmp;
MRI *AParc;
MRI *Dist;
MRIS *lhwhite, *rhwhite;
MRIS *lhpial, *rhpial;
MHT *lhwhite_hash, *rhwhite_hash;
MHT *lhpial_hash, *rhpial_hash;
VERTEX vtx;
int  lhwvtx, lhpvtx, rhwvtx, rhpvtx;
MATRIX *Vox2RAS, *CRS, *RAS;
float dlhw, dlhp, drhw, drhp;
float dminctx = 6.0;


char tmpstr[2000];
char annotfile[1000];

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, err, asegid, c, r, s, nctx, annot;
  int annotid, IsCortex=0, hemi=0, segval=0;
  float dmin=0.0;

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

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

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
  sprintf(annotfile,"%s/%s/label/lh.aparc.annot",SUBJECTS_DIR,subject);
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
  sprintf(annotfile,"%s/%s/label/rh.aparc.annot",SUBJECTS_DIR,subject);
  printf("\nLoading rh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(rhwhite, annotfile);
  if(err){
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  if(lhwhite->ct) printf("Have color table for annotation\n");
  //print_annotation_table(stdout);

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

	// If it's not labeled as cortex or wm in the aseg, skip
	asegid = MRIgetVoxVal(ASeg,c,r,s,0);
	if(asegid != 3 && asegid != 42 && asegid != 2 && asegid != 41) continue;
	if(asegid == 3 || asegid == 42) IsCortex = 1;
	else                            IsCortex = 0;

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

	if( IsCortex && hemi == 1) segval = annotid+1000;
	if( IsCortex && hemi == 2) segval = annotid+2000;
	if(!IsCortex && hemi == 1) segval = annotid+3000;
	if(!IsCortex && hemi == 2) segval = annotid+4000;

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
    else if (!strcmp(option, "--s")){
      if(nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaseg")){
      if(nargc < 1) argnerr(option,1);
      OutASegFile = pargv[0];
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
  printf("   --oaseg  file : output aseg+aparc volume file\n");
  printf("   --oaparc file : output aparc-only volume file\n");
  printf("\n");
  printf("\n");
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
  printf("Help! Currently does nothing.\n");
  printf(" \n"
	 "mri_aparc2aseg --s fbirn-anat-101 \n"
	 "  --oaseg aseg_aparc.mgz --oaparc aparc.mgz\n");
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
    printf("ERROR: must specify an output aseg file\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
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
