/**
 * @brief Prints out information about a surface file
 *
 */
/*
 * Original Author: Yasunari Tosa
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

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>
#include <sys/utsname.h>

 
#include "fio.h"
#include "mri.h"
#include "utils.h"
#include "gcsa.h"
#include "colortab.h"
#include "diag.h"
#include "transform.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "version.h"
#include "proto.h"
#include "error.h"
#include "gifti.h"
#include "surfgrad.h"
#include "mrisurf_metricProperties.h"

const char *Progname = "mris_info";


static int  parse_commandline(int argc, char **argv);
static void print_usage(void);
static void usage_exit(void);
static void print_help(void);
static void argnerr(char *option, int n);
static void print_version(void);

// copied from mrisurf.c
#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)

using namespace std;
char *surffile=NULL, *outfile=NULL, *curvfile=NULL, *annotfile=NULL;
char *SUBJECTS_DIR=NULL, *subject=NULL, *hemi=NULL, *surfname=NULL;
int debug = 0;
char tmpstr[2000];
struct utsname uts;
int talairach_flag = 0 ;
MATRIX *XFM=NULL;
int rescale = 0;
double scale=0;
int diag_vno=-1;
int vnox = -1;
int DoAreaStats = 0;
int DoEdgeStats = 0;
int EdgeMetricId = 0;
char *edgefile = NULL;
MRI *mask=NULL;
LABEL *label=NULL;
int DoQuality = 0;
int vmatlab = -1;
char *vmatlabfile=NULL;
int  edgenox = -1;
int CountIntersections = 0;
char *patchname=NULL;
const char *mtxfmt = NULL;
int gifti_disp_image = 1;

static bool doTkrRASConvert = false;

int MRISsaveMarkedAsPointSet(char *fname, MRIS *surf);
int MRISedgeVertices2Pointset(MRIS *surf, const MRI *mask, const int metricid, const double thresh, const char *fname);
int MRISprintEdgeInfo(FILE *fp, const MRIS *surf, int edgeno);

/*------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  double InterVertexDistAvg, InterVertexDistStdDev,avgvtxarea,avgfacearea;
  char ext[STRLEN] ;
  vector<string> type;
  FILE *fp;
  type.push_back("MRIS_BINARY_QUADRANGLE_FILE");
  type.push_back("MRIS_ASCII_TRIANGLE_FILE");
  type.push_back("MRIS_GEO_TRIANGLE_FILE");
  type.push_back("MRIS_TRIANGULAR_SURFACE=MRIS_ICO_SURFACE");
  type.push_back("MRIS_ICO_FILE");
  type.push_back("MRIS_VTK_FILE");

  if (argc < 2) usage_exit();
  parse_commandline(argc, argv);
  if (surffile == NULL) {
    printf("ERROR: must specify a surface file\n");
    exit(1);
  }

  // set environment variable FS_GII to 0, overwrite the value
  setenv("FS_GII", "0", 1);
  
  // Check whether it's a gcs file. If so, just print ctab
  if (!stricmp(FileNameExtension(surffile, ext), (char*)"gcs")) {
    GCSA *gcsa = GCSAread(surffile) ;
    if (!gcsa) {
      cerr << "could not open " << surffile << endl;
      return -1;
    }
    printf("GCSA file %s opened\n", surffile) ;
    if (gcsa->ct != NULL) {
      CTABprintASCII(gcsa->ct,stdout) ;
    }
    return(0) ;
  }

  // Check whether it's a .annot file. If so, just print ctab
  if (!stricmp(FileNameExtension(surffile, ext), (char*)"annot")) {
    COLOR_TABLE* ctab =  NULL;
    int return_code = MRISreadCTABFromAnnotationIfPresent(surffile, &ctab);
    if (NO_ERROR != return_code) {
      fprintf(stderr,"ERROR: could not open %s\n", surffile);
      return -1;
    }
    if (ctab != NULL) {
      CTABprintASCII(ctab,stdout) ;
    }
    return(0) ;
  }

  // Check whether it's a Gifti file. If so, just dump gifti struct
  if (!stricmp(FileNameExtension(surffile, ext), (char*)"gii") && gifti_disp_image) {

    gifti_image* image = gifti_read_image (surffile, 1);
    if (NULL == image)
    {
      fprintf (stderr,"gifti_read_image() returned NULL\n");
      return 1;
    }
    int valid = gifti_valid_gifti_image (image, 1);
    if (valid == 0)
    {
      fprintf (stderr,"GIFTI file %s is invalid!\n", surffile);
      gifti_free_image (image);
      return 1;
    }

    gifti_disp_gifti_image(NULL,image,1);
    
    return(0) ;
  }

  // Open as a surface
  MRIS *mris = MRISread(surffile, doTkrRASConvert);
  if (!mris) {
    fprintf(stderr, "ERROR: cannot open either %s or %s.gii\n", surffile, surffile);
    return -1;
  }
  if(patchname){
    // do we need to read surffile if --patch is specified ???
    int err = MRISreadPatch(mris, patchname, doTkrRASConvert);
    if(err) exit(1);
  }
  MRIScomputeMetricProperties(mris) ;
  if(DoEdgeStats){
    MRISedges(mris);
    MRIScorners(mris);
    MRISfaceMetric(mris,0);
    MRISedgeMetric(mris,0);
    MRIScornerMetric(mris,0);
  }

  if(label){
    // Create a mask from the label
    int n;
    mask = MRIalloc(mris->nvertices,1,1,MRI_INT);
    for (n = 0; n < label->n_points; n++){
      MRIsetVoxVal(mask,label->lv[n].vno,0,0,0, 1);
    }
  }

  if(DoQuality){
    MRISprettyPrintSurfQualityStats(stdout, mris);
    exit(0);
  }

  if(CountIntersections){
    mrisMarkIntersections(mris,0);
    int n, nintersections=0;
    for(n=0; n < mris->nvertices; n++){
      if(mris->vertices[n].marked) nintersections++;
    }
    printf("Found %d intersections\n",nintersections);
    exit(0);
  }

  if(DoAreaStats){
    double *stats = MRIStriangleAreaStats(mris, mask, NULL);
    printf("%d %g %g %g %g\n",(int)stats[0],stats[1],stats[2],stats[3],stats[4]);
    free(stats);
    exit(0);
  }
  if(DoEdgeStats){
    if(EdgeMetricId >= 0){
      double *stats = MRISedgeStats(mris, EdgeMetricId, mask, NULL);
      printf("%d %g %g %g %g\n",(int)stats[0],stats[1],stats[2],stats[3],stats[4]);
      free(stats);
    }
    else {
      int k;
      for(k=0; k < 3; k++){
	double *stats = MRISedgeStats(mris, k, mask, NULL);
	printf("%d %d %g %g %g %g\n",k,(int)stats[0],stats[1],stats[2],stats[3],stats[4]);
	free(stats);
      }
    }
    exit(0);
  }
  if(edgefile){
    MRISedgeWrite(edgefile, mris);
    //MRISedgeVertices2Pointset(mris, mask, 2, 145, "tmp.pointset");
    exit(0);
  }

  if (diag_vno >= 0)
  {
    printf("v %d: (%2.1f, %2.1f, %2.1f)\n", 
           diag_vno,
           mris->vertices[diag_vno].x, 
           mris->vertices[diag_vno].y, 
           mris->vertices[diag_vno].z) ;
    return(0);
  }
  if(vnox >= 0){
    MRISprintVertexInfo(stdout, mris, vnox);
    exit(0);
  }
  if(vmatlab > -1){
    fp = fopen(vmatlabfile,"w");
    MatlabPlotVertexNbhd(fp, mris, vmatlab, 3, 'k', .1);
    fclose(fp);
  }
  if(edgenox >= 0){
    MRISprintEdgeInfo(stdout, mris, edgenox);
    exit(0);
  }

  // attempt to load curvature info, which has the side-effect of checking
  // that the number of vertices is the same (exiting with error if not)
  if (curvfile) {
    printf("\n");
    if (MRISreadCurvatureFile(mris,curvfile)) {
      printf("\n");
      exit(1);
    }
    printf((char*)"%s has the same number of vertices as %s\n",
          curvfile,surffile);
  }

  // read number of vertices in annotation file, and check
  // that the number of vertices is the same (exiting with error if not)
  if (annotfile) {
    printf("\n");  
    /* Open the file. */
    FILE  *fp = fopen(annotfile,"r");
    if (fp==NULL)
      ErrorExit(ERROR_NOFILE, 
                 "ERROR: could not read annot file %s",
                 annotfile) ;

    /* First int is the number of elements. */
    int numVertices = freadInt(fp) ;
    fclose(fp);

    if (numVertices != mris->nvertices ) {
      printf("\n");
      printf("ERROR: %s has %d vertices, while %s has %d vertices!\n",
             annotfile, numVertices, surffile, mris->nvertices);
      exit(1);
    }
    printf((char*)"%s has the same number of vertices as %s\n",
          annotfile,surffile);

    // also dump the colortable
    COLOR_TABLE* ctab =  NULL;
    int return_code = MRISreadCTABFromAnnotationIfPresent(annotfile, &ctab);
    if (NO_ERROR != return_code) {
      fprintf(stderr,"ERROR: could not open %s\n", annotfile);
      return -1;
    }
    if (ctab != NULL) {
      CTABprintASCII(ctab,stdout) ;
    }
  }

  if (rescale) {
    if (mris->group_avg_surface_area == 0) {
      printf("ERROR: cannot rescale a non-group surface\n");
      exit(1);
    }
    scale = MRISrescaleMetricProperties(mris);
    //scale = sqrt((double)mris->group_avg_surface_area/mris->total_area);
    //printf("scale = %lf\n",scale);
    //MRISscale(mris,scale);
  }

  if (talairach_flag) {
    if (! subject) {
      printf("ERROR: need --s with --t\n");
      exit(1);
    }
    XFM = DevolveXFM(subject, NULL, NULL);
    if (XFM == NULL) exit(1);
    printf("Applying talairach transform\n");
    MatrixPrint(stdout,XFM);
    MRISmatrixMultiply(mris,XFM);
  }

  cout << "SURFACE INFO ======================================== " << endl;
  cout << "type        : " << type[mris->type].c_str() << endl;
  if (mris->type == MRIS_BINARY_QUADRANGLE_FILE) {
    FILE *fp = fopen(surffile, "rb") ;
    int magic;
    fread3(&magic, fp) ;
    if (magic == QUAD_FILE_MAGIC_NUMBER)
      cout << "              QUAD_FILE_MAGIC_NUMBER" << endl;
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER)
      cout << "              NEW_QUAD_FILE_MAGIC_NUMBER" << endl;
    fclose(fp);
  }

  InterVertexDistAvg    = mris->avg_vertex_dist;
  InterVertexDistStdDev = mris->std_vertex_dist;
  avgvtxarea = mris->avg_vertex_area;
  avgfacearea = mris->total_area/mris->nfaces;

  cout << "num vertices: " << mris->nvertices << endl;
  cout << "num faces   : " << mris->nfaces << endl;
  cout << "num strips  : " << mris->nstrips << endl;
  cout << "surface area: " << mris->total_area << endl;
  cout << "vg.valid: "     << mris->vg.valid << endl;
  printf("AvgFaceArea      %lf\n",avgfacearea);
  printf("AvgVtxArea       %lf\n",avgvtxarea);
  printf("AvgVtxDist       %lf\n",InterVertexDistAvg);
  printf("StdVtxDist       %lf\n",InterVertexDistStdDev);

  if (mris->group_avg_surface_area > 0) {
    cout << "group avg surface area: " << mris->group_avg_surface_area << endl;
  }
  cout << "ctr                     : (" << mris->xctr 
       << ", " << mris->yctr 
       << ", " << mris->zctr 
       << ")" << endl;
  cout << "original vertex space   : " 
       << (mris->orig_xyzspace ? "scannerRAS" : "surfaceRAS") << endl;
  cout << "vertex locs             : " 
       << (mris->useRealRAS ? "scannerRAS" : "surfaceRAS") 
       << ((mris->useRealRAS != mris->orig_xyzspace) ? "(converted)" : "") << endl;
  if (mris->lta) {
    cout << "talairch.xfm: " << endl;
    MatrixPrint(stdout, mris->lta->xforms[0].m_L);
    cout << "surfaceRAS to talaraiched surfaceRAS: " << endl;
    MatrixPrint(stdout, mris->SRASToTalSRAS_);
    cout << "talairached surfaceRAS to surfaceRAS: " << endl;
    MatrixPrint(stdout, mris->TalSRASToSRAS_);
  }
  printf("Volume Geometry (vg)\n");
  mris->vg.vgprint();

  if (mris->vg.valid)
  {
    printf("Volume Geometry vox2ras\n");
    if(mtxfmt==NULL) MatrixPrint(stdout,vg_i_to_r(&mris->vg));
    else             MatrixPrintFmt(stdout, mtxfmt, vg_i_to_r(&mris->vg));
    printf("Volume Geometry vox2ras-tkr\n");
    if(mtxfmt==NULL) MatrixPrint(stdout,TkrVox2RASfromVolGeom(&mris->vg));
    else MatrixPrintFmt(stdout,mtxfmt,TkrVox2RASfromVolGeom(&mris->vg));
  }
  else
  {
    printf("Volume Geometry vox2ras:     Not available\n");
    printf("Volume Geometry vox2ras-tkr: Not available\n");
  }

  {
    int i ;
    for (i = 0 ; i < mris->ncmds ; i++)
      printf("cmd[%d]: %s\n", i, mris->cmdlines[i]) ;
  }
  if (argc > 2) {
    int vno = atoi(argv[2]) ;
    VERTEX *v ;
    v= &mris->vertices[vno] ;
    printf("mris[%d] = (%2.1f, %2.1f, %2.1f)\n", vno, v->x, v->y, v->z) ;
  }

  uname(&uts);

  // Open an output file to capture values
  if (outfile != NULL) {
    fp = fopen(outfile,"w");
    if (fp == NULL) {
      printf("ERROR: cannot open %s\n",outfile);
      exit(1);
    }
  } else fp = stdout;

  fprintf(fp,"mris_info\n");
  fprintf(fp,"creationtime %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"surfacefile %s\n",surffile);
  if (SUBJECTS_DIR) fprintf(fp,"SUBJECTS_DIR  %s\n",SUBJECTS_DIR);
  if (subject)      fprintf(fp,"subject       %s\n",subject);
  if (hemi)         fprintf(fp,"hemi          %s\n",hemi);
  if (surfname)     fprintf(fp,"surfname      %s\n",surfname);
  fprintf(fp,"hemicode    %d\n",mris->hemisphere);
  fprintf(fp,"talairach_flag  %d\n",talairach_flag);
  fprintf(fp,"rescale     %lf\n",scale);
  fprintf(fp,"nvertices   %d\n",mris->nvertices);
  fprintf(fp,"nfaces      %d\n",mris->nfaces);
  fprintf(fp,"total_area  %f\n",mris->total_area);
  if (mris->group_avg_surface_area > 0)
    fprintf(fp,"group_avg_surf_area  %f\n",mris->group_avg_surface_area);
  printf("group_avg_vtxarea_loaded %d\n",mris->group_avg_vtxarea_loaded);
  fprintf(fp,"avgfacearea  %lf\n",avgfacearea);
  fprintf(fp,"avgvtxarea  %lf\n",avgvtxarea);
  fprintf(fp,"avgvtxdist  %lf\n",InterVertexDistAvg);
  fprintf(fp,"stdvtxdist  %lf\n",InterVertexDistStdDev);
  fprintf(fp,"vtx0xyz   %f %f %f\n",mris->vertices[0].x,mris->vertices[0].y,
          mris->vertices[0].z);
  fprintf(fp,"num_intersecting_faces %d\n",mrisMarkIntersections(mris,0));

  if (outfile != NULL) fclose(fp);

  MRISfree(&mris);
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc-1;
  pargv = &argv[1];
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--tal"))   talairach_flag = 1;
    else if (!strcasecmp(option, "--t"))   talairach_flag = 1;
    else if (!strcasecmp(option, "--r"))   rescale = 1;
    else if (!strcasecmp(option, "--nogifti-disp-image")) gifti_disp_image = 0;
    else if ( !strcmp(option, "--o") ) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--v") ) {
      if (nargc < 1) argnerr(option,1);
      diag_vno = atoi(pargv[0]);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--vx") ) {
      if (nargc < 1) argnerr(option,1);
      vnox = atoi(pargv[0]);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--v-matlab") ) {
      if (nargc < 2) argnerr(option,2);
      vmatlab = atoi(pargv[0]);
      vmatlabfile = pargv[1];
      nargsused = 2;
    } 
    else if ( !strcmp(option, "--ex") ) {
      if (nargc < 1) argnerr(option,1);
      edgenox = atoi(pargv[0]);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--patch") ) {
      if (nargc < 1) argnerr(option,1);
      patchname = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--quality") ) {
      DoQuality = 1;
    } 
    else if ( !strcmp(option, "--intersections") ) {
      CountIntersections = 1;
    } 
    else if ( !strcmp(option, "--area-stats") ) {
      DoAreaStats = 1;
    } 
    else if ( !strcmp(option, "--edge-stats") ) {
      if (nargc < 1) argnerr(option,1);
      EdgeMetricId = atoi(pargv[0]);
      DoEdgeStats = 1;
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--edge-file") ) {
      if (nargc < 1) argnerr(option,1);
      edgefile = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--label") ) {
      if(mask){
	printf("ERROR: cannot spec --mask and  --label\n");
	exit(1);
      }
      if (nargc < 1) argnerr(option,1);
      label = LabelRead(NULL, pargv[0]);
      if(label == NULL) exit(1);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--mask") ) {
      // Only compute edge and area stats from vertices inside the mask
      if(label){
	printf("ERROR: cannot spec --mask and  --label\n");
	exit(1);
      }
      if (nargc < 1) argnerr(option,1);
      mask = MRIread(pargv[0]);
      if(mask == NULL) exit(1);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mtx-fmt")) {
      mtxfmt = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--c") ) {
      if (nargc < 1) argnerr(option,1);
      curvfile = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--a") ) {
      if (nargc < 1) argnerr(option,1);
      annotfile = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--cog") ) {
      if(nargc < 2) argnerr(option,2);
      MRIS *surftmp = MRISread(pargv[0]);
      if(surftmp == NULL) exit(1);
      double xCOG, yCOG, zCOG;
      MRIScalculateCenterCOG2(surftmp, &xCOG, &yCOG, &zCOG);
      FILE *fp=stdout;
      if(strcmp(pargv[1],"nofile")!=0) fp = fopen(pargv[2],"w");
      fprintf(fp,"%9.6f %9.6f %9.6f\n",xCOG, yCOG, zCOG);
      if(fp != stdout) fclose(fp);
      exit(0);
    } 
    else if ( !strcmp(option, "--cog-zero") ) {
      if(nargc < 2) argnerr(option,2);
      MRIS *surftmp = MRISread(pargv[0]);
      if(surftmp == NULL) exit(1);
      double xCOG, yCOG, zCOG;
      MRIScalculateCenterCOG2(surftmp, &xCOG, &yCOG, &zCOG);
      fprintf(stdout,"%9.6f %9.6f %9.6f\n",xCOG, yCOG, zCOG);
      for(int vno=0; vno < surftmp->nvertices; vno++){
	VERTEX *v = &(surftmp->vertices[vno]);
	v->x -= xCOG;
	v->y -= yCOG;
	v->z -= zCOG;
      }
      int err = MRISwrite(surftmp,pargv[1]);
      exit(err);
    } 
    else if ( !strcmp(option, "--s") ) {
      if (nargc < 3) argnerr(option,3);
      subject  = pargv[0];
      hemi     = pargv[1];
      surfname = pargv[2];
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      if (SUBJECTS_DIR==NULL) {
        printf("ERROR: SUBJECTS_DIR not defined in environment\n");
        exit(1);
      }
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
      surffile = strcpyalloc(tmpstr);
      nargsused = 3;
    } else {
      if (NULL == surffile) surffile = option;
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
  printf("USAGE: %s [options] <surfacefile>\n",Progname) ;
  printf("\nOptions:\n");
  printf("  --o outfile : save some data to outfile\n");
  printf("  --s subject hemi surfname : instead of surfacefile\n");
  printf("  --t : apply talairach xfm before reporting info\n");
  printf("  --r : rescale group surface so metrics same as avg of individuals\n");
  printf("  --patch patchfile : load patch before reporting\n");
  printf("  --v vnum : print out vertex information for vertex vnum\n") ;
  printf("  --vx vnum : print out extended vertex information for vertex vnum\n") ;
  printf("     NbrInfo: nbrno, nbrvno, nbrdist, nbrarea, faceno, facearea\n");
  printf("  --c curvfile : check if the specified curvature file has the\n");
  printf("                 same number of vertices as the surface, and\n");
  printf("                 exit with error if not. This is a QA check.\n");
  printf("  --a annotfile : check if the specified annotation file has the\n");
  printf("                  same number of vertices as the surface, and\n");
  printf("                  exit with error if not. This is a QA check.\n");
  printf("                  Also, the colortable is dumped.\n");
  printf("  --area-stats : compute stats on triangle area (n, mean, std, min, max)\n");
  printf("  --edge-stats id : compute stats on edge metric (n, mean, std, min, max)\n");
  printf("                    id=0=length, id=1=dot, id=2=angle, id<0= all\n");
  printf("  --ex edgeno : printout extended into about edge\n");
  printf("  --v-matlab vtxno mfile : write matlab file to plot vertex neighborhood\n");
  printf("  --quality: print out surface quality stats\n");
  printf("  --intersections : print the number of vertices that belong to a face that intersects another face\n");
  printf("  --mask mask.mgz : only compute edge and area stats using vertices in mask\n");
  printf("  --label labelfile : only compute edge and area stats using vertices in label\n");
  printf("  --edge-file file : print edge info for all edges into file\n");
  printf("  --mtx-fmt format : set format for matrix printing (eg, %%12.8f)\n");
  printf("\n");
  printf("  --version   : print version and exits\n");
  printf("  --help      : no clue what this does\n");
  printf("  --nogifti-disp-image: no dump of GIFTI struct, read .gii as surface instead\n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
  printf("Prints out information about a surface file.\n");
  printf("\n");
  exit(1);
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
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/*!
  \fn int MRISsaveMarkedAsPointSet(char *fname, MRIS *surf)
  \brief Outputs a file that can be loaded into freeview with -c
 */
int MRISsaveMarkedAsPointSet(char *fname, MRIS *surf)
{
  int n,vtxno;
  VERTEX *v;
  FILE *fp;

  fp = fopen(fname,"w");
  n = 0;
  for(vtxno=0; vtxno < surf->nvertices; vtxno++){
    v = &(surf->vertices[vtxno]);
    if(! v->marked) continue;
    fprintf(fp,"%g %g %g\n",v->x,v->y,v->z);
    n++;
  }
  fprintf(fp,"info\n");
  fprintf(fp,"numpoints %d\n",n);
  fprintf(fp,"useRealRAS 0\n");
  fclose(fp);

  return(0);
}

int MRISedgeVertices2Pointset(MRIS *surf, const MRI *mask, const int metricid, const double thresh, const char *fname)
{
  int edgeno, nthv, npoints;
  MRI_EDGE *e;
  double metric;
  FILE *fp;

  if(surf->edges == NULL) MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,0);

  fp = fopen(fname,"w");
  npoints = 0;
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    e = &(surf->edges[edgeno]);
    int skip = 0;
    for(nthv=0; nthv < 4; nthv++){
      int vno = e->vtxno[nthv];
      VERTEX  * const v = &(surf->vertices[vno]);
      if(v->ripflag) skip = 1;
      if(mask && MRIgetVoxVal(mask,vno,0,0,0) < 0.5) skip = 1;
    }
    if(skip) continue;
    switch(metricid){
    case 0: metric = e->len; break;
    case 1: metric = e->dot; break;
    case 2: metric = e->angle; break;
    default:
      printf("ERROR: MRISmarkEdgeVertices) metricid %d unrecognized\n",metricid);
      return(1);
    }
    if(metric < thresh) continue;
    for(nthv=0; nthv < 2; nthv++){
      int vno = e->vtxno[nthv];
      VERTEX  * const v = &(surf->vertices[vno]);
      fprintf(fp,"%g %g %g\n",v->x,v->y,v->z);
      npoints++;
    }
  }
  fprintf(fp,"info\n");
  fprintf(fp,"numpoints %d\n",npoints);
  fprintf(fp,"useRealRAS 0\n");
  fclose(fp);
  return(0);
}


int MRISprintEdgeInfo(FILE *fp, const MRIS *surf, int edgeno)
{
  MRI_EDGE *e = &(surf->edges[edgeno]);
  int k,j;
  VERTEX *v;
  FACE *f;

  fprintf(fp,"edgeno %d\n",edgeno);

  fprintf(fp,"vertices ");
  for(k=0; k < 4; k++)  fprintf(fp,"%5d ",e->vtxno[k]);
  fprintf(fp,"\n");


  for(k=0; k < 4; k++){
    v = &(surf->vertices[e->vtxno[k]]);
    fprintf(fp,"v%d = [%10.8f %10.8f %10.8f]'; \n",k,v->x,v->y,v->z);
  }

  fprintf(fp,"faces ");
  for(k=0; k < 2; k++)  fprintf(fp,"%5d ",e->faceno[k]);
  fprintf(fp,"\n");

  for(k=0; k < 2; k++){
    f = &(surf->faces[e->faceno[k]]);
    if(f->norm == NULL) continue;
    fprintf(fp,"f%d = [%10.8f %10.8f %10.8f]'; \n",k,f->norm->rptr[1][1],
	    f->norm->rptr[2][1],f->norm->rptr[3][1]);
  }

  for(j=0; j < 2; j++){
    fprintf(fp,"corner%d ",j);
    for(k=0; k < 4; k++)  fprintf(fp,"%d ",(int)(e->corner[k][j]));
    fprintf(fp,"\n");
  }

  fprintf(fp,"dot %g\n",e->dot);
  fprintf(fp,"angle %g\n",e->angle);
  fprintf(fp,"J %g\n",e->J);

  fprintf(fp,"u = [ ");
  for(k=0; k < 3; k++)  fprintf(fp,"%10.8f ",e->u[k]);
  fprintf(fp,"];\n");

  if(e->gradU){
    fprintf(fp,"gradU = [ ");
    DMatrixPrint(fp,e->gradU);
    fprintf(fp,"];\n");
  }

  fflush(fp);
  return(0);
}


