/**
 * @file  mris_diff.c
 * @brief Compare two surfaces.
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/12/07 20:58:33 $
 *    $Revision: 1.19 $
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



// things to do:
// --test-surf-vtx vtxno val field
// --test-surf-face faceno val field
// --test-aparc vtxno val
// --test-curv vtxno val
// --log
// --thresh
// --debug

/*
  BEGINHELP

  Determines whether two surfaces or surface-based files differ.
  See below for what 'differ' means.

  The basic usage is something like:

  mris_diff surf1 surf2

  mris_diff --s1 subj1 --s2 subj2 --hemi lh --surf white
  same as:
  mris_diff SD/subj1/surf/hemi.surf SD/subj2/surf/hemi.surf

  Reads in surf1 and surf2, checks:
  1. nvertices
  2. nfaces
  3. vtx->x
  4. vtx->y
  5. vtx->z
  6. vtx->nx
  7. vtx->ny
  8. vtx->nz
  9. number of neighbors
  10. neighbor identity (6?)
  11. ripflag
  12. face->nx
  13. face->ny
  14. face->nz
  15. face->area
  16. face->ripflag
  17. face vertex identities (3)


  mris_diff --s1 subj1 --s2 subj2 --hemi lh --curv curv
  SD/subj1/surf/hemi.curv SD/subj2/surf/hemi.curv

  mris_diff --s1 subj1 --s2 subj2 --hemi lh --annot aparc
  SD/subj1/label/hemi.aparc SD/subj2/label/hemi.aparc.annot

  ENDHELP
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
#include "mrisurf.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "annotation.h"
#include "cmdargs.h"
#include "timer.h"
#include "matfile.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_diff.c,v 1.19 2011/12/07 20:58:33 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
static int debug=0;
static int checkoptsonly=0;
static struct utsname uts;

static char *subject1=NULL, *subject2=NULL, *hemi=NULL;
static char *SUBJECTS_DIR=NULL, *SUBJECTS_DIR1=NULL, *SUBJECTS_DIR2=NULL;
static char *curvname=NULL, *aparcname=NULL,*aparc2name=NULL, *surfname=NULL;
static char *surf1path=NULL, *surf2path=NULL;
static char *out_fname ;
static char tmpstr[2000];
static char *xyzRMSFile=NULL;
static char *angleRMSFile=NULL;

static MRIS *surf1, *surf2;

static int CheckSurf=0;
static int CheckXYZ=1;
static int CheckNXYZ=1;
static int ComputeNormalDist=0;
static int CheckCurv=0;
static int CheckAParc=0;
static double thresh=0;
static double seed=1234;

static int error_count=0;
static int MAX_NUM_ERRORS=10; // in loops, stop after this many errors found
// set by cmd-line parm --maxerrs

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, nthvtx, nnbrs1, nnbrs2, nthnbr, nbrvtxno1, nbrvtxno2;
  int nthface, annot1, annot2;
  VERTEX *vtx1, *vtx2;
  FACE *face1, *face2;
  float diff, maxdiff, rms;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
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

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("INFO: SUBJECTS_DIR not defined in environment\n");
    //exit(1);
  }
  if (SUBJECTS_DIR1 == NULL) SUBJECTS_DIR1 = SUBJECTS_DIR;
  if (SUBJECTS_DIR2 == NULL) SUBJECTS_DIR2 = SUBJECTS_DIR;

  if (surf1path == NULL && surfname == NULL) surfname = "orig";

  if (surf1path == NULL) {
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR1,subject1,hemi,surfname);
    surf1path = strcpyalloc(tmpstr);
  }
  if (surf2path == NULL) {
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR2,subject2,hemi,surfname);
    surf2path = strcpyalloc(tmpstr);
  }
  dump_options(stdout);

  //read-in each surface.  notice that the random number generator is
  //seeded with the same value prior to each read.  this is because in
  //the routine MRIScomputeNormals, if it finds a zero-length vertex
  //normal, is adds a random value to the x,y,z and recomputes the normal.
  //so if comparing identical surfaces, the seed must be the same so that
  //any zero-length vertex normals appear the same.
  setRandomSeed(seed) ;
  surf1 = MRISread(surf1path);
  if (surf1 == NULL) {
    printf("ERROR: could not read %s\n",surf1path);
    exit(1);
  }
  setRandomSeed(seed) ;
  surf2 = MRISread(surf2path);
  if (surf2 == NULL) {
    printf("ERROR: could not read %s\n",surf2path);
    exit(1);
  }
  printf("Number of vertices %d %d\n",surf1->nvertices,surf2->nvertices);
  printf("Number of faces    %d %d\n",surf1->nfaces,surf2->nfaces);

  //Number of Vertices ----------------------------------------
  if (surf1->nvertices != surf2->nvertices) {
    printf("Surfaces differ in number of vertices %d %d\n",
           surf1->nvertices,surf2->nvertices);
    exit(101);
  }
  //Number of Faces ------------------------------------------
  if (surf1->nfaces != surf2->nfaces) {
    printf("Surfaces differ in number of faces %d %d\n",
           surf1->nfaces,surf2->nfaces);
    exit(101);
  }

  //surf1->faces[10000].area = 100;
  //surf1->vertices[10000].x = 100;

  if (ComputeNormalDist) {
    double dist, dx, dy, dz, dot ;
    MRI    *mri_dist ;

    mri_dist = MRIalloc(surf1->nvertices,1,1,MRI_FLOAT) ;
    MRIScomputeMetricProperties(surf1) ;
    MRIScomputeMetricProperties(surf2) ;
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      vtx1 = &(surf1->vertices[nthvtx]);
      vtx2 = &(surf2->vertices[nthvtx]);
      dx = vtx2->x-vtx1->x ;
      dy = vtx2->y-vtx1->y ;
      dz = vtx2->z-vtx1->z ;
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      dot = dx*vtx1->nx + dy*vtx1->ny + dz*vtx1->nz ;
      dist = dist * dot / fabs(dot) ;
      MRIsetVoxVal(mri_dist, nthvtx, 0, 0, 0, dist) ;
    }
    MRIwrite(mri_dist, out_fname) ;
    MRIfree(&mri_dist) ;
    exit(0);
  }

  if(xyzRMSFile){
    printf("Computing xyz RMS\n");
    MRI *xyzRMS;
    xyzRMS = MRIalloc(surf1->nvertices,1,1,MRI_FLOAT) ;    
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      vtx1 = &(surf1->vertices[nthvtx]);
      vtx2 = &(surf2->vertices[nthvtx]);
      rms = sqrt(pow(vtx1->x-vtx2->x,2) + pow(vtx1->y-vtx2->y,2) + pow(vtx1->z-vtx2->z,2));
      MRIsetVoxVal(xyzRMS,nthvtx,0,0,0,rms);
    }
    MRIwrite(xyzRMS,xyzRMSFile);
    exit(0);
  }

  if(angleRMSFile){
    printf("Computing angle RMS\n");
    MRI *angleRMS;
    double dot, radius1, radius2;
    angleRMS = MRIalloc(surf1->nvertices,1,1,MRI_FLOAT) ;    
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      vtx1 = &(surf1->vertices[nthvtx]);
      vtx2 = &(surf2->vertices[nthvtx]);
      radius1 = sqrt(vtx1->x*vtx1->x + vtx1->y*vtx1->y + vtx1->z*vtx1->z);
      radius2 = sqrt(vtx2->x*vtx2->x + vtx2->y*vtx2->y + vtx2->z*vtx2->z);
      dot = (vtx1->x*vtx2->x + vtx1->y*vtx2->y + vtx1->z*vtx2->z)/(radius1*radius2);
      //printf("%6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f  %6.2f %6.2f  %5.4f\n",
      // vtx1->x,vtx1->y,vtx1->z, vtx2->x,vtx2->y,vtx2->z, radius1, radius2, acos(dot)*180/M_PI);
      MRIsetVoxVal(angleRMS,nthvtx,0,0,0,acos(dot)*180/M_PI);
    }
    MRIwrite(angleRMS,angleRMSFile);
    exit(0);
  }

  maxdiff=0;
  //------------------------------------------------------------
  if (CheckSurf) {
    printf("Comparing surfaces\n");

    // Loop over vertices ---------------------------------------
    error_count=0;
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      vtx1 = &(surf1->vertices[nthvtx]);
      vtx2 = &(surf2->vertices[nthvtx]);
      if (vtx1->ripflag != vtx2->ripflag) {
        printf("Vertex %d differs in ripflag %c %c\n",
               nthvtx,vtx1->ripflag,vtx2->ripflag);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
      if (CheckXYZ) {
        diff=fabs(vtx1->x - vtx2->x);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Vertex %d differs in x %g %g\t(%g)\n",
                 nthvtx,vtx1->x,vtx2->x,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
        diff=fabs(vtx1->y - vtx2->y);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Vertex %d differs in y %g %g\t(%g)\n",
                 nthvtx,vtx1->y,vtx2->y,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
        diff=fabs(vtx1->z - vtx2->z);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Vertex %d differs in z %g %g\t(%g)\n",
                 nthvtx,vtx1->z,vtx2->z,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
      }
      if (CheckNXYZ) {
        diff=fabs(vtx1->nx - vtx2->nx);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Vertex %d differs in nx %g %g\t(%g)\n",
                 nthvtx,vtx1->nx,vtx2->nx,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
        diff=fabs(vtx1->ny - vtx2->ny);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Vertex %d differs in ny %g %g\t(%g)\n",
                 nthvtx,vtx1->ny,vtx2->ny,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
        diff=fabs(vtx1->nz - vtx2->nz);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Vertex %d differs in nz %g %g\t(%g)\n",
                 nthvtx,vtx1->nz,vtx2->nz,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
      }
      nnbrs1 = surf1->vertices[nthvtx].vnum;
      nnbrs2 = surf2->vertices[nthvtx].vnum;
      if (nnbrs1 != nnbrs2) {
        printf("Vertex %d has a different number of neighbors %d %d\n",
               nthvtx,nnbrs1,nnbrs2);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
      for (nthnbr=0; nthnbr < nnbrs1; nthnbr++) {
        nbrvtxno1 = surf1->vertices[nthvtx].v[nthnbr];
        nbrvtxno2 = surf2->vertices[nthvtx].v[nthnbr];
        if (nbrvtxno1 != nbrvtxno2) {
          printf("Vertex %d differs in the identity of the "
                 "%dth neighbor %d %d\n",nthvtx,nthnbr,nbrvtxno1,nbrvtxno2);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
      }
      if (error_count>=MAX_NUM_ERRORS) break;
    }// loop over vertices
    if (maxdiff>0) printf("maxdiff=%g\n",maxdiff);
    if (error_count > 0) {
      printf("Exiting after finding %d errors\n",error_count);
      if (error_count>=MAX_NUM_ERRORS) {
        printf("Exceeded MAX_NUM_ERRORS loop guard\n");
      }
      exit(103);
    }

    // Loop over faces ----------------------------------------
    error_count=0;
    for (nthface=0; nthface < surf1->nfaces; nthface++) {
      face1 = &(surf1->faces[nthface]);
      face2 = &(surf2->faces[nthface]);
      if (CheckNXYZ) {
        diff=fabs(face1->nx - face2->nx);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Face %d differs in nx %g %g\t(%g)\n",
                 nthface,face1->nx,face2->nx,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
        diff=fabs(face1->ny - face2->ny);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Face %d differs in ny %g %g\t(%g)\n",
                 nthface,face1->ny,face2->ny,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
        diff=fabs(face1->nz - face2->nz);
        if (diff>maxdiff) maxdiff=diff;
        if (diff>thresh) {
          printf("Face %d differs in nz %g %g\t(%g)\n",
                 nthface,face1->nz,face2->nz,diff);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
      }
      diff=fabs(face1->area - face2->area);
      if (diff>maxdiff) maxdiff=diff;
      if (diff>thresh) {
        printf("Face %d differs in area %g %g\t(%g)\n",
               nthface,face1->area,face2->area,diff);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
      if (face1->ripflag != face2->ripflag) {
        printf("Face %d differs in ripflag %c %c\n",
               nthface,face1->ripflag,face2->ripflag);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
      for (nthvtx = 0; nthvtx < 3; nthvtx++) {
        if (face1->v[nthvtx] != face2->v[nthvtx]) {
          printf("Face %d differs in identity of %dth vertex %d %d\n",
                 nthface,nthvtx,face1->ripflag,face2->ripflag);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
      } // end loop over nthface vertex
      if (error_count>=MAX_NUM_ERRORS) break;
    } // end loop over faces
    if (maxdiff>0) printf("maxdiff=%g\n",maxdiff);
    if (error_count > 0) {
      printf("Exiting after finding %d errors\n",error_count);
      if (error_count>=MAX_NUM_ERRORS) {
        printf("Exceeded MAX_NUM_ERRORS loop guard\n");
      }
      exit(103);
    }

    printf("Surfaces are the same\n");
    exit(0);
  } // end check surf

  // -----------------------------------------------------------------
  if (CheckCurv) {
    printf("Checking curv file %s\n",curvname);
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR1,subject1,hemi,curvname);
    printf("Loading curv file %s\n",tmpstr);
    if (MRISreadCurvatureFile(surf1, tmpstr) != 0) {
      printf("ERROR: reading curvature file %s\n",tmpstr);
      exit(1);
    }
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR2,subject2,hemi,curvname);
    printf("Loading curv file %s\n",tmpstr);
    if (MRISreadCurvatureFile(surf2, tmpstr) != 0) {
      printf("ERROR: reading curvature file %s\n",tmpstr);
      exit(1);
    }
    error_count=0;
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      vtx1 = &(surf1->vertices[nthvtx]);
      vtx2 = &(surf2->vertices[nthvtx]);
      diff=fabs(vtx1->curv - vtx2->curv);
      if (diff>maxdiff) maxdiff=diff;
      if (diff > thresh) {
        printf("curv files differ at vertex %d %g %g\t(%g)\n",
               nthvtx,vtx1->curv,vtx2->curv,diff);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
    } // end loop over vertices
    if (maxdiff>0) printf("maxdiff=%g\n",maxdiff);
    if (error_count > 0) {
      printf("Exiting after finding %d errors\n",error_count);
      if (error_count>=MAX_NUM_ERRORS) {
        printf("Exceeded MAX_NUM_ERRORS loop guard\n");
      }
      exit(103);
    }
    printf("Curv files are the same\n");
    exit(0);
  } // end check curv

  // ---------------------------------------------------------
  if (CheckAParc) {
    printf("Checking AParc %s\n",aparcname);
    sprintf(tmpstr,"%s/%s/label/%s.%s.annot",
            SUBJECTS_DIR1,subject1,hemi,aparcname);
    printf("Loading aparc file %s\n",tmpstr);
    fflush(stdout);
    if (MRISreadAnnotation(surf1, tmpstr)) {
      printf("ERROR: MRISreadAnnotation() failed %s\n",tmpstr);
      exit(1);
    }
    if (aparc2name) aparcname = aparc2name;
    sprintf(tmpstr,"%s/%s/label/%s.%s.annot",
            SUBJECTS_DIR2,subject2,hemi,aparcname);
    printf("Loading aparc file %s\n",tmpstr);
    fflush(stdout);
    if (MRISreadAnnotation(surf2, tmpstr)) {
      printf("ERROR: MRISreadAnnotation() failed %s\n",tmpstr);
      exit(1);
    }
    error_count=0;
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      annot1 = surf1->vertices[nthvtx].annotation;
      annot2 = surf2->vertices[nthvtx].annotation;
      if (annot1 != annot2) {
        printf("aparc files differ at vertex %d: 1:%s 2:%s\n",
               nthvtx,
               CTABgetAnnotationName(surf1->ct,annot1),
               CTABgetAnnotationName(surf2->ct,annot2));
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
    } // end loop over vertices
    if (error_count > 0) {
      printf("Exiting after finding %d errors\n",error_count);
      if (error_count>=MAX_NUM_ERRORS) {
        printf("Exceeded MAX_NUM_ERRORS loop guard\n");
      }
      exit(103);
    }
    printf("\n"
           "AParc files are the same\n"
           "------------------------\n");
    exit(0);
  }

  return 0;
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc = argc;
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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--no-check-xyz")) CheckXYZ = 0;
    else if (!strcasecmp(option, "--no-check-nxyz")) CheckNXYZ = 0;
    else if (!strcasecmp(option, "--ndist")) {
      ComputeNormalDist = 1;
      out_fname = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--xyz-rms")) {
      if (nargc < 1) CMDargNErr(option,1);
      xyzRMSFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--angle-rms")) {
      if (nargc < 1) CMDargNErr(option,1);
      angleRMSFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--s1")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject1 = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--s2")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject2 = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--sd1")) {
      if (nargc < 1) CMDargNErr(option,1);
      SUBJECTS_DIR1 = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--sd2")) {
      if (nargc < 1) CMDargNErr(option,1);
      SUBJECTS_DIR2 = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--hemi")) {
      if (nargc < 1) CMDargNErr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      CheckSurf=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--curv")) {
      if (nargc < 1) CMDargNErr(option,1);
      curvname = pargv[0];
      CheckCurv=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--aparc")) {
      if (nargc < 1) CMDargNErr(option,1);
      aparcname = pargv[0];
      CheckAParc=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--aparc2")) {
      if (nargc < 1) CMDargNErr(option,1);
      aparc2name = pargv[0];
      CheckAParc=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&thresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--maxerrs")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&MAX_NUM_ERRORS);
      nargsused = 1;
    } else {
      if (surf1path == NULL) {
        surf1path = option;
        CheckSurf=1;
      } else if (surf2path == NULL) surf2path = option;
      else {
        fprintf(stderr,"ERROR: Option %s unknown\n",option);
        if (CMDsingleDash(option))
          fprintf(stderr,"       Did you really mean -%s ?\n",option);
        exit(-1);
      }
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
  printf("USAGE: %s surf1 surf2\n",Progname) ;
  printf("OR: %s required:\n",Progname);
  printf("   --s1 subj1 \n");
  printf("   --s2 subj2 \n");
  printf("   --sd1 subj1_directory (default is SUBJECTS_DIR)\n");
  printf("   --sd2 subj2_directory (default is SUBJECTS_DIR)\n");
  printf("   --hemi hemi (rh or lh)\n");
  printf("   and one of:\n");
  printf("   --surf surf\n");
  printf("   --curv curv\n");
  printf("   --aparc aparc\n");
  printf("   --aparc2 aparc2   optional different name to compare to aparc\n");
  printf("\n");
  printf("other options:\n");
  printf("   --thresh N    threshold (default=0)\n");
  printf("   --maxerrs N   stop looping after N errors (default=%d)\n",
         MAX_NUM_ERRORS);
  printf("\n");
  printf("   --no-check-xyz  : do not check vertex xyz\n");
  printf("   --no-check-nxyz : do not check vertex normals\n");
  printf("   --xyz-rms xyzrmsfile : compute and save rms diff between xyz\n");
  printf("   --angle-rms anglermsfile : compute angle on sphere between xyz\n");
  printf("\n");
  printf("   --debug       turn on debugging\n");
  printf("   --checkopts   don't run anything, just check options and exit\n");
  printf("   --help        print out information on how to use program\n");
  printf("   --version     print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
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
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (surf1path == NULL && subject1==NULL) {
    printf("ERROR: surface 1 not specified\n");
    exit(1);
  }
  if (surf2path == NULL && subject2==NULL) {
    printf("ERROR: surface 2 not specified\n");
    exit(1);
  }
  if ( (surf1path == NULL && surf2path != NULL) ||
       (surf2path == NULL && surf1path != NULL) ) {
    printf("ERROR: must specify absolute path to both or neither.\n");
    exit(1);
  }
  if (surf1path != NULL && surfname != NULL) {
    printf("ERROR: cannot specify both absolute path and --surf .\n");
    exit(1);
  }
  if (surfname != NULL && curvname != NULL) {
    printf("ERROR: cannot specify both surf and curv.\n");
    exit(1);
  }
  if (surfname != NULL && aparcname != NULL) {
    printf("ERROR: cannot specify both surf and aparc.\n");
    exit(1);
  }
  if (curvname != NULL && aparcname != NULL) {
    printf("ERROR: cannot specify both curv and aparc.\n");
    exit(1);
  }
  if (surf1path == NULL && (subject1 == NULL || subject2 == NULL) ) {
    printf("ERROR: need subject names for relative path.\n");
    exit(1);
  }
  if (surf1path == NULL && surfname == NULL && curvname == NULL
      && aparcname == NULL ) {
    printf("ERROR: need --surf or --curv or --aparc with relative path.\n");
    exit(1);
  }
  if (surf1path == NULL && hemi == NULL) {
    printf("ERROR: need hemi for relative path.\n");
    exit(1);
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"SUBJECTS_DIR    %s\n",getenv("SUBJECTS_DIR"));
  if (SUBJECTS_DIR1) fprintf(fp,"SUBJECTS_DIR1   %s\n",SUBJECTS_DIR1);
  if (SUBJECTS_DIR2) fprintf(fp,"SUBJECTS_DIR2   %s\n",SUBJECTS_DIR2);
  fprintf(fp,             "cwd       %s\n",cwd);
  fprintf(fp,             "cmdline   %s\n",cmdline);
  fprintf(fp,             "timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,             "sysname   %s\n",uts.sysname);
  fprintf(fp,             "hostname  %s\n",uts.nodename);
  fprintf(fp,             "machine   %s\n",uts.machine);
  fprintf(fp,             "user      %s\n",VERuser());
  fprintf(fp,             "surf1path %s\n",surf1path);
  fprintf(fp,             "surf2path %s\n",surf2path);
  if (subject1) fprintf(fp,"subject1  %s\n",subject1);
  if (subject2) fprintf(fp,"subject2  %s\n",subject2);
  if (hemi)     fprintf(fp,"hemi      %s\n",hemi);
  if (surfname) fprintf(fp,"surfname  %s\n",surfname);
  if (curvname) fprintf(fp,"curvname  %s\n",curvname);
  if (aparcname)fprintf(fp,"aparcname %s\n",aparcname);
  if (aparc2name)fprintf(fp,"aparc2name %s\n",aparc2name);
  fprintf(fp,"\n");
  fprintf(fp,"\n");

  return;
}
