/**
 * @brief Compare two surfaces.
 *
 */
/*
 * Original Author: Doug Greve
 * Modifications: Bevin R Brett
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



// things to do:
// --test-surf-vtx vtxno val field
// --test-surf-face faceno val field
// --test-aparc vtxno val
// --test-curv vtxno val
// --log
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

    mris_diff --worst-bucket bucket_file --okayBucketMax int file1 file2
    
    mris_diff --grid {[xyz]} spacingfloat grid_file ...
    
  ENDHELP
*/


#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <string>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisurf_metricProperties.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "label.h"
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

const char *Progname = NULL;
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
static const char *xyzRMSFile=NULL;
static const char *angleRMSFile=NULL;
static const char *worstBucketFile=NULL;
static int   okayBucketMax=1;
static int   gridx=0, gridy=0, gridz=0;
static float gridspacing=0;
static char* gridFile=NULL;
int UseScannerRAS = 0;

static MRIS *surf1, *surf2;

// this flag doesn't apply to files passed in with command line options
static bool doTkrRASConvert = false;

static int CheckSurf=0;
static int CheckXYZ=1;
static int CheckNXYZ=1;
static int ComputeNormalDist=0;
static int CheckCurv=0;
static int CheckAParc=0;

static int renumberedSpecified=0;

static long seed=1234;

static int error_count=0;
static int MAX_NUM_ERRORS=10; // in loops, stop after this many errors found
// set by cmd-line parm --maxerrs

MRI *MRISminDist(MRIS *srcsurf, MRIS *trgsurf);

// Because slight numerical differences can cause triangles, especially those on folds, 
// to end up at significantly different orientations (seen as big variations in nx,ny,nz)
// it is important to understand how well the whole surface fits, rather than just looking for the
// few bad matches.  So histograms of the various properties are used...
//
static std::vector<char> vnoToWorstBucket;

#define HistogramSize 20
struct HistogramOfFit {
  bool contributesToWorstBucket;
  HistogramOfFit(bool contributesToWorstBucket = false) : contributesToWorstBucket(contributesToWorstBucket) {}
  double maxV;
  double maxDiff;
  unsigned int v[HistogramSize];
};

static void initHistogramOfFit(HistogramOfFit* histogramOfFit) {
  histogramOfFit->maxV = 0.0;
  histogramOfFit->maxDiff = 0.0;
  int i; for (i = 0; i < HistogramSize; i++) histogramOfFit->v[i] = 0;
}

static int populationHistogramOfFit(HistogramOfFit* histogramOfFit) {
  int population = 0;
  int i; for (i = 0; i < HistogramSize; i++) population += histogramOfFit->v[i];
  return population;
}

static int headHistogramOfFit(HistogramOfFit* histogramOfFit) {
  int i; 
  for (i = HistogramSize; i > 0; i--)
    if (histogramOfFit->v[i - 1] > 0) 
      return i;
  return 0;
}

static void insertHistogramOfFit(int vnoOrNegative, HistogramOfFit* histogramOfFit, double diff, double v) {
  if (histogramOfFit->maxV < v) histogramOfFit->maxV = v;
  if (histogramOfFit->maxDiff < diff) histogramOfFit->maxDiff = diff;
  double fit = 0.01; int i = 0;
  while (fit < diff) { fit *= 3; i++; }
  if (i >= HistogramSize) i = HistogramSize-1;
  histogramOfFit->v[i]++;
  if (histogramOfFit->contributesToWorstBucket && vnoOrNegative >= 0) {
    auto & e = vnoToWorstBucket[vnoOrNegative];
    e = std::max(e,char(i));
  }
}

static int printfHistogramOfFit(HistogramOfFit* histogramOfFit, double const* requiredFit) {
  int countOfBad = 0;
  const int pop = populationHistogramOfFit(histogramOfFit);
  double fit = 0.01; 
  int popSoFar     = 0;
  int requiredFitI = 0;
  const int head = headHistogramOfFit(histogramOfFit);
  int i = 0;
  while (i < head) { 
    const char* comment = "";
    popSoFar += histogramOfFit->v[i];
    double fractionSoFar = (double)popSoFar / (double)pop;
    if (fractionSoFar < requiredFit[requiredFitI]) {
      countOfBad++;
      comment = " *** too few";
    }
    printf("    %8.2g %9d %4.2g%c of the required %g %s\n", 
      (i+1<head)?fit:histogramOfFit->maxDiff, histogramOfFit->v[i], 
      fractionSoFar*100.0, '%', requiredFit[requiredFitI], comment);
    if (requiredFit[requiredFitI+1] >= 0.0) requiredFitI++;
    fit *= 3; i++;
  }
  return countOfBad;
}

static HistogramOfFit vertexXyzHistogram;
static HistogramOfFit vertexRelativeXyzHistogram(true);
static HistogramOfFit vertexNxnynzHistogram;
static HistogramOfFit faceNxnynzHistogram;
static HistogramOfFit faceAreaHistogram;
static HistogramOfFit vertexCurvHistogram;

static void compare(int vnoOrNegative, HistogramOfFit* histogramOfFit, double lhs, double rhs) {
  double absLhs = fabs(lhs);
  double absRhs = fabs(rhs);
  double diffAbs = fabs(lhs - rhs);
  insertHistogramOfFit(vnoOrNegative, histogramOfFit, diffAbs, absLhs>absRhs?absLhs:absRhs);
}

static void initHistograms() {
  initHistogramOfFit(&vertexXyzHistogram);
  initHistogramOfFit(&vertexRelativeXyzHistogram);
  initHistogramOfFit(&vertexNxnynzHistogram);
  initHistogramOfFit(&faceNxnynzHistogram);
  initHistogramOfFit(&faceAreaHistogram);
  initHistogramOfFit(&vertexCurvHistogram);
}

static void printOneHistogram(
  HistogramOfFit* histogramOfFit, 
  const char*     name,
  double const*   requiredFit,
  const char**    badHistogram) {
  printf("%s  largest:%g\n", name, histogramOfFit->maxV);
  if (populationHistogramOfFit(histogramOfFit) == 0) { printf(" empty\n"); return; }
  if (printfHistogramOfFit(histogramOfFit, requiredFit) > 0) {
    *badHistogram = name;
  }
  printf("\n");
}

static const char* printHistograms() {
  const char* badHistogram = NULL;
  
  // The following numbers are tunable guesses
  // 	90%   should be within 0.01
  // 	90%   should be within 0.03
  // 	95%   should be within 0.09
  //	99%   should be within 1
  //
  const double vertexRequiredFit[9] = {0.0, 0.0, 0.0, 0.05, 0.1, 0.5, 0.95, 0.99, -1};
  const double relVtxRequiredFit[9] = {0.5, 0.6, 0.90, 0.95, 0.99, -1};
  const double otherRequiredFit [9] = {0.2, 0.6, 0.95, 0.99, -1};
  
  printOneHistogram(&vertexXyzHistogram        , "vertex xyz"     , vertexRequiredFit, &badHistogram);
  printOneHistogram(&vertexRelativeXyzHistogram, "vertex rel xyz" , relVtxRequiredFit, &badHistogram);
  printOneHistogram(&vertexNxnynzHistogram     , "vertex nxnynz"  , otherRequiredFit,  &badHistogram);
  printOneHistogram(&faceNxnynzHistogram       , "face nxnynz"    , otherRequiredFit,  &badHistogram);
  printOneHistogram(&faceAreaHistogram         , "face area"      , otherRequiredFit,  &badHistogram);
  printOneHistogram(&vertexCurvHistogram       , "vertex curv"    , otherRequiredFit,  &badHistogram);
  return badHistogram;
}


static bool compareVertexPositions(MRIS * const lhs, MRIS * const rhs, std::vector<int> & lhsVno2rhsVno) {

  std::fill(lhsVno2rhsVno.begin(), lhsVno2rhsVno.end(), -1);
  
  initHistogramOfFit(&vertexXyzHistogram);

  static const float maxDistortion = 0.001;
  
  static auto closeEnoughF = [](float lhs, float rhs)->bool {
    return fabs(lhs - rhs) < maxDistortion;
  };

  static auto closeEnoughV = [](VERTEX const & lhs, VERTEX const & rhs) {
    return 
       closeEnoughF(lhs.x, rhs.x)
    && closeEnoughF(lhs.y, rhs.y)
    && closeEnoughF(lhs.z, rhs.z);
  };

  size_t  matched = 0, missing = 0;
  
  typedef std::vector<size_t> List;

  List lhsList(lhs->nvertices), rhsList(rhs->nvertices);

  // sort along the x dimension
  auto initList = [&](MRIS * const mris, List & list) {
    for (size_t i = 0; i < list.size(); i++) list[i] = i;
    std::sort(list.begin(), list.end(), [&](size_t li, size_t ri)->bool { return mris->vertices[li].x < mris->vertices[ri].x; });
  };
  initList(lhs, lhsList);
  initList(rhs, rhsList);

  // slide along the lhs, which is in order
  // and slide 2*maxDistortion window along the rhs around the lhs entry
  // 
  size_t rLo = 0, rHi = 0;
  for (size_t li = 0; li < lhsList.size(); li++) {
      auto const lhsVno = lhsList[li];
      auto const & lv = lhs->vertices[lhsVno];
      while (rLo < rhsList.size() && rhs->vertices[rhsList[rLo]].x < lv.x - maxDistortion) rLo++;
      while (rHi < rhsList.size() && rhs->vertices[rhsList[rHi]].x < lv.x + maxDistortion) rHi++;
      if (rLo == rhsList.size()) break;       // no more candidates
      // the candidates are now [rLo..rhi)
      size_t found = 0;
      VERTEX* rv = nullptr;
      for (auto ri = rLo; ri < rHi; ri++) {
        if (closeEnoughV(lv, rhs->vertices[rhsList[ri]])) {
          found++;
          rv = &rhs->vertices[rhsList[ri]];
        }
      }
      if (found == 0) {
        missing++;
      } else if (found != 1) {
        std::cout << "compareVertexPositions more than 1 candidate matched, in fact found: " << found << std::endl;
      } else {
        lhsVno2rhsVno[lhsVno] = rv - rhs->vertices;
        matched++;
        if (rv) {
          compare(lhsVno,&vertexXyzHistogram, lv.x, rv->x);
          compare(lhsVno,&vertexXyzHistogram, lv.y, rv->y);
          compare(lhsVno,&vertexXyzHistogram, lv.z, rv->z);
        }
      }
  }

  std::cout << "compareVertexPositions matched:" << matched << " missing:" << missing << " out of " << lhs->nvertices << std::endl;

  const char*  badHistogram         = NULL;
  const double vertexRequiredFit[9] = {0.99, -1};
  printOneHistogram(&vertexXyzHistogram, "vertex xyz", vertexRequiredFit, &badHistogram);
  
  if (matched < 0.99*lhs->nvertices) badHistogram = "too many not matched";
  
  return !badHistogram;
}

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, nthvtx, nnbrs1, nnbrs2, nthnbr, nbrvtxno1, nbrvtxno2;
  int nthface, annot1, annot2;
  FACE *face1, *face2;
  float maxdiff, rms;

  nargs = handleVersionOption(argc, argv, "mris_diff");
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

  printf("%s\n", cmdline);
  
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("INFO: SUBJECTS_DIR not defined in environment\n");
    //exit(1);
  }
  if (SUBJECTS_DIR1 == NULL) SUBJECTS_DIR1 = SUBJECTS_DIR;
  if (SUBJECTS_DIR2 == NULL) SUBJECTS_DIR2 = SUBJECTS_DIR;

  // set environment variable FS_GII to 0, overwrite the value
  setenv("FS_GII", "0", 1);
  
  if (surf1path == NULL && surfname == NULL) surfname = const_cast<char*>("orig"); // This is.... nasty

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
  surf1 = MRISread(surf1path, doTkrRASConvert);
  if (surf1 == NULL) {
    printf("ERROR: could not read %s\n",surf1path);
    exit(1);
  }
  setRandomSeed(seed) ;
  surf2 = MRISread(surf2path, doTkrRASConvert);
  if (surf2 == NULL) {
    printf("ERROR: could not read %s\n",surf2path);
    exit(1);
  }
  if(UseScannerRAS){
    printf("Converting to scanner RAS\n");
    MRIStkr2Scanner(surf1);
    MRIStkr2Scanner(surf2);
  }

  if (surf1->vg.valid != surf2->vg.valid)
    printf("WARN: Surface validity of the geometry differs.\n");
  if (surf1->useRealRAS != surf2->useRealRAS)
    printf("WARN: Surface coordinates are in different space.\n");

  printf("Number of vertices %d %d\n",surf1->nvertices,surf2->nvertices);
  printf("Number of faces    %d %d\n",surf1->nfaces,surf2->nfaces);

  //Number of Vertices ----------------------------------------
  if (surf1->nvertices > surf2->nvertices) {
    printf("Swapping surf1 and surf2 to make the surf1 be the one with the least vertices\n");
    std::swap(surf1,surf2);
  }

  if (surf1->nfaces > surf2->nfaces) {
    printf("surf1 has more faces than surf2\n");
    exit(100);
  }
  
  std::vector<int> surf1Vno_to_surf2Vno(surf1->nvertices);

  if (surf1->nvertices != surf2->nvertices) {
    printf("Surfaces differ in number of vertices %d %d%s\n",
         surf1->nvertices,surf2->nvertices,
         renumberedSpecified?"":" Consider using --renumbered\n");
  }

  if (!renumberedSpecified) {

    if (surf1->nvertices != surf2->nvertices) {
      exit(101);
    }

    for (int i = 0; i < surf1->nvertices; i++) surf1Vno_to_surf2Vno[i] = i;   

  } else {

    bool closeEnough =
      compareVertexPositions(surf1,surf2,surf1Vno_to_surf2Vno);

    if (!closeEnough) {
      exit(101);
    }
    
    printf("Surfaces have enough vertices in about the same locations to try to match up the faces\n");
  }
  

  //Number of Faces ------------------------------------------

  // for every face1 in surf1, find a vno1 on it, find the corresponding vno2 in surf2, find the corresponding surf2 face 

  // Even if there are the same number of faces, they need to be matched this way

  if (surf1->nfaces != surf2->nfaces) {
    printf("Surfaces differ in number of faces %d %d%s\n",
      surf1->nfaces,surf2->nfaces,
      renumberedSpecified?"":" Consider using --renumbered\n");
  }

  std::vector<int> surf1Fno_to_surf2Fno(surf1->nfaces);

  if (!renumberedSpecified) {

    if (surf1->nfaces != surf2->nfaces) {
      exit(101);
    }

    for (int fno1 = 0; fno1 < surf1->nfaces; fno1++) {
      surf1Fno_to_surf2Fno[fno1] = fno1;
    }

  } else {

    int matchedFaces = 0;

    for (int fno1 = 0; fno1 < surf1->nfaces; fno1++) {
      surf1Fno_to_surf2Fno[fno1] = -1;

      FACE const * f1 = &surf1->faces[fno1];

      int vnos2[3];
      for (int i = 0; i < 3; i++) vnos2[i] = surf1Vno_to_surf2Vno[f1->v[i]]; 
      std::sort(vnos2+0,vnos2+3);

      if (vnos2[0] < 0) continue;

      size_t found = 0;
      VERTEX_TOPOLOGY const * v2 = &surf2->vertices_topology[vnos2[0]];
      for (int fi2 = 0; fi2 < v2->num; fi2++) {
        FACE const * candidateF2 = &surf2->faces[v2->f[fi2]];
        int candidateVnos[3];
        for (int i = 0; i < 3; i++) candidateVnos[i] = candidateF2->v[i]; 
        std::sort(candidateVnos+0,candidateVnos+3);
        if (candidateVnos[0] == vnos2[0]
        &&  candidateVnos[1] == vnos2[1]
        &&  candidateVnos[2] == vnos2[2]) {
          cheapAssert(!found);
          found = 1;
          surf1Fno_to_surf2Fno[fno1] = v2->f[fi2];
        }
      } 

      if (found) matchedFaces++;
    }

    if (matchedFaces != surf1->nfaces) {
      printf("Matched %d of %d faces\n", matchedFaces, surf1->nfaces);
      if (matchedFaces < 0.99 * surf1->nfaces) {
        printf("Not enough faces matched to try to compare the surfaces\n");
        exit(101);
      }
    }
  }
  
  if (ComputeNormalDist) {
    MRI* mri_dist = MRIalloc(surf1->nvertices,1,1,MRI_FLOAT) ;
    
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      int nthvtx2 = surf1Vno_to_surf2Vno[nthvtx];
      if (nthvtx2 < 0) continue;
    
      VERTEX const * const vtx1 = &(surf1->vertices[nthvtx ]);
      VERTEX const * const vtx2 = &(surf2->vertices[nthvtx2]);
      double dx = vtx2->x - vtx1->x ;
      double dy = vtx2->y - vtx1->y ;
      double dz = vtx2->z - vtx1->z ;
      double dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      double dot  = dx*vtx1->nx + dy*vtx1->ny + dz*vtx1->nz ;
      dist = dist * dot / fabs(dot) ;
      MRIsetVoxVal(mri_dist, nthvtx, 0, 0, 0, dist) ;
    }
    MRIwrite(mri_dist, out_fname) ;
    MRIfree(&mri_dist) ;
    exit(0);
  }

  if(xyzRMSFile){
    printf("Computing xyz RMS\n");
    MRI *xyzRMS = MRIalloc(surf1->nvertices,1,1,MRI_FLOAT) ;    
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      int nthvtx2 = surf1Vno_to_surf2Vno[nthvtx];
      if (nthvtx2 < 0) continue;

      VERTEX const * const vtx1 = &(surf1->vertices[nthvtx ]);
      VERTEX const * const vtx2 = &(surf2->vertices[nthvtx2]);
      rms = sqrt(pow(vtx1->x - vtx2->x,2) + pow(vtx1->y - vtx2->y,2) + pow(vtx1->z - vtx2->z,2));
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
      int nthvtx2 = surf1Vno_to_surf2Vno[nthvtx];
      if (nthvtx2 < 0) continue;

      VERTEX const * const vtx1 = &(surf1->vertices[nthvtx ]);
      VERTEX const * const vtx2 = &(surf2->vertices[nthvtx2]);
      
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

  vnoToWorstBucket.resize(surf1->nvertices);
  std::fill(vnoToWorstBucket.begin(),vnoToWorstBucket.end(),0);
  
  maxdiff=0;
  //------------------------------------------------------------
  if (CheckSurf) {
    printf("Comparing surfaces\n");

    initHistograms();
    // Loop over vertices ---------------------------------------
    error_count=0;
    
    int vertices_with_bad_neighbours_count = 0;
    
    for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
      int nthvtx2 = surf1Vno_to_surf2Vno[nthvtx];
      if (nthvtx2 < 0) continue;

      VERTEX_TOPOLOGY const * const vtx1t = &(surf1->vertices_topology[nthvtx ]);
      VERTEX          const * const vtx1  = &(surf1->vertices         [nthvtx ]);
      VERTEX_TOPOLOGY const * const vtx2t = &(surf2->vertices_topology[nthvtx2]);
      VERTEX          const * const vtx2  = &(surf2->vertices         [nthvtx2]);
      
      if (vtx1->ripflag != vtx2->ripflag) {
        printf("Vertex %d differs in ripflag %c %c\n",
               nthvtx,vtx1->ripflag,vtx2->ripflag);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
      
      bool has_bad_neighbours = false;
      
      if (CheckXYZ) {
        compare(nthvtx,&vertexXyzHistogram, vtx1->x, vtx2->x);
        compare(nthvtx,&vertexXyzHistogram, vtx1->y, vtx2->y);
        compare(nthvtx,&vertexXyzHistogram, vtx1->z, vtx2->z);

	// The problem with comparing xyz is that a whole "continent" of
	// faces can drift in the same internal shape and they all
	// appear to be bad, whereas in reality the internals are as good
	// as elsewhere.  So, in addition to this dubious compare, compare
	// the distances between the non-ripped vertices of the faces that come
	// together at a vertex.
	//
	if (vtx1t->num != vtx2t->num) {
          
          has_bad_neighbours = true;
          
	} else {
	  int fn;
	  for (fn = 0; fn < vtx1t->num; fn++) {
	    FACE* f1 = &(surf1->faces[vtx1t->f[fn]]);
	    FACE* f2 = &(surf2->faces[vtx2t->f[fn]]);
	    if (f1->ripflag || f2->ripflag) continue;
	    int vn;
	    for (vn = 0; vn < VERTICES_PER_FACE; vn++) {
	      VERTEX* v1 = &(surf1->vertices[f1->v[vn]]);
	      VERTEX* v2 = &(surf2->vertices[f2->v[vn]]);
	      if (v1->ripflag || v2->ripflag) continue;
      	      double dx1 = v1->x - vtx1->x ;
              double dy1 = v1->y - vtx1->y ;
              double dz1 = v1->z - vtx1->z ;
              double dist1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
      	      double dx2 = v2->x - vtx2->x ;
              double dy2 = v2->y - vtx2->y ;
              double dz2 = v2->z - vtx2->z ;
              double dist2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
	      compare(f1->v[vn],&vertexRelativeXyzHistogram, dist1, dist2);
	    }
	  }
	}
      }
      if (CheckNXYZ) {
        compare(nthvtx,&vertexNxnynzHistogram, vtx1->nx, vtx2->nx);
        compare(nthvtx,&vertexNxnynzHistogram, vtx1->ny, vtx2->ny);
        compare(nthvtx,&vertexNxnynzHistogram, vtx1->nz, vtx2->nz);
      }
      
      nnbrs1 = vtx1t->vnum;
      nnbrs2 = vtx2t->vnum;

      int nthnbr2=0;      
      for (nthnbr=0; nthnbr < nnbrs1; nthnbr++) {
      
        nbrvtxno1 = vtx1t->v[nthnbr];
        if (surf1Vno_to_surf2Vno[vtx1t->v[nthnbr]] < 0) continue;    // this neighbor no longer exists
 
        if (nnbrs2 == nthnbr2) {
          if (0) printf("Surf2->vertices[%d] doesn't have enough neighbors\n", nthvtx);
          has_bad_neighbours = true;
          break;
        }
        
        nbrvtxno2 = vtx2t->v[nthnbr2];

        if (surf1Vno_to_surf2Vno[nbrvtxno1] != nbrvtxno2) {
          if (0) printf("Surf1->Vertices[%d] and Surf2->Vertices[%d] differs in the identity of the "
                 "v[%d:%d], namely Surf2 vertices[%d:%d]\n",
                 nthvtx, nthvtx2,
                 nthnbr,nthnbr2,
                 surf1Vno_to_surf2Vno[nbrvtxno1], nbrvtxno2);
          has_bad_neighbours = true;
        }
        
        nthnbr2++;
      }

      if (nthnbr2 != nnbrs2) {
        if (0) printf("Surf2->vertices[%d] has extra neighbors\n",
               nthvtx);
        has_bad_neighbours = true;
      }
      
      if (has_bad_neighbours) {
        vertices_with_bad_neighbours_count++;
        if (vertices_with_bad_neighbours_count > surf1->nvertices*0.01) {
          error_count = MAX_NUM_ERRORS;
          break;
        }
      }
      
    }// loop over vertices

    if (vertices_with_bad_neighbours_count) {
      
      printf("%d vertices have differing neighbour info\n",
        vertices_with_bad_neighbours_count);
        
      if (surf1->nvertices == surf2->nvertices) error_count++;
    }

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
    int faces_with_no_equiv_count = 0;
    for (nthface=0; nthface < surf1->nfaces; nthface++) {
      auto nthface2 = surf1Fno_to_surf2Fno[nthface];
      
      if (nthface2 < 0) {
        faces_with_no_equiv_count++;
        continue;
      }
      
      face1 = &(surf1->faces[nthface ]); FaceNormCacheEntry const * fNorm1 = getFaceNorm(surf1, nthface );
      face2 = &(surf2->faces[nthface2]); FaceNormCacheEntry const * fNorm2 = getFaceNorm(surf2, nthface2);
      
      if (CheckNXYZ) {
        compare(-1,&faceNxnynzHistogram, fNorm1->nx, fNorm2->nx);
        compare(-1,&faceNxnynzHistogram, fNorm1->ny, fNorm2->ny);
        compare(-1,&faceNxnynzHistogram, fNorm1->nz, fNorm2->nz);
      }
      compare(-1,&faceAreaHistogram, face1->area, face2->area);
      if (face1->ripflag != face2->ripflag) {
        printf("Face %d:%d differs in ripflag %c %c\n",
               nthface,nthface2,face1->ripflag,face2->ripflag);
        if (++error_count>=MAX_NUM_ERRORS) break;
      }
      for (nthvtx = 0; nthvtx < 3; nthvtx++) {
        if (surf1Vno_to_surf2Vno[face1->v[nthvtx]] != face2->v[nthvtx]) {
          printf("Face %d:%d differs in identity of %dth vertex %d %d\n",
                 nthface,nthface2,nthvtx,face1->ripflag,face2->ripflag);
          if (++error_count>=MAX_NUM_ERRORS) break;
        }
      } // end loop over nthface vertex
      if (error_count>=MAX_NUM_ERRORS) break;
    } // end loop over faces
    if (maxdiff>0) printf("maxdiff=%g\n",maxdiff);
    
    if (faces_with_no_equiv_count > 0) {
      printf("%d surf1 faces of %d have no equivalent in surf2\n", faces_with_no_equiv_count, surf1->nfaces);
    }
    
    if (faces_with_no_equiv_count > surf1->nfaces * 0.01) {
      error_count++;
    }
    
    if (error_count > 0) {
      printf("Exiting after finding %d errors\n",error_count);
      if (error_count>=MAX_NUM_ERRORS) {
        printf("Exceeded MAX_NUM_ERRORS loop guard\n");
      }
      exit(103);
    }

    const char* badHistogram = printHistograms();
    if (badHistogram) {
      printf("Too many differences in %s (and maybe others)\n", badHistogram);
      exit(103);
    }

    if(worstBucketFile){
      printf("Writing worstBucket\n");
      for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
        auto & v = surf1->vertices[nthvtx];
        bool interesting = (vnoToWorstBucket[nthvtx] > okayBucketMax);
        
        // yellow interesting, grey uninteresting by default
        if (interesting) {
          v.stat   = 1;                     
          v.marked = 1;
        }
        
      }
      LABEL* area = LabelFromMarkedSurface(surf1);
      LabelWrite(area,worstBucketFile);
    }

    if(gridFile && gridspacing > 0.0){
      printf("Writing gridFile\n");
      size_t interestingCount = 0;
      for (nthvtx=0; nthvtx < surf1->nvertices; nthvtx++) {
        auto & vt = surf1->vertices_topology[nthvtx];
        auto & v  = surf1->vertices         [nthvtx];
        int interesting = 0;
        
        for (int ni = 0; ni < vt.vnum; ni++) {
          auto spans = [&](const char* which, float d0,float d1)->bool { 
            bool result = int(d0/gridspacing) != int(d1/gridspacing); 
            return result;
          };
          auto & v2 = surf1->vertices[vt.v[ni]];
          interesting |= (gridx && spans("x",v.x,v2.x)) ? 1 : 0;
          interesting |= (gridy && spans("y",v.y,v2.y)) ? 2 : 0;
          interesting |= (gridz && spans("z",v.z,v2.z)) ? 4 : 0;
          if (false) {
            static size_t count,limit = 1;
            if (count++ > limit) { if (limit < 100) limit++; else limit *= 2;
              printf("%6ld d0:(%f,%f,%f) d1:(%f,%f,%f) result:%d\n", count, v.x,v.y,v.z,v2.x,v2.y,v2.z,interesting); 
            }
          }
        }
        
        if (interesting) {
          v.stat   = 1;                     
          v.marked = 1;
          interestingCount++; 
        }
      }
      LABEL* area = LabelFromMarkedSurface(surf1);
      printf("gridfile populated with interestingCount:%ld out of %d\n", interestingCount, surf1->nvertices);
      if (area) LabelWrite(area,gridFile);
      else printf("No interesting points, so gridfile not written\n");
    }

    exit(0);
  } // end check surf

  // -----------------------------------------------------------------
  if (CheckCurv) {
    initHistograms();
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
      int nthvtx2 = surf1Vno_to_surf2Vno[nthvtx];
      if (nthvtx2 < 0) continue;
      VERTEX const * const vtx1 = &(surf1->vertices[nthvtx ]);
      VERTEX const * const vtx2 = &(surf2->vertices[nthvtx2]);
      compare(nthvtx,&vertexCurvHistogram, vtx1->curv, vtx2->curv);
    } // end loop over vertices
    if (maxdiff>0) printf("maxdiff=%g\n",maxdiff);
    if (error_count > 0) {
      printf("Exiting after finding %d errors\n",error_count);
      if (error_count>=MAX_NUM_ERRORS) {
        printf("Exceeded MAX_NUM_ERRORS loop guard\n");
      }
      exit(103);
    }
    const char* badHistogram = printHistograms();
    if (badHistogram) {
      printf("Too many differences in %s (and maybe others)\n", badHistogram);
      exit(103);
    }
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
      int nthvtx2 = surf1Vno_to_surf2Vno[nthvtx];
      if (nthvtx2 < 0) continue;
      annot1 = surf1->vertices[nthvtx ].annotation;
      annot2 = surf2->vertices[nthvtx2].annotation;
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
    else if (!strcasecmp(option, "--renumbered")) renumberedSpecified = 1;
    else if (!strcasecmp(option, "--scanner-ras")) UseScannerRAS = 1;
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
    else if (!strcasecmp(option, "--worst-bucket")) {
      if (nargc < 1) CMDargNErr(option,1);
      worstBucketFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--okayBucketMax")) {
      if (nargc < 1) CMDargNErr(option,1);
      long int val;
      sscanf(pargv[0],"%ld",&val);
      okayBucketMax = int(val);
      nargsused = 1;
    }     
    else if (!strcasecmp(option, "--grid")) {
      if (nargc < 3) CMDargNErr(option,0);
      char const * p = pargv[0];
      while (int c = *p++) { 
        switch (c) {
          case 0              :  break;
          case 'x' : case 'X' : gridx = 1; continue;   
          case 'y' : case 'Y' : gridy = 1; continue;   
          case 'z' : case 'Z' : gridz = 1; continue;   
          default:   fprintf(stderr, "--grid only supports {[xyz]}, not %s\n",pargv[1]); continue;
        }
        break;
      } 
      sscanf(pargv[1],"%f",&gridspacing);
      gridFile = pargv[2];
      nargsused = 3;
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
      // ignore --thresh for now
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--maxerrs")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&MAX_NUM_ERRORS);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--gdiag_no")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      printf("Gdiag_no %d\n",Gdiag_no );
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--min-dist")) {
      if(nargc < 4) CMDargNErr(option,4);
      surf1 = MRISread(pargv[0]);
      if(surf1==NULL) exit(1);
      surf2 = MRISread(pargv[1]);
      if(surf2==NULL) exit(1);
      if(UseScannerRAS){
	printf("Converting to scanner RAS\n");
	MRIStkr2Scanner(surf1);
	MRIStkr2Scanner(surf2);
      }

      // mindist will be on surf2
      int UseExact;
      sscanf(pargv[2],"%d",&UseExact);
      printf("Use Exact = %d\n",UseExact);
      MRI *mindist;
      if(UseExact){
        MRISdistanceBetweenSurfacesExact(surf2, surf1);
        mindist = MRIcopyMRIS(NULL, surf2, 0, "curv");
      }
      else 
        mindist = MRISminDist(surf1, surf2);
      if(mindist==NULL) exit(1);
      printf("Writing mindist to %s\n",pargv[3]);
      MRIwrite(mindist,pargv[3]);
      MRISfree(&surf1);
      MRISfree(&surf2);
      MRIfree(&mindist);
      printf("mris_diff done\n");
      exit(0);
      nargsused = 4;
    } 
    else if(!strcasecmp(option, "--simple") || !strcasecmp(option, "--si")) {
      if (nargc < 2) CMDargNErr(option,2);
      MRIS *surf1tmp = MRISread(pargv[0]);
      if(surf1tmp==NULL) exit(1);
      MRIS *surf2tmp = MRISread(pargv[1]);
      if(surf2tmp==NULL) exit(1);
      printf("Checking for differences between %s and %s\n",pargv[0],pargv[1]);
      printf("Gdiag_no %d\n",Gdiag_no );
      int res = MRISdiffSimple(surf1tmp, surf2tmp, 0, .00000001, 0);
      if(nargc > 2){
	printf("Writing RMS diff to %s\n",pargv[2]);
	const char **field;
	field = (const char **)calloc(sizeof(char*),1);
	field[0] = strcpyalloc("val");
	MRISwriteField(surf1tmp,field,1,pargv[2]);
      }
      exit(res);
    }
    else if (!strcasecmp(option, "--simple-patch")) {
      // surf surf2 patch2
      if (nargc < 3) CMDargNErr(option,3);
      MRIS *surf1tmp = MRISread(pargv[0]);
      if(surf1tmp==NULL) exit(1);
      int err = MRISreadPatch(surf1tmp,pargv[1]);
      if(err) exit(1);
      MRIS *surf2tmp = MRISread(pargv[0]); // re-read
      if(surf2tmp==NULL) exit(1);
      err = MRISreadPatch(surf2tmp,pargv[2]);
      if(err) exit(1);
      printf("Checking for differences in patches between %s and %s\n",pargv[0],pargv[1]);
      int res = MRISdiffSimple(surf1tmp, surf2tmp, 0, .00000001, 0);
      if(nargc > 2){
	printf("Writing RMS diff to %s\n",pargv[2]);
	const char **field;
	field = (const char **)calloc(sizeof(char*),1);
	field[0] = strcpyalloc("val");
	MRISwriteField(surf1tmp,field,1,pargv[2]);
      }
      exit(res);
    }
    else {
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
  printf("   --simple surf1 surf2 <rmsdiff.mgz>: just report whether the surfaces are different\n");
  printf("   --simple-patch surf patch1 patch2 : just report whether the patches are different\n");
  printf("   --thresh N    threshold (default=0) [note: not currently implemented!] \n");
  printf("   --maxerrs N   stop looping after N errors (default=%d)\n",
         MAX_NUM_ERRORS);
  printf("   --renumbered  the vertices or faces may have been renumbered and a few deleted\n");
  printf("   --worst-bucket worstbucketfile : compute the worst histogram bucket each vertex is in\n");
  printf("   --grid {[xyz]} spacingfloat grid_file : label the vertices of edges that span a grid\n");
  printf("\n");
  printf("   --no-check-xyz  : do not check vertex xyz\n");
  printf("   --no-check-nxyz : do not check vertex normals\n");
  printf("   --xyz-rms xyzrmsfile : compute and save rms diff between xyz\n");
  printf("   --angle-rms anglermsfile : compute angle on sphere between xyz\n");
  printf("   --seed seed : set random seed for degenerate normals\n");
  printf("   --min-dist surf1 surf2 exactflag mindist : compute vertex-by-vert RMS distance between surfs\n");
  printf("     surfs do not need to have the same number of vertices. Output on surf2\n");
  printf("   --scanner-ras : convert each surface to scanner RAS (applies to --min-dist as well)\n");
  printf("\n");
  printf("   --debug       turn on debugging\n");
  printf("   --gdiag_no Gdiag_no\n");
  printf("   --checkopts   don't run anything, just check options and exit\n");
  printf("   --help        print out information on how to use program\n");
  printf("   --version     print out version and exit\n");
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
  fprintf(fp,"%s\n", getVersion().c_str());
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
  fprintf(fp,"seed %ld\n",seed);
  fprintf(fp,"\n");
  fprintf(fp,"\n");

  return;
}

/*!
  \fn MRI *MRISminDist(MRIS *srcsurf, MRIS *trgsurf)
  \brief Computes the RMS distance between two surfaces when they do
  not have the same number of vertices. Output is in the target
  surface space. Algorithm is not perfect as it just computes the
  distance between vertices. So if there is a case where the vertex of
  one surface is very close to the face of the other but not close to
  a vertex, the distance could look largeish when it is really quite
  small.
 */
MRI *MRISminDist(MRIS *srcsurf, MRIS *trgsurf)
{
  int svtx = 0, tvtx;
  VERTEX *vtrg,*vsrc;
  float dmin;
  MHT *srchash = NULL, *trghash = NULL;
  MRI *mindist;

  mindist = MRIallocSequence(trgsurf->nvertices, 1, 1, MRI_FLOAT, 1);

  srchash = MHTcreateVertexTable_Resolution(srcsurf, CURRENT_VERTICES, 16);
  trghash = MHTcreateVertexTable_Resolution(trgsurf, CURRENT_VERTICES, 16);

  /* Go through the forward loop (finding closest srcvtx to each trgvtx).
  This maps each target vertex to a source vertex */
  for(tvtx = 0; tvtx < trgsurf->nvertices; tvtx++) {
    // Compute the source vertex that corresponds to this target vertex
    vtrg = &(trgsurf->vertices[tvtx]);
    svtx = MHTfindClosestVertexNo2(srchash, srcsurf, trgsurf,vtrg, &dmin);
    MRIsetVoxVal(mindist,tvtx,0,0,0,dmin);
  }
  // Go through the reverse loop
  for(svtx = 0; svtx < srcsurf->nvertices; svtx++) {
    vsrc = &(srcsurf->vertices[svtx]);
    tvtx = MHTfindClosestVertexNo2(trghash, trgsurf, srcsurf, vsrc, &dmin);
    if(dmin > MRIgetVoxVal(mindist,tvtx,0,0,0)) MRIsetVoxVal(mindist,tvtx,0,0,0,dmin);
  }
  MHTfree(&srchash);
  MHTfree(&trghash);
  return(mindist);
}

