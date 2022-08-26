/**
 * @brief optimally correcting the topology of triangulated surface
 *
 * "Genetic Algorithm for the Topology Correction of Cortical Surfaces",
 * F. Segonne, E. Grimson, B. Fischl
 * (2005) IPMI pp393--405.
 */
/*
 * Original Author: Florent Segonne
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

 
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "timer.h"
#include "topo_parms.h"

#include "mris_topology.h"
#include "patchdisk.h"

int main(int argc, char *argv[]) ;

static int asc = 0;
static int noint = 1 ;
static int sphere = 1 ;

static int  get_option(int argc, char *argv[]) ;
static void initTopoFixerParameters();
static void freeTopoFixerParameters();

const char *Progname ;

static const char *brain_name    = "brain" ;
static const char *wm_name       = "wm" ;
static const char *input_name   = "input" ;
static const char *sphere_name = "qsphere" ;
//static const char *defect_name  = "defects" ;
const char *orig_name = "orig" ;
const char *out_name = "orig_corrected" ;

static char sdir[STRLEN] = "";
static TOPOFIX_PARMS parms;
static int MGZ = 1; // set to 1 for MGZ

static double pct_over = 1.1;

static void initTopoFixerParameters() {
  parms.verbose = 1; //minimal mode
  parms.smooth = 0;    //smoothing
  parms.match = 1;     //matching using local intensity estimates
  parms.mode = 0 ;
  parms.minimal_mode=0;
  parms.nminattempts = 10;
  parms.l_mri = 0.0;
  parms.l_curv = 4.0 ;
  parms.l_qcurv = 0.0;
  parms.l_unmri = 1.0;
  parms.volume_resolution = 3;
  parms.nattempts_percent=0.15;
  parms.minimal_loop_percent = 0.4;
  parms.no_self_intersections = 1;
  parms.contrast = -2; //contrast (regular==1;inversion== -1 ; detect = -2)
  //does not write out information
  parms.write = 0;

  //creation of the patching surfaces
  PatchDisk *disk = new PatchDisk[4];
  for (int n = 0 ; n < 4 ; n++)
    disk[n].Create(n);
  parms.patchdisk = (void*)disk;

  //mri
  parms.mri = NULL;
  parms.mri_wm = NULL;
  parms.mri_gray_white = NULL;
  parms.mri_k1_k2 = NULL;
  //histo
  parms.h_k1 = NULL;
  parms.h_k2 = NULL;
  parms.h_gray = NULL;
  parms.h_white = NULL;
  parms.h_dot = NULL;
  parms.h_border = NULL;
  parms.h_grad = NULL;
  ;
  parms.transformation_matrix=NULL;
  parms.defect_list=NULL;
}

static void freeTopoFixerParameters() {
  PatchDisk *disk = (PatchDisk*)parms.patchdisk;
  delete [] disk;
  //mri
  if (parms.mri) MRIfree(&parms.mri);
  if (parms.mri_wm) MRIfree(&parms.mri_wm);
  if (parms.mri_gray_white) MRIfree(&parms.mri_gray_white);
  if (parms.mri_k1_k2) MRIfree(&parms.mri_k1_k2);
  //histo
  if (parms.h_k1) HISTOfree(&parms.h_k1);
  if (parms.h_k2) HISTOfree(&parms.h_k2);
  if (parms.h_gray) HISTOfree(&parms.h_gray);
  if (parms.h_white) HISTOfree(&parms.h_white);
  if (parms.h_dot) HISTOfree(&parms.h_dot);
  if (parms.h_border) HISTOfree(&parms.h_border);
  if (parms.h_grad) HISTOfree(&parms.h_grad);
  if (parms.transformation_matrix)  MatrixFree(&parms.transformation_matrix);
  //defect_list
  if (parms.defect_list) {
    DEFECT_LIST *dl = (DEFECT_LIST*)parms.defect_list;
    for (int i = 0 ; i < dl->ndefects ; i++) {
      if (dl->defects[i].vertices)
        free(dl->defects[i].vertices) ;
      if (dl->defects[i].status)
        free(dl->defects[i].status) ;
      if (dl->defects[i].border)
        free(dl->defects[i].border) ;
      if (dl->defects[i].edges)
        free(dl->defects[i].edges);
    }
    free(dl) ;
  }
}


int main(int argc, char *argv[]) {

  char          *hemi, *sname, *cp, fname[STRLEN] ;
  int           nargs ;
  MRI_SURFACE   *mris, *mris_corrected ;
  // MRI           *mri, *mri_wm ;
  int           msec, nvert, nfaces, nedges, eno ,is_valid;
  Timer then ;

  std::string cmdline = getAllInfo(argc, argv, "mris_topo_fixer");

  nargs = handleVersionOption(argc, argv, "mris_topo_fixer");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  //parameters for mris_topo_fixer
  initTopoFixerParameters();

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  /*if(argc < 2)  usage_exit() ;

  print_parameters();

  printf("%s\n",getVersion().c_str());
  printf("  %s\n",getVersion().c_str());
  fflush(stdout); */

  then.reset() ;
  sname = argv[1] ;
  hemi = argv[2] ;
  if (strlen(sdir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  int req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.%s", sdir, sname, hemi,orig_name) ;
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("reading input surface %s...\n", fname) ;
  mris = MRISreadOverAlloc(fname,pct_over) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read input surface %s",
              Progname, fname) ;
  MRISaddCommandLine(mris, cmdline) ;
  strcpy(mris->subject_name, sname) ;

  eno = MRIScomputeEulerNumber(mris, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "Before topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;
  if (eno == 2 ) {
    fprintf(stderr,"The Euler Number of this surface is 2"
            "\nNothing to do\nProgram exiting\n");
    exit(0);
  }

  is_valid = MRISisSurfaceValid(mris,0,1);
  if (is_valid == 0 ) {
    fprintf
    (stderr,
     "The original surface is not a valid tessellated manifold!!!\n");
    fprintf(stderr,"\nAbort !!!\n");
    MRISfree(&mris);
    exit(-1);
  } else
    fprintf(stderr,"   The original surface is a valid manifold\n");

	//checking if we have only one single component
	fprintf(stderr,"   Counting the number of connected components\n");
	int ncpts;
	int did_extract_main_component = 0;
	MRISextractMainComponent(mris,1, 0, &ncpts);
	if(ncpts != 1){
		fprintf
      (stderr,
       "The original surface has more than one component (%d components)!!\n",
       ncpts);
		fprintf
      (stderr,
       "   ->   Extracting the largest component of the initial surface\n");
    /*before extracting the largest component, 
      we will load the spherical coordinates*/
		MRISsaveVertexPositions(mris,ORIGINAL_VERTICES);
		req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.%s", sdir, sname, hemi, sphere_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
		if (MRISreadVertexPositions(mris, (char*)sphere_name) != NO_ERROR)
				ErrorExit(ERROR_NOFILE, "%s: could not read original surface %s",
									Progname, orig_name) ;
		MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
		MRISrestoreVertexPositions(mris,ORIGINAL_VERTICES);
		MRIS* mris_temp = MRISextractMainComponent(mris,0 , 0, &ncpts);
    MRISfree(&mris);
		mris=mris_temp;
		did_extract_main_component = 1;
		//fprintf(stderr,"\nAbort !!!\n");
    //MRISfree(&mris);
    //exit(-1);
	}else
		fprintf(stderr,"   The original surface has one component\n");

  if (parms.no_self_intersections) {
    int self_intersect =  IsMRISselfIntersecting(mris);
    if (self_intersect) {
      fprintf(stderr,"The original surface self-intersects!!!\n");
      fprintf(stderr,"\nAbort !!!\n");
      MRISfree(&mris);
      exit(-1);
    } else
      fprintf(stderr,"   The original surface does not self-intersect\n");
  }

	fprintf(stderr,"\n");
  MRISsaveVertexPositions(mris,ORIGINAL_VERTICES);
  MRISaverageVertexPositions(mris, 2) ;
  if (parms.no_self_intersections) {
    int self_intersect =  IsMRISselfIntersecting(mris);
    if (self_intersect) {
      fprintf
      (stderr,
       "Smoothing led to self-intersection: restoring original vertices\n");
      MRISrestoreVertexPositions(mris,ORIGINAL_VERTICES);
    };
  }
  MRISsaveVertexPositions(mris,ORIGINAL_VERTICES);

  req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.%s", sdir, sname, hemi, input_name) ;
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  MRISwrite(mris,fname);

  //number of loops
  int nloops = (2-eno)/2;
  fprintf(stderr,"The surface has %d loops (X=%d)\n\n",nloops,eno);

  //  if(curv){
  //    //read the topological defects
  //    sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, defect_name) ;
  //    MRISreadCurvature(mris,fname);
  //    for(int n = 0 ; n < mris->nvertices ; n++)
  //      mris->vertices[n].marked2 = int(mris->vertices[n].curv);
  //  }else {
  if (sphere) {
		if(did_extract_main_component==0){
			//read the spherical coordinates
			req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.%s", sdir, sname, hemi, sphere_name) ;
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
			if (MRISreadVertexPositions(mris, (char*)sphere_name) != NO_ERROR)
				ErrorExit(ERROR_NOFILE, "%s: could not read original surface %s",
									Progname, orig_name) ;
		}else //surface already extracted
			MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  } else
    MRISmapOntoSphere(mris);

  //fprintf(stderr,"smoothing spherical vertices...\n");
  //MRISsmoothOnSphere(mris,300);
  //MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  /* centering the surface using CANONICAL_VERTICES */
  MRIScenterSphere(mris);
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  fprintf(stderr,"identify defects...\n");
  MRISidentifyDefects(mris,&parms);
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  //  }

  //read the mri volume
  if (MGZ) {
    req = snprintf(fname, STRLEN, "%s/%s/mri/%s.mgz", sdir, sname, brain_name);
  } else {
    req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, brain_name);
  }
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("reading brain volume from %s...\n", brain_name) ;
  parms.mri = MRIread(fname) ;
  if (!parms.mri)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read brain volume from %s", Progname, fname) ;

  if (MGZ) {
    req = snprintf(fname, STRLEN, "%s/%s/mri/%s.mgz", sdir, sname, wm_name);
  }
  else {
    req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, wm_name);
  }
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("reading wm segmentation from %s...\n", wm_name) ;
  parms.mri_wm = MRIread(fname) ;
  if (!parms.mri_wm)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read wm volume from %s", Progname, fname) ;

  DEFECT_LIST *dl=(DEFECT_LIST*)parms.defect_list;
  for (int n = 0 ; n < mris->nvertices ; n++) {
    mris->vertices[n].curv=0;
    mris->vertices[n].marked2=0;
  }
  for (int i = 0 ; i < dl->ndefects ; i++) {
    DEFECT *defect = &dl->defects[i];
    for (int n = 0 ; n < defect->nvertices ; n++) {
      mris->vertices[defect->vertices[n]].curv=i+1;
      mris->vertices[defect->vertices[n]].marked2=i+1;
    }
    for (int n = 0 ; n < defect->nborder ; n++) {
      mris->vertices[defect->border[n]].curv=i+1;
      mris->vertices[defect->border[n]].marked2=i+1;
    }
  }

  req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.defects", sdir, sname, hemi ) ;
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  MRISwriteCurvature(mris,fname);

  mris_corrected=MRISduplicateOver(mris);
  mris_corrected->type=MRIS_TRIANGULAR_SURFACE;
  //mris becomes useless!
  MRISfree(&mris) ;

  MRISsaveVertexPositions(mris_corrected,ORIGINAL_VERTICES);
  MRIScomputeMetricProperties(mris_corrected);

  MRISinitTopoFixParameters(mris_corrected,&parms);

  int def=-1;
  if (argc > 3) def = atoi(argv[3]);

  dl=(DEFECT_LIST*)parms.defect_list;
  for (int i = 0 ; i < dl->ndefects ; i++) {
    //    fprintf
    // (stderr, " defect %d has %d vertices\n",i,dl->defects[i].nvertices);
    if (def > 0 && def != i+1) continue;
    DEFECT *defect = &dl->defects[i];
    for (int n = 0 ; n < mris_corrected->nvertices ; n++)
      mris_corrected->vertices[n].marked2=0;
    for (int n = 0 ; n < defect->nvertices ; n++)
      mris_corrected->vertices[defect->vertices[n]].marked2=i+1;
    for (int n = 0 ; n < defect->nborder ; n++)
      mris_corrected->vertices[defect->border[n]].marked2=i+1;
    parms.defect_number = i+1;
    MRIScorrectDefect(mris_corrected,i+1,parms);
    fprintf
    (stderr,
     "AFTER CORRECTION, EULER IS %d \n",
     MRISgetEuler(mris_corrected));
    if (0) {//parms.no_self_intersections){
      int self_intersect =  IsMRISselfIntersecting(mris_corrected);
      if (self_intersect) {
        fprintf(stderr,"\nThe final surface self-intersects !\n");
        exit(-1);
      }
    }
    //check if the surface is still a valid one!
    if (parms.verbose >= VERBOSE_MODE_MEDIUM)
      is_valid = MRISisSurfaceValid(mris_corrected,0,1);
    else
      is_valid = MRISisSurfaceValid(mris_corrected,0,0);
    if (is_valid == 0 ) {
      if (parms.verbose < VERBOSE_MODE_MEDIUM)
        MRISisSurfaceValid(mris_corrected,0,1);
      fprintf
      (stderr,
       "The original surface is not a valid tessellated "
       "manifold anymore!!!\\n");
      fprintf(stderr,"\nAbort !!!\n");
      MRISfree(&mris_corrected);
      ErrorExit
      (ERROR_BADPARM,
       "Defect #%d: The topology correction generated an invalid surface\n",
       def) ;
    };
  }

  //checking if we have some self-intersections (should not happen)
  if (parms.no_self_intersections) {
    int self_intersect =  IsMRISselfIntersecting(mris_corrected);
    if (self_intersect) {
      fprintf(stderr,"\nThe final surface self-intersects !\n");
      if (noint) MRISremoveIntersections(mris_corrected,0) ;
    } else
      fprintf(stderr,"\nThe final surface does not self-intersect !\n");
  }

  eno = MRIScomputeEulerNumber(mris_corrected, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "after topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;

  /* compute the orientation changes */
  MRISmarkOrientationChanges(mris_corrected);

  if (asc) {
    req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.%s.asc", sdir, sname, hemi,out_name) ;
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISwrite(mris_corrected,fname);
  } else {
    req = snprintf(fname, STRLEN,  "%s/%s/surf/%s.%s", sdir, sname, hemi,out_name) ;
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISwrite(mris_corrected,fname);
  }

  fprintf(stderr,"\n\n");

  msec = then.milliseconds() ;
  fprintf(stderr,"topology fixing took %2.1f minutes\n",
          (float)msec/(60*1000.0f));

  freeTopoFixerParameters();

  return(0) ;  /* for ansi */
}


static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, (char*)"-help"))
    exit(1);//print_help() ;
  else if (!stricmp(option, (char*)"mgz")) {
    printf("INFO: assuming .mgz format\n");
    MGZ = 1;
  } else if (!stricmp(option , (char*)"fast")) {
    parms.mode=1;
    fprintf(stderr,"fast mode on\n");
  } else if (!stricmp(option , (char*)"verbose")) {
    Gdiag = DIAG_VERBOSE;
  } else if (!stricmp(option ,(char*) "write")) {
    parms.write=1;
  } else if (!stricmp(option, (char*)"verbose_low")) {
    parms.verbose=VERBOSE_MODE_LOW;
    fprintf(stderr,"verbose mode on (default+low mode)\n");
    nargs = 0 ;
  } else if (!stricmp(option, (char*)"warnings")) {
    parms.verbose=VERBOSE_MODE_MEDIUM;
    fprintf(stderr,"verbose mode on (medium mode): printing warnings\n");
    nargs = 0 ;
  } else if (!stricmp(option, (char*)"errors")) {
    parms.verbose=VERBOSE_MODE_HIGH;
    fprintf(stderr,"verbose mode on (high mode): "
            "exiting when warnings appear\n");
    nargs = 0 ;
  } else if (!stricmp(option, (char*)"mri")) {
    parms.l_mri = atof(argv[2]) ;
    fprintf(stderr,"setting l_mri = %2.2f\n", parms.l_mri) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"curv")) {
    parms.l_curv = atof(argv[2]) ;
    fprintf(stderr,"setting l_curv = %2.2f\n", parms.l_curv) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"qcurv")) {
    parms.l_qcurv = atof(argv[2]) ;
    fprintf(stderr,"setting l_qcurv = %2.2f\n", parms.l_qcurv) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"unmri")) {
    parms.l_unmri = atof(argv[2]) ;
    //parms.l_unmri = 0.0;
    fprintf(stderr,"setting l_unmri = %2.2f\n", parms.l_unmri) ;
    //      if(parms.volume_resolution<1)
    //parms.volume_resolution = 2;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"pct")) {
    parms.nattempts_percent = atof(argv[2]) ;
    fprintf
    (stderr,
     "setting pct of attempts = %2.2f\n",
     parms.nattempts_percent) ;
    nargs = 1 ;
  } else if (!stricmp(option,(char*) "inverted_contrast")) {
    parms.contrast = -1 ;
    fprintf(stderr,"contrast inversion\n") ;
  } else if (!stricmp(option, (char*)"detect_contrast")) {
    parms.contrast = -2 ;
    fprintf(stderr,"automatically detecting contrast inversion\n") ;
  } else if (!stricmp(option, (char*)"usual_contrast")) {
    parms.contrast = 1 ;
    fprintf(stderr,"contrast in right direction\n") ;
  } else if (!stricmp(option, (char*)"nmin")) {
    parms.nminattempts = __MAX(1,atoi(argv[2])) ;
    fprintf
    (stderr,
     "setting minimal number of attempts = %d\n",
     parms.nminattempts) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"asc")) {
    asc=1;
    fprintf(stderr,"writting out surface in asci mode\n") ;
  } else if (!stricmp(option, (char*)"int")) {
    noint = 0 ;
    nargs = 0 ;
    parms.no_self_intersections = 0 ;
  } else if (!stricmp(option,(char*) "no_intersection")) {
    parms.no_self_intersections = 1;
    fprintf(stderr,"avoiding self-intersecting patches\n") ;
  } else if (!stricmp(option, (char*)"minimal")) {
    parms.minimal_mode=1;
    parms.nminattempts = 1;
    parms.nattempts_percent=0.0f;
    fprintf(stderr,"cuting minimal loop only\n") ;
  } else if (!stricmp(option, (char*)"loop_pct")) {
    parms.minimal_loop_percent = atof(argv[2]) ;
    fprintf
    (stderr,
     "setting loop_pct = %2.2f\n",
     parms.minimal_loop_percent) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"smooth")) {
    parms.smooth = atoi(argv[2]) ;
    fprintf(stderr,"smoothing defect with mode %d\n", parms.smooth) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"match")) {
    if (atoi(argv[2])) parms.match = 1 ;
    else  parms.match = 0;
    fprintf(stderr,"matching mode : %d \n", parms.match) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"seed")) {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr, "setting seed for random number genererator to %d\n", atoi(argv[2])) ;
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"out_name")) {
    out_name = argv[2];
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"orig_name")) {
    orig_name = argv[2];
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      //      print_usage() ;
      exit(1) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    }

  return(nargs) ;
}

