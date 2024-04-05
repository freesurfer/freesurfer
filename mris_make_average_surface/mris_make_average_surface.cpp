/**
 * @brief create single surface from averages of a subject list
 *
 * This program will average the orig surfaces from the given subject
 * list into a single surface using Talairach coords and the spherical
 * transform.  The cooridinates of the vertices are the average of the
 * talairach coordinates (as defined by mri/transforms/talairach.xfm)
 * of the vertices from the input subjects.  The results will be saved
 * in a the specified subject's directory.  This default behavior can
 * be changed with option flags.
 */
/*
 * Original Author: Bruce Fischl
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

  This program will average the orig surfaces from the given subject
  list into a single surface using Talairach coords and the spherical
  transform.  The cooridinates of the vertices are the average of the
  talairach coordinates (as defined by mri/transforms/talairach.xfm)
  of the vertices from the input subjects.  The results will be saved
  in a the specified subject's directory.  This default behavior can
  be changed with option flags.

  The user must supply at least 4 arguments, plus subject list:

  mris_make_average_surface [options] hemi outsurfname cansurfname outsubject
    subj1 subj2 subj3 ...

  hemi - hemisphere,  lh or rh
  outsurfname - output surface name (eg, avg_orig)
  cansurfname - registration surface (eg, sphere.reg)
  outsubject  - name of subject to store the results in

  OPTIONS

  -help

  Print help and exit.

  -version

  Print version and exit

  -sdir sdir

  Use sdir instead of SUBJECTS_DIR

  -sdir-out sdirout

  Save results in sdirout/outsubject instead of SUBJECTS_IDR/outsubject.

  -nonorm

  Do not normalize area

  -i icoorder

  Use given icosahedron order (default is 7)

  -x xfmname

  Use transforms/xfmname instead of talairach.xfm

  -s surfname

  Use surfname instead of orig

  -v diagno

  Set Gdiag_no = diagno

ENDHELP
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/stat.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mri.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "macros.h"
#include "icosahedron.h"
#include "transform.h"
#include "version.h"
#include "fio.h"
#include "gca.h"
#include "gcamorph.h"
#include "mrisurf_metricProperties.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int normalize_area = 1 ;

static const char *orig_name = "orig" ;
static const char *xform_name = "talairach.xfm" ;

static int ico_no = 6 ;
int UseSurf2Surf = 1; // use surf2surf instead of parametric surface

const char *Progname ;
static char *sdir = NULL, *sdirout = NULL;

char *TargTempVolPath = NULL;
LTA *DestLTA=NULL;
int  RemoveIntersections = 0;
int Conform=0;

int main(int argc, char *argv[]) {
  char         **av, *avg_surf_name, *canon_surf_name, fname[STRLEN],
  *mdir, ico_fname[STRLEN], *hemi, *out_sname ;
  int          ac, nargs, i, vno, n, err ;
  MRI_SURFACE  *mris_ico, *surf ;
  MRI_SP       *mrisp_total ;
  LTA          *lta ;
  VOL_GEOM     vg;
  float        average_surface_area = 0.0 ;
  MATRIX *XFM=NULL;
  GCA_MORPH *gcam=NULL;
  MRI *mritemplate;

  memset((void *) &vg, 0, sizeof (VOL_GEOM));

  nargs = handleVersionOption(argc, argv, "mris_make_average_surface");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  mdir = getenv("FREESURFER_HOME") ;
  if (!mdir)
    ErrorExit(ERROR_BADPARM, 
              "%s: no FREESURFER_HOME in environment.\n",Progname);
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if (sdir == NULL) {
    sdir =  getenv("SUBJECTS_DIR");
    if (!sdir)
      ErrorExit(ERROR_BADPARM, 
                "%s: no SUBJECTS_DIR in environment.\n",Progname);
  }
  if (sdirout == NULL) sdirout = sdir;
  if (argc < 6) usage_exit() ;

  hemi = argv[1] ;
  avg_surf_name = argv[2] ;
  canon_surf_name = argv[3] ;
  out_sname = argv[4] ;

  printf("---------------------------------------------------\n");
  printf("hemi            = %s\n",hemi);
  printf("avg_surf_name   = %s\n",avg_surf_name);
  printf("canon_surf_name = %s\n",canon_surf_name);
  printf("out_sname       = %s\n",out_sname);
  printf("xform           = %s\n",xform_name);
  printf("---------------------------------------------------\n");
  printf("\n\n");
  fflush(stdout);

  // create the output directory now
  sprintf(fname, "%s/%s/surf", sdirout,out_sname);
  err = fio_mkdirp(fname,0777);
  if (err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",fname);
    perror(NULL);
    exit(1);
  }

  if(TargTempVolPath == NULL)
    sprintf(fname, "%s/average/mni305.cor.mgz", mdir);
  else
    sprintf(fname, "%s",TargTempVolPath);
  printf("Adding volume geometry from %s\n",fname);
  mritemplate = MRIread(fname);

  if(UseSurf2Surf){
    printf("Using surf2surf instead of parametric surface.\n");
    printf(" To use the old method include -no-surf2surf\n");

    AVERAGE_SURFACE_PARAMS *asp;
    asp = MRISaverageSurfaceParamAlloc(argc-5);
    asp->icoorder = ico_no;
    asp->hemi = hemi;
    asp->surfname = avg_surf_name;
    asp->surfregname = canon_surf_name;
    asp->xform_name = xform_name;
    if(DestLTA && Conform){
      printf("ERROR: cannot have -d (destination LTA) and -c (conform)\n");
      exit(1);
    }
    asp->DestLTA = DestLTA;
    asp->Conform = Conform;
    n=0;
    for (i = 5 ; i < argc ; i++) {
      asp->subjectlist[n] = strcpyalloc(argv[i]);
      n++;
    }
    surf = MakeAverageSurf(asp);
    MRISaverageSurfaceParamFree(&asp);

    initVolGeom(&surf->vg);
    getVolGeom(mritemplate, &surf->vg);
    MRIfree(&mritemplate);
    
    if(RemoveIntersections){
      printf("Remvoing intersections\n");
      int FillHoles = 1;
      mrisMarkIntersections(surf,FillHoles);
      int nintersections=0;
      for(int n=0; n < surf->nvertices; n++) if(surf->vertices[n].marked) nintersections++;
      printf("Found %d intersections\n",nintersections);
      MRISremoveIntersections(surf,FillHoles) ;
    }
    sprintf(fname, "%s/%s/surf/%s.%s", sdirout,out_sname, hemi, avg_surf_name) ;
    printf("writing average %s surface to %s\n", avg_surf_name, fname);
    MRISwrite(surf,fname) ;
    printf("#VMPC# mris_make_average_surface VmPeak  %d\n",GetVmPeak());
    printf("mris_make_average_surface done\n");
    exit(0);
    //====================================================================
  }

  printf("Using parametric surface instead of surf2surf.\n");
  printf(" To use the surf2surf specify -surf2surf\n");


#define SCALE 1
  mrisp_total = MRISPalloc(SCALE, 3) ;
  for (n = 0, i = 5 ; i < argc ; i++) {
    MRI *mri;
    MRI_SURFACE *mris;
    MRI_SP *mrisp;

    printf("\n---------------------------------------------------\n");
    printf("#@# processing subject %d/%d %s...\n", i-4,argc-5,argv[i]) ;
    fflush(stdout);

    // read sphere.reg
    sprintf(fname, "%s/%s/surf/%s.%s", sdir, argv[i], hemi, canon_surf_name) ;
    printf("  Reading %s\n",fname);
    fflush(stdout);
    mris = MRISread(fname) ;
    if (!mris) {
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
      exit(1);
    }
    // get "pial" surface vertex into ->origx, origy, origz
    if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
      ErrorExit(ERROR_BADFILE,"%s: could not read orig file for %s.\n",
                Progname, argv[1]);
    // read transform
    if (0) {
      sprintf(fname, "%s/%s/mri/transforms/%s", sdir, argv[i], xform_name) ;
      lta = LTAreadEx(fname) ;
      if (!lta)
        ErrorExit(ERROR_BADPARM, 
                  "%s: could not read transform from %s", Progname, fname) ;
    }

    // read T1 volume
    sprintf(fname, "%s/%s/mri/T1.mgz", sdir, argv[i]) ;
    if (fio_FileExistsReadable(fname)) mri = MRIreadHeader(fname, MRI_VOLUME_TYPE_UNKNOWN);
    else {
      sprintf(fname, "%s/%s/mri/T1", sdir, argv[i]) ;
      mri = MRIreadHeader(fname, MRI_UCHAR); // MRI_CORONAL_SLICE_DIRECTORY) ;
    }
    printf("  Read %s\n",fname);
    fflush(stdout);

    if (!mri)
      ErrorExit(ERROR_BADPARM, 
                "%s: could not read reference MRI volume from %s",
                Progname, fname) ;

    // save current vertex position into ->cx
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    // get the vertex position from ->origx, ... 
    // (get the "pial" vertex position)
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    printf("  Surface area: %2.1f cm^2\n", mris->total_area/100) ;
    fflush(stdout);
    average_surface_area += mris->total_area ;

    // this means that we transform "pial" surface

    if (xform_name)
    {
      if (!strcmp(xform_name,"talairach.xfm")) {
        printf("  Applying linear transform\n");
        fflush(stdout);
        XFM = DevolveXFMWithSubjectsDir(argv[i], NULL, "talairach.xfm", sdir);
        if (XFM == NULL) exit(1);
        MRISmatrixMultiply(mris, XFM);
        MatrixFree(&XFM);
      } else if (!strcmp(xform_name,"talairach.m3z")) {
        printf("  Applying GCA Morph\n");
        fflush(stdout);
        sprintf(fname, "%s/%s/mri/transforms/talairach.m3z", sdir, argv[i]) ;
        gcam = GCAMreadAndInvert(fname);
        if (gcam == NULL) exit(1);
        GCAMmorphSurf(mris, gcam);
        GCAMfree(&gcam);
      } else {
        printf("ERROR: don't know what to do with %s\n",xform_name);
        exit(1);
      }
    }

    // save transformed position in ->orig 
    // (store "pial" vertices position in orig)
    MRIScomputeMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    // get the vertex position from ->cx 
    // (note that this is not transformed)  sphere.reg vertices
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    // mris contains sphere.reg in vertex and pial vertices in orig
    // map to a theta-phi space and accumulate values
    mrisp = MRIScoordsToParameterization(mris, NULL, SCALE, ORIGINAL_VERTICES) ;
    MRISPaccumulate(mrisp, mrisp_total, 0) ;
    MRISPaccumulate(mrisp, mrisp_total, 1) ;
    MRISPaccumulate(mrisp, mrisp_total, 2) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris) ;
    MRIfree(&mri) ;
    //LTAfree(&lta) ;
    fflush(stdout);
    n++ ;
  }
  printf("Finished loading all data\n");
  average_surface_area /= (float)n ;
  printf("Avg surf area = %g cm\n",average_surface_area/100.0);
  fflush(stdout);

  // mrisp_total lost info on the modified surface
  sprintf(ico_fname, "%s/lib/bem/ic%d.tri", mdir, ico_no) ;
  printf("Reading icosahedron from %s...\n", ico_fname) ;
  mris_ico = ICOread(ico_fname) ;
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not read icosahedron file %s\n",
              Progname,ico_fname) ;
  MRISscaleBrain(mris_ico, mris_ico,
                 DEFAULT_RADIUS/MRISaverageRadius(mris_ico)) ;
  // save current ico position to ->cx, cy, cz
  MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES) ;
  // using mrisp_total to calculate position into ->origx, origy, origz 
  // (orig is the "pial" vertices)
  MRIScoordsFromParameterization(mrisp_total, mris_ico, ORIGINAL_VERTICES) ;
  // copy geometry info
  memcpy((void *) &mris_ico->vg, (void *) &vg, sizeof (VOL_GEOM));

  if (Gdiag_no >= 0 && Gdiag_no < mris_ico->nvertices) {
    int n ;

    VERTEX_TOPOLOGY const * const vt = &mris_ico->vertices_topology[Gdiag_no] ;
    VERTEX          const * const v  = &mris_ico->vertices         [Gdiag_no] ;
    printf( "v %d: x = (%2.2f, %2.2f, %2.2f)\n",
            Gdiag_no, v->origx, v->origy, v->origz) ;
    for (n = 0 ; n < vt->vnum ; n++) {
      VERTEX const * const vn = &mris_ico->vertices[vt->v[n]] ;
      printf( "v %d: x = (%2.2f, %2.2f, %2.2f)\n",
              vt->v[n], vn->origx, vn->origy, vn->origz) ;
    }
  }
  // write *h.sphere.reg
  sprintf(fname, "%s/%s/surf/%s.%s", 
          sdirout, out_sname, hemi, canon_surf_name) ;
  if (Gdiag & DIAG_SHOW)
    printf("writing average canonical surface to %s\n", fname);
  MRISwrite(mris_ico, fname) ;

  // get "pial vertices" from orig
  MRISrestoreVertexPositions(mris_ico, ORIG_VERTICES);
  
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++) {   // MRISscale looks suitable but sets other things also
    VERTEX * const v = &mris_ico->vertices[vno];
    MRISsetXYZ(mris_ico,vno,
    // n = number of subjects
      v->x / (float)n,
      v->y / (float)n,
      v->z / (float)n);
  }
  if (normalize_area) {
    MRIScomputeMetricProperties(mris_ico) ;
    
    printf("setting group surface area to be %2.1f cm^2 (scale=%2.2f)\n",
           average_surface_area/100.0,
           sqrt(average_surface_area/mris_ico->total_area)) ;

    mris_ico->group_avg_surface_area = average_surface_area ;
    
    MRIScomputeMetricProperties(mris_ico) ; // surely this is unnecessary
  }

  initVolGeom(&mris_ico->vg);
  getVolGeom(mritemplate, &mris_ico->vg);
  MRIfree(&mritemplate);

  // This catches cases where vertex 0 and vertex 40969 have the same coordinate
  if(mris_ico->nvertices == 163842)
    MRISfixAverageSurf7(mris_ico);

  sprintf(fname, "%s/%s/surf/%s.%s", sdirout,out_sname, hemi, avg_surf_name) ;
  printf("writing average %s surface to %s\n", avg_surf_name, fname);
  MRISwrite(mris_ico,  fname) ;

  if (0) {
    char path[STRLEN] ;
    LTA  *lta ;

    FileNamePath(fname, path) ;
    lta = LTAalloc(1, NULL) ;
    // write to a different location
    sprintf(fname, "%s/../mri/transforms/%s", path,xform_name) ;
    LTAwriteEx(lta, fname) ;
    LTAfree(&lta) ;
  }

  MRISfree(&mris_ico) ;
  MRISPfree(&mrisp_total) ;

  printf("#VMPC# mris_make_average_surface VmPeak  %d\n",GetVmPeak());
  printf("mris_make_average_surface done\n");

  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "help") || !stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "version") || !stricmp(option, "version"))
    print_version() ;
  else if (!stricmp(option, "sdir")) {
    sdir = argv[2];
    nargs = 1 ;
  } else if (!stricmp(option, "identity")) {
    xform_name = NULL ;
    printf("using identity transform\n") ;
  } else if (!stricmp(option, "sdir-out")) {
    sdirout = argv[2];
    nargs = 1 ;
  } 
  else if (!stricmp(option, "nonorm")) {
    normalize_area = 0 ;
    printf("not normalizing surface area\n") ;
  } 
  else if (!stricmp(option, "surf2surf")) {
    UseSurf2Surf = 1;
  } 
  else if (!stricmp(option, "no-surf2surf")) {
    UseSurf2Surf = 0;
  } 
  else if (!stricmp(option, "threads")) {
    int threads=1;
    sscanf(argv[2],"%d",&threads);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #endif
    nargs = 1 ;
  } 
  else if (!stricmp(option, "simple")) {
    // Stand-alone function to compute an average of surfaces
    // -simple averagesurf surf1 surf2 ...
    if(argc < 4){
      printf("ERROR: -simple averagesurf surf1 surf2 ...\n");
      exit(1);
    }
    std::vector<MRIS*> surfs;
    for(int n = 3; n < argc; n++){
      printf("Reading surf %s\n",argv[n]);
      MRIS *surf = MRISread(argv[n]);
      if(surf==NULL) exit(1);
      surfs.push_back(surf);
    }
    MRIS *avgsurf = MRISaverageSurfaces(surfs, NULL);
    if(avgsurf==NULL) exit(1);
    int err = MRISwrite(avgsurf,argv[2]);
    exit(err);
  } 
  else if (!stricmp(option, "conform")) {
    Conform = 1;
    printf("Mapping surfaces to conformed output space\n");
  }
  else if (!stricmp(option, "noconform")) {
    Conform = 0;
    printf("Not conforming output space\n");
  }
  else switch (toupper(*option)) {
    case 'I':
      ico_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'X':
      xform_name = argv[2] ;
      nargs = 1 ;
      printf("using xform %s...\n", xform_name) ;
      break ;
    case 'D':
      DestLTA = LTAread(argv[2]);
      printf("using dest LTA %s...\n", argv[2]) ;
      nargs = 1 ;
      break ;
    case 'R':
      RemoveIntersections = 1;
      printf("Removing intersections\n");
      break ;
    case 'T':
      TargTempVolPath = argv[2] ;
      nargs = 1 ;
      printf("using target %s...\n", TargTempVolPath);
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'S':
    case 'O':
      orig_name = argv[2] ;
      printf("reading vertex positions from %s...\n", orig_name) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    default:
      printf( "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  printf("mris_make_average_surface [options] hemi "
         "outsurfname cansurfname outsubject \n");
  printf("       subj1 subj2 subj3 ...\n");
  printf("  Run with -help for more info\n");
}

static void
print_help(void) {
  print_usage() ;
  printf("\n");
  printf("  This program will average the orig surfaces "
         "from the given subject\n");
  printf("  list into a single surface using "
         "Talairach coords and the spherical\n");
  printf("  transform.  The cooridinates of the vertices are the "
         "average of the\n");
  printf("  talairach coordinates (as defined by "
         "mri/transforms/talairach.xfm)\n");
  printf("  of the vertices from the input subjects.  "
         "The results will be saved\n");
  printf("  in a the specified subject's directory.  "
         "This default behavior can\n");
  printf("  be changed with option flags.\n");
  printf("\n");
  printf("  The user must supply at least 4 arguments, plus subject list:\n");
  printf("\n");
  printf("  mris_make_average_surface [options] "
         "hemi outsurfname cansurfname outsubject \n");
  printf("    subj1 subj2 subj3 ...\n");
  printf("\n");
  printf("  hemi - hemisphere,  lh or rh\n");
  printf("  outsurfname - output surface name (eg, avg_orig)\n");
  printf("  cansurfname - registration surface (eg, sphere.reg)\n");
  printf("  outsubject  - name of subject to store the results in\n");
  printf("\n");
  printf("  OPTIONS\n");
  printf("  -nonorm : Do not normalize area\n");
  printf("  -d DestLTA : apply LTA to the average surface\n");
  printf("  -i icoorder : Use given icosahedron order (default is 7)\n");
  printf("  -r : remove intersections in the output\n");
  printf("  -s surfname Use surfname instead of orig\n");
  printf("  -t templatename : Volume to use as geometry template for output surfaces\n");
  printf("  -x xfmname : Use transforms/xfmname instead of talairach.xfm\n");
  printf("  -conform : map output to the conformed space (-noconform)\n");
  printf("\n");
  printf("  -help Print help and exit.\n");
  printf("  -version\n");
  printf("  -sdir sdir : Use sdir instead of SUBJECTS_DIR\n");
  printf("  -sdir-out sdirout : Save results in sdirout/outsubject instead of SUBJECTS_IDR/outsubject.\n");
  printf("\n");
  printf("-surf2surf, -no-surf2surf\n");
  printf("  Use (don't use) surf2surf transform instead of parametric surface.\n");
  printf("  The parametric surface often creates large faces near the poles.\n");
  printf("\n");
  printf("-simple avgsurf surf1 surf2 ...\n");
  printf("  Stand-alone option to compute an average surface from the list.\n");
  printf("  All surfaces must have same number of vertices \n");
  printf("\n");
  printf("  -v diagno\n");
  printf("\n");
  printf("  Set Gdiag_no = diagno\n");
  printf("  \n");

  exit(1) ;
}

static void
print_version(void) {
  printf( "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

