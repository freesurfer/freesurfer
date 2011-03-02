/**
 * @file  mris_show.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:34 $
 *    $Revision: 1.40 $
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <GL/gl.h>
#include "glut.h"

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "oglutil.h"
#include "label.h"
#include "version.h"

static char vcid[] = "$Id: mris_show.c,v 1.40 2011/03/02 00:04:34 nicks Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

#define BIG_FOV                 300   /* for unfolded hemispheres */
#define SMALL_FOV               200

#define FRAME_SIZE         600

#define DELTA_ANGLE  16.0f
#define RIGHT_HEMISPHERE_ANGLE   90.0
#define LEFT_HEMISPHERE_ANGLE   -90.0

#define ORIG_SURFACE_LIST   1
#define ELLIPSOID_LIST      1
#define MRISP_LIST          2
#define NLISTS              2+1  /* 1 based */

#ifndef GLUT_ESCAPE_KEY
#define GLUT_ESCAPE_KEY  27
#endif

#define MAX_MARKED       20000

static int which_norm = NORM_MEAN;

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
/* static void rotate_brain(float angle, char c) ;*/

static void keyboard_handler(unsigned char key, int x, int y) ;
static void special_key_handler(int key, int x, int y) ;
static void reshape_handler(int width, int height) ;
static void mouse_handler(int button, int state, int x, int y) ;
#if 0
static void controlLights(int value) ;
#endif

static void display_handler(void) ;
static void home(MRI_SURFACE *mris) ;

int MRIScheckSphericalOrientation(MRI_SURFACE *mris) ;
int MRISinvertSurface(MRI_SURFACE *mris) ;
int MRISfillSurface(MRI_SURFACE *mris) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;
static char *surf_fname ;
static MRI_SURFACE  *mris ;
static int nmarked = 0 ;
static int marked_vertices[MAX_MARKED+1] = {
      -1
    } ;
static long frame_xdim = FRAME_SIZE;
static long frame_ydim = FRAME_SIZE;

static int talairach_flag = 0 ;
static int ellipsoid_flag = 0 ;

static float x_angle = 0.0f ;
static float y_angle = 0.0f ;
static float z_angle = 0.0f ;

static float delta_angle = 2*DELTA_ANGLE ;


static int current_list = ORIG_SURFACE_LIST ;
static int normalize_flag = 0 ;
static INTEGRATION_PARMS  parms ;


static void findAreaExtremes(MRI_SURFACE *mris) ;
static char curvature_fname[100] = "" ;

static int compiled[NLISTS] ;

static float cslope = 1.75f ;
static int patch_flag = 0 ;
static int compile_flags = 0 ;
static int mean_curvature_flag = 0 ;
static int gaussian_curvature_flag = 0 ;
static int fit_flag = 0 ;
static int fov = -1 ;
static int noscale = 1 ;
static MRI_SURFACE_PARAMETERIZATION *mrisp = NULL ;
static int navgs = 0 ;
static int nbrs = 2 ;
static int nonmax_flag = 0 ;
static char *coord_fname = NULL ;
static int param_no = -1 ;
static int normalize_param = 0 ;


static float total_alpha = 0.0f ;
static float total_beta = 0.0f ;
static float total_gamma = 0.0f ;

static char wname[200] ;
static char title[200] ;

static double starting_mse = 0.0f ;

static float scale = 1.0f ;

static float lighting_offset = 0.25f ;
static float light = 0.0f ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, *out_fname, fname[100], hemi[10],
  *cp, path[100], name[100] ;
  int          ac, nargs, i ;
  float        angle ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_show.c,v 1.40 2011/03/02 00:04:34 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  memset(compiled, 0, sizeof(compiled)) ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.projection = PROJECT_ELLIPSOID ;
  parms.tol = TOL ;
  parms.n_averages = N_AVERAGES ;
  parms.l_angle = L_ANGLE ;
  parms.l_area = L_AREA ;
  parms.niterations = NITERATIONS ;
  parms.write_iterations = WRITE_ITERATIONS ;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit() ;

  surf_fname = in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (patch_flag && !strstr(surf_fname, ".geo")) {
    /* read in orig surface before reading in patch */
    FileNamePath(surf_fname, path) ;
    FileNameOnly(surf_fname, name) ;
    cp = strchr(name, '.') ;
    if (cp) {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    } else
      strcpy(hemi, "lh") ;
    sprintf(fname, "%s/%s.orig", path, hemi) ;
    mris = MRISfastRead(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    FileNameOnly(surf_fname, name) ;
    if (MRISreadPatch(mris, name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, surf_fname) ;
  } else {
    mris = MRISfastRead(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }
  if (!FEQUAL(scale,1)) {
    MRIScenter(mris, mris) ;
    MRISscaleBrain(mris, mris,scale) ;
  }
  if (coord_fname) {
    if (MRISreadCanonicalCoordinates(mris, coord_fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: could not read canonical coordinate system",
                Progname) ;
  } else
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  /*  MRIScomputeMetricProperties(mris) ;*/
  if (mean_curvature_flag || gaussian_curvature_flag) {
    MRISsetNeighborhoodSize(mris, nbrs) ;  /* so curvatures are reasonable */
    MRIScomputeSecondFundamentalForm(mris) ;
    if (mean_curvature_flag)
      MRISuseMeanCurvature(mris) ;
    if (gaussian_curvature_flag)
      MRISuseGaussianCurvature(mris) ;
  }
  if (curvature_fname[0])
    MRISreadCurvatureFile(mris, curvature_fname) ;
  if (mrisp) {
    MRISnormalizeCurvature(mris, which_norm) ;
    MRISstoreMeanCurvature(mris) ;
    starting_mse = MRIScomputeCorrelationError(mris, mrisp, param_no) ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    if (normalize_param)
      MRISnormalizeFromParameterization(mrisp, mris, param_no) ;
    else
      MRISfromParameterization(mrisp, mris, param_no) ;
    MRISnormalizeCurvature(mris, which_norm) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    current_list = MRISP_LIST ;
  }
  if (ellipsoid_flag) {
    MRISupdateEllipsoidSurface(mris) ;
    findAreaExtremes(mris) ;
  }
  if (normalize_flag)
    MRISnormalizeCurvature(mris, which_norm) ;
  MRISaverageCurvatures(mris, navgs) ;
  if (nonmax_flag)
    /*MRISnonmaxSuppress(mris)*/ ;
  if (talairach_flag)
    MRIStalairachTransform(mris, mris) ;
  MRIScenter(mris, mris) ;
  if (noscale)
    OGLUnoscale() ;
  if (fov > 0)
    OGLUsetFOV(fov) ;
  else   /* try and pick an appropriate one */
  {
    if (mris->patch) {
      int   ngood, vno ;
      float pct_vertices ;

      for (vno = ngood = 0 ; vno < mris->nvertices ; vno++)
        if (!mris->vertices[vno].ripflag)
          ngood++ ;
      pct_vertices = (float)ngood / (float)mris->nvertices ;
      if (pct_vertices > 0.75f)
        OGLUsetFOV(BIG_FOV) ;
      else if (pct_vertices < 0.4f)
        OGLUsetFOV(SMALL_FOV) ;
    }
  }
#if 0
  fprintf(stderr, "(%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f), ctr "
          "(%2.1f, %2.1f, %2.1f)\n",
          mris->xlo, mris->ylo, mris->zlo, mris->xhi, mris->yhi, mris->zhi,
          mris->xctr, mris->yctr, mris->zctr);
#endif

  glutInit(&argc, argv) ;  /* pass -geometry and such to  GLUT */
  glutInitWindowSize(frame_xdim, frame_ydim) ;

  /*
     use double buffering, RBF color stuff, and a depth buffer for
     hidden surface removal
     */
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH) ;
  sprintf(wname, "%scortical surface reconstruction '%s'",
          talairach_flag ? "talairached " : "", in_fname) ;
  strcpy(title, wname) ;
  glutCreateWindow(title) ;          /* open the window */
  glutHideWindow() ;
  glutDisplayFunc(display_handler) ; /* specify a drawing callback function */
  OGLUinit(mris, frame_xdim, frame_ydim) ; /* specify lighting and such */
  if (fit_flag)
    oglu_fov = 0.0 ;
  if (!FZERO(light))
    OGLUsetLightingModel(-1.0f, -1.0f, -1.0f, -1.0f, light) ;

  if (compile_flags & TP_FLAG) {
    MRISsetNeighborhoodSize(mris, nbrs) ;
    MRIScomputeSecondFundamentalForm(mris) ;
  }

  /* now compile the surface tessellation */
  glNewList(current_list, GL_COMPILE) ;
  compiled[current_list] = 1 ;
  for (i = 0 ; i < nmarked ; i++) {
    if (marked_vertices[i] >= mris->nvertices) {
      nmarked = i ;
      marked_vertices[i] = -1 ;
    } else
      mris->vertices[marked_vertices[i]].marked = 1 ;
  }
  OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
  glEndList() ;

  glMatrixMode(GL_MODELVIEW);
  if (mris->patch) {
    angle = 90.0f ;
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  } else {
    if (mris->hemisphere == RIGHT_HEMISPHERE)
      angle = RIGHT_HEMISPHERE_ANGLE ;
    else
      angle = LEFT_HEMISPHERE_ANGLE ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  }


  /* handle reshaping, keyboard and mouse events */
  glutKeyboardFunc(keyboard_handler) ;
  glutSpecialFunc(special_key_handler) ;
  glutMouseFunc(mouse_handler) ;
  glutReshapeFunc(reshape_handler) ;

  /* create a menu and attach it to the window */
#if 0
  glutCreateMenu(controlLights) ;
  glutAddMenuEntry("toggle light 1", 1) ;
  glutAddMenuEntry("toggle light 2", 2) ;
  glutAddMenuEntry("toggle light 3", 3) ;
  glutAddMenuEntry("toggle light 4", 4) ;
  glutAttachMenu(GLUT_RIGHT_BUTTON) ;
#endif

  glutShowWindow() ;

  glRotatef(x_angle, 1.0f, 0.0f, 0.0f) ;
  glRotatef(y_angle, 0.0f, 1.0f, 0.0f) ;
  glRotatef(z_angle, 0.0f, 0.0f, 1.0f) ;
  glutPostRedisplay() ;

  glutMainLoop() ;               /* enter event handling loop */
  MRISfree(&mris) ;

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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "light")) {
    light = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting lighting to %2.2f\n", light) ;
  } else if (!stricmp(option, "nparam")) {
    char *cp ;
    cp = strchr(argv[2], '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    } else
      param_no = 0 ;
    mrisp = MRISPread(argv[2]) ;
    normalize_param = 1 ;
    if (!mrisp)
      ErrorExit(ERROR_NOFILE, "%s: could not read parameterization file %s",
                Progname, argv[2]) ;
    fprintf(stderr, "reading parameterized curvature from %s\n", argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "param")) {
    char *cp ;
    cp = strchr(argv[2], '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    } else
      param_no = 0 ;
    mrisp = MRISPread(argv[2]) ;
    if (!mrisp)
      ErrorExit(ERROR_NOFILE, "%s: could not read parameterization file %s",
                Progname, argv[2]) ;
    fprintf(stderr, "reading parameterized curvature from %s\n", argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "bw"))
    compile_flags |= BW_FLAG ;
  else if (!stricmp(option, "coord")) {
    coord_fname = argv[2] ;
    fprintf(stderr, "reading canonical coordinate system from %s\n", argv[2]) ;
    nargs = 1 ;
    compile_flags |= COORD_FLAG ;
  } else if (!stricmp(option, "canon")) {
    coord_fname = argv[2] ;
    fprintf(stderr, "reading canonical coordinate system from %s\n", argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "cslope")) {
    sscanf(argv[2], "%f", &cslope) ;
    nargs = 1 ;
    fprintf(stderr, "using color slope compression = %2.4f\n", cslope) ;
  } else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using %d neighbors for curvature calculation\n", nbrs) ;
  } else if (!stricmp(option, "tp"))
    compile_flags |= TP_FLAG ;
  else if (!stricmp(option, "mesh"))
    compile_flags |= MESH_FLAG ;
  else if (!stricmp(option, "noscale"))
    noscale = 1 ;
  else if (!stricmp(option, "scale"))
    noscale = 0 ;
  else if (!stricmp(option, "rescale")) {
    scale = atof(argv[2]) ;
    fprintf(stderr, "scaling brain by %2.2f\n", scale) ;
    nargs = 1 ;
  } else if (!stricmp(option, "nonmax"))
    nonmax_flag = 1 ;
  else if (!stricmp(option, "fov")) {
    fov = atoi(argv[2]) ;
    fprintf(stderr, "using fov = %d\n", fov) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "averaging curvature pattern %d times.\n", navgs) ;
      nargs = 1 ;
      break ;
    case 'M':
      mean_curvature_flag = 1 ;
      break ;
    case 'G':
      gaussian_curvature_flag = 1 ;
      break ;
    case 'P':
      patch_flag = 1 ;
      compile_flags |= PATCH_FLAG ;
      nargs = 0 ;
      break ;
    case 'C':
      strcpy(curvature_fname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'L': {
      LABEL *area ;
      int   i ;

      fprintf(stderr, "reading label file %s\n",argv[2]) ;
      nargs = 1 ;
      area = LabelRead(NULL, argv[2]) ;
      if (nmarked+area->n_points >= MAX_MARKED)
        area->n_points = MAX_MARKED - nmarked ;
      for (i = 0 ; i < area->n_points ; i++)
        marked_vertices[nmarked++] = area->lv[i].vno ;
      LabelFree(&area) ;
    }
    case 'V':
      if (nmarked >= MAX_MARKED)
        ErrorExit(ERROR_NOMEMORY, "%s: too many vertices marked (%d)",
                  Progname, nmarked) ;
      if (*argv[2] == '#' || *argv[2] == ':') {
        int  skip, vno, vindex ;

        sscanf(argv[2], "#%d", &skip) ;
        for (vindex = vno = 0 ; vindex < MAX_MARKED ; vindex++, vno += skip)
          marked_vertices[vindex] = vno ;
        nmarked = MAX_MARKED ;
      } else
        sscanf(argv[2], "%d", &marked_vertices[nmarked++]) ;
      marked_vertices[nmarked] = -1 ;  /* terminate list */
      nargs = 1 ;
      break ;
    case 'X':
      sscanf(argv[2], "%f", &x_angle) ;
      x_angle = RADIANS(x_angle) ;
      nargs = 1 ;
      break ;
    case 'F':
      fit_flag = 1 ;
      break ;
    case 'Y':
      sscanf(argv[2], "%f", &y_angle) ;
      y_angle = RADIANS(y_angle) ;
      nargs = 1 ;
      break ;
    case 'Z':
      sscanf(argv[2], "%f", &z_angle) ;
      z_angle = RADIANS(z_angle) ;
      nargs = 1 ;
      break ;
    case 'S':
      ellipsoid_flag = 1 ;
      break ;
    case 'T':
      talairach_flag = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'N':
      normalize_flag = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
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
  fprintf(stderr,
          "usage: %s [options] <input image file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will display an MRI surface reconstruction.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


#define SHOW_AXES  0

#define COMPILE_SURFACE 1
static void
display_handler(void) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
#if COMPILE_SURFACE
  glCallList(current_list) ;      /* render sphere display list */
#else
  OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
#endif

#if SHOW_AXES
  glCallList(2) ;      /* show axes */
#endif
  glutSwapBuffers() ;  /* swap rendered buffer into display */
}

static void
keyboard_handler(unsigned char key, int x, int y) {
  int   redraw = 1 ;
  char  wname[100] ;

  glMatrixMode(GL_MODELVIEW);
  switch (key) {
  case 'E':
  case 'e':
    if (mrisp) {
      double mse ;
      fprintf(stderr, "computing correlation error...") ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      MRISrotate(mris, mris, RADIANS(total_alpha),
                 RADIANS(total_beta), RADIANS(total_gamma)) ;
      mse = MRIScomputeCorrelationError(mris, mrisp, param_no) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      fprintf(stderr, "done.\n") ;
      sprintf(title, "(%2.2f, %2.2f, %2.2f), mse: %2.3f (%2.3f)",
              total_alpha, total_beta, total_gamma, mse, starting_mse) ;
      glutSetWindowTitle(title) ;
    }
    break ;
  case 'W':  /* write it out */
    if (mris->hemisphere == RIGHT_HEMISPHERE)
      sprintf(wname, "rh.surf") ;
    else
      sprintf(wname, "lh.surf") ;
    fprintf(stderr, "writing surface to %s\n", wname) ;
    MRISwrite(mris, wname) ;
    break ;
  case '7':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(-10.0f, 10.0f, 10.0f) ;
    break ;
  case '9':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(10.0f, 10.0f, 0.0f) ;
    break ;
  case '1':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(-10.0f, -10.0f, 0.0f) ;
    break ;
  case '3':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(10.0f, -10.0f, 0.0f) ;
    break ;
  case '8':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(0.0f, 10.0f, 0.0f) ;
    break ;
  case '2':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(0.0f, -10.0f, 0.0f) ;
    break ;
  case '4':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(-10.0f, 0.0f, 0.0f) ;
    break ;
  case '6':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(10.0f, 0.0f, 0.0f) ;
    break ;
  case 'b':
  case 'B':   /* brighter */
    lighting_offset *= 1.5f ;
    OGLUsetLightingModel(-1.0f, -1.0f, -1.0f, -1.0f, lighting_offset) ;
    break ;
  case 'd':   /* dimmer */
  case 'D':
    lighting_offset /= 1.5f ;
    OGLUsetLightingModel(-1.0f, -1.0f, -1.0f, -1.0f, lighting_offset) ;
    break ;
  case '+':
    delta_angle *= 2.0f ;
    fprintf(stderr, "delta angle = %2.2f\n", delta_angle) ;
    redraw = 0 ;
    break ;
  case '/':
    cslope /= 1.5f ;
    fprintf(stderr, "cslope = %2.3f\n", cslope) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    break ;
  case '*':
    cslope *= 1.5f ;
    fprintf(stderr, "cslope = %2.3f\n", cslope) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    break ;
  case '-':
    delta_angle *= 0.5f ;
    fprintf(stderr, "delta angle = %2.2f\n", delta_angle) ;
    redraw = 0 ;
    break ;
  case GLUT_ESCAPE_KEY:
  case 'q':
    exit(0) ;
    break ;
  case 'h':
    home(mris) ;
    break ;
#if 0
  case 'x':   /* local x axis is reversed */
    glRotatef(-delta_angle, 1.0f, 0.0f, 0.0f) ;
    break ;
  case 'X':
    glRotatef(delta_angle, 1.0f, 0.0f, 0.0f) ;
    break ;
  case 'y':   /* local y is actually z coordinate */
    glRotatef(delta_angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case 'Y':
    glRotatef(-delta_angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case 'z':   /* local z is actually y axis */
    glRotatef(delta_angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case 'Z':
    glRotatef(-delta_angle, 0.0f, 1.0f, 0.0f) ;
    break ;
#endif
  case 'r':
    redraw = 1 ;
    break ;
  case 'T':
  case 't':
#if 1
    current_list = ORIG_SURFACE_LIST ;
    fprintf(stderr, "rotating brain by (%2.2f, %2.2f, %2.2f)\n",
            total_alpha, total_beta, total_gamma) ;
    MRISrotate(mris, mris, RADIANS(total_alpha),
               RADIANS(total_beta), RADIANS(total_gamma)) ;
    home(mris) ;   /* reset the view to original state */
    glNewList(current_list, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    compiled[current_list] = 1 ;
#else
    sprintf(wname, "talairached cortical surface reconstruction '%s'",
            surf_fname);
    MRIStalairachTransform(mris, mris) ;
    MRIScenter(mris, mris) ;
#if COMPILE_SURFACE
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
#endif
    glutSetWindowTitle(wname) ;
#endif
    break ;
  case 'O':
  case 'o':
    if (!mrisp)
      break ;
    switch (current_list) {
    default:
    case ORIG_SURFACE_LIST:
      current_list = MRISP_LIST ;
      break ;
    case MRISP_LIST:
      current_list = ORIG_SURFACE_LIST ;
      break ;
    }
    switch (current_list) {
    default:
    case ORIG_SURFACE_LIST:
      if (!compiled[current_list]) {
        MRISuseMeanCurvature(mris) ;
        glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
        OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
        glEndList() ;
      }
      break ;
    case MRISP_LIST:
      if (!compiled[current_list]) {
        MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
        if (normalize_param)
          MRISnormalizeFromParameterization(mrisp, mris, param_no) ;
        else
          MRISfromParameterization(mrisp, mris, 0) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        glNewList(current_list, GL_COMPILE) ;
        OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
        glEndList() ;
      }
      break ;
    }
    compiled[current_list] = 1 ;
    break ;
  case 'm':
  case 'M':
    if (compile_flags & MESH_FLAG)
      compile_flags &= ~MESH_FLAG ;
    else
      compile_flags |= MESH_FLAG ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    break ;
  case 'F':
  case 'f':
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    MRISfillSurface(mris) ;
    glEndList() ;
    break ;
  case 'i':
  case 'I':
    MRISinvertSurface(mris) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    break ;
  default:
    fprintf(stderr, "unknown normal key=%d\n", key) ;
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
}

#if 0
static void
rotate_brain(float angle, char c) {
  int i,j,k;
  float m[4][4], m1[4][4], m2[4][4]; /* Matrix m,m1,m2; */
  float sa,ca;

  glGetFloatv(GL_MODELVIEW_MATRIX,(float *)(m)) ;

  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  }
  angle = angle*PI/180;
  sa = sin(angle);
  ca = cos(angle);
  if (c=='y') {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } else if (c=='x') {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  } else if (c=='z') {
    m1[0][0] = m1[1][1] = ca;
    m1[1][0] = -(m1[0][1] = sa);
  } else {
    printf("surfer: ### Illegal axis %c\n",c);
    return;
  }
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += m[i][k]*m1[k][j];
    }
  }

  glLoadMatrixf((float *)m) ;
}
#endif
static void
reshape_handler(int width, int height) {
  int new_dim ;

  new_dim = MIN(width, height) ;
  glViewport(0,0,frame_xdim = new_dim, frame_ydim = new_dim);
}

static void
mouse_handler(int button, int state, int x, int y) {
  int kb_state, redraw = 1, width, height ;

  kb_state = glutGetModifiers() ;
  switch (button) {
  case GLUT_LEFT_BUTTON:
    width = glutGet(GLUT_WINDOW_WIDTH) ;
    height = glutGet(GLUT_WINDOW_HEIGHT) ;
    x -= width/2 ;
    y -= height/2 ;  /* put origin to center of window */
    if (kb_state & GLUT_ACTIVE_CTRL && state == GLUT_UP)/* zoom around click */
    {
      glMatrixMode(GL_PROJECTION) ;
      glTranslatef(-SCALE_FACTOR*(float)x, SCALE_FACTOR*(float)y, 0.0f) ;
      glScalef(2.0f, 2.0f, 2.0f) ;
      /*      glTranslatef(0.0f, SCALE_FACTOR*y, SCALE_FACTOR*x) ;*/
    } else {
      redraw = 0 ;
      fprintf(stderr, "(%d, %d)      \r", x, y) ;
    }
    break ;
  case GLUT_RIGHT_BUTTON:
    if (kb_state & GLUT_ACTIVE_CTRL && state == GLUT_UP)
      home(mris) ;
    break ;
  case GLUT_MIDDLE_BUTTON:
  default:
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
}

static void
special_key_handler(int key, int x, int y) {
  int   redraw = 1, kb_state ;
  float angle ;

  kb_state = glutGetModifiers() ;

  glMatrixMode(GL_MODELVIEW);
  if (kb_state & GLUT_ACTIVE_SHIFT)
    angle = 180.0f ;
  else
    angle = delta_angle ;

  switch (key) {
  case GLUT_KEY_F1:
    break ;
  case GLUT_KEY_LEFT:
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
    total_alpha = total_alpha - angle ;
    break ;
  case GLUT_KEY_RIGHT:
    glRotatef(-angle, 0.0f, 1.0f, 0.0f) ;
    total_alpha = total_alpha + angle ;
    break ;
  case GLUT_KEY_DOWN:
    glRotatef(-angle, 1.0f, 0.0f, 0.0f) ;
    total_gamma = total_gamma + angle ;
    break ;
  case GLUT_KEY_UP:
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    total_gamma = total_gamma - angle ;
    break ;
#define GLUT_KEY_DELETE (GLUT_KEY_INSERT+1)
  case GLUT_KEY_DELETE:
  case GLUT_KEY_PAGE_UP:
    glRotatef(angle, 0.0f, 0.0f, 1.0f) ;
    total_beta = total_beta - angle ;
    break ;
  case GLUT_KEY_PAGE_DOWN:
    glRotatef(-angle, 0.0f, 0.0f, 1.0f) ;
    total_beta = total_beta + angle ;
    break ;
  case GLUT_KEY_INSERT:
    break ;
  case GLUT_KEY_HOME:
    break ;
#if 0
  case GLUT_KEY_ESCAPE:
    exit(0) ;
    break ;
#endif
  default:
    fprintf(stderr, "unknown special key=%d\n", key) ;
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
  if (!mrisp) {
    sprintf(title, "%s (%2.2f, %2.2f, %2.2f)",
            wname, total_alpha, total_beta, total_gamma) ;
    glutSetWindowTitle(title) ;
  } else {
    sprintf(title, "%s (%2.2f, %2.2f, %2.2f)",
            wname, total_alpha, total_beta, total_gamma) ;
    glutSetWindowTitle(title) ;
  }
}


#if 0
#define NUM_LIGHTS 4
static void
controlLights(int value) {
  static int light_enabled[NUM_LIGHTS] = {
                                           1, 1, 1, 1
                                         } ;

  light_enabled[value] = !light_enabled[value] ;

  switch (value) {
  case 1:
    if (light_enabled[value])
      glEnable(GL_LIGHT0);
    else
      glDisable(GL_LIGHT0) ;
    break ;
  case 2:
    if (light_enabled[value])
      glEnable(GL_LIGHT1);
    else
      glDisable(GL_LIGHT1) ;
    break ;
  case 3:
    if (light_enabled[value])
      glEnable(GL_LIGHT2) ;
    else
      glDisable(GL_LIGHT2) ;
    break ;
  case 4:
    if (light_enabled[value])
      glEnable(GL_LIGHT3);
    else
      glDisable(GL_LIGHT3) ;
    break ;
  default:
    break ;
  }

  glutPostRedisplay() ;
}
#endif

static void
home(MRI_SURFACE *mris) {
  float angle ;

  total_alpha = total_beta = total_gamma = 0.0f ;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (mris->patch) {
    angle = 90.0f ;
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  } else {
    if (mris->hemisphere == RIGHT_HEMISPHERE)
      angle = RIGHT_HEMISPHERE_ANGLE ;
    else
      angle = LEFT_HEMISPHERE_ANGLE ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  }
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-oglu_fov, oglu_fov, -oglu_fov, oglu_fov,
          -10.0f*oglu_fov, 10.0f*oglu_fov);
}
static void
findAreaExtremes(MRI_SURFACE *mris) {
  int     fno, fmax, fmin ;
  FACE    *face ;
  float   min_area, max_area;

  fmax = fmin = 0 ;
  min_area = 10000.0f ;
  max_area = -10000.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->area > max_area) {
      max_area = face->area ;
      fmax = fno ;
      if (face->area < min_area) {
        min_area = face->area ;
        fmin = fno ;
      }
    }
  }
  fprintf(stderr, "min_area = %2.3f at f %d\n", min_area, fmin) ;
  fprintf(stderr, "max_area = %2.3f at f %d\n", max_area, fmax) ;
}

int
MRISinvertSurface(MRI_SURFACE *mris) {
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->nx *= -1 ;
    v->ny *= -1 ;
    v->nz *= -1 ;
  }
  return(NO_ERROR) ;
}
int
MRIScheckSphericalOrientation(MRI_SURFACE *mris) {
  int     vno ;
  VERTEX  *v ;
  float   dot ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dot = v->nx * v->x + v->ny * v->y + v->nz * v->z ;
    if (dot < 0) {
      fprintf(stderr, "vertex %d: x = (%2.1f,%2.1f,%2.1f), "
              "n = (%2.1f,%2.1f,%2.1f), dot = %2.1f.\n",
              vno, v->x, v->y, v->z, v->nx, v->ny, v->nz, dot) ;
    }
  }
  return(NO_ERROR) ;
}


static int mrisFillFace(MRI_SURFACE *mris, int fno) ;
static void load_brain_coords(float x,float y, float z, float *v) ;
int
MRISfillSurface(MRI_SURFACE *mris) {
  int    fno ;
  float  vn[3], v0[3], v1[3], v2[3] ;
  VERTEX *v_0, *v_1, *v_2 ;
  FACE   *f ;

  /* draw surface in red */
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    f = &mris->faces[fno] ;
    v_0 = &mris->vertices[f->v[0]] ;
    v_1 = &mris->vertices[f->v[1]] ;
    v_2 = &mris->vertices[f->v[2]] ;
    load_brain_coords(f->nx,f->ny,f->nz,vn);
    load_brain_coords(v_0->x,v_0->y,v_0->z,v0);
    load_brain_coords(v_1->x,v_1->y,v_1->z,v1);
    load_brain_coords(v_2->x,v_2->y,v_2->z,v2);
    glBegin(GL_LINES) ;
    glColor3ub(0,255,0) ;
    glNormal3fv(vn);
    glVertex3fv(v0);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v0);
    glEnd() ;
  }

  for (fno = 0 ; fno < 1 /*mris->nfaces*/ ; fno++)
    mrisFillFace(mris, fno) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        Description
        each face has 2 triangles defined by it:

       V0    b     V2
        o----------o
        |        /
        |      /
      a |    /
        |  /
        |/
        o
       V1      b        V2
------------------------------------------------------*/
#if 0
#define SAMPLE_DIST   0.25
#else
#define SAMPLE_DIST   25
#endif

static int
mrisFillFace(MRI_SURFACE *mris, int fno) {
  Real   x, y, z, xa, ya, za, xc, yc, zc, t0, t1, adx, ady, adz, dx, dy, dz,
  cdx, cdy, cdz, alen, clen, delta_t0, delta_t1, len ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *face ;
  float  vn[3], vstart[3], vend[3] ;
  int    i ;

  face = &mris->faces[fno] ;
  if (face->ripflag)
    return(NO_ERROR) ;

  for (i = 0 ; i < 3 ; i++) {
    switch (i) {
    default:
    case 0:
      v0 = &mris->vertices[face->v[0]] ;
      v1 = &mris->vertices[face->v[1]] ;
      v2 = &mris->vertices[face->v[2]] ;
      break ;
    case 1:
      v0 = &mris->vertices[face->v[1]] ;
      v1 = &mris->vertices[face->v[2]] ;
      v2 = &mris->vertices[face->v[0]] ;
      break ;
    case 2:
      v0 = &mris->vertices[face->v[2]] ;
      v1 = &mris->vertices[face->v[0]] ;
      v2 = &mris->vertices[face->v[1]] ;
      break ;
    }

    adx = v1->x - v0->x ;
    ady = v1->y - v0->y ;
    adz = v1->z - v0->z ;
    alen = sqrt(SQR(adx)+SQR(ady)+SQR(adz)) ;
    cdx = v2->x - v0->x ;
    cdy = v2->y - v0->y ;
    cdz = v2->z - v0->z ;
    clen = sqrt(SQR(cdx)+SQR(cdy)+SQR(cdz)) ;

    /*
      sample along legs of the triangle making sure the maximum spacing
      between samples (along the longer leg) is SAMPLE_DIST.
    */

    /*move along v0->v1 and v3->v2 lines and draw in crossing line to fill face*/
    /* t0 parameterizes lines from v0->v1 and v0->v2 */
    if (FZERO(alen) && FZERO(clen))
      delta_t0 = 0.99 ;
    else
      delta_t0 = (alen > clen) ? (SAMPLE_DIST / alen) : (SAMPLE_DIST / clen ) ;
    if (FZERO(delta_t0))
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "mrisFillFace: face %d has infinite leg (%d, %d)\n",
                   fno, alen, clen)) ;

    if (delta_t0 >= 1.0)
      delta_t0 = 0.99 ;

    /* delta_t0 is % of alen or clen (whichever is bigger) of SAMPLE_DIST */
    for (t0 = 0 ; t0 <= 1.0f ; t0 += delta_t0) {
      /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
      xa = v0->x + t0*adx ;
      ya = v0->y + t0*ady ;
      za = v0->z + t0*adz ;
      xc = v0->x + t0*cdx ;
      yc = v0->y + t0*cdy ;
      zc = v0->z + t0*cdz ;
      dx = xc-xa ;
      dy = yc-ya ;
      dz = zc-za ;
      len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
      if (FZERO(len))
        delta_t1 = 0.99 ;
      else {
        delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
        if (delta_t1 >= 1.0f)
          delta_t1 = 0.99 ;
      }

      /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
      load_brain_coords(face->nx,face->ny,face->nz,vn);
      load_brain_coords(xa, ya, za, vstart);
      for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1) {
        /* compute a point on the line connecting a and c */
        x = xa + t1*dx ;
        y = ya + t1*dy ;
        z = za + t1*dz ;
#if 1
        glBegin(GL_POINTS) ;
        glColor3ub(255,255,255) ;
        glNormal3fv(vn);
        load_brain_coords(x, y, z, vend);
        glVertex3fv(vstart);
        glVertex3fv(vend);
        glEnd() ;
        memmove(vstart, vend, 3*sizeof(float)) ;
#else
        // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
        MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
        xv = nint(x) ;
        yv = nint(y) ;
        zv = nint(z) ;  /* voxel coordinate */
        MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
#endif
      }
      glBegin(GL_POINTS) ;
      glColor3ub(255,255,255) ;
      glNormal3fv(vn);
      load_brain_coords(x, y, z, vstart);
      load_brain_coords(xc, yc, zc, vend);
      glVertex3fv(vstart);
      glVertex3fv(vend);
      glEnd() ;
    }
    t0 = 1.0f ;
    /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
    xa = v0->x + t0*adx ;
    ya = v0->y + t0*ady ;
    za = v0->z + t0*adz ;
    xc = v0->x + t0*cdx ;
    yc = v0->y + t0*cdy ;
    zc = v0->z + t0*cdz ;
    dx = xc-xa ;
    dy = yc-ya ;
    dz = zc-za ;
    len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
    if (FZERO(len))
      delta_t1 = 0.99 ;
    else {
      delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
      if (delta_t1 >= 1.0f)
        delta_t1 = 0.99 ;
    }

    /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
    load_brain_coords(face->nx,face->ny,face->nz,vn);
    load_brain_coords(xa, ya, za, vstart);
    for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1) {
      /* compute a point on the line connecting a and c */
      x = xa + t1*dx ;
      y = ya + t1*dy ;
      z = za + t1*dz ;
#if 1
      glBegin(GL_POINTS) ;
      glColor3ub(255,255,255) ;
      glNormal3fv(vn);
      load_brain_coords(x, y, z, vend);
      glVertex3fv(vstart);
      glVertex3fv(vend);
      glEnd() ;
      memmove(vstart, vend, 3*sizeof(float)) ;
#else
      // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      xv = nint(x) ;
      yv = nint(y) ;
      zv = nint(z) ;  /* voxel coordinate */
      MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
#endif
    }
    glBegin(GL_POINTS) ;
    glColor3ub(255,255,255) ;
    glNormal3fv(vn);
    load_brain_coords(x, y, z, vstart);
    load_brain_coords(xc, yc, zc, vend);
    glVertex3fv(vstart);
    glVertex3fv(vend);
    glEnd() ;
  }

  return(NO_ERROR) ;
}
static void
load_brain_coords(float x,float y, float z, float *v) {
  v[0] = -x;
  v[1] = z;
  v[2] = y;
}
