#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <glut.h>
#include <gl.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"

static char vcid[] = "$Id: mris_show.c,v 1.9 1997/09/25 22:33:12 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static void load_brain_coords(float x,float y, float z, float v[]) ;
static int  MRISGLcompile(MRI_SURFACE *mris) ;
/* static void rotate_brain(float angle, char c) ;*/

static void keyboard_handler(unsigned char key, int x, int y) ;
static void special_key_handler(int key, int x, int y) ;
static void reshape_handler(int width, int height) ;
static void mouse_handler(int button, int state, int x, int y) ;
static void controlLights(int value) ;

static void display_handler(void) ;
static void home(MRI_SURFACE *mris) ;
static void gfxinit(MRI_SURFACE *mris) ;
static void set_lighting_model(float lite0, float lite1, float lite2, 
             float lite3, float newoffset) ;

char *Progname ;
char *surf_fname ;
static MRI_SURFACE  *mris, *mris_orig, *mris_ellipsoid = NULL ;
static int marked_vertex = -1 ;
static long frame_xdim = 600;
static long frame_ydim = 600;

static int show_temporal_region_flag = 0 ;
static int talairach_flag = 0 ;
static int ellipsoid_flag = 0 ;

static float x_angle = 0.0f ;
static float y_angle = 0.0f ;
static float z_angle = 0.0f ;

#define DELTA_ANGLE  18.0f
static float delta_angle = 2*DELTA_ANGLE ;

#define RIGHT_HEMISPHERE_ANGLE   90.0
#define LEFT_HEMISPHERE_ANGLE   -90.0

#define ORIG_SURFACE_LIST   1
#define ELLIPSOID_LIST      1

static int current_list = ORIG_SURFACE_LIST ;

#define SCALE_FACTOR   0.55f
#define FOV            (256.0f*SCALE_FACTOR)

static INTEGRATION_PARMS  parms ;

#ifndef GLUT_ESCAPE_KEY
#define GLUT_ESCAPE_KEY  27
#endif

static void findAreaExtremes(MRI_SURFACE *mris) ;
static char curvature_fname[100] = "" ;


static float cslope = 5.0f ;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, wname[200] ;
  int          ac, nargs ;
  float        angle ;

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
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit() ;

  surf_fname = in_fname = argv[1] ;
  out_fname = argv[2] ;
  mris_orig = mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  if (curvature_fname[0])
    MRISreadCurvatureFile(mris, curvature_fname) ;
  if (ellipsoid_flag)
  {
    MRISupdateEllipsoidSurface(mris) ;
    findAreaExtremes(mris) ;
  }
  if (talairach_flag)
    MRIStalairachTransform(mris_orig, mris_orig) ;
  MRIScenter(mris_orig, mris_orig) ;
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
  glutCreateWindow(wname) ;          /* open the window */
  glutHideWindow() ;
  glutDisplayFunc(display_handler) ; /* specify a drawing callback function */
  gfxinit(mris) ;                    /* specify lighting and such */

  if (mris->hemisphere == RIGHT_HEMISPHERE)
    angle = RIGHT_HEMISPHERE_ANGLE ;
  else
    angle = LEFT_HEMISPHERE_ANGLE ;
  glMatrixMode(GL_MODELVIEW);
  glRotatef(angle, 0.0f, 1.0f, 0.0f) ;

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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "cslope"))
  {
    sscanf(argv[2], "%f", &cslope) ;
    nargs = 1 ;
    fprintf(stderr, "using color slope compression = %2.4f\n", cslope) ;
  }
  else switch (toupper(*option))
  {
  case 'C':
    strcpy(curvature_fname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'V':
    sscanf(argv[2], "%d", &marked_vertex) ;
    nargs = 1 ;
    break ;
  case 'X':
    sscanf(argv[2], "%f", &x_angle) ;
    x_angle = RADIANS(x_angle) ;
    nargs = 1 ;
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
  case 'E':
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
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input image file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program will display an MRI surface reconstruction.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int meshr = 100;   /* mesh color */
static int meshg = 100;
static int meshb = 100;

float pup = 0.07;  /* dist point above surface */

#define POLE_DISTANCE   1.0f   /* color everything within 1 mm of pole */

int
MRISGLcompile(MRI_SURFACE *mris)
{
  int          k,n, red, green, blue, error ;
  face_type    *f;
  vertex_type  *v ;
  float        v1[3], xf, yf, zf, xt, yt, zt, xo, yo, zo, xd, yd, zd,td,fd,od,
               min_curv, max_curv, offset ;
  Real         x, y, z, xtal,ytal,ztal ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "compiling surface tesselation...") ;

  if (mris->v_frontal_pole)
  {
    xf = mris->v_frontal_pole->x ;
    yf = mris->v_frontal_pole->y ;
    zf = mris->v_frontal_pole->z ;
  }
  if (mris->v_temporal_pole)
  {
    xt = mris->v_temporal_pole->x ;
    yt = mris->v_temporal_pole->y ;
    zt = mris->v_temporal_pole->z ;
  }
  else
    xt = yt = zt = 0.0f ;

  if (mris->v_occipital_pole)
  {
    xo = mris->v_occipital_pole->x ;
    yo = mris->v_occipital_pole->y ;
    zo = mris->v_occipital_pole->z ;
  }
  min_curv = mris->min_curv ;
  max_curv = mris->max_curv ;
  if (-min_curv > max_curv)
    max_curv = -min_curv ;

  ytal = 10000.0f ;  /* disable temporal region stuff */
  for (k=0;k<mris->nfaces;k++) if (!mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    glBegin(GL_QUADS) ;
    for (n=0;n<4;n++)
    {
      v = &mris->vertices[f->v[n]];

#if 0
      xd = v->x - xf ; yd = v->y - yf ; zd = v->z - zf ;
      fd = xd*xd + yd*yd + zd*zd ;  /* distance to frontal pole */
      if (mris->v_temporal_pole)
      {
        xd = v->x - xt ; yd = v->y - yt ; zd = v->z - zt ;
        td = xd*xd + yd*yd + zd*zd ;  /* distance to temporal pole */
      }
      else
        td = 1000.0f ;
      xd = v->x - xo ; yd = v->y - yo ; zd = v->z - zo ;
      od = xd*xd + yd*yd + zd*zd ;  /* distance to occipital pole */

      if (show_temporal_region_flag)
      {
        x = (Real)v->x ; y = (Real)v->y ; z = (Real)v->z ;
        transform_point(mris->linear_transform, -x, z, y, &xtal, &ytal,&ztal);
      }
#endif
      if (k == marked_vertex)
        glColor3ub(0,0,240) ;      /* paint the marked vertex blue */
#if 0
      else if (ytal < MAX_TALAIRACH_Y)
        glColor3ub(0,0,meshb) ;      /* paint the temporal pole blue */
      else if (td < POLE_DISTANCE)
        glColor3ub(0,0,255) ;          /* paint the temporal pole blue */
      else if (fd < POLE_DISTANCE)
        glColor3ub(255,255,0);         /* paint the frontal pole yellow */
      else if (od < POLE_DISTANCE)
        glColor3ub(255,255,255);       /* paint the occipital pole white */
#endif
      else   /* color it depending on curvature */
      {
#define MIN_GRAY  50
#define MAX_COLOR ((float)(255 - MIN_GRAY))

        red = green = blue = MIN_GRAY ;
        if (FZERO(max_curv))  /* no curvature info */
          red = green = blue = meshr ;    /* paint it gray */

/* curvatures seem to be sign-inverted??? */
        if (v->curv > 0)  /* color it red */
        {
          offset = MAX_COLOR*tanh(cslope * (v->curv/max_curv)) ;
          red = MIN_GRAY + nint(offset);
        }
        else              /* color it green */
        {
          offset = MAX_COLOR*tanh(cslope*(-v->curv/max_curv)) ;
          green = MIN_GRAY + nint(offset);
        }
        glColor3ub(red,green,blue);  /* specify the RGB color */
      }
      load_brain_coords(v->nx,v->ny,v->nz,v1);
      glNormal3fv(v1);                /* specify the normal for lighting */
      load_brain_coords(v->x,v->y,v->z,v1);
      glVertex3fv(v1);                /* specify the position of the vertex*/
    }
    glEnd() ;   /* done specifying this polygon */
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  error = glGetError() ;
  if (error != GL_NO_ERROR)
  {
    const char *errstr ;

    errstr = gluErrorString(error) ;
    fprintf(stderr, "GL error %d: %s\n", error, errstr) ;
    exit(1) ;
  }
  return(NO_ERROR) ;
}

static void
load_brain_coords(float x,float y, float z, float v[])
{
  v[0] = -x;
  v[1] = z;
  v[2] = y;
}

#define SHOW_AXES  0

#define COMPILE_SURFACE 1
static void
display_handler(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
#if COMPILE_SURFACE
  glCallList(current_list) ;      /* render sphere display list */
#else
  MRISGLcompile(mris) ;
#endif

#if SHOW_AXES
  glCallList(2) ;      /* show axes */
#endif
  glutSwapBuffers() ;  /* swap rendered buffer into display */
}

static void
gfxinit(MRI_SURFACE *mris)
{
  int   error ;

  glViewport(0,0,frame_xdim,frame_ydim);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);

  glOrtho(-FOV, FOV, -FOV, FOV, -10.0f*FOV, 10.0f*FOV);

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  set_lighting_model(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f) ;

  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);

  /* now specify where the model is in space */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

#if 0
  /* MRI coords (vs. viewer:translate_brain_) */
  glTranslatef(mris->xctr, -mris->zctr,-mris->yctr); 
#endif

#if COMPILE_SURFACE
  glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
  MRISGLcompile(mris) ;
  glEndList() ;
#endif

  error = glGetError() ;
  if (error != GL_NO_ERROR)
  {
    const char *errstr ;

    errstr = gluErrorString(error) ;
    fprintf(stderr, "GL error %d: %s\n", error, errstr) ;
    exit(1) ;
  }
}

static void
keyboard_handler(unsigned char key, int x, int y)
{
  int   redraw = 1, niter = 1 ;
  char  wname[100] ;
  float angle ;

  glMatrixMode(GL_MODELVIEW);
  switch (key)
  {
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
  case '+':
    delta_angle *= 2.0f ;
    redraw = 0 ;
    break ;
  case '/':
    cslope /= 1.5f ;
    fprintf(stderr, "cslope = %2.3f\n", cslope) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    MRISGLcompile(mris) ;
    glEndList() ;
    break ;
  case '*':
    cslope *= 1.5f ;
    fprintf(stderr, "cslope = %2.3f\n", cslope) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    MRISGLcompile(mris) ;
    glEndList() ;
    break ;
  case '-':
    delta_angle *= 0.5f ;
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
  case 'D':
  case 'd':
  case 'r':
    redraw = 1 ;
    break ;
  case 'T':
  case 't':
    sprintf(wname, "talairached cortical surface reconstruction '%s'", 
            surf_fname);
    MRIStalairachTransform(mris_orig, mris_orig) ;
    MRIScenter(mris_orig, mris_orig) ;
#if COMPILE_SURFACE
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    MRISGLcompile(mris_orig) ;
    glEndList() ;
#endif
    glutSetWindowTitle(wname) ;
    break ;
  case 'P':
  case 'p':
    mris_ellipsoid = MRISprojectOntoEllipsoid(mris, NULL, 0.0f, 0.0f, 0.0f);
    MRISfree(&mris_orig) ;
    mris = mris_ellipsoid ;
#if COMPILE_SURFACE
    glDeleteLists(ORIG_SURFACE_LIST, 1) ;
    glNewList(ELLIPSOID_LIST, GL_COMPILE) ;
    MRISGLcompile(mris) ;
    glEndList() ;
#endif
    current_list = ELLIPSOID_LIST ;
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    if (mris_ellipsoid->hemisphere == RIGHT_HEMISPHERE)
      angle = RIGHT_HEMISPHERE_ANGLE ;
    else
      angle = LEFT_HEMISPHERE_ANGLE ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case 'u':
    niter = parms.niterations ;
  case 'U':
    MRISunfold(mris, &parms) ;

#if COMPILE_SURFACE
    glDeleteLists(ELLIPSOID_LIST, 1) ;
    glNewList(ELLIPSOID_LIST, GL_COMPILE) ;
    MRISGLcompile(mris) ;
    glEndList() ;
#endif
    break ;
  default:
    fprintf(stderr, "unknown normal key=%d\n", key) ;
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
}

#define LIGHT0_BR  0.4 /* was 0.2 */
#define LIGHT1_BR  0.0 
#define LIGHT2_BR  0.6 /* was 0.3 */
#define LIGHT3_BR  0.2 /* was 0.1 */
#define OFFSET 0.25   /* was 0.15 */
static void
set_lighting_model(float lite0, float lite1, float lite2, float lite3, 
                   float newoffset)
{
  float offset = OFFSET ;
  static GLfloat mat0_ambient[] =  { 0.0, 0.0, 0.0, 1.0 };
  static GLfloat mat0_diffuse[] =  { OFFSET, OFFSET, OFFSET, 1.0 };
  static GLfloat mat0_emission[] = { 0.0, 0.0, 0.0, 1.0 };
  static GLfloat light0_diffuse[] = { LIGHT0_BR, LIGHT0_BR, LIGHT0_BR, 1.0 };
  static GLfloat light1_diffuse[] = { LIGHT1_BR, LIGHT1_BR, LIGHT1_BR, 1.0 };
  static GLfloat light2_diffuse[] = { LIGHT2_BR, LIGHT2_BR, LIGHT2_BR, 1.0 };
  static GLfloat light3_diffuse[] = { LIGHT3_BR, LIGHT3_BR, LIGHT3_BR, 1.0 };
  static GLfloat light0_position[] = { 0.0, 0.0, 1.0, 0.0 };
  static GLfloat light1_position[] = { 0.0, 0.0,-1.0, 0.0 };
  static GLfloat light2_position[] = { 0.6, 0.6, 1.6, 0.0 };
  static GLfloat light3_position[] = {-1.0, 0.0, 0.0, 0.0 };
  static GLfloat lmodel_ambient[] =  { 0.0, 0.0, 0.0, 0.0 };

  if (lite0 < 0.0)     lite0 = light0_diffuse[0];
  if (lite1 < 0.0)     lite1 = light1_diffuse[0];
  if (lite2 < 0.0)     lite2 = light2_diffuse[0];
  if (lite3 < 0.0)     lite3 = light3_diffuse[0];
  newoffset = offset;

  /* material: change DIFFUSE,EMISSION (purpler: EMISSION=0.05*newoffset) */
  mat0_diffuse[0] = mat0_diffuse[1] = mat0_diffuse[2] = newoffset;
  mat0_emission[0] = mat0_emission[1] = mat0_emission[2] = 0.1*newoffset;

  /* lights: change DIFFUSE */
  light0_diffuse[0] = light0_diffuse[1] = light0_diffuse[2] = lite0;
  light1_diffuse[0] = light1_diffuse[1] = light1_diffuse[2] = lite1;
  light2_diffuse[0] = light2_diffuse[1] = light2_diffuse[2] = lite2;
  light3_diffuse[0] = light3_diffuse[1] = light3_diffuse[2] = lite3;

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glLoadIdentity();

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat0_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat0_diffuse);
    glMaterialfv(GL_FRONT, GL_EMISSION, mat0_emission);

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
    glLightfv(GL_LIGHT3, GL_DIFFUSE, light3_diffuse);

    /* might need to move */
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
    glLightfv(GL_LIGHT3, GL_POSITION, light3_position);

    /* turn off default global 0.2 ambient (inf. viewpnt, not 2-sided OK) */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHT3);
  glPopMatrix();
}
#if 0
static void
rotate_brain(float angle, char c)
{
  int i,j,k;
  float m[4][4], m1[4][4], m2[4][4]; /* Matrix m,m1,m2; */
  float sa,ca;

  glGetFloatv(GL_MODELVIEW_MATRIX,(float *)(m)) ;

  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  }
  angle = angle*PI/180;
  sa = sin(angle);
  ca = cos(angle);
  if (c=='y')
  {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } 
  else if (c=='x') 
  {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  } 
  else if (c=='z') 
  {
    m1[0][0] = m1[1][1] = ca;
    m1[1][0] = -(m1[0][1] = sa);
  } 
  else 
  {
    printf("surfer: ### Illegal axis %c\n",c);return;
  }
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++) 
    {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += m[i][k]*m1[k][j];
    }
  }

  glLoadMatrixf((float *)m) ;
}
#endif
static void
reshape_handler(int width, int height)
{
  int new_dim ;

  new_dim = MIN(width, height) ;
  glViewport(0,0,frame_xdim = new_dim, frame_ydim = new_dim);
}

static void
mouse_handler(int button, int state, int x, int y)
{
  int kb_state, redraw = 1, width, height ;

  kb_state = glutGetModifiers() ;
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    width = glutGet(GLUT_WINDOW_WIDTH) ;
    height = glutGet(GLUT_WINDOW_HEIGHT) ;
    x -= width/2 ; y -= height/2 ;  /* put origin to center of window */
    if (kb_state & GLUT_ACTIVE_CTRL && state == GLUT_UP)/* zoom around click */
    {
      glMatrixMode(GL_PROJECTION) ;
      glTranslatef(-SCALE_FACTOR*(float)x, SCALE_FACTOR*(float)y, 0.0f) ;
      glScalef(2.0f, 2.0f, 2.0f) ;
      /*      glTranslatef(0.0f, SCALE_FACTOR*y, SCALE_FACTOR*x) ;*/
    }
    else
    {
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
special_key_handler(int key, int x, int y)
{
  int   redraw = 1, kb_state ;
  float angle ;

  kb_state = glutGetModifiers() ;

  glMatrixMode(GL_MODELVIEW);
  if (kb_state & GLUT_ACTIVE_SHIFT)
    angle = 180.0f ;
  else
    angle = delta_angle ;

  switch (key)
  {
  case GLUT_KEY_F1:
    break ;
  case GLUT_KEY_LEFT:
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case GLUT_KEY_RIGHT:
    glRotatef(-angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case GLUT_KEY_DOWN:
    glRotatef(-angle, 1.0f, 0.0f, 0.0f) ;
    break ;
  case GLUT_KEY_UP:
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    break ;
#define GLUT_KEY_DELETE (GLUT_KEY_INSERT+1)
  case GLUT_KEY_DELETE:
  case GLUT_KEY_PAGE_UP:
    glRotatef(angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case GLUT_KEY_PAGE_DOWN:
    glRotatef(-angle, 0.0f, 0.0f, 1.0f) ;
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
}


#define NUM_LIGHTS 4
static void
controlLights(int value)
{
  static int light_enabled[NUM_LIGHTS] = {1, 1, 1, 1} ;

  light_enabled[value] = !light_enabled[value] ;

  switch (value)
  {
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

#if SHOW_AXES
  /* build axes */
  glNewList(2, GL_COMPILE) ;

  /* x axis */
  glBegin(GL_LINE) ;
  glColor3ub(meshr,0,0) ;
  load_brain_coords(mris->xlo-30,mris->yctr,mris->zctr,v1);
  glVertex3fv(v1);               
  load_brain_coords(mris->xhi+30,mris->yctr,mris->zctr,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  glEnd() ;

  /* y axis */
  glBegin(GL_LINE) ;
  load_brain_coords(mris->xctr,mris->ylo-30,mris->zctr,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  load_brain_coords(mris->xctr,mris->yhi+30,mris->zctr,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  glEnd() ;

  /* z axis */
  glBegin(GL_LINE) ;
  load_brain_coords(mris->xctr,mris->yctr,mris->zlo-30,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  load_brain_coords(mris->xctr,mris->yctr,mris->zhi+30,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  glEnd() ;

  glEndList() ;
#endif

static void
home(MRI_SURFACE *mris)
{
  float angle ;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (mris->hemisphere == RIGHT_HEMISPHERE)
    angle = RIGHT_HEMISPHERE_ANGLE ;
  else
    angle = LEFT_HEMISPHERE_ANGLE ;
  glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-FOV, FOV, -FOV, FOV, -10.0f*FOV, 10.0f*FOV);
}
static void
findAreaExtremes(MRI_SURFACE *mris)
{
  int     fno, tno, fmax, tmax, fmin, tmin ;
  FACE    *face ;
  float   min_area, max_area;

  min_area = 10000.0f ; max_area = -10000.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      if (face->area[tno] > max_area)
      {
        max_area = face->area[tno] ;
        fmax = fno ;
        tmax = tno ;
      }
      if (face->area[tno] < min_area)
      {
        min_area = face->area[tno] ;
        fmin = fno ;
        tmin = tno ;
      }
    }
  }
  fprintf(stderr, "min_area = %2.3f at f %d, t %d\n", min_area, fmin, tmin) ;
  fprintf(stderr, "max_area = %2.3f at f %d, t %d\n", max_area, fmax, tmax) ;
}

#if 0
static void
set_color(float val, float curv, int mode)
{
  short r,g,b;
  float f,fr,fg,fb,tmpoffset;

  if (curv<0)  tmpoffset = cvfact*offset;
  else         tmpoffset = offset;

  if (mode==GREEN_RED_CURV)
  {
    f = tanh(cslope*(curv-cmid));
    if (f>0) {
      r = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
      g = 255 * (offset/blufact*(1 - fabs(f)));
    }
    else {
      r = 255 * (offset/blufact*(1 - fabs(f)));
      g = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
    }
    b = 255 * (offset*blufact*(1 - fabs(f)));
  }

  if (mode==GREEN_RED_VAL)   /* single val positive or signed */
  {
    if (colscale==HEAT_SCALE)
    {
      set_stat_color(val,&fr,&fg,&fb,tmpoffset);
      r=fr; g=fg; b=fb;
    }
    else
    if (colscale==CYAN_TO_RED || colscale==BLU_GRE_RED || colscale==JUST_GRAY)
    {
      if (val<fthresh) 
      {
        r = g = 255 * (tmpoffset/blufact);
        b =     255 * (tmpoffset*blufact);
      }
      else 
      {
        if (fslope!=0)
          f = (tanh(fslope*fmid)+tanh(fslope*(val-fmid)))/(2-tanh(fslope*fmid));
        else
          f = (val<0)?0:((val>1)?1:val);
        set_positive_color(f,&fr,&fg,&fb,tmpoffset);
        r=fr; g=fg; b=fb;
      }
    }
    else
    {
      if (fabs(val)<fthresh) 
      {
        r = g = 255 * (tmpoffset/blufact);
        b =     255 * (tmpoffset*blufact);
      }
      else 
      {
        if (fslope!=0)
        {
          if (fmid==0)
            f = tanh(fslope*(val));
          else
          {
            if (val<0)
              f = -(tanh(fslope*fmid) + tanh(fslope*(-val-fmid)))/
                   (2-tanh(fslope*fmid));
            else
              f = (tanh(fslope*fmid) + tanh(fslope*( val-fmid)))/
                  (2-tanh(fslope*fmid));
          }
        }
        else
          f = (val<-1)?-1:((val>1)?1:val);
        if (revphaseflag)
          f = -f;
        set_signed_color(f,&fr,&fg,&fb,tmpoffset);
        r=fr; g=fg; b=fb;
      }
    }
  }

  if (mode==FIELDSIGN_POS || mode==FIELDSIGN_NEG) {
    if (val<fthresh) {
      r = g = 255 * (tmpoffset/blufact);
      b =     255 * (tmpoffset*blufact);
    }
    else {
      f = (1.0 + tanh(fslope*(val-fmid)))/2.0;
      if (mode==FIELDSIGN_POS) {
        b = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
        r = g = 255* (tmpoffset*(1 - fabs(f)));
      }
      else {
        b = 255 * (tmpoffset*(1 - fabs(f)));
        r = g = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
      }
    }
  }

  if (mode==BORDER)  /* AMD 5/27/95 */
  {
    r = 255;
    g = 255;
    b = 0;
  }

  if (mode==MARKED)
  {
    r = 255;
    g = 255;
    b = 255;
  }

  r = (r<0)?0:(r>255)?255:r;
  g = (g<0)?0:(g>255)?255:g;
  b = (b<0)?0:(b>255)?255:b;
  glColor3ub(r,g,b);
}
#endif
