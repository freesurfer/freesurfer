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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <strings.h>

#include "glut.h"
#include "error.h"
#include "const.h"
#include "macros.h"
#include "diag.h"
#include "version.h"

#define EPS 1.0e-6

/* Global variables */
const char *Progname ;
static double flickerFreq;
static double stimPeriod;
static double minEcc;
static double maxEcc;
static int numRings, numSpokes, OUTWARD;
static int counter = 0;
static int currentFrame = 0 ;

static double timeTick;
static double timePoint = 0;

static GLfloat On = 0.0;
static GLint width = 500, height = 500;

static int numFrames = 1024 ;

static float ***redFrames = NULL ;
static float ***greenFrames = NULL ;
static float ***blueFrames = NULL ;

/* m-sequence stuff */
static int use_mseq = 0 ;
static int mseq_order = 16 ;
static int mseq_tap = 45 ;   /* could also be 57, 63, 83, 189, 215, 303, 317, 335, 349 and others */

/* for order=15 could use 3, 17, 23, 45, 53, 95, 119, 129, 135, 147, 165, 195, 207, 221 etc... */
/* for order=17 could use 9, 15, 33, 45, 51, 63, 65, 85, 105, 123, 141, 153, 163, 175, 187, 197, 245 */


static int get_option(int argc, char *argv[]) ;
static unsigned char *generate_msequence(int order, unsigned int tap, int *len) ;

void display_expanding_rings(void) {

  GLUquadricObj *qobj = NULL;
  double innerRadiusOfRing;
  double outerRadiusOfRing;
  double start_angle;
  double angular_step;
  static int current_list = 0 ;

  int ring, spoke, bit;
  double scaleFactor;

  double t, T, eccRange;
  double innerRadiusOfAnnulus;
  double outerRadiusOfAnnulus;

  putenv(const_cast<char*>("__GL_SYNC_TO_VBLANK=23")) ;  /* will force syncing with vertical retrace */
  glClear(GL_COLOR_BUFFER_BIT);
  current_list = !current_list ;
  glNewList(current_list, GL_COMPILE) ;

  angular_step = 360.0/numSpokes;

  scaleFactor = pow( (maxEcc/minEcc),
                     (1.0/( (double) numRings) ) );

  innerRadiusOfRing = minEcc;
  outerRadiusOfRing = minEcc * scaleFactor;

  /* draw radial checkerboard */
  for (ring = 0; ring < numRings; ring++) {
    if ( ring % 2 )
      bit = 0;
    else
      bit = 1;

    for (spoke = 0; spoke < numSpokes; spoke++) {

      start_angle = spoke * angular_step;

      if ( (spoke + bit) % 2 )
        glColor3f(On, On, On);
      else
        glColor3f(1-On, 1-On, 1-On);

      qobj = gluNewQuadric();
      gluPartialDisk(qobj,
                     (GLdouble) innerRadiusOfRing,  /* inner radius */
                     (GLdouble) outerRadiusOfRing,  /* outer radius */
                     (GLint) angular_step, /* #vert. on disk perimeter */
                     1,   /* # vertices on disk ray (minus 1) */
                     (GLdouble) start_angle,  /* start angle for wedge */
                     (GLdouble) angular_step); /* wedge width in degrees */

    }

    innerRadiusOfRing *= scaleFactor;
    outerRadiusOfRing *= scaleFactor;

  }


  T = stimPeriod;
  t = timePoint - counter*T;
  if (t > T) {
    outerRadiusOfAnnulus = maxEcc;
    innerRadiusOfAnnulus = minEcc;
    counter++;
  }

  eccRange = log(maxEcc) - log(minEcc);

  if (OUTWARD == 1) { /*** --------- EXP. OUTWARD ---------- ***/
    if ( t < T/2 ) {
      outerRadiusOfAnnulus = maxEcc;
      innerRadiusOfAnnulus = (log(minEcc) + log(maxEcc))/2 +
                             eccRange * (1/T) * t;
      innerRadiusOfAnnulus = exp(innerRadiusOfAnnulus);

      /* Threshold, and add 1 (pixel?) so that mask overcovers stim */
      if ( outerRadiusOfAnnulus >= maxEcc)
        outerRadiusOfAnnulus = maxEcc + 1;

      /* Threshold, and subtract 1 so that mask overcovers stim */
      if ( innerRadiusOfAnnulus <= minEcc)
        innerRadiusOfAnnulus = minEcc - 1;

      /*** GL calls to draw annulus ***/
      glColor3f(0.5, 0.5, 0.5);
      gluDisk(qobj,
              (GLdouble) innerRadiusOfAnnulus,
              (GLdouble) outerRadiusOfAnnulus,
              360,
              1);

      outerRadiusOfAnnulus = log(minEcc) +
                             eccRange * (1/T) * t;
      outerRadiusOfAnnulus = exp(outerRadiusOfAnnulus);
      innerRadiusOfAnnulus = minEcc;

      /* Threshold, and add 1 (pixel?) so that mask overcovers stim */
      if ( outerRadiusOfAnnulus >= maxEcc)
        outerRadiusOfAnnulus = maxEcc + 1;

      /* Threshold, and subtract 1 so that mask overcovers stim */
      if ( innerRadiusOfAnnulus <= minEcc)
        innerRadiusOfAnnulus = minEcc - 1;

      /*** GL calls to draw annulus ***/
      glColor3f(0.5, 0.5, 0.5);
      gluDisk(qobj,
              (GLdouble) innerRadiusOfAnnulus,
              (GLdouble) outerRadiusOfAnnulus,
              360,
              1);

    } else {
      innerRadiusOfAnnulus = log(minEcc) +
                             eccRange * (1/T) * (t - T/2);
      innerRadiusOfAnnulus = exp(innerRadiusOfAnnulus);

      outerRadiusOfAnnulus = (log(minEcc) + log(maxEcc))/2 +
                             eccRange * (1/T) * (t - T/2);
      outerRadiusOfAnnulus = exp(outerRadiusOfAnnulus);

      /* Threshold, and add 1 (pixel?) so that mask overcovers stim */
      if ( outerRadiusOfAnnulus >= maxEcc)
        outerRadiusOfAnnulus = maxEcc + 1;

      /* Threshold, and subtract 1 so that mask overcovers stim */
      if ( innerRadiusOfAnnulus <= minEcc)
        innerRadiusOfAnnulus = minEcc - 1;

      /*** GL calls to draw annulus ***/
      glColor3f(0.5, 0.5, 0.5);
      gluDisk(qobj,
              (GLdouble) innerRadiusOfAnnulus,
              (GLdouble) outerRadiusOfAnnulus,
              360,
              1);


    }
  } else { /*** --------- EXP. INWARD ---------- ***/
    if ( t < T/2 ) {
      outerRadiusOfAnnulus = (log(minEcc) + log(maxEcc))/2 -
                             eccRange * (1/T) * t;
      outerRadiusOfAnnulus = exp(outerRadiusOfAnnulus);
      innerRadiusOfAnnulus = minEcc;

      /* Threshold, and add 1 (pixel?) so that mask overcovers stim */
      if ( outerRadiusOfAnnulus >= maxEcc)
        outerRadiusOfAnnulus = maxEcc + 1;

      /* Threshold, and subtract 1 so that mask overcovers stim */
      if ( innerRadiusOfAnnulus <= minEcc)
        innerRadiusOfAnnulus = minEcc - 1;

      /*** GL calls to draw annulus ***/
      glColor3f(0.5, 0.5, 0.5);
      gluDisk(qobj,
              (GLdouble) innerRadiusOfAnnulus,
              (GLdouble) outerRadiusOfAnnulus,
              360,
              1);

      outerRadiusOfAnnulus = maxEcc;
      innerRadiusOfAnnulus = log(maxEcc) -
                             eccRange * (1/T) * t;
      innerRadiusOfAnnulus = exp(innerRadiusOfAnnulus);

      /* Threshold, and add 1 (pixel?) so that mask overcovers stim */
      if ( outerRadiusOfAnnulus >= maxEcc)
        outerRadiusOfAnnulus = maxEcc + 1;

      /* Threshold, and subtract 1 so that mask overcovers stim */
      if ( innerRadiusOfAnnulus <= minEcc)
        innerRadiusOfAnnulus = minEcc - 1;

      /*** GL calls to draw annulus ***/
      glColor3f(0.5, 0.5, 0.5);
      gluDisk(qobj,
              (GLdouble) innerRadiusOfAnnulus,
              (GLdouble) outerRadiusOfAnnulus,
              360,
              1);

    } else {
      innerRadiusOfAnnulus = (log(minEcc) + log(maxEcc))/2 -
                             eccRange * (1/T) * (t - T/2);
      innerRadiusOfAnnulus = exp(innerRadiusOfAnnulus);

      outerRadiusOfAnnulus = log(maxEcc) -
                             eccRange * (1/T) * (t - T/2);
      outerRadiusOfAnnulus = exp(outerRadiusOfAnnulus);

      /* Threshold, and add 1 (pixel?) so that mask overcovers stim */
      if ( outerRadiusOfAnnulus >= maxEcc)
        outerRadiusOfAnnulus = maxEcc + 1;

      /* Threshold, and subtract 1 so that mask overcovers stim */
      if ( innerRadiusOfAnnulus <= minEcc)
        innerRadiusOfAnnulus = minEcc - 1;

      /*** GL calls to draw annulus ***/
      glColor3f(0.5, 0.5, 0.5);
      gluDisk(qobj,
              (GLdouble) innerRadiusOfAnnulus,
              (GLdouble) outerRadiusOfAnnulus,
              360,
              1);


    }
  }

  /*  glTranslatef( -Width/2 , -Height/2 , 0 ); */

  glEndList() ;
  glCallList(current_list) ;
  /*  glFlush();*/
  glutSwapBuffers();
}

void display_all(void) {

  GLUquadricObj *qobj = NULL;
  double innerRadiusOfRing;
  double outerRadiusOfRing;
  double start_angle;
  double angular_step;
  static int current_list = 0 ;
  int ring, spoke ;
  double scaleFactor;


  glClear(GL_COLOR_BUFFER_BIT);
  current_list = !current_list ;
  glNewList(current_list, GL_COMPILE) ;

  angular_step = 360.0/numSpokes;

  scaleFactor = pow( (maxEcc/minEcc),
                     (1.0/( (double) numRings) ) );

  innerRadiusOfRing = minEcc;
  outerRadiusOfRing = minEcc * scaleFactor;


  /* draw radial checkerboard */
  for (ring = 0; ring < numRings; ring++) {
    for (spoke = 0; spoke < numSpokes; spoke++) {
      start_angle = spoke * angular_step;

      glColor3f(redFrames[ring][spoke][currentFrame], greenFrames[ring][spoke][currentFrame],
                blueFrames[ring][spoke][currentFrame]) ;

      qobj = gluNewQuadric();
      gluPartialDisk(qobj,
                     (GLdouble) innerRadiusOfRing,  /* inner radius */
                     (GLdouble) outerRadiusOfRing,  /* outer radius */
                     (GLint) angular_step, /* #vert. on disk perimeter */
                     1,   /* # vertices on disk ray (minus 1) */
                     (GLdouble) start_angle,  /* start angle for wedge */
                     (GLdouble) angular_step); /* wedge width in degrees */

    }

    innerRadiusOfRing *= scaleFactor;
    outerRadiusOfRing *= scaleFactor;

  }

  /*  glTranslatef( -Width/2 , -Height/2 , 0 ); */

  glEndList() ;
  glCallList(current_list) ;
  /*  glFlush();*/
  glutSwapBuffers();
}

void idle(void) {
  float time ;

  time = glutGet( GLUT_ELAPSED_TIME ); /* cumulative time */

  if (time > timePoint) {
    /* printf("time %f\n", time); */
    timePoint += timeTick;
    if (++currentFrame >= numFrames)
      currentFrame = 0 ;

    glutPostRedisplay();
  }

}

void init(void) {
  int  ring, spoke, frame, bit, start_bit  ;

  glClearColor (0.5, 0.5, 0.5, 0.5);
  glShadeModel (GL_FLAT);
  timeTick = 1000/flickerFreq;

  redFrames = (float ***)calloc(numRings, sizeof(float **)) ;
  blueFrames = (float ***)calloc(numRings, sizeof(float **)) ;
  greenFrames = (float ***)calloc(numRings, sizeof(float **)) ;
  if (!redFrames || !blueFrames || !greenFrames)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate ring list of %d x %d x %d buffer",
              Progname, numRings, numSpokes, numFrames) ;
  for (ring = 0 ; ring < numRings ; ring++) {
    redFrames[ring] = (float **)calloc(numSpokes, sizeof(float *)) ;
    blueFrames[ring] = (float **)calloc(numSpokes, sizeof(float *)) ;
    greenFrames[ring] = (float **)calloc(numSpokes, sizeof(float *)) ;
    if (!redFrames[ring] || !blueFrames[ring] || !greenFrames[ring])
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate spoke list [%d] of %d x %d x %d buffer",
                Progname, ring, numRings, numSpokes, numFrames) ;
    for (spoke = 0 ; spoke < numSpokes ; spoke++) {
      redFrames[ring][spoke] = (float *)calloc(numFrames, sizeof(float)) ;
      blueFrames[ring][spoke] = (float *)calloc(numFrames, sizeof(float)) ;
      greenFrames[ring][spoke] = (float *)calloc(numFrames, sizeof(float)) ;
      if (!redFrames[ring][spoke] || !blueFrames[ring][spoke] || !greenFrames[ring][spoke])
        ErrorExit(ERROR_NOMEMORY, "%s: could not allocate frame list [%d,%d] of %d x %d x %d buffer",
                  Progname, ring, spoke, numRings, numSpokes, numFrames) ;
    }
  }

  if (use_mseq) {
    unsigned char *mseq ;
    int      nbr_delay, delay, index, len ;

    nbr_delay = rint(pow(2,mseq_order) / (numRings*numSpokes)) ;

    printf("generating m-sequence, order=%d, tap=%d...", mseq_order, mseq_tap) ;
    mseq = generate_msequence(mseq_order, mseq_tap, &len) ;
    printf("finished\n") ;
    for (frame = 0 ; frame < numFrames ; frame++) {
      delay = 0 ;
      for (ring = 0 ; ring < numRings ; ring++) {
        for (spoke = 0 ; spoke < numSpokes ; spoke++, delay += nbr_delay) {
          index = (frame + delay) % len ;
          if (mseq[index]) {
            redFrames[ring][spoke][frame] = 1 ;
            blueFrames[ring][spoke][frame] = 1 ;
            greenFrames[ring][spoke][frame] = 1 ;
          } else {
            redFrames[ring][spoke][frame] = 0 ;
            blueFrames[ring][spoke][frame] = 0 ;
            greenFrames[ring][spoke][frame] = 0 ;
          }

        }
      }
    }
    free(mseq) ;
  } else {
    for (frame = 0 ; frame < numFrames ; frame++) {
      start_bit = frame % 2 ;
      for (ring = 0 ; ring < numRings ; ring++) {
        if ( ring % 2 )
          bit = start_bit ;
        else
          bit = !start_bit ;
        for (spoke = 0 ; spoke < numSpokes ; spoke++) {
          if ( (spoke + bit) % 2 ) {
            redFrames[ring][spoke][frame] = On ;
            blueFrames[ring][spoke][frame] = On ;
            greenFrames[ring][spoke][frame] = On ;
          } else {
            redFrames[ring][spoke][frame] = 1-On ;
            blueFrames[ring][spoke][frame] = 1-On ;
            greenFrames[ring][spoke][frame] = 1-On ;
          }

        }
      }
    }
  }

}

void reshape(int w, int h) {
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-50.0, 50.0, -50.0, 50.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/*#include <term.h>*/

int main(int argc, char** argv) {
  int     nargs ;

  {
    char c ;

    while (feof(stdin))
      printf(".") ;
    c = fgetc(stdin) ;
    printf("c = %c\n", c) ;
    exit(0) ;
  }

  nargs = handleVersionOption(argc, argv, "stim_polar");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  setenv("__GL_SYNC_TO_VBLANK", "23", 1) ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  glutInit(&argc, argv);
  if (argc != 8) {
    printf("\nUsage: eccen <flickerFreq> <stimPeriod (ms)> <minEcc> <maxEcc> "
           "<numRings> <numSpokes> <OUTWARD>\n\n");
    exit(1);
  }

  flickerFreq = atof(argv[1]);
  stimPeriod = atof(argv[2]);
  minEcc = atof(argv[3]);
  maxEcc = atof(argv[4]);
  numRings = atoi(argv[5]);
  numSpokes = atoi(argv[6]);
  OUTWARD = atoi(argv[7]);

  if (minEcc < EPS) {
    printf("\nArgument 3 (minEcc) must be larger than %.1e\n\n", EPS);
    exit(1);
  }


  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize (width, height);
  glutInitWindowPosition (100, 100);
  glutCreateWindow (argv[0]);
  init ();
#if 0
  glNewList(0, GL_COMPILE) ;
  glNewList(1, GL_COMPILE) ;
#endif

  glutDisplayFunc(display_all);
  glutReshapeFunc(reshape);
  glutIdleFunc(idle);
  glutMainLoop();
  return 0;
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
  if (!stricmp(option, "MSEQ")) {
    use_mseq = 1 ;
  } else switch (toupper(*option)) {
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static unsigned char *
generate_msequence(int order, unsigned int tap, int *plen) {
  unsigned int    len, q, R, i ;
  unsigned char   *mseq, *used ;

  len = pow(2, order)-1 ;
  q = pow(2, order-1) ;    /* high bit */
  mseq = (unsigned char *)calloc(len, sizeof(unsigned char)) ;
  used = (unsigned char *)calloc(len, sizeof(unsigned char)) ;

  R = 1 ;  /* shift register */
  for (i = 0 ; i < len ; i++) {
    if (used[R])
      ErrorExit(ERROR_BADPARM, "\n%s: tap %x does not generate an m-sequence!\n",
                Progname, tap) ;
    used[R] = 1 ;
    mseq[i] = (R & 0x01) ;
    if ((R & q) != q) {
      R = R << 1 ;
    } else {
      R = R << 1 ;
      R = R ^ tap ;
    }
    R = R & len ;
  }

  *plen = len ;

  free(used) ;
  return(mseq) ;
}

