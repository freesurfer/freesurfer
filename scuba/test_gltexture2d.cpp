/**
 * @file  test_gltexture2d.cpp
 * @brief tests
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.4 $
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


#include <string>
#include "glut.h"
#include <GL/glx.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <map>
#include <sstream>
#include "VolumeCollection.h"
#include "Scuba-impl.h"

using namespace std;

char* Progname = "test_gltexture";

GLuint* textureID = NULL;

GLenum error;

#define CheckError()  \
  error = glGetError(); \
  while( error != GL_NO_ERROR ) { \
    cerr << __LINE__ << " error: " << gluErrorString( error ) << endl; \
    error = glGetError(); \
  } \


#define USE_3D_TEXTURE 1

int nSlice = 0;
int width, height, depth;

static GLubyte image[16][16][16][3];

void LoadTextures () {

  VolumeCollection vol;
  string fn( "/home/kteich/freesurfer/subjects/bert/mri/T1" );
  vol.SetFileName( fn );
#if 0
  MRI* mri = vol.GetMRI();

  nSlice = mri->depth / 2;

  GLubyte* texture = NULL;
  texture = (GLubyte*) malloc( mri->width * mri->height * mri->depth *
                               sizeof(GLubyte) );

  for ( int nZ = 0; nZ < mri->depth; nZ++ ) {
    for ( int nY = 0; nY < mri->height; nY++ ) {
      for ( int nX = 0; nX < mri->width; nX++ ) {
        texture[nZ*mri->height*mri->width + nY*mri->width + nX] =
          MRIvox(mri,nZ,nY,nX);
      }
    }
  }
#endif

  glClearColor( 0.0, 0.0, 0.0, 0.0 );
  glShadeModel( GL_FLAT );
  glEnable( GL_DEPTH_TEST );

#if USE_3D_TEXTURE

  for (int s = 0; s < 16; s++)
    for (int t = 0; t < 16; t++)
      for (int r = 0; r < 16; r++) {
        image[r][t][s][0] = (GLubyte) (s * 17);
        image[r][t][s][1] = (GLubyte) (t * 17);
        image[r][t][s][2] = (GLubyte) (r * 17);
      }


  width = height = depth = 16;
  nSlice = depth / 2;

  glPixelStorei( GL_UNPACK_ROW_LENGTH, width );
  CheckError();
  glPixelStorei( GL_UNPACK_IMAGE_HEIGHT, height );
  CheckError();
  glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
  CheckError();

  glTexParameterf( GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP );
  CheckError();
  glTexParameterf( GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP );
  CheckError();
  glTexParameterf( GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP );
  CheckError();
  glTexParameterf( GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  CheckError();
  glTexParameterf( GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  CheckError();

  glTexImage3D( GL_TEXTURE_3D, 0, GL_RGB,
                width, height, depth, 0,
                GL_RGB, GL_UNSIGNED_BYTE,
                image );

#if 0
  glTexImage3D( GL_TEXTURE_3D, 0, GL_RGB,
                width, height, depth, 0,
                GL_LUMINANCE, GL_UNSIGNED_BYTE,
                &(texture[0]) );
  CheckError();
#endif

  glEnable( GL_TEXTURE_3D  );
  CheckError();


#else

  textureID = (GLuint*) malloc( mri->depth * sizeof(GLuint) );
  glGenTextures( mri->depth, textureID );
  CheckError();

  for ( int nZ = 0; nZ < mri->depth; nZ++ ) {


    glBindTexture( GL_TEXTURE_2D, textureID[nZ] );
    CheckError();
    glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
    CheckError();
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB,
                  width, height, 0,
                  GL_LUMINANCE, GL_UNSIGNED_BYTE,
                  &(texture[nZ*mri->width*mri->height]) );
    CheckError();


    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
    CheckError();
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
    CheckError();
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    CheckError();
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    CheckError();

    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
    CheckError();
  }

  glEnable( GL_TEXTURE_2D );

#endif

  glClearColor (0.0, 0.0, 0.0, 0.0);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);

  glShadeModel( GL_FLAT );
  CheckError();

  // free( texture );
}





int windowWidth = 500;
int windowHeight = 500;
GLfloat rotation = 0.0f;

void DrawWindow () {

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  CheckError();

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  glTranslatef( 0.0f, 0.0f, -5.0f);
  glRotatef( rotation, 0.0f, 1.0f, 0.0f );


#if USE_3D_TEXTURE

  glPushAttrib( GL_ALL_ATTRIB_BITS );
  glEnable( GL_TEXTURE_3D );

#if 0
  glEnable( GL_ALPHA_TEST );
  glAlphaFunc( GL_GREATER, 0.5 );
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
#endif

  glMatrixMode( GL_TEXTURE );
  glLoadIdentity();

  glBegin( GL_QUADS );

  glTexCoord3f( 0.0f, 0.0f, (float)nSlice/(float)depth );
  glVertex3f  ( -1.0f, -1.0f, 0 );

  glTexCoord3f( 0.0f, 1.0f, (float)nSlice/(float)depth );
  glVertex3f  ( -1.0f, 1.0f, 0 );

  glTexCoord3f( 1.0f, 1.0f, (float)nSlice/(float)depth );
  glVertex3f  ( 1.0f, 1.0f, 0 );

  glTexCoord3f( 1.0f, 0.0f, (float)nSlice/(float)depth );
  glVertex3f  ( 1.0f, -1.0f, 0 );

  glEnd();
  glGetError();  // get rid of error

  glPopAttrib();
  CheckError();

#else

  glBindTexture( GL_TEXTURE_2D, textureID[nSlice] );
  CheckError();

  glBegin( GL_QUADS );

  glColor3f   ( 1.0f, 0.0f, 0.0f );
  glTexCoord2f( 0.0f, 0.0f );
  glVertex3f  ( -1.0f, -1.0f, 0.0f );

  glColor3f   ( 0.0f, 1.0f, 0.0f );
  glTexCoord2f( 0.0f, 1.0f );
  glVertex3f  ( -1.0f, 1.0f, 0.0f );

  glColor3f   ( 0.0f, 0.0f, 1.0f );
  glTexCoord2f( 1.0f, 1.0f );
  glVertex3f  ( 1.0f, 1.0f, 0.0f );

  glColor3f   ( 1.0f, 0.0f, 1.0f );
  glTexCoord2f( 1.0f, 0.0f );
  glVertex3f  ( 1.0f, -1.0f, 0.0f );

  glEnd();
  glGetError();  // get rid of error

#endif

  glutSwapBuffers();
  CheckError();
}

void ReshapeWindow ( int iWidth, int iHeight ) {

  windowWidth = iWidth;
  windowHeight = iHeight;

  glViewport( 0, 0, windowWidth, windowHeight );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 45.0f, (GLfloat) windowWidth / windowHeight, 0.1f, 100.0f );

}

void HandleKeyboard ( unsigned char iKey, int iX, int iY ) {


  switch ( iKey ) {
  case '.':
    nSlice++;
    if ( nSlice >= depth ) {
      nSlice = 0;
    }
    break;
  case ',':
    nSlice--;
    if ( nSlice < 0 ) {
      nSlice = depth - 1;
    }
    break;
  }
}

void Timer ( int value ) {

  rotation++;
  if ( rotation > 360.0f ) rotation = 0.0f;
  glutPostRedisplay();

  glutTimerFunc( 10, Timer, value );
}

void BuildWindow () {

  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowSize( windowWidth, windowHeight );
  glutCreateWindow( "gltexture" );
  glutDisplayFunc( DrawWindow );
  glutReshapeFunc( ReshapeWindow );
  glutTimerFunc( 10, Timer, 0 );
  glutKeyboardFunc( HandleKeyboard );
}


int main ( int argc, char** argv ) {

  glutInit( &argc, argv );

  BuildWindow();
  LoadTextures();

  glutMainLoop();
}
