/**
 * @file  test_gltexture.cpp
 * @brief test gltexturing
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.6 $
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

const char* Progname = "test_gltexture";

GLuint* textureID = NULL;

GLenum error;

#define CheckError()  \
  error = glGetError(); \
  while( error != GL_NO_ERROR ) { \
    cerr << __LINE__ << " error: " << gluErrorString( error ) << endl; \
    error = glGetError(); \
  } \

int nSlice = 0;
int nTexture = 0;
int width, height, depth;
int tSize[3];
GLuint textures[8];
GLfloat m[16];

float intensity = 1.0;
float density = 0.0;

void LoadTextures () {

  VolumeCollection vol;
  string fn( "/home/kteich/freesurfer/subjects/test.mgh" );
  vol.SetFileName( fn );
  MRI* mri = vol.GetMRI();

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  Matrix44& t = vol.GetWorldToIndexTransform().ExtractRotation();
  cerr << t << endl;
  m[0] = t(0,0);
  m[4] = t(1,0);
  m[8] = t(2,0);
  m[12] = t(3,0);
  m[1] = t(0,1);
  m[5] = t(1,1);
  m[9] = t(2,1);
  m[13] = t(3,1);
  m[2] = t(0,2);
  m[6] = t(1,2);
  m[10] = t(2,2);
  m[14] = t(3,2);
  m[3] = t(0,3);
  m[7] = t(1,3);
  m[11] = t(2,3);
  m[15] = t(3,3);

  nSlice = depth / 2;

  glClearColor( 0.0, 0.0, 0.0, 0.0 );
  glShadeModel( GL_FLAT );




  tSize[0] = width/2;
  tSize[1] = height/2;
  tSize[2] = depth/2;

  cerr << "Data size: " << width << "x" << height << "x" << depth << endl;
  cerr << "Texture size: " << tSize[0] << "x" << tSize[1] << "x" << tSize[2] << endl;

#if 0
  // check memory
  glTexImage3D( GL_PROXY_TEXTURE_3D, 0, GL_RGB,
                tSize[0], tSize[1], tSize[2], 0,
                GL_RGB, GL_UNSIGNED_BYTE,
                NULL );
  CheckError();

  GLint testWidth;
  glGetTexLevelParameteriv( GL_PROXY_TEXTURE_3D, 0,
                            GL_TEXTURE_WIDTH, &testWidth );
  CheckError();
  if ( testWidth != tSize[0] ) {
    fprintf( stderr, "test width was %d\n", testWidth );
    exit( 0 );
  }
#endif


  // gen texture names
  glGenTextures( 8, textures );


  GLubyte* buffer = NULL;
  buffer = (GLubyte*) malloc( tSize[0] * tSize[1] * tSize[2] *
                              sizeof(GLubyte) );

  // make textures
  int i = 0;
  for ( int nTexZ = 0; nTexZ < 2; nTexZ++ ) {
    for ( int nTexY = 0; nTexY < 2; nTexY++ ) {
      for ( int nTexX = 0; nTexX < 2; nTexX++ ) {

        for ( int nBufZ = 0; nBufZ < tSize[2]; nBufZ++ ) {
          for ( int nBufY = 0; nBufY < tSize[1]; nBufY++ ) {
            for ( int nBufX = 0; nBufX < tSize[0]; nBufX++ ) {

              int nMRIX = nBufX + (nTexX * tSize[0]);
              int nMRIY = nBufY + (nTexY * tSize[1]);
              int nMRIZ = nBufZ + (nTexZ * tSize[2]);

              buffer[nBufZ*tSize[1]*tSize[0] + nBufY*tSize[0] + nBufX] =
                MRIvox(mri,nMRIZ,nMRIY,nMRIX);
            }
          }
        }

        glBindTexture( GL_TEXTURE_3D, textures[i] );

        glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
        CheckError();

        glTexImage3D( GL_TEXTURE_3D, 0, GL_ALPHA,
                      tSize[0], tSize[1], tSize[2], 0,
                      GL_INTENSITY, GL_UNSIGNED_BYTE,
                      buffer );
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

        i++;
      }
    }
  }


  glEnable( GL_TEXTURE_3D  );
  CheckError();

  //  glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
  //  CheckError();

  glShadeModel( GL_FLAT );
  CheckError();

  free( buffer );
}

int windowWidth = 500;
int windowHeight = 500;
GLfloat rotation = 0.0f;

void DrawWindow () {

  GLfloat xPlane[] = { 1, 0, 0, 0 };
  GLfloat yPlane[] = { 0, 1, 0, 0 };
  GLfloat zPlane[] = { 0, 0, 1, 0 };
  GLfloat qPlane[] = { 0, 0, 0, 1 };

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  CheckError();

  glPushAttrib( GL_ALL_ATTRIB_BITS );
  glEnable( GL_TEXTURE_3D );
  glEnable( GL_TEXTURE_GEN_S );
  CheckError();
  glEnable( GL_TEXTURE_GEN_T );
  CheckError();
  glEnable( GL_TEXTURE_GEN_R );
  CheckError();
  glEnable( GL_TEXTURE_GEN_Q );
  CheckError();

  glEnable( GL_ALPHA_TEST );
  glAlphaFunc( GL_GREATER, intensity*density );
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

#if 0
  int i = 0;
  for ( int nTexZ = 0; nTexZ < 2; nTexZ++ ) {
    for ( int nTexY = 0; nTexY < 2; nTexY++ ) {
      for ( int nTexX = 0; nTexX < 2; nTexX++ ) {
#endif

        for ( int i = nTexture; i < nTexture+4; i++ ) {

          int t = i % 8;

          glBindTexture( GL_TEXTURE_3D, textures[t] );

          glMatrixMode( GL_TEXTURE );
          glLoadIdentity();

          switch ( t ) {
          case 0:
            break;
          case 1:
            break;
          case 2:
            break;
          case 3:
            break;
          case 4:
            break;
          case 5:
            break;
          case 6:
            break;
          case 7:
            break;
          }

          int textureSlice = nSlice % 128;

          glTranslatef( 0, 0, (float)textureSlice/(float)tSize[2] );

          //   glMultMatrixf( m );

          glRotatef( rotation, 0, 1, 0 );

          glTexGeni( GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
          CheckError();
          glTexGeni( GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
          CheckError();
          glTexGeni( GL_R, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
          CheckError();
          glTexGeni( GL_Q, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
          CheckError();
          glTexGenfv( GL_S, GL_OBJECT_PLANE, xPlane );
          CheckError();
          glTexGenfv( GL_T, GL_OBJECT_PLANE, yPlane );
          CheckError();
          glTexGenfv( GL_R, GL_OBJECT_PLANE, zPlane );
          CheckError();
          glTexGenfv( GL_Q, GL_OBJECT_PLANE, qPlane );
          CheckError();

          glMatrixMode( GL_MODELVIEW );
          glLoadIdentity();
          glTranslatef( 0.0f, 0.0f, -5.0f);
          switch ( t ) {
          case 0:
            glTranslatef( -1, -1,  0 );
            break;
          case 1:
            glTranslatef(  0, -1,  0 );
            break;
          case 2:
            glTranslatef( -1,  0,  0 );
            break;
          case 3:
            glTranslatef(  0,  0,  0 );
            break;
          case 4:
            glTranslatef( -1, -1,  0 );
            break;
          case 5:
            glTranslatef(  0, -1,  0 );
            break;
          case 6:
            glTranslatef( -1,  0,  0 );
            break;
          case 7:
            glTranslatef(  0,  0,  0 );
            break;
          }

          glColor4f( 1, 1, 1, intensity );
          glBegin( GL_QUADS );

          glVertex3f( 0, 0, 0 );
          glVertex3f( 0, 1, 0 );
          glVertex3f( 1, 1, 0 );
          glVertex3f( 1, 0, 0 );

          glEnd();
          glGetError();  // get rid of error

          glRasterPos2f( 0, 0 );
          glColor3f( 1, 1, 0 );
          char sLabel[60];
          sprintf( sLabel, "%d", t );
          for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
            glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
          }
        }

#if 0
        i++;

      }
    }
  }
#endif


  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  glTranslatef( 0.0f, 0.0f, -5.0f);

  glDisable( GL_ALPHA_TEST );
  glDisable( GL_BLEND );
  char sLabel[60];
  sprintf( sLabel, "tex %d-%d, slice %d, rotation %.2f i %.1f d %.1f",
           nTexture, (nTexture+4)%8, nSlice, rotation, intensity, density );
  glRasterPos2f( -1.4, -1.4 );
  glColor3f( 0, 1, 0 );
  for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }

  glPopAttrib();
  CheckError();

  glutSwapBuffers();
  CheckError();
}

void ReshapeWindow ( int iWidth, int iHeight ) {

  windowWidth = iWidth;
  windowHeight = iHeight;

  glViewport( 0, 0, windowWidth, windowHeight );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho( -1.5, 1.5, -1.5, 1.5, 0.1, 100 );
}

void HandleKeyboard ( unsigned char iKey, int iX, int iY ) {


  switch ( iKey ) {
  case '.':
    nSlice++;
    if ( nSlice >= depth ) nSlice = 0;
    break;
  case ',':
    nSlice--;
    if ( nSlice < 0 ) nSlice = depth - 1;
    break;
  case ']':
    rotation++;
    if ( rotation > 360.0f ) rotation = 0.0f;
    break;
  case '[':
    rotation--;
    if ( rotation < 0 ) rotation = 359;
    break;
  case 'T':
    nTexture++;
    if ( nTexture >= 8 ) nTexture = 0;
    break;
  case 't':
    nTexture--;
    if ( nTexture < 0 ) nTexture = 8;
    break;
  case 'd':
    density -= 0.1;
    if ( density < 0 ) density = 1;
    break;
  case 'D':
    density += 0.1;
    if ( density > 1 ) density = 0;
    break;
  case 'i':
    intensity -= 0.1;
    if ( intensity < 0 ) intensity = 1;
    break;
  case 'I':
    intensity += 0.1;
    if ( intensity > 1 ) intensity = 0;
    break;
  }
  glutPostRedisplay();
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
  //glutTimerFunc( 10, Timer, 0 );
  glutKeyboardFunc( HandleKeyboard );
}


int main ( int argc, char** argv ) {

  glutInit( &argc, argv );

  BuildWindow();
  LoadTextures();

  glutMainLoop();
}
