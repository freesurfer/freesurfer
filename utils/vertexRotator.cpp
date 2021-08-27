/**
 * @brief rapidly rotate vertexs around the x=0 y=0 z=0 axes
 */
/*
 * Original Author: Bevin Brett based on code in mrisurf.c
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
#include "vertexRotator.h"

//  The command line for the self-test is
//
//      gcc -I../include -DTEST vertexRotator.c -lm

// Note: it is expected that the caller will do any needed parallelism
// because they can tile along the vertices dimension to optimize for cache behaviour

// This is the special case of two of the angles being zero
// The invariant dimension is optional
// 
// It can be faster when some of the angles are zero
// The inliner is relied on to make those fixes
// Only one implementation is needed because the axes are permuted below accordingly
//
static void rotateVertices1axis_wkr(
  float*         xv_out, float*       yv_out, 
  float const *  xv_inp, float const* yv_inp, 
  float*         zv_out,
  float const*   zv_inp,
  size_t         nvertices,
  float          alpha)       // rotate around x axis
{
  float const sa = sin(alpha);
  float const ca = cos(alpha);

  for (unsigned vno = 0; vno < nvertices; vno++) {
    float x  = xv_inp[vno];
    float y  = yv_inp[vno];
    float xp =  x * ca + y * sa;
    float yp = -x * sa + y * ca;
    xv_out[vno] = xp;
    yv_out[vno] = yp;
  }

  if (zv_out) {
    for (unsigned int vno = 0; vno < nvertices; vno++) {
      zv_out[vno] = zv_inp[vno];
    }
  }
}


// This is the general case
// It can be faster when some of the angles are zero
// The inliner is relied on to make those fixes
//
static inline void rotateVertices_wkr(
  float*         xv_out, float*       yv_out, float      * zv_out,
  float const *  xv_inp, float const* yv_inp, float const* zv_inp,
  size_t         nvertices,
  float          alpha,       // rotate around z axis - last rotation
  float          beta,        // rotate around y axis - middle rotation
  float          gamma)       // rotate around x axis - first rotation
{
  float const sa = (alpha == 0.0) ? 0.0 : sin(alpha);
  float const sb = (beta  == 0.0) ? 0.0 : sin(beta);
  float const sg = (gamma == 0.0) ? 0.0 : sin(gamma);
  float const ca = (alpha == 0.0) ? 1.0 : cos(alpha);
  float const cb = (beta  == 0.0) ? 1.0 : cos(beta);
  float const cg = (gamma == 0.0) ? 1.0 : cos(gamma);

  float const cacb = ca * cb;
  float const cacgsb = ca * cg * sb;
  float const sasg = sa * sg;
  float const cgsa = cg * sa;
  float const casbsg = ca * sb * sg;
  float const cbsa = cb * sa;
  float const cgsasb = cg * sa * sb;
  float const casg = ca * sg;
  float const cacg = ca * cg;
  float const sasbsg = sa * sb * sg;
  float const cbcg = cb * cg;
  float const cbsg = cb * sg;
  
  for (unsigned int vno = 0; vno < nvertices; vno++) {

    float x = xv_inp[vno];
    float y = yv_inp[vno];
    float z = zv_inp[vno];
                                                                            // when sa == sb == 0       when sa == sg == 0      when sb == sg == 0
    float xp =  x * cacb + z * (-cacgsb - sasg) + y * (cgsa - casbsg);      //             x            x cb - z sb             x ca  + y sa
    float yp = -x * cbsa + z * ( cgsasb - casg) + y * (cacg + sasbsg);      // - z sg + y cg            y                       -x sa + y ca
    float zp =  z * cbcg + x * sb               + y * cbsg;                 //   z cg + y sg            x sb - z cb             z
    
    xv_out[vno] = xp;
    yv_out[vno] = yp;
    zv_out[vno] = zp;
  }
}



void rotateVertices(
    float*         xv_out, float*       yv_out, float      * zv_out,
    float const *  xv_inp, float const* yv_inp, float const* zv_inp,
    size_t         nvertices,
    float          alpha,       // alpha is the last rotation
    float          beta,        // beta  is the second rotation
    float          gamma)       // gamma is the first rotation
{
    // One axis rotation
    if (alpha == 0.0 && beta  == 0.0) { rotateVertices1axis_wkr( yv_out, zv_out, yv_inp, zv_inp, xv_out, xv_inp, nvertices, -gamma); return; }
    if (alpha == 0.0 && gamma == 0.0) { rotateVertices1axis_wkr( xv_out, zv_out, xv_inp, zv_inp, yv_out, yv_inp, nvertices, -beta ); return; }
    if (beta  == 0.0 && gamma == 0.0) { rotateVertices1axis_wkr( xv_out, yv_out, xv_inp, yv_inp, zv_out, zv_inp, nvertices,  alpha); return; }
    
    // Three axis rotation
    if (alpha != 0.0 && beta != 0.0 && gamma != 0.0) 
                                      { rotateVertices_wkr( xv_out, yv_out, zv_out, xv_inp, yv_inp, zv_inp, nvertices, alpha, beta, gamma); return; }

    // Two axis rotation
    if (alpha != 0.0 && beta  != 0.0) { rotateVertices_wkr( xv_out, yv_out, zv_out, xv_inp, yv_inp, zv_inp, nvertices, alpha, beta, 0.0);   return; }
    if (alpha != 0.0 && gamma != 0.0) { rotateVertices_wkr( xv_out, yv_out, zv_out, xv_inp, yv_inp, zv_inp, nvertices, alpha, 0.0,  gamma); return; }
    if (beta  != 0.0 && gamma != 0.0) { rotateVertices_wkr( xv_out, yv_out, zv_out, xv_inp, yv_inp, zv_inp, nvertices, 0.0,   beta, gamma); return; }
}


#if defined(TEST) 
#include <stdio.h>

static const char* final = "PASSED";
static int errorCount;
static void check_wkr(float t, float o, float xi, float yi, float zi, const char* msg) {
    if (fabsf(t - o) < 0.01) return;
    printf("FAILED %s, t:%g o:%g (x:%g,y:%g,z:%g)\n", msg, t,o,xi,yi,zi);
    final = "FAILED";
    errorCount++;
}

#define check(t,o) check_wkr(t,o,xi,yi,zi,#t "!=" #o)

static void test(float xi, float yi, float zi, float xo, float yo, float zo, float alpha, float beta, float gamma) {
    float xt,yt,zt;
    rotateVertices(&xt,&yt,&zt,&xi,&yi,&zi,1,alpha,beta,gamma);
    int saved_errorCount = errorCount;
    check(xt,xo); check(yt,yo); check(zt,zo);
    if (saved_errorCount != errorCount) printf(" when alpha:%g beta:%g gamma:%g\n", alpha,beta,gamma);
}

int main() {
    double const pi = 3.141592653589793;
    printf("vectorRotator test\n");

    test(1,0,0, 1,0,0, 0,0,0);
    test(0,1,0, 0,1,0, 0,0,0);
    test(0,0,1, 0,0,1, 0,0,0);

    // Test singles
    //
    float tiny = 0.00001;
    for (tiny = 0.00001; ; tiny = 0) { 
        printf("around z\n");
        test(1,0,0,  0,-1,0,  pi/2,tiny,tiny);
        test(0,1,0, +1, 0,0,  pi/2,tiny,tiny);
        test(0,0,1,  0, 0,1,  pi/2,tiny,tiny);
        test(1,0,0,  0,+1,0, -pi/2,tiny,tiny);
        test(0,1,0, -1, 0,0, -pi/2,tiny,tiny);

        printf("around y\n");
        test(1,0,0,  0,0, 1, tiny, pi/2,tiny);
        test(0,1,0,  0,1, 0, tiny, pi/2,tiny);
        test(0,0,1, -1,0, 0, tiny, pi/2,tiny);
        test(1,0,0,  0,0,-1, tiny,-pi/2,tiny);
        test(0,0,1,  1,0, 0, tiny,-pi/2,tiny);

        printf("around x\n");
        test(1,0,0, 1, 0, 0, tiny,tiny, pi/2);
        test(0,1,0, 0, 0,+1, tiny,tiny, pi/2);
        test(0,0,1, 0,-1, 0, tiny,tiny, pi/2);
        test(0,1,0, 0, 0,-1, tiny,tiny,-pi/2);
        test(0,0,1, 0, 1, 0, tiny,tiny,-pi/2);

        if (tiny == 0.0f) break;
    }

    // Test composition
    //
    float xi=1,yi=2,zi=3;
    float x1,y1,z1;
    float x2a,y2a,z2a;
    float x2b,y2b,z2b;
    
    //  beta then alpha
    printf("beta then alpha\n");
    rotateVertices(&x1,&y1,&z1,    &xi,&yi,&zi, 1,    0, pi/2,0);
    rotateVertices(&x2a,&y2a,&z2a, &x1,&y1,&z1, 1, pi/2,    0,0);
    rotateVertices(&x2b,&y2b,&z2b, &xi,&yi,&zi, 1, pi/2, pi/2,0);
    check(x2a, x2b); check(y2a, y2b); check(z2a, z2b);
            
    //  gamma then alpha
    printf("gamma then alpha\n");
    rotateVertices(&x1,&y1,&z1,    &xi,&yi,&zi, 1,    0,0, pi/2);
    rotateVertices(&x2a,&y2a,&z2a, &x1,&y1,&z1, 1, pi/2,0,    0);
    rotateVertices(&x2b,&y2b,&z2b, &xi,&yi,&zi, 1, pi/2,0, pi/2);
    check(x2a, x2b); check(y2a, y2b); check(z2a, z2b);
            
    //  gamma then beta
    printf("gamma then beta\n");
    rotateVertices(&x1,&y1,&z1,    &xi,&yi,&zi, 1, 0, 0,   pi/2);
    rotateVertices(&x2a,&y2a,&z2a, &x1,&y1,&z1, 1, 0, pi/2,   0);
    rotateVertices(&x2b,&y2b,&z2b, &xi,&yi,&zi, 1, 0, pi/2,pi/2);
    check(x2a, x2b); check(y2a, y2b); check(z2a, z2b);
            
    printf("%s\n",final);
}
#endif
