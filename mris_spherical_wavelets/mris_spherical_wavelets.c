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


///////////////////////////////////////////
// mris_spherical_wavelets.c
//
// written by Peng Yu
// date: 11/10/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "icosahedron.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "matrix.h"


int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char            *Progname ;
static MRI_SURFACE *center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst);


static int      ANALYSIS=0;
static int      SYNTHESIS=0;
static int      COMPARE=0;
static int      CURV=0;
static char     *fname;

int
main(int argc, char *argv[]) {
  int           nargs, msec, order, i, number, vno, nnum, m, k, b1, b2, cno, flag=0, fno;
  Timer then ;
  MRIS          *mris_in, *mris_out, *mris_high;
  MRI_SP        *mrisp ;
  VERTEX        *vm_out, *vm_high, *v;
  float         s_jkm, area;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input surface> <orig surface> <finest order> <output surface>", Progname);

  then.reset() ;

  order = atoi (argv[3]);
  fprintf(stdout, "Set %s as the finest scale level\n", argv[3]);
  if (order > 7)
    ErrorExit(ERROR_BADPARM, "the highest order is 7\n");

  /*Spherical Wavelet Analysis*/

  if (ANALYSIS&&!CURV) {
    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);
    MRISreadOriginalProperties(mris_in, argv[2]) ;
    fprintf(stdout, "Reading original surface from %s orig area is %f\n", argv[2],mris_in->orig_area);

    mris_out = ReadIcoByOrder(order, 100);
    for (m = 0; m<mris_out->nvertices; m++)
      mris_out->vertices[m].nsize=1;
    mrisp = MRISPalloc(1, 3);
#if 1
    MRIScoordsToParameterization(mris_in, mrisp, 1, ORIGINAL_VERTICES) ;
    MRISPblur(mrisp, mrisp, 1, 0);
    MRISPblur(mrisp, mrisp, 1, 1);
    MRISPblur(mrisp, mrisp, 1, 2);
    MRIScoordsFromParameterization(mrisp, mris_out) ;
#else
    MRISreadOriginalProperties(mris_out, argv[2]) ;
#endif
#if 1 /*just to test if the parameterization is correct */
    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    MRISupdateSurface(mris_out);
    fprintf(stderr, "original area becomes %f\n", mris_out->total_area);
    center_brain(mris_out, mris_out);
    MRISscaleBrain(mris_out, mris_out, sqrt(100000.0f/mris_out->total_area)) ;
    MRISupdateSurface(mris_out);
    for (fno=0; fno<mris_out->nfaces; fno++)
      area += mris_out->faces[fno].area;
    fprintf(stderr, "original area becomes %f\n", area);
    //MRISwrite(mris_out, "/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.sampled") ;
    MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
#endif

    /* Initialize Ij,k*/
    for (vno = 0 ; vno<mris_out->nvertices; vno++) {
      vm_out = &mris_out->vertices[vno];
      vm_out->val = 1;
    }

    /*Iteratively compute Ij,k*/
    for (i=order;i>0;i--) {
      mris_high = ReadIcoByOrder(i, 100); //higher order surface
      for (m = 0; m<mris_high->nvertices; m++)
        mris_high->vertices[m].nsize=1;
      MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
      number = IcoNVtxsFromOrder(i-1); //the start of m vertices
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        flag=0;
        for (nnum=0; nnum<vm_high->vnum; nnum++)
          if ( vm_high->v[nnum]<number ) //A(j,m)
          {
            k = vm_high->v[nnum];
            v = &mris_out->vertices[k];
            v->val += 0.5*vm_out->val ;
          }
        for (; nnum<vm_high->v2num; nnum++)
          if ( vm_high->v[nnum]<number ) //B(j,m)
          {
            k = vm_high->v[nnum];
            if (flag==0) b1=k;
            else b2=k;
            flag++;
            v = &mris_out->vertices[k];
            v->val += 0.125*vm_out->val ;
          }
        for (; nnum<vm_high->v3num; nnum++)
          if ( vm_high->v[nnum]<number ) //C(j,m)
          {
            k = vm_high->v[nnum];
            flag=0; //C has to be a second-order neighbor of B
            for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
              if (mris_high->vertices[b1].v[cno]==k) flag=1;
            for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
              if (mris_high->vertices[b2].v[cno]==k) flag=1;
            if (flag) {
              v = &mris_out->vertices[k];
              v->val -= 0.0625*vm_out->val ;
            }
          }
      }
    }


    /*Analysis Stage I:*/
    for (i=order;i>0;i--) {
      mris_high = ReadIcoByOrder(i, 100); //higher order surface
      for (m = 0; m<mris_high->nvertices; m++)
        mris_high->vertices[m].nsize=1;
      MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;

      number = IcoNVtxsFromOrder(i-1); //the start of m vertices
      /* compute Yj,m for each m vertices */
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        flag=0;
        for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
          {
            k = vm_high->v[nnum] ;
            v = &mris_out->vertices[k];
            vm_out->origx -= 0.5*v->origx;
            vm_out->origy -= 0.5*v->origy;
            vm_out->origz -= 0.5*v->origz;
          }
        for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
          {
            k = vm_high->v[nnum] ;
            if (flag==0) b1=k;
            else b2=k;
            flag++;
            v = &mris_out->vertices[k];
            vm_out->origx -= 0.125*v->origx;
            vm_out->origy -= 0.125*v->origy;
            vm_out->origz -= 0.125*v->origz;
          }
        for (; nnum<vm_high->v3num; nnum++)
          if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
          {
            k = vm_high->v[nnum] ;
            flag=0; //C has to be a second-order neighbor of B
            for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
              if (mris_high->vertices[b1].v[cno]==k) flag=1;
            for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
              if (mris_high->vertices[b2].v[cno]==k) flag=1;
            if (flag) {
              v = &mris_out->vertices[k];
              vm_out->origx += 0.0625*v->origx;
              vm_out->origy += 0.0625*v->origy;
              vm_out->origz += 0.0625*v->origz;
            }
          }
      }


      /*Analysis Stage II: */
      /*Compute Lamda(j,k) using the Yita(j,m)*/
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        for (nnum=0; nnum<vm_high->vnum; nnum++)
          if ( vm_high->v[nnum]<number ) //A(j,m)
          {
            k = vm_high->v[nnum];
            v = &mris_out->vertices[k];
            s_jkm = vm_out->val/2/v->val;
            v->origx += s_jkm*vm_out->origx;
            v->origy += s_jkm*vm_out->origy;
            v->origz += s_jkm*vm_out->origz;
          }

      }

    }

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
#if 0
    for (m=0;m<mris_out->nvertices;m++)
      if (mris_out->vertices[m].z>6)
        fprintf(stdout, "%d %f %f %f\n", m,mris_out->vertices[m].x, mris_out->vertices[m].y, mris_out->vertices[m].z);
    //mris_high = ReadIcoByOrder(0, 100);
    //for (m=0;m<mris_high->nvertices;m++)
    //{mris_high->vertices[m].x=mris_out->vertices[m].x;
    //mris_high->vertices[m].y=mris_out->vertices[m].y;
    //mris_high->vertices[m].z=mris_out->vertices[m].z;
    //}
    //MRISwrite(mris_high, "/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.sampled") ;
#endif
    fprintf(stdout, "Writing wavelets coefficient of original surface to %s\n", argv[4]);
    MRISwrite(mris_out,argv[4] ) ;
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris_in) ;
    /*End of Analysis*/
  } else if (ANALYSIS&&CURV) {
    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);

    MRISreadCurvatureFile(mris_in, argv[2]) ;
    fprintf(stdout, "Reading input from %s\n", argv[2]);

    mris_out = ReadIcoByOrder(order, 100);
    for (m = 0; m<mris_out->nvertices; m++)
      mris_out->vertices[m].nsize=1;
    //mrisp = MRISPalloc(1, 3);
    mrisp = MRIStoParameterization(mris_in, NULL, 1, 0) ;
    //MRISPblur(mrisp, mrisp, 1, 0);
    MRISfromParameterization(mrisp, mris_out, 0) ;
    //MRISwriteCurvature(mris_out,"/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.thickness.sampled");
    /* Initialize Ij,k*/
    for (vno = 0 ; vno<mris_out->nvertices; vno++) {
      vm_out = &mris_out->vertices[vno];
      vm_out->val = 1;
    }

    /*Iteratively compute Ij,k*/
    for (i=order;i>0;i--) {
      mris_high = ReadIcoByOrder(i, 100); //higher order surface
      for (m = 0; m<mris_high->nvertices; m++)
        mris_high->vertices[m].nsize=1;
      MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
      number = IcoNVtxsFromOrder(i-1); //the start of m vertices
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        flag=0;
        for (nnum=0; nnum<vm_high->vnum; nnum++)
          if ( vm_high->v[nnum]<number ) //A(j,m)
          {
            k = vm_high->v[nnum];
            v = &mris_out->vertices[k];
            v->val += 0.5*vm_out->val ;
          }
        for (; nnum<vm_high->v2num; nnum++)
          if ( vm_high->v[nnum]<number ) //B(j,m)
          {
            k = vm_high->v[nnum];
            if (flag==0) b1=k;
            else b2=k;
            flag++;
            v = &mris_out->vertices[k];
            v->val += 0.125*vm_out->val ;
          }
        for (; nnum<vm_high->v3num; nnum++)
          if ( vm_high->v[nnum]<number ) //C(j,m)
          {
            k = vm_high->v[nnum];
            flag=0; //C has to be a second-order neighbor of B
            for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
              if (mris_high->vertices[b1].v[cno]==k) flag=1;
            for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
              if (mris_high->vertices[b2].v[cno]==k) flag=1;
            if (flag) {
              v = &mris_out->vertices[k];
              v->val -= 0.0625*vm_out->val ;
            }
          }
      }
    }


    /*Analysis Stage I:*/
    for (i=order;i>0;i--) {
      mris_high = ReadIcoByOrder(i, 100); //higher order surface
      for (m = 0; m<mris_high->nvertices; m++)
        mris_high->vertices[m].nsize=1;
      MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;

      number = IcoNVtxsFromOrder(i-1); //the start of m vertices
      /* compute Yj,m for each m vertices */
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        flag=0;
        for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
          {
            k = vm_high->v[nnum] ;
            v = &mris_out->vertices[k];
            vm_out->curv -= 0.5*v->curv;
          }
        for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
          {
            k = vm_high->v[nnum] ;
            if (flag==0) b1=k;
            else b2=k;
            flag++;
            v = &mris_out->vertices[k];
            vm_out->curv -= 0.125*v->curv;
          }
        for (; nnum<vm_high->v3num; nnum++)
          if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
          {
            k = vm_high->v[nnum] ;
            flag=0; //C has to be a second-order neighbor of B
            for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
              if (mris_high->vertices[b1].v[cno]==k) flag=1;
            for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
              if (mris_high->vertices[b2].v[cno]==k) flag=1;
            if (flag) {
              v = &mris_out->vertices[k];
              vm_out->curv += 0.0625*v->curv;
            }
          }
      }


      /*Analysis Stage II: */
      /*Compute Lamda(j,k) using the Yita(j,m)*/
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        for (nnum=0; nnum<vm_high->vnum; nnum++)
          if ( vm_high->v[nnum]<number ) //A(j,m)
          {
            k = vm_high->v[nnum];
            v = &mris_out->vertices[k];
            s_jkm = vm_out->val/2/v->val;
            v->curv += s_jkm*vm_out->curv;
          }

      }
    }

    fprintf(stdout, "Writing wavelets coefficient of original surface to %s\n", argv[4]);
    MRISwriteCurvature(mris_out,argv[4] ) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris_in) ;
    /*End of Analysis*/
  } else if (SYNTHESIS) /*Spherical Wavelet Synthesis*/
  {
    mris_out = ReadIcoByOrder(order, 100); //higher order surface
    fprintf(stdout, "Creating a %d order spherical surface\n", order);
    MRISreadOriginalProperties(mris_out, argv[1]) ;
    fprintf(stdout, "Reading wavelet coefficients from %s\n", argv[1]);
    for (m = 0; m<mris_out->nvertices; m++)
      mris_out->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_out, 3) ;

    if (COMPARE) {
      mris_in = MRISread(fname);
      for (i=1; i<IcoNVtxsFromOrder(order-1); i++) {
        if (mris_out->vertices[i].origx==0)
          area =  fabs(mris_out->vertices[i].origx-mris_in->vertices[i].x);
        else area = fabs((mris_out->vertices[i].origx-mris_in->vertices[i].x)/mris_out->vertices[i].origx);
        if ( area>5 ) {
          mris_out->vertices[i].origx = mris_in->vertices[i].x ;
          fprintf(stdout, "%d %f\n", i, area);
        }
        if (mris_out->vertices[i].origy==0)
          area =  fabs(mris_out->vertices[i].origy-mris_in->vertices[i].y);
        else area = fabs((mris_out->vertices[i].origy-mris_in->vertices[i].y)/mris_out->vertices[i].origy);
        if ( area>5 ) {
          mris_out->vertices[i].origy = mris_in->vertices[i].y ;
          fprintf(stdout, "%d %f\n", i, area);
        }
        if (mris_out->vertices[i].origz==0)
          area =  fabs(mris_out->vertices[i].origz-mris_in->vertices[i].z);
        else area = fabs((mris_out->vertices[i].origz-mris_in->vertices[i].z)/mris_out->vertices[i].origz);
        if ( area>5 ) {
          mris_out->vertices[i].origz = mris_in->vertices[i].z ;
          fprintf(stdout, "%d %f\n", i, area);
        }
      }
      MRISfree(&mris_in);
    }

    fprintf(stdout, "Recover the surface using %s order coefficients\n",argv[2]);
    number = IcoNVtxsFromOrder(atoi(argv[2]));
    for (m = number; m<mris_out->nvertices; m++) {
      mris_out->vertices[m].origx = 0;
      mris_out->vertices[m].origy = 0;
      mris_out->vertices[m].origz = 0;
    }

    /*Initialize Ij,k*/
    for (vno = 0; vno<mris_out->nvertices; vno++) {
      vm_out = &mris_out->vertices[vno];
      vm_out->val = 1;
    }

    /*Iteratively compute Ij,k*/
    for (i=order;i>0;i--) {
      mris_high = ReadIcoByOrder(i, 100); //higher order surface
      for (m = 0; m<mris_high->nvertices; m++)
        mris_high->vertices[m].nsize=1;
      MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
      number = IcoNVtxsFromOrder(i-1); //the start of m vertices
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        flag=0;
        for (nnum=0; nnum<vm_high->vnum; nnum++)
          if ( vm_high->v[nnum]<number ) //A(j,m)
          {
            k = vm_high->v[nnum];
            v = &mris_out->vertices[k];
            v->val += 0.5*vm_out->val ;
          }
        for (; nnum<vm_high->v2num; nnum++)
          if ( vm_high->v[nnum]<number ) //B(j,m)
          {
            k = vm_high->v[nnum];
            if (flag==0) b1=k;
            else b2=k;
            flag++;
            v = &mris_out->vertices[k];
            v->val += 0.125*vm_out->val ;
          }
        for (; nnum<vm_high->v3num; nnum++)
          if ( vm_high->v[nnum]<number ) //C(j,m)
          {
            k = vm_high->v[nnum];
            flag=0; //C has to be a second-order neighbor of B
            for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
              if (mris_high->vertices[b1].v[cno]==k) flag=1;
            for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
              if (mris_high->vertices[b2].v[cno]==k) flag=1;
            if (flag) {
              v = &mris_out->vertices[k];
              v->val -= 0.0625*vm_out->val ;
            }
          }
      }
    }


    for (i=1;i<=order;i++) {
      mris_high = ReadIcoByOrder(i, 100); //higher order surface
      for (m = 0; m<mris_high->nvertices; m++)
        mris_high->vertices[m].nsize=1;
      MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
      number = IcoNVtxsFromOrder(i-1); //the start of m vertices

      /* Synthesis Stage I */
      /* Compute Lamda(j+1,k) using the Yita(j,m) */
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        for (nnum=0; nnum<vm_high->vnum; nnum++)
          if ( vm_high->v[nnum]<number ) //A(j,m)
          {
            k = vm_high->v[nnum];
            v = &mris_out->vertices[k];
            s_jkm = vm_out->val/2/v->val;
            v->origx -= s_jkm*vm_out->origx;
            v->origy -= s_jkm*vm_out->origy;
            v->origz -= s_jkm*vm_out->origz;
          }
      }

      /* compute Lamda(j+1,m) for each m vertices */
      for (m = number; m<mris_high->nvertices; m++) {
        vm_out = &mris_out->vertices[m];
        vm_high = &mris_high->vertices[m];
        flag=0;
        for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
          {
            k = vm_high->v[nnum] ;
            v = &mris_out->vertices[k];
            vm_out->origx += 0.5*v->origx;
            vm_out->origy += 0.5*v->origy;
            vm_out->origz += 0.5*v->origz;
          }
        for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
          {
            k = vm_high->v[nnum] ;
            if (flag==0) b1=k;
            else b2=k;
            flag++;
            v = &mris_out->vertices[k];
            vm_out->origx += 0.125*v->origx;
            vm_out->origy += 0.125*v->origy;
            vm_out->origz += 0.125*v->origz;
          }
        for (; nnum<vm_high->v3num; nnum++) //third order neighborhood
          if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
          {
            k = vm_high->v[nnum] ;
            flag=0; //C has to be a second-order neighbor of B
            for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
              if (mris_high->vertices[b1].v[cno]==k) flag=1;
            for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
              if (mris_high->vertices[b2].v[cno]==k) flag=1;
            if (flag) {
              v = &mris_out->vertices[k];
              vm_out->origx -= 0.0625*v->origx;
              vm_out->origy -= 0.0625*v->origy;
              vm_out->origz -= 0.0625*v->origz;
            }
          }
      }
    }

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    fprintf(stdout, "Writing recovered surface to %s\n", argv[4]);
    MRISwrite(mris_out, argv[4]) ;
#if 0
    mris_high = ReadIcoByOrder(4, 100);
    for (m=0;m<mris_high->nvertices;m++) {
      mris_high->vertices[m].x=mris_out->vertices[m].x;
      mris_high->vertices[m].y=mris_out->vertices[m].y;
      mris_high->vertices[m].z=mris_out->vertices[m].z;
    }
    MRISwrite(mris_high, "/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.wavelet.recon") ;
#endif
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    /*End of Synthesis*/
  }

  MRISfree(&mris_out);
  MRISfree(&mris_high) ;
  msec = then.milliseconds() ;
  fprintf(stdout, "spherical wavelet took %2.1f minutes\n", (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
}

/*----------------------------------------------------------------------

 Parameters:

 Description:
 ----------------------------------------------------------------------*/
MRI_SURFACE *
center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) {
  int         fno, vno ;
  FACE        *face ;
  VERTEX      *vdst;
  float       x, y, z, x0, y0, z0 ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = y0 = z0 = 0 ;   /* silly compiler warning */

  for (fno = 0 ; fno < mris_src->nfaces ; fno++) {
    face = &mris_dst->faces[fno] ;
    if (face->ripflag)
      continue ;
    x = mris_dst->vertices[face->v[0]].x;
    y = mris_dst->vertices[face->v[0]].y;
    z = mris_dst->vertices[face->v[0]].z;
    x += mris_dst->vertices[face->v[1]].x;
    y += mris_dst->vertices[face->v[1]].y;
    z += mris_dst->vertices[face->v[1]].z;
    x += mris_dst->vertices[face->v[2]].x;
    y += mris_dst->vertices[face->v[2]].y;
    z += mris_dst->vertices[face->v[2]].z;
    x /= face->area;
    y/= face->area;
    z/= face->area;
    x0 += x;
    y0 += y;
    z0 += z;
  }
  x0 /= mris_dst->total_area ;
  y0 /= mris_dst->total_area  ;
  z0 /= mris_dst->total_area ;

  for (vno = 0 ; vno < mris_src->nvertices ; vno++) {
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    vdst->x -= x0 ;
    vdst->y -= y0 ;
    vdst->z -= z0 ;
  }

  mris_dst->xctr = mris_dst->yctr = mris_dst->zctr = 0 ;
  return(mris_dst) ;
}


static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "A")) {
    ANALYSIS = 1 ;
    fprintf(stdout,"spherical wavelet analysis\n");
  } else if (!stricmp(option, "S")) {
    SYNTHESIS = 1 ;
    fprintf(stdout,"Spherical wavelt synthesis\n");
  } else if (!stricmp(option, "C")) {
    COMPARE = 1;
    fname = argv[2] ;
    fprintf(stdout,"Reconstruct the surface using only differences\n");
    nargs = 1;
  } else if (!stricmp(option, "CURV")) {
    CURV = 1;
    fprintf(stdout,"Decompose scalar function\n");
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      fprintf(stdout,
              "usage: %s <input volumes> <output volume>\n",
              Progname) ;
      exit(1) ;
      break ;
    default:
      fprintf(stdout, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}






