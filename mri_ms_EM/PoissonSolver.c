/**
 * @file  PoissonSolver.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:23 $
 *    $Revision: 1.3 $
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


/* Solve Poisson Equation using Multigrid techniques.
 * The equation to solve is given by
 * div( g(x) \nabla u(x)) - h(x) u = f(x), which
 * is the GGVF equation
 * For convenience, in the program, g(x) is called weight,
 * h is called hu
 */
/* All the volumes here are assumed to be FLOAT */
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "PoissonSolver.h"
#define JACOBI 0

static int maxlevel;
static int n1, n2;

void rstrct(MRI *out, MRI *in, int Xo, int Yo, int Zo,
            int Xi, int Yi, int Zi) {
  int inew, jnew, knew, iold, jold, kold;

  for (knew = 0; knew < Zo; knew++)
    for (inew = 0; inew < Yo; inew++)
      for (jnew =0; jnew < Xo; jnew++) {
        kold = knew<<1;
        iold = inew <<1;
        jold = jnew<<1;
        if (kold > (Zi-2)) kold = Zi-2;
        if (iold > (Yi-2)) iold = Yi-2;
        if (jold > (Xi-2)) jold = Xi-2;

        MRIFvox(out, jnew, inew, knew) = 0.125*(MRIFvox(in, jold, iold, kold)+ MRIFvox(in,jold+1, iold, kold) + MRIFvox(in, jold, iold+1, kold) + MRIFvox(in, jold+1, iold+1, kold) + MRIFvox(in, jold, iold, kold+1) + MRIFvox(in, jold + 1, iold, kold+1) + MRIFvox(in, jold, iold+1, kold+1) + MRIFvox(in, jold+1, iold+1, kold+1));

      }


  return;
}

void interp(MRI *out, MRI *in, int Xo, int Yo, int Zo,
            int Xi, int Yi, int Zi) {

  int ic,jc,kc, i,j,k;
  int HX, HY, HZ;

  for (kc=0; kc < Zo; kc++)
    for (ic=0; ic < Yo; ic++)
      for (jc=0; jc < Xo; jc++) {
        k = kc>>1;
        j = jc >> 1;
        i = ic >>1;
        HX = (j == (Xi-1));
        HY = (i == (Yi-1));
        HZ = (k==(Zi-1));
        MRIFvox(out, jc, ic, kc) = MRIFvox(in, j, i, k); /* 0.5*(in[k][i][j] + in[k+1-HZ][i+1-HY][j+1-HX]); */

      }
  return;
}

void addin(MRI *uf, MRI *uc, MRI *res, int Xo, int Yo, int Zo, int Xi, int Yi, int Zi) {
  int i,j,k;

  interp(res, uc, Xo, Yo,Zo, Xi,Yi,Zi);

  for (k=0; k < Zo; k++) {
    for (i=0; i<Yo; i++) {
      for (j=0; j< Xo; j++) {
        MRIFvox(uf, j, i, k)  += MRIFvox(res, j, i, k);
      }
    }
  }

  return;
}

void resid(MRI *res, MRI *u, MRI *rhs, MRI *hu, int XN, int YN, int ZN, int lev) {
  int i,j,k;
  int h;
  float cmu;

  int prek, nextk, prei, nexti, prej, nextj;

  h = 1<<(maxlevel-lev);
  h = h*h;
  cmu = (float)1.0/(float)h;

  /* If change to half-point symmetric extension, no longer converges */
  for (k=0;k<ZN; k++) {
    prek = (k>0) ? (k-1) : 0;
    nextk = (k == (ZN-1)) ? (ZN-1) : (k+1);
    for (i=0; i<YN;i++) {
      prei = (i>0) ? (i-1) : 0;
      nexti = (i == (YN-1)) ? (YN-1) : (i+1);

      for (j=0;j<XN;j++) {
        prej = (j>0) ? (j-1) : 0;
        nextj = (j == (XN-1)) ? (XN-1) : (j+1);

        MRIFvox(res, j, i, k) = MRIFvox(hu, j, i, k) * MRIFvox(u, j, i, k) + MRIFvox(rhs, j, i, k)  -
                                cmu*
                                (MRIFvox(u, j, i, prek) +
                                 MRIFvox(u, j, i, nextk) +
                                 MRIFvox(u, j, prei, k) +
                                 MRIFvox(u, j, nexti, k) +
                                 MRIFvox(u, prej, i, k) +
                                 MRIFvox(u, nextj, i, k) -
                                 MRIFvox(u, j, i, k)*6);

        // if( j== (XN >> 1) && (i == (YN >> 1)) && (k == (ZN >> 1)))
        //  printf("residue = %g\n", MRIFvox(res, j, i, k));
      }
    }
  }

  return;

}

void copymem(MRI *out, MRI *in, int XN, int YN, int ZN) {
  int k,i,j;

  for (k=0; k<ZN; k++)
    for (i=0;i <YN; i++)
      for (j=0; j<XN; j++) {
        MRIFvox(out, j, i, k) = MRIFvox(in, j, i, k);
      }



  return;
}

void octant(MRI *u, MRI *rhs, MRI *hu, float hsquare, int XN, int YN, int ZN, int kc, int ic, int jc) {

  int k, i, j;
  int prek, nextk, prei,nexti, prej, nextj;
  float tmpv;

  /* hsquare = 2*hsquare; For the first case that g is inside div() */

  for (k=kc; k<ZN; k += 2 ) {
    prek = (k>0) ? (k-1) : 0;
    nextk = (k == (ZN-1)) ? (ZN-1) : (k+1);
    for (i=ic;i<YN; i+= 2) {
      prei = (i>0) ? (i-1) : 0;
      nexti = (i == (YN-1)) ? (YN-1) : (i+1);

      for (j=jc;j<XN; j+= 2) {
        prej = (j>0) ? (j-1) : 0;
        nextj = (j == (XN-1)) ? (XN-1) : (j+1);

        tmpv = MRIFvox(hu, j, i, k)*hsquare;

        MRIFvox(u, j, i, k) = (( MRIFvox(u, j, i, prek) +
                                 MRIFvox(u, j, i, nextk) +
                                 MRIFvox(u, j, prei, k) +
                                 MRIFvox(u, j, nexti, k) +
                                 MRIFvox(u, prej, i, k) +
                                 MRIFvox(u, nextj, i, k))
                               - MRIFvox(rhs, j, i, k)*hsquare)/(6 + tmpv);

        // if( j== (XN >> 1) && (i == (YN >> 1)) && (k == (ZN >> 1)))
        //  printf("u = %g\n", MRIFvox(u, j, i, k));
      }
    }
  }

  return;
}

void Jacobi(MRI *u, MRI *rhs, MRI *hu, int lev, int XN, int YN, int ZN, int iter) {

  MRI *tmpvol;

  int k, i, j, index, h;
  int prek, nextk, prei,nexti, prej, nextj;

  float alpha = 0.9; /* 0.66666666666; alpha has to be large to converge in 3D */
  float hsquare;

  float tmpv;

  h = 1 << (maxlevel - lev);
  hsquare = 2*h*h;

  tmpvol = MRIclone(u, NULL);

  for (index =1; index <= iter; index++) {
    for (k=0; k<ZN; k ++ ) {
      prek = (k>0) ? (k-1) : 0;
      nextk = (k == (ZN-1)) ? (ZN-1) : (k+1);
      for (i=0;i<YN; i++) {
        prei = (i>0) ? (i-1) : 0;
        nexti = (i == (YN-1)) ? (YN-1) : (i+1);

        for (j=0;j<XN; j++) {
          prej = (j>0) ? (j-1) : 0;
          nextj = (j == (XN-1)) ? (XN-1) : (j+1);

          tmpv = MRIFvox(hu, j, i,k)*hsquare;

          MRIFvox(tmpvol, j, i, k) = (1-alpha)*MRIFvox(u, j, i, k) +
                                     alpha*( MRIFvox(u, j, i, prek) +
                                             MRIFvox(u, j, i, nextk) +
                                             MRIFvox(u, j, prei, k) +
                                             MRIFvox(u, j, nexti, k) +
                                             MRIFvox(u, prej, i, k) +
                                             MRIFvox(u, nextj, i, k)
                                             - MRIFvox(rhs, j, i, k)*hsquare)/(12 + tmpv);

        }
      }
    }

    MRIcopy(tmpvol, u);

    //    for(k=0;k<ZN;k++)
    //  for(i=0;i<YN;i++)
    // for(j=0;j<XN;j++){
    //   MRIFvox(u, j, i, k) = MRIFvox(tmpvol, j, i, k);
    // }

    /*
    memcpy((char *)&u[0][0][0], (char *)&tmpvol[0][0][0], size);
    */

    /*
    sum = 0.0;
    for(k=0;k<ZN;k++)
    for(i=0;i<YN;i++)
    for(j=0;j<XN;j++){
    sum += u[k][i][j] ;
    }

    sum /= (float)(XN*YN*ZN);

    for(k=0;k<ZN;k++)
    for(i=0;i<YN;i++)
    for(j=0;j<XN;j++){
    u[k][i][j] -= sum;
    }

    */
  }

  MRIfree(&tmpvol);

  return;

}

void Gauss_Seidel(MRI *u, MRI *rhs, MRI *hu, int lev, int XN, int YN, int ZN, int iter) {
  int h, i;
  float hsquare;

  h = 1 << (maxlevel - lev);
  hsquare = h*h;

  for (i=1; i<= iter;i++) {
    /* Process black first like follows gives slightly better convergence */

    octant(u,rhs,hu,hsquare, XN, YN, ZN, 0, 0, 0);
    octant(u,rhs,hu,hsquare, XN, YN, ZN, 0, 1, 1);
    octant(u,rhs,hu,hsquare, XN, YN, ZN, 1, 0, 1);
    octant(u,rhs,hu,hsquare, XN, YN, ZN, 1, 1, 0);

    octant(u,rhs,hu,hsquare, XN, YN, ZN, 1, 0, 0);
    octant(u,rhs,hu,hsquare, XN, YN, ZN, 0, 0, 1);
    octant(u,rhs,hu,hsquare, XN, YN, ZN, 0, 1, 0);
    octant(u,rhs,hu,hsquare, XN, YN, ZN, 1, 1, 1);
  }

  return;
}



void slvsml(MRI *u, MRI *rhs, MRI *hu, int XN, int YN, int ZN) {
  int i,j,k;
  //  float sum;

  for (k=0;k<ZN;k++)
    for (i=0;i<YN;i++)
      for (j=0;j<XN;j++) {
        MRIFvox(u, j, i, k) = 0;
      }




  if (JACOBI)
    Jacobi(u,rhs,hu,1,XN,YN,ZN,2);
  else
    Gauss_Seidel(u,rhs,hu,1,XN,YN,ZN,2);

  return;
}

void multigrid(MRI *u, MRI *f, MRI *hu, int XN, int YN, int ZN) {
  /* u is the solution for
   * div( \nabla u(x)) - hu(x)*u(x) - f(x) = 0;
   */

  int m, n, l, mf;
  int iters;
  int *mo, *no, *lo;

  int cx = 0, cy= 0, cz =0;
  int j, jcycle, jj;
  int msize, nsize, lsize;

  float maxerr=10, tmpv;

  MRI **ires, **irho, **irhs, **iu, **ih;

  MRI *uold;

  /* Computer maximum level */
  m = 0;
  mf = ZN;
  while (mf >>= 1) m++;
  maxlevel = m;
  m = 0;
  mf = YN;
  while (mf >>= 1) m++;
  if (maxlevel > m) maxlevel = m;
  m = 0;
  mf = XN;
  while (mf >>= 1) m++;
  if (maxlevel > m) maxlevel = m;

  printf("Simple scheme: maxlevel = %d \n", maxlevel);

  /* Set pre-smoothing steps */
  n1 = 1; //2;
  /* Set post-smoothing steps */
  n2 = 2; //2

  /* Allocate memory */
  ires = (MRI **)malloc((maxlevel+1)*sizeof(MRI *));
  irho = (MRI **)malloc((maxlevel+1)*sizeof(MRI *));
  irhs = (MRI **)malloc((maxlevel+1)*sizeof(MRI *));
  iu = (MRI **)malloc((maxlevel+1)*sizeof(MRI *));
  ih = (MRI **)malloc((maxlevel+1)*sizeof(MRI *));
  mo = (int *)malloc((maxlevel+1)*sizeof(int));
  no = (int *)malloc((maxlevel+1)*sizeof(int));
  lo = (int *)malloc((maxlevel+1)*sizeof(int));

  uold = MRIclone(u, NULL);
  msize = ZN;
  nsize = YN;
  lsize = XN;


  for (j=maxlevel;j>=1;j--) {
    ires[j] = MRIalloc(lsize, nsize, msize, MRI_FLOAT);
    irho[j] = MRIalloc(lsize, nsize, msize, MRI_FLOAT);
    irhs[j] = MRIalloc(lsize, nsize, msize, MRI_FLOAT);
    iu[j] = MRIalloc(lsize, nsize, msize, MRI_FLOAT);
    ih[j] = MRIalloc(lsize, nsize, msize, MRI_FLOAT);

    mo[j] = msize;
    no[j] = nsize;
    lo[j] = lsize;
    msize = (msize+1)>>1;
    nsize = (nsize + 1)>>1;
    lsize = (lsize+1)>>1;
  }

  printf("n1=%d, n2=%d, x1=%d,y1=%d,z1=%d\n",n1,n2,lo[1],no[1],mo[1]);

  for (m=0;m<ZN;m++)
    for (n=0;n<YN;n++)
      for (l=0;l<XN;l++) {
        MRIFvox(ih[maxlevel], l, n, m) = MRIFvox(hu, l, n,  m);
        MRIFvox(uold, l, n, m) = 0.0;
      }

  for (j = maxlevel-1; j >= 1; j--) {
    rstrct(ih[j],ih[j+1],lo[j], no[j],mo[j],lo[j+1], no[j+1], mo[j+1]);
  }


  /* Now start multigrid processing */
  for (iters = 1; iters <= 10; iters++) {

    /* Several iterations of fmgv */
    /* Compute the initial residue */
    resid(irho[maxlevel], u, f, ih[maxlevel], XN,YN,ZN, maxlevel);
#if 1
    maxerr = 0.0;
    for (m=0;m<ZN;m++)
      for (n=0;n<YN;n++)
        for (l=0;l<XN;l++) {
          tmpv = MRIFvox(irho[maxlevel], l, n, m);
          if (tmpv < 0) tmpv = 0-tmpv;
          if (maxerr < tmpv) {
            maxerr = tmpv;
            cx = l;
            cy = n;
            cz = m;
          }
        }


    printf("Maximum Residue (at (%d, %d, %d)) is %g \n", cx, cy, cz, maxerr);
    //    if(maxerr < 0.01 && iters > 1) break;
    if (maxerr < 0.005) break;
#endif

    for (j = maxlevel-1; j >= 1; j--) {
      rstrct(irho[j],irho[j+1],lo[j],no[j],mo[j],lo[j+1],no[j+1], mo[j+1]);
    }

    /*Now solve the PDE at level 1, the coarsest level */
    slvsml(iu[1], irho[1], ih[1], lo[1], no[1], mo[1]); /* Initial solution on coarest grid */


    /* Now start Full-Multigrid V or W cycle */
    for (j=2; j<= maxlevel;j++) {

      interp(iu[j],iu[j-1], lo[j], no[j], mo[j], lo[j-1], no[j-1],mo[j-1]); /* Get the initial guess of the solution of the original eq */

      copymem(irhs[j], irho[j],lo[j], no[j],mo[j]); /* Set up right hand side */
      /* Now begin mgv.m 12-9-03 */
      for (jcycle = 1; jcycle <= 1; jcycle ++) {
        for (jj=j; jj>=2; jj--) { /* Down stroke of the V */
          if (JACOBI)
            Jacobi(iu[jj], irhs[jj],  ih[jj], jj, lo[jj],no[jj],mo[jj], n1);
          else
            Gauss_Seidel(iu[jj], irhs[jj],  ih[jj], jj, lo[jj],no[jj],mo[jj], n1);

          resid(ires[jj], iu[jj], irhs[jj],  ih[jj],lo[jj],no[jj],mo[jj],jj); /* Defect */

          /* Does the following initialization only need to be done at jj=2
           * or doesn't need to be done at all ?? No! It need to be done
           * at all level, otherwise, the addin need to be changed!
           * Why not get rid of += in addin, and get rid of all these
           * initialization?? Am I missing something?? Indeed, no change
           * should be made here!
           */
          for (m=0; m < mo[jj-1]; m++) /*loop over image */
            for (n=0;n< no[jj-1];n++)
              for (l=0;l<lo[jj-1];l++) {
                MRIFvox(iu[jj-1], l, n, m) = 0; /* Initial value for errors at levels from j-1 to 1 which are iu[j-1] ~ iu[1] and are zeros */
              }

          rstrct(irhs[jj-1], ires[jj], lo[jj-1], no[jj-1],mo[jj-1],lo[jj],no[jj],mo[jj]); /* mf, nf is the size of the first parameter */

        }


        slvsml(iu[1], irhs[1], ih[1],lo[1],no[1],mo[1]); /* Bottom of the V */
        /* Now iu[1] stores the solution of the error at the coarest level */

        for (jj=2; jj <=j; jj++) { /*Upard stroke of V */
          addin(iu[jj], iu[jj-1], ires[jj],lo[jj],no[jj],mo[jj],lo[jj-1],no[jj-1],mo[jj-1]);
          /*ires[jj] is used for temporary storage inside addint */

          /* Post-smooyhing */
          if (JACOBI)
            Jacobi(iu[jj], irhs[jj], ih[jj], jj, lo[jj],no[jj],mo[jj], n1);
          else
            Gauss_Seidel(iu[jj], irhs[jj], ih[jj], jj,lo[jj],no[jj],mo[jj], n2);
        }

      }

    }

    /*
    maxerr = 0.0;
    for(m=0;m<=(ZN-1);m++)
      for(n=0;n<=(YN-1);n++)
    for(l=0;l<=(XN-1);l++){
    tmpv = iu[maxlevel][m][n][l];
    if(tmpv < 0) tmpv = 0-tmpv;
    if(maxerr < tmpv) maxerr = tmpv;
    }

    if(maxerr < 0.01){
      printf("Converged after %d iterations\n", iters);
      break;
    }else{
      printf("Error after %d multigrid iters is %g \n", iters, maxerr);
    }
    */

    /* Update solution */
    maxerr = 0;
    for (m=0; m < ZN; m++)
      for (n=0;n< YN;n++)
        for (l=0;l<XN; l++) {
          tmpv =  MRIFvox(iu[maxlevel], l, n, m);
          MRIFvox(u, l, n, m) += tmpv;
          if (tmpv < 0) tmpv = 0-tmpv;
          if (maxerr < tmpv) maxerr = tmpv;
        }
#if 0
    if (maxerr < 0.01) {
      printf("Converged after %d iterations\n", iters);
      break;
    } else {
      printf("Error after %d multigrid iters is %g \n", iters, maxerr);
    }
#endif
  } /* endof for (iters) */

  MRIfree(&uold);

  for (j=maxlevel;j>=1;j--) {
    msize = mo[j];
    nsize = no[j];
    lsize = lo[j];
    MRIfree(&ires[j]);
    MRIfree(&irhs[j]);
    MRIfree(&irho[j]);
    MRIfree(&iu[j]);
    MRIfree(&ih[j]);
  }

  free(ih);
  free(iu);
  free(ires);
  free(irho);
  free(irhs);


  free(mo);
  free(no);
  free(lo);

  return;
}

