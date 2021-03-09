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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

static int gNb(int *tab) {
  int nb,n;
  nb=0;
  n=0;
  while (tab[n]>=0) {
    switch (tab[n]) {
    case 0:
      switch (tab[n+1]) {
      case 0:
        if (tab[n+2]==0)
          nb+=1;
        else
          nb+=16;
        break;
      case 1:
        if (tab[n+2]==0)
          nb+=4;
        else
          nb+=64;
        break;
      }
      break;
    case 1:
      switch (tab[n+1]) {
      case 0:
        if (tab[n+2]==0)
          nb+=2;
        else
          nb+=32;
        break;
      case 1:
        if (tab[n+2]==0)
          nb+=8;
        else
          nb+=128;
        break;
      }
      break;
    }
    n=n+3;
  }
  return nb;
}


static int Reverse(int *f,int *v,int *fd, int *vd) {
  int n,p,test,k;
  int tb[8][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0},
                {0,0,1},{1,0,1},{0,1,1},{1,1,1}};
  n=0;
  while (f[n]>=0) {
    fd[n]=f[n+2];
    fd[n+1]=f[n+1];
    fd[n+2]=f[n];
    n+=3;
  }

  n=0;
  for (k=0;k<8;k++) {
    p=0;
    test=0;
    while (v[p]>=0) {
      if (v[p]==tb[k][0]&&v[p+1]==tb[k][1]&&v[p+2]==tb[k][2])
        test=1;
      p+=3;
    }
    if (!test) {
      vd[n]=tb[k][0];
      vd[n+1]=tb[k][1];
      vd[n+2]=tb[k][2];
      n+=3;
    }
  }
  return gNb(vd);
}


static const int RfOx[12]= {
                             2,6,10,7,1,3,9,11,0,4,8,5
                           };
static int ROx(int *f,int *v,int *fd, int *vd) {
  int n;
  n=0;
  while (f[n]>=0) {
    fd[n]=RfOx[f[n]];
    fd[n+1]=RfOx[f[n+1]];
    fd[n+2]=RfOx[f[n+2]];
    n+=3;
  }
  n=0;
  while (v[n]>=0) {
    if (v[n+1]==0&&v[n+2]==0) {
      vd[n+1]=1;
      vd[n+2]=0;
    } else if (v[n+1]==1&&v[n+2]==0) {
      vd[n+1]=1;
      vd[n+2]=1;
    } else if (v[n+1]==1&&v[n+2]==1) {
      vd[n+1]=0;
      vd[n+2]=1;
    } else {
      vd[n+1]=0;
      vd[n+2]=0;
    }
    vd[n]=v[n];
    n+=3;
  }
  return gNb(vd);
}
static const int RfOy[12]= {
                             4,9,6,1,8,0,10,2,5,11,7,3
                           };
static int ROy(int *f,int *v,int *fd, int *vd) {
  int n;
  n=0;
  while (f[n]>=0) {
    fd[n]=RfOy[f[n]];
    fd[n+1]=RfOy[f[n+1]];
    fd[n+2]=RfOy[f[n+2]];
    n+=3;
  }
  n=0;
  while (v[n]>=0) {
    if (v[n]==0&&v[n+2]==0) {
      vd[n]=0;
      vd[n+2]=1;
    } else if (v[n]==0&&v[n+2]==1) {
      vd[n]=1;
      vd[n+2]=1;
    } else if (v[n]==1&&v[n+2]==1) {
      vd[n]=1;
      vd[n+2]=0;
    } else {
      vd[n]=0;
      vd[n+2]=0;
    }
    vd[n+1]=v[n+1];
    n+=3;
  }
  return gNb(vd);
}

static const int RfOz[12]= {
                             3,0,1,2,5,7,4,6,11,8,9,10
                           };
static int ROz(int *f,int *v,int *fd, int *vd) {
  int n;
  n=0;
  while (f[n]>=0) {
    fd[n]=RfOz[f[n]];
    fd[n+1]=RfOz[f[n+1]];
    fd[n+2]=RfOz[f[n+2]];
    n+=3;
  }
  n=0;
  while (v[n]>=0) {
    if (v[n]==0&&v[n+1]==0) {
      vd[n]=1;
      vd[n+1]=0;
    } else if (v[n]==1&&v[n+1]==0) {
      vd[n]=1;
      vd[n+1]=1;
    } else if (v[n]==1&&v[n+1]==1) {
      vd[n]=0;
      vd[n+1]=1;
    } else {
      vd[n]=0;
      vd[n+1]=0;
    }
    vd[n+2]=v[n+2];
    n+=3;
  }
  return gNb(vd);
}

static const int SmfOx[12]= {
                              0,3,2,1,5,4,7,6,8,11,9,10
                            };
static const int SmfOy[12]= {
                              2,1,0,3,6,7,4,5,10,9,8,11
                            };
static const int SmfOz[12]= {
                              8,9,10,11,4,5,6,7,0,1,2,3
                            };
#if 0
static int SmOx(int *f,int *v,int *fd, int *vd) {
  int n;
  n=0;
  while (f[n]>=0) //take into account the permutation
  {
    fd[n]=SmfOx[f[n+2]];
    fd[n+1]=SmfOx[f[n+1]];
    fd[n+2]=SmfOx[f[n]];
    n+=3;
  }
  n=0;
  while (v[n]>=0) {
    if (v[n]==0)
      vd[n]=1;
    else
      vd[n]=0;
    vd[n+1]=v[n+1];
    vd[n+2]=v[n+2];
    n+=3;
  }
  return gNb(vd);
}

static int SmOy(int *f,int *v,int *fd, int *vd) {
  int n;
  n=0;
  while (f[n]>=0) //take into account the permutation
  {
    fd[n]=SmfOy[f[n+2]];
    fd[n+1]=SmfOy[f[n+1]];
    fd[n+2]=SmfOy[f[n]];
    n+=3;
  }
  n=0;
  while (v[n]>=0) {
    if (v[n+1]==0)
      vd[n+1]=1;
    else
      vd[n+1]=0;
    vd[n]=v[n];
    vd[n+2]=v[n+2];
    n+=3;
  }
  return gNb(vd);
}
static int SmOz(int *f,int *v,int *fd, int *vd) {
  int n;
  n=0;
  while (f[n]>=0) //take into account the permutation
  {
    fd[n]=SmfOz[f[n+2]];
    fd[n+1]=SmfOz[f[n+1]];
    fd[n+2]=SmfOz[f[n]];
    n+=3;
  }
  n=0;
  while (v[n]>=0) {
    if (v[n+2]==0)
      vd[n+2]=1;
    else
      vd[n+2]=0;
    vd[n]=v[n];
    vd[n+1]=v[n+1];
    n+=3;
  }
  return gNb(vd);
}
#endif


static void reinit(int *tab,int nb) {
  int n;
  for (n=0;n<nb;n++)
    tab[n]=-1;
}


static void cp(int *s,int *d, int nb) {
  int n;
  for (n=0;n<nb;n++)
    d[n]=s[n];
}

static void copy(int *s,int *d) {
  int n;
  n=0;
  if (d[n]>=0)
    return;
  while (s[n]>=0) {
    d[n]=s[n];
    n++;
  }
}

#define NB_TRIANGLES 6
#define NB_VERTICES (3*NB_TRIANGLES +1)
#define NV NB_VERTICES

#define NB_POINTS 25
#define NP NB_POINTS

static int Case6[15][NV],Case18[15][NV],Case26[15][NV];
static int Pt[15][NP];

static void initCases(void) {
  int k;

  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //// 6-connectivity /////////////////////////////////
  /////////////////////////////////////////////////////

  for (k=0;k<15;k++)
    reinit(Case6[k],NV);

  //case 0
  //nothing to do

  //case 1: 1 face
  k=1;
  Case6[k][0]=4;
  Case6[k][1]=0;
  Case6[k][2]=1;


  //case 2: 2 faces
  k=2;
  Case6[k][0]=4;
  Case6[k][1]=5;
  Case6[k][2]=1;
  Case6[k][3]=1;
  Case6[k][4]=5;
  Case6[k][5]=3;

  //case 3: 2 faces (6-connectivity)
  k=3;
  Case6[k][0]=4;
  Case6[k][1]=0;
  Case6[k][2]=1;
  Case6[k][3]=5;
  Case6[k][4]=8;
  Case6[k][5]=11;


  //case 4: 2 faces (6,18-connectivity)
  k=4;
  Case6[k][0]=4;
  Case6[k][1]=0;
  Case6[k][2]=1;
  Case6[k][3]=7;
  Case6[k][4]=11;
  Case6[k][5]=10;


  //case 5: 3 faces
  k=5;
  Case6[k][0]=5;
  Case6[k][1]=7;
  Case6[k][2]=6;
  Case6[k][3]=5;
  Case6[k][4]=6;
  Case6[k][5]=1;
  Case6[k][6]=5;
  Case6[k][7]=1;
  Case6[k][8]=0;

  //case 6: 3 faces (6-connectivity)
  k=6;
  Case6[k][0]=4;
  Case6[k][1]=5;
  Case6[k][2]=1;
  Case6[k][3]=1;
  Case6[k][4]=5;
  Case6[k][5]=3;
  Case6[k][6]=7;
  Case6[k][7]=11;
  Case6[k][8]=10;

  //case 7: 3 faces (6-connectivity)
  k=7;
  Case6[k][0]=7;
  Case6[k][1]=11;
  Case6[k][2]=10;
  Case6[k][3]=4;
  Case6[k][4]=9;
  Case6[k][5]=8;
  Case6[k][6]=0;
  Case6[k][7]=5;
  Case6[k][8]=3;

  //case 8: 2 faces
  k=8;
  Case6[k][0]=4;
  Case6[k][1]=5;
  Case6[k][2]=6;
  Case6[k][3]=6;
  Case6[k][4]=5;
  Case6[k][5]=7;

  //case 9: 4 faces
  k=9;
  Case6[k][0]=4;
  Case6[k][1]=0;
  Case6[k][2]=9;
  Case6[k][3]=0;
  Case6[k][4]=10;
  Case6[k][5]=9;
  Case6[k][6]=0;
  Case6[k][7]=3;
  Case6[k][8]=10;
  Case6[k][9]=10;
  Case6[k][10]=3;
  Case6[k][11]=7;

  //case 10: 4 faces (6-connectivity)
  k=10;
  Case6[k][0]=0;
  Case6[k][1]=9;
  Case6[k][2]=8;
  Case6[k][3]=0;
  Case6[k][4]=1;
  Case6[k][5]=9;
  Case6[k][6]=3;
  Case6[k][7]=11;
  Case6[k][8]=10;
  Case6[k][9]=3;
  Case6[k][10]=10;
  Case6[k][11]=2;


  //case 11: 4 faces
  k=11;
  Case6[k][0]=4;
  Case6[k][1]=0;
  Case6[k][2]=6;
  Case6[k][3]=0;
  Case6[k][4]=3;
  Case6[k][5]=6;
  Case6[k][6]=6;
  Case6[k][7]=3;
  Case6[k][8]=10;
  Case6[k][9]=10;
  Case6[k][10]=3;
  Case6[k][11]=11;


  //case 12: 4 faces (6-connectivity)
  k=12;
  Case6[k][0]=4;
  Case6[k][1]=9;
  Case6[k][2]=8;
  Case6[k][3]=0;
  Case6[k][4]=6;
  Case6[k][5]=1;
  Case6[k][6]=0;
  Case6[k][7]=5;
  Case6[k][8]=6;
  Case6[k][9]=6;
  Case6[k][10]=5;
  Case6[k][11]=7;

  //case 13: 4 faces (6-connectivity)
  k=13;
  Case6[k][0]=0;
  Case6[k][1]=1;
  Case6[k][2]=4;
  Case6[k][3]=3;
  Case6[k][4]=7;
  Case6[k][5]=2;
  Case6[k][6]=5;
  Case6[k][7]=8;
  Case6[k][8]=11;
  Case6[k][9]=9;
  Case6[k][10]=6;
  Case6[k][11]=10;


  //case 14: 4 faces
  k=14;
  Case6[k][0]=0;
  Case6[k][1]=9;
  Case6[k][2]=1;
  Case6[k][3]=0;
  Case6[k][4]=5;
  Case6[k][5]=9;
  Case6[k][6]=5;
  Case6[k][7]=10;
  Case6[k][8]=9;
  Case6[k][9]=5;
  Case6[k][10]=7;
  Case6[k][11]=10;


  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //// 18-connectivity /////////////////////////////////
  /////////////////////////////////////////////////////

  for (k=0;k<15;k++)
    reinit(Case18[k],NV);

  //case 0
  //nothing to do

  //case 1: 1 face
  k=1;
  Case18[k][0]=4;
  Case18[k][1]=0;
  Case18[k][2]=1;


  //case 2: 2 faces
  k=2;
  Case18[k][0]=4;
  Case18[k][1]=5;
  Case18[k][2]=1;
  Case18[k][3]=1;
  Case18[k][4]=5;
  Case18[k][5]=3;

  //case 3: 4 faces (18-26 connectivity)
  k=3;
  Case18[k][0]=1;
  Case18[k][1]=4;
  Case18[k][2]=8;
  Case18[k][3]=1;
  Case18[k][4]=8;
  Case18[k][5]=11;
  Case18[k][6]=0;
  Case18[k][7]=1;
  Case18[k][8]=11;
  Case18[k][9]=5;
  Case18[k][10]=0;
  Case18[k][11]=11;

  //case 4: 2 faces (6,18-connectivity)
  k=4;
  Case18[k][0]=4;
  Case18[k][1]=0;
  Case18[k][2]=1;
  Case18[k][3]=7;
  Case18[k][4]=11;
  Case18[k][5]=10;


  //case 5: 3 faces
  k=5;
  Case18[k][0]=5;
  Case18[k][1]=7;
  Case18[k][2]=6;
  Case18[k][3]=5;
  Case18[k][4]=6;
  Case18[k][5]=1;
  Case18[k][6]=5;
  Case18[k][7]=1;
  Case18[k][8]=0;

  //case 6: 5 faces (18-26-connectivity)
  k=6;
  Case18[k][0]=4;
  Case18[k][1]=10;
  Case18[k][2]=1;
  Case18[k][3]=4;
  Case18[k][4]=5;
  Case18[k][5]=10;
  Case18[k][6]=5;
  Case18[k][7]=11;
  Case18[k][8]=10;
  Case18[k][9]=1;
  Case18[k][10]=10;
  Case18[k][11]=3;
  Case18[k][12]=3;
  Case18[k][13]=10;
  Case18[k][14]=7;

  //case 7: 15 faces (18-26-connectivity)
  k=7;
  Case18[k][0]=5;
  Case18[k][1]=11;
  Case18[k][2]=8;
  Case18[k][3]=0;
  Case18[k][4]=4;
  Case18[k][5]=9;
  Case18[k][6]=0;
  Case18[k][7]=9;
  Case18[k][8]=10;
  Case18[k][9]=0;
  Case18[k][10]=10;
  Case18[k][11]=3;
  Case18[k][12]=3;
  Case18[k][13]=10;
  Case18[k][14]=7;


  //case 8: 2 faces
  k=8;
  Case18[k][0]=4;
  Case18[k][1]=5;
  Case18[k][2]=6;
  Case18[k][3]=6;
  Case18[k][4]=5;
  Case18[k][5]=7;

  //case 9: 4 faces
  k=9;
  Case18[k][0]=4;
  Case18[k][1]=0;
  Case18[k][2]=9;
  Case18[k][3]=0;
  Case18[k][4]=10;
  Case18[k][5]=9;
  Case18[k][6]=0;
  Case18[k][7]=3;
  Case18[k][8]=10;
  Case18[k][9]=10;
  Case18[k][10]=3;
  Case18[k][11]=7;

  //case 10: 4 faces (18-26-connectivity)
  k=10;
  Case18[k][0]=1;
  Case18[k][1]=9;
  Case18[k][2]=10;
  Case18[k][3]=1;
  Case18[k][4]=10;
  Case18[k][5]=2;
  Case18[k][6]=0;
  Case18[k][7]=11;
  Case18[k][8]=8;
  Case18[k][9]=0;
  Case18[k][10]=3;
  Case18[k][11]=11;


  //case 11: 4 faces
  k=11;
  Case18[k][0]=4;
  Case18[k][1]=0;
  Case18[k][2]=6;
  Case18[k][3]=0;
  Case18[k][4]=3;
  Case18[k][5]=6;
  Case18[k][6]=6;
  Case18[k][7]=3;
  Case18[k][8]=10;
  Case18[k][9]=10;
  Case18[k][10]=3;
  Case18[k][11]=11;


  //case 12: 4 faces (18-26-connectivity)
  k=12;
  Case18[k][0]=0;
  Case18[k][1]=4;
  Case18[k][2]=1;
  Case18[k][3]=6;
  Case18[k][4]=9;
  Case18[k][5]=8;
  Case18[k][6]=6;
  Case18[k][7]=8;
  Case18[k][8]=5;
  Case18[k][9]=6;
  Case18[k][10]=5;
  Case18[k][11]=7;

  //case 13: 4 faces (18-26-connectivity)
  k=13;
  Case18[k][0]=7;
  Case18[k][1]=10;
  Case18[k][2]=11;
  Case18[k][3]=4;
  Case18[k][4]=8;
  Case18[k][5]=9;
  Case18[k][6]=0;
  Case18[k][7]=3;
  Case18[k][8]=5;
  Case18[k][9]=1;
  Case18[k][10]=6;
  Case18[k][11]=2;

  //case 14: 4 faces
  k=14;
  Case18[k][0]=0;
  Case18[k][1]=9;
  Case18[k][2]=1;
  Case18[k][3]=0;
  Case18[k][4]=5;
  Case18[k][5]=9;
  Case18[k][6]=5;
  Case18[k][7]=10;
  Case18[k][8]=9;
  Case18[k][9]=5;
  Case18[k][10]=7;
  Case18[k][11]=10;

  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //// 26-connectivity /////////////////////////////////
  /////////////////////////////////////////////////////

  for (k=0;k<15;k++)
    reinit(Case26[k],NV);

  //case 0
  //nothing to do

  //case 1: 1 face
  k=1;
  Case26[k][0]=4;
  Case26[k][1]=0;
  Case26[k][2]=1;


  //case 2: 2 faces
  k=2;
  Case26[k][0]=4;
  Case26[k][1]=5;
  Case26[k][2]=1;
  Case26[k][3]=1;
  Case26[k][4]=5;
  Case26[k][5]=3;

  //case 3: 4 faces (18-26 connectivity)
  k=3;
  Case26[k][0]=1;
  Case26[k][1]=4;
  Case26[k][2]=8;
  Case26[k][3]=1;
  Case26[k][4]=8;
  Case26[k][5]=11;
  Case26[k][6]=0;
  Case26[k][7]=1;
  Case26[k][8]=11;
  Case26[k][9]=5;
  Case26[k][10]=0;
  Case26[k][11]=11;

  //case 4: 6 faces (26-connectivity)
  k=4;
  Case26[k][0]=1;
  Case26[k][1]=4;
  Case26[k][2]=10;
  Case26[k][3]=4;
  Case26[k][4]=11;
  Case26[k][5]=10;
  Case26[k][6]=0;
  Case26[k][7]=7;
  Case26[k][8]=11;
  Case26[k][9]=4;
  Case26[k][10]=0;
  Case26[k][11]=11;
  Case26[k][12]=0;
  Case26[k][13]=1;
  Case26[k][14]=7;
  Case26[k][15]=1;
  Case26[k][16]=10;
  Case26[k][17]=7;


  //case 5: 3 faces
  k=5;
  Case26[k][0]=5;
  Case26[k][1]=7;
  Case26[k][2]=6;
  Case26[k][3]=5;
  Case26[k][4]=6;
  Case26[k][5]=1;
  Case26[k][6]=5;
  Case26[k][7]=1;
  Case26[k][8]=0;

  //case 6: 5 faces (18-26-connectivity)
  k=6;
  Case26[k][0]=4;
  Case26[k][1]=10;
  Case26[k][2]=1;
  Case26[k][3]=4;
  Case26[k][4]=5;
  Case26[k][5]=10;
  Case26[k][6]=5;
  Case26[k][7]=11;
  Case26[k][8]=10;
  Case26[k][9]=1;
  Case26[k][10]=10;
  Case26[k][11]=3;
  Case26[k][12]=3;
  Case26[k][13]=10;
  Case26[k][14]=7;

  //case 7: 15 faces (18-26-connectivity)
  k=7;
  Case26[k][0]=5;
  Case26[k][1]=11;
  Case26[k][2]=8;
  Case26[k][3]=0;
  Case26[k][4]=4;
  Case26[k][5]=9;
  Case26[k][6]=0;
  Case26[k][7]=9;
  Case26[k][8]=10;
  Case26[k][9]=0;
  Case26[k][10]=10;
  Case26[k][11]=3;
  Case26[k][12]=3;
  Case26[k][13]=10;
  Case26[k][14]=7;


  //case 8: 2 faces
  k=8;
  Case26[k][0]=4;
  Case26[k][1]=5;
  Case26[k][2]=6;
  Case26[k][3]=6;
  Case26[k][4]=5;
  Case26[k][5]=7;

  //case 9: 4 faces
  k=9;
  Case26[k][0]=4;
  Case26[k][1]=0;
  Case26[k][2]=9;
  Case26[k][3]=0;
  Case26[k][4]=10;
  Case26[k][5]=9;
  Case26[k][6]=0;
  Case26[k][7]=3;
  Case26[k][8]=10;
  Case26[k][9]=10;
  Case26[k][10]=3;
  Case26[k][11]=7;

  //case 10: 4 faces (18-26-connectivity)
  k=10;
  Case26[k][0]=1;
  Case26[k][1]=9;
  Case26[k][2]=10;
  Case26[k][3]=1;
  Case26[k][4]=10;
  Case26[k][5]=2;
  Case26[k][6]=0;
  Case26[k][7]=11;
  Case26[k][8]=8;
  Case26[k][9]=0;
  Case26[k][10]=3;
  Case26[k][11]=11;


  //case 11: 4 faces
  k=11;
  Case26[k][0]=4;
  Case26[k][1]=0;
  Case26[k][2]=6;
  Case26[k][3]=0;
  Case26[k][4]=3;
  Case26[k][5]=6;
  Case26[k][6]=6;
  Case26[k][7]=3;
  Case26[k][8]=10;
  Case26[k][9]=10;
  Case26[k][10]=3;
  Case26[k][11]=11;


  //case 12: 4 faces (18-26-connectivity)
  k=12;
  Case26[k][0]=0;
  Case26[k][1]=4;
  Case26[k][2]=1;
  Case26[k][3]=6;
  Case26[k][4]=9;
  Case26[k][5]=8;
  Case26[k][6]=6;
  Case26[k][7]=8;
  Case26[k][8]=5;
  Case26[k][9]=6;
  Case26[k][10]=5;
  Case26[k][11]=7;

  //case 13: 4 faces (18-26-connectivity)
  k=13;
  Case26[k][0]=7;
  Case26[k][1]=10;
  Case26[k][2]=11;
  Case26[k][3]=4;
  Case26[k][4]=8;
  Case26[k][5]=9;
  Case26[k][6]=0;
  Case26[k][7]=3;
  Case26[k][8]=5;
  Case26[k][9]=1;
  Case26[k][10]=6;
  Case26[k][11]=2;

  //case 14: 4 faces
  k=14;
  Case26[k][0]=0;
  Case26[k][1]=9;
  Case26[k][2]=1;
  Case26[k][3]=0;
  Case26[k][4]=5;
  Case26[k][5]=9;
  Case26[k][6]=5;
  Case26[k][7]=10;
  Case26[k][8]=9;
  Case26[k][9]=5;
  Case26[k][10]=7;
  Case26[k][11]=10;

}

static void initPoints(void) {
  int k;

  for (k=0;k<15;k++)
    reinit(Pt[k],NP);

  //case 0
  //nothing to do

  //case 1: 1 face
  k=1;
  //1 point
  Pt[k][0]=Pt[k][1]=Pt[k][2]=0;


  //case 2: 2 faces
  k=2;
  //2 points
  Pt[k][0]=Pt[k][1]=Pt[k][2]=0;
  Pt[k][3]=1;
  Pt[k][4]=Pt[k][5]=0;


  //case 3: 2 faces (6-connectivity)
  k=3;
  //2 points
  Pt[k][0]=Pt[k][1]=Pt[k][2]=0;
  Pt[k][4]=0;
  Pt[k][3]=Pt[k][5]=1;

  //case 4: 2 faces (6,18-connectivity)
  k=4;
  //2 points
  Pt[k][0]=Pt[k][1]=Pt[k][2]=0;
  Pt[k][3]=Pt[k][4]=Pt[k][5]=1;

  //case 5: 3 faces (6-connectivity)
  k=5;
  //3 points
  Pt[k][0]=1;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=1;
  Pt[k][4]=1;
  Pt[k][5]=0;
  Pt[k][6]=0;
  Pt[k][7]=1;
  Pt[k][8]=0;

  //case 6: 3 faces (6-connectivity)
  k=6;
  //3 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=1;
  Pt[k][4]=0;
  Pt[k][5]=0;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=1;

  //case 7: 3 faces (6-connectivity)
  k=7;
  //3 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=1;
  Pt[k][3]=1;
  Pt[k][4]=0;
  Pt[k][5]=0;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=1;

  //case 8: 2 faces
  k=8;
  //4 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=1;
  Pt[k][4]=0;
  Pt[k][5]=0;
  Pt[k][6]=0;
  Pt[k][7]=1;
  Pt[k][8]=0;
  Pt[k][9]=1;
  Pt[k][10]=1;
  Pt[k][11]=0;

  //case 9: 4 faces
  k=9;
  //4 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=0;
  Pt[k][4]=1;
  Pt[k][5]=0;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=0;
  Pt[k][9]=0;
  Pt[k][10]=1;
  Pt[k][11]=1;

  //case 10: 4 faces
  k=10;
  //4 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=0;
  Pt[k][4]=0;
  Pt[k][5]=1;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=1;
  Pt[k][9]=1;
  Pt[k][10]=1;
  Pt[k][11]=0;

  //case 11: 4 faces
  k=11;
  //4 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=0;
  Pt[k][4]=1;
  Pt[k][5]=0;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=0;
  Pt[k][9]=1;
  Pt[k][10]=1;
  Pt[k][11]=1;

  //case 12: 4 faces
  k=12;
  //4 points
  Pt[k][0]=1;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=0;
  Pt[k][4]=1;
  Pt[k][5]=0;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=0;
  Pt[k][9]=0;
  Pt[k][10]=0;
  Pt[k][11]=1;

  //case 13: 4 faces
  k=13;
  //4 points
  Pt[k][0]=0;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=0;
  Pt[k][4]=1;
  Pt[k][5]=1;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=0;
  Pt[k][9]=1;
  Pt[k][10]=0;
  Pt[k][11]=1;

  //case 14: 4 faces
  k=14;
  //4 points
  Pt[k][0]=1;
  Pt[k][1]=0;
  Pt[k][2]=0;
  Pt[k][3]=0;
  Pt[k][4]=1;
  Pt[k][5]=0;
  Pt[k][6]=1;
  Pt[k][7]=1;
  Pt[k][8]=0;
  Pt[k][9]=0;
  Pt[k][10]=1;
  Pt[k][11]=1;
}


int main(int argc, char *argv[]) {
  int Tess[256][NV];
  int Case[NV],Cased[NV],Ctmp[NV];
  int pt[NP],Ptd[NP],Ptmp[NP];
  int k,l,count,a,b,c,d,e,f;//,g,h,i;
  int con;

  FILE *file=NULL;


  for (con=1;con<5;con++) {
    initPoints();
    initCases();

    //now build all the different cases
    //have to take into account 3 sym. + 1 revers.

    for (k=0;k<256;k++)
      reinit(Tess[k],NV);

    for (k=0;k<15;k++) {

      reinit(Case,NV);
      reinit(pt,NP);

      switch (con) {
      case 1:
        cp(Case6[k],Case,NV);
        break;
      case 2:
        cp(Case18[k],Case,NV);
        break;
      case 3:
        cp(Case6[k],Case,NV);
        break;
      case 4:
        cp(Case26[k],Case,NV);
        break;
      default:
        cp(Case6[k],Case,NV);
        break;
      }
      cp(Pt[k],pt,NP);

      //Id
      copy(Case,Tess[gNb(pt)]);

      //Rotations
      for (a=0;a<4;a++)
        for (b=0;b<4;b++)
          for (c=0;c<4;c++)
            for (d=0;d<4;d++)
              for (e=3;e<4;e++)
                for (f=3;f<4;f++) {
                  reinit(Cased,NV);
                  reinit(Ptd,NP);
                  reinit(Ctmp,NV);
                  reinit(Ptmp,NP);
                  switch (a) {
                  case 0:
                    ROx(Case,pt,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Case,pt,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Case,pt,Cased,Ptd);
                    break;
                  case 3:
                    cp(Case,Cased,NV);
                    cp(pt,Ptd,NP);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (b) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (c) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (d) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (e) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (f) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }

                  copy(Cased,Tess[gNb(Ptd)]);
                }
    }
    count=0;
    for (k=0;k<256;k++)
      if (Tess[k][0]>=0)
        count++;
    fprintf(stderr,"\n%d cases ",count);

    for (k=0;k<15;k++) {
      reinit(Case,NV);
      reinit(pt,NP);

      switch (con) {
      case 1:
        cp(Case18[k],Case,NV);
        break;
      case 2:
        cp(Case6[k],Case,NV);
        break;
      case 3:
        cp(Case26[k],Case,NV);
        break;
      case 4:
        cp(Case6[k],Case,NV);
        break;
      default:
        cp(Case18[k],Case,NV);
        break;
      }
      cp(Pt[k],pt,NP);

      reinit(Ctmp,NV);
      reinit(Ptmp,NP);
      Reverse(Case,pt,Ctmp,Ptmp);

      cp(Ctmp,Case,NV);
      cp(Ptmp,pt,NP);

      //Id
      copy(Case,Tess[gNb(pt)]);


      //Rotations
      for (a=0;a<4;a++)
        for (b=0;b<4;b++)
          for (c=0;c<4;c++)
            for (d=0;d<4;d++)
              for (e=3;e<4;e++)
                for (f=3;f<4;f++) {
                  reinit(Cased,NV);
                  reinit(Ptd,NP);
                  reinit(Ctmp,NV);
                  reinit(Ptmp,NP);
                  switch (a) {
                  case 0:
                    ROx(Case,pt,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Case,pt,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Case,pt,Cased,Ptd);
                    break;
                  case 3:
                    cp(Case,Cased,NV);
                    cp(pt,Ptd,NP);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (b) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (c) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (d) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (e) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  cp(Cased,Ctmp,NV);
                  cp(Ptd,Ptmp,NP);
                  switch (f) {
                  case 0:
                    ROx(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 1:
                    ROy(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  case 2:
                    ROz(Ctmp,Ptmp,Cased,Ptd);
                    break;
                  }
                  copy(Cased,Tess[gNb(Ptd)]);
                }
    }


    count=0;
    for (k=0;k<256;k++) {
      l=0;
      if (Tess[k][l]>=0) {
        count++;
        /*fprintf(stderr,"\n %d -> ",k);
        while(Tess[k][l]>=0)
        fprintf(stderr,"%d ",Tess[k][l++]);*/
      }
    }
    fprintf(stderr,"\n%d cases for con=%d\n",count,con);
    if (count<254)
      ;

    for (k=0;k<256;k++) {
      /*      fprintf(stderr,"\n{");
       for(l=0;l<NV-1;l++)
       fprintf(stderr,"%d,",Tess[k][l]);
       fprintf(stderr,"%d},",Tess[k][l]);*/
    }
    fprintf(stderr,"\n");

    {
      FILE *f;
      char fname[100];
      sprintf(fname,"./MC%d",con);
      f=fopen(fname,"w");
      for (k=0;k<256;k++) {
        fprintf(f,"\n{");
        for (l=0;l<NV-1;l++)
          fprintf(f,"%d,",Tess[k][l]);
        fprintf(f,"%d},",Tess[k][l]);
      }
      fclose(f);


      if (con==1)
        file=fopen("../mri_test/MC.h","w");

      switch (con) {
      case 1:
        fprintf(file,"int MC6p[256][19]={");
        break;
      case 2:
        fprintf(file,"int MC18[256][19]={");
        break;
      case 3:
        fprintf(file,"int MC6[256][19]={");
        break;
      case 4:
        fprintf(file,"int MC26[256][19]={");
        break;
      default:
        fprintf(file,"int MC6p[256][19]={");
        break;
      }
      for (k=0;k<255;k++) {
        fprintf(file,"\n{");
        for (l=0;l<NV-1;l++)
          fprintf(file,"%d,",Tess[k][l]);
        fprintf(file,"%d},",Tess[k][l]);
      }
      fprintf(file,"\n{");
      for (l=0;l<NV-2;l++)
        fprintf(file,"%d,",Tess[255][l]);
      fprintf(file,"%d}\n};\n",Tess[255][NV-2]);
      if (con==4)
        fclose(file);
    }



  }
  return 0;
}








