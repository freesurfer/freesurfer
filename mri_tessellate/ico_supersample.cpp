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


/* isobol.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "trilib.h"

char name[2048], line[80];
int  dhk[4000000][3], npnt, ndhk;
float pnt[2000000][3], x[3], y[3];
FILE *outfile;


void icosa(void) {
  int index[3], i, j, k, l;
  int i1, i2, i3, j1, npnt0, ndhk0;
  float h[3][3], eps=0.001, r, ra, rho, fpi, d;
  char line[80];
  int ico[3][20]=
    {{1, 1, 1, 1, 1, 4, 4, 4, 5, 5, 6, 6, 2, 2, 3, 9,10,11, 7, 8},
     {5, 6, 2, 3, 4, 9,10, 5,11, 6, 7, 2, 8, 3, 9,10,11, 7, 8, 9},
     {4, 5, 6, 2, 3, 3, 9,10,10,11,11, 7, 7, 8, 8,12,12,12,12,12}};

  /* contruct icosahedron */

  for (i=0; i<20; i++)
    for (j=0; j<3; j++)
      dhk[i][j]=ico[j][i]-1;
  r=1;
  pnt[0][0]=0;
  pnt[0][1]=0;
  pnt[0][2]=r;
  pnt[11][0]=0;
  pnt[11][1]=0;
  pnt[11][2]=-r;
  rho=0.4*sqrt(5.0)*r;
  fpi=0.4*M_PI;
  for (i=1; i<6; i++) {
    pnt[i][0]=rho*cos((i-2)*fpi);
    pnt[i][1]=rho*sin((i-2)*fpi);
    pnt[i][2]=rho/2;
  }
  for (i=6; i<11; i++) {
    pnt[i][0]=rho*cos((i-2.5)*fpi);
    pnt[i][1]=rho*sin((i-2.5)*fpi);
    pnt[i][2]=-rho/2;
  }
  npnt=12;
  ndhk=20;

  /* refine icosahedron */

  for ( ; ; ) {
    printf("Currently %d points and %d triangels.\n",npnt,ndhk);
    printf("Refine? (y/n) ");
    fgets(line,3,stdin);
    if (tolower(line[0])!='y') return;

    /* vorm nieuwe hoekpunten en projecteer deze op de bol  */
    /* controleer of deze al op lijst voorkomen.    */
    /* zo ja: noteer index zo nee: neem op in lijst en noteer index */

    npnt0=npnt;
    ndhk0=ndhk;
    for (i=0; i<ndhk0; i++) {
      for (j=0; j<3; j++) {
        j1=(j+1)%3;
        i1=dhk[i][j];
        i2=dhk[i][j1];
        for (k=0; k<3; k++)
          h[j][k]=(pnt[i1][k]+pnt[i2][k])/2;
        ra=sqrt(h[j][0]*h[j][0]+h[j][1]*h[j][1]+h[j][2]*h[j][2]);
        for (k=0; k<3; k++)
          h[j][k]=h[j][k]*r/ra;
        for (l=npnt0; l<npnt; l++) {
          for (k=0; k<3; k++) {
            d=fabs(h[j][k]-pnt[l][k]);
            if (d>=eps) break;
          }
          if (k==3) {
            index[j]=l;
            goto nextj;
          }
        }
        index[j]=npnt;
        for (k=0; k<3; k++)
          pnt[npnt][k]=h[j][k];
        npnt=npnt+1;
nextj:
        ;
      }

      /* vul hoekpunten van nieuwe vlakken in,gooi oude driehoekje weg */

      i2=dhk[i][1];
      i3=dhk[i][2];
      dhk[i][1]=index[0];
      dhk[i][2]=index[2];
      dhk[ndhk][0]=index[0];
      dhk[ndhk][1]=i2;
      dhk[ndhk][2]=index[1];
      ndhk++;
      dhk[ndhk][0]=index[2];
      dhk[ndhk][1]=index[0];
      dhk[ndhk][2]=index[1];
      ndhk++;
      dhk[ndhk][0]=index[2];
      dhk[ndhk][1]=index[1];
      dhk[ndhk][2]=i3;
      ndhk++;
    }
  }
}

int main(int argc, char *argv[]) {
  int i,k;
  float radius, a, b, c, rlab, r[3]={0.0, 0.0, 0.0}, ri[3];
  char line[280];

  icosa();

  printf("Radius of sphere: ");
  fgets(line, 50, stdin);
  if (line[0]=='\0') exit(0);
  sscanf(line,"%f",&radius);

  printf("Projection point (default origin): ");
  fgets(line,50,stdin);
  if (line[0]!='\0')
    if (sscanf(line,"%f%f%f", r, r+1, r+2)!=3) {
      printf("Enter three floating point values separated by spaces!\n\n\a");
      exit(1);
    }
  c=r[0]*r[0]+r[1]*r[1]+r[2]*r[2] - radius*radius;
  if (c>=0) {
    printf("Projection point must be within sphere\n");
    exit(1);
  }

  printf("number of points:%4d\n", npnt);
  for (i=0; i<npnt; i++) {
    printf("   now at point :%4d\r", i+1);
    a=b=0;
    for (k=0; k<3; k++) {
      ri[k]=pnt[i][k];
      a += ri[k]*ri[k];
      b += ri[k]*r[k];
    }
    b = 2*b;
    rlab=(-b+sqrt(b*b-4*a*c))/(2*a);
    if (rlab<0) rlab=(-b-sqrt(b*b-4*a*c))/(2*a);
    for (k=0; k<3; k++)
      pnt[i][k]=r[k]+rlab*ri[k];
  }

  printf("\nName output file: ");
  fgets(name,2000,stdin);
  if (name[0]=='\0') exit(0);
  outfile=fopen(name, "w");
  if (outfile==NULL) {
    printf("Error opening %s\n\n\a", name);
    exit(1);
  }
  fprintf(outfile,"%4d\n",npnt);
  for (i=0; i<npnt; i++)
    fprintf(outfile,"%4d %7.4f %7.4f %7.4f\n",
            i+1,pnt[i][0],pnt[i][1],pnt[i][2]);
  fprintf(outfile,"%4d\n",ndhk);
  for (i=0; i<ndhk; i++)
    fprintf(outfile,"%4d %3d %3d %3d\n",
            i+1,dhk[i][0]+1,dhk[i][1]+1,dhk[i][2]+1);
  fclose(outfile);
  return(0) ;
}

