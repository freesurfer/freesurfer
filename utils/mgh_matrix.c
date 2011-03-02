/**
 * @file  mgh_matrix.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.7 $
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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "typedefs.h"
#include "proto.h"
#include "mgh_matrix.h"
#include "error.h"

FLOATTYPE * MGH_vector(int n)
{
  FLOATTYPE *h;

  h = (FLOATTYPE *)calloc(n,sizeof(FLOATTYPE));
  if (h==NULL)
  {
    printf("out of memory in vector(%d)\n",n);
    exit(0);
  }
  return h;
}

int* MGH_ivector(int n)
{
  int *h;

  h = (int *)calloc(n,sizeof(int));
  if (h==NULL)
  {
    printf("out of memory in ivector(%d)\n",n);
    exit(0);
  }
  return h;
}

FLOATTYPE ** MGH_matrix(int n,int m)
{
  FLOATTYPE **h;
  int i;

  h = (FLOATTYPE **)calloc(n,sizeof(FLOATTYPE *));
  if (h==NULL)
  {
    printf("out of memory in matrix(%d,%d)\n",n,m);
    exit(0);
  }
  for (i=0;i<n;i++) h[i]=MGH_vector(m);
  return h;
}

void print_matrix(FLOATTYPE **a,int m,int n)
{
  int i,j;

  if (n>0)
  {
    for (i=0;i<m;i++)
    {
      for (j=0;j<n;j++)
      {
        printf("%10.3e ",a[i][j]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

void read_matrix(FILE *fptr,FLOATTYPE **a,int m,int n)
{
  int i,j;
  float f;

  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
    {
      fscanf(fptr,"%f",&f);
      a[i][j] = f;
    }
}

void print_vector(FLOATTYPE *v,int n)
{
  int i;

  for (i=0;i<n;i++)
  {
    printf("%10.3e ",v[i]);
  }
  printf("\n\n");
}

void row_vector(FLOATTYPE **a,FLOATTYPE *v,int i,int n)
{
  int j;

  for (j=0;j<n;j++)
  {
    v[j] = a[i][j];
  }
}

void vector_to_matrix(FLOATTYPE *v,FLOATTYPE **a,int m,int n)
{
  int i,j,indx=0;

  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
    {
      a[i][j] = v[indx++];
    }
}

void scale_matrix(FLOATTYPE **a,FLOATTYPE s,int n,int m)
{
  int i,j;

  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      a[i][j] *= s;
}

void normalize_matrix(FLOATTYPE **a,int n,int m)
{
  FLOATTYPE sum=0;
  int i,j;

  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      sum += a[i][j];
  /*  scale_matrix(a,n/sum,n,m); */
}

void matrix_copy(FLOATTYPE **a,FLOATTYPE **b,int n,int m)
{
  int i,j;

  for (i=0;i<n;i++) for (j=0;j<m;j++) b[i][j] = a[i][j];
}

void matrix_copy2(FLOATTYPE **a,FLOATTYPE **b,int n,int m,int sno,int smo,
                  int tno,int tmo)
{
  int i,j;

  for (i=0;i<n;i++) for (j=0;j<m;j++) b[i+tno][j+tmo] = a[i+sno][j+smo];
}

void matrix_transpose(FLOATTYPE **a,FLOATTYPE **at,int n,int m)
{
  int i,j;

  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      at[i][j] = a[j][i];
}

void matrix_add(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c, int n,int m)
{
  int i,j;

  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      c[i][j] = a[i][j]+b[i][j];
}

void matrix_multiply(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c,int n,int m)
{
  FLOATTYPE sum;
  int i,j,k;

  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
    {
      sum = 0.0;
      for (k=0;k<m;k++)
        sum += a[i][k]*b[k][j];
      c[i][j] = sum;
    }
  }
}

void matrix_multiply2(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c,int n,int m,
                      int l)
{
  FLOATTYPE sum;
  int i,j,k;

  for (i=0;i<n;i++)
  {
    for (j=0;j<l;j++)
    {
      sum = 0.0;
      for (k=0;k<m;k++)
        sum += a[i][k]*b[k][j];
      c[i][j] = sum;
    }
  }
}

void matrix_angles(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c,int n,int m)
{
  FLOATTYPE sum,asum,bsum;
  int i,j,k;

  for (i=0;i<n;i++)
  {
    for (j=0;j<m;j++)
    {
      sum = asum = bsum = 0.0;
      for (k=0;k<m;k++)
      {
        sum += a[i][k]*b[k][j];
        asum += a[i][k]*a[i][k];
        bsum += b[k][j]*b[k][j];
      }
      c[i][j] = (sum==0.0)?0.0:sum/sqrt(asum*bsum);
    }
  }
}

void vector_subtract(FLOATTYPE *a,FLOATTYPE *b,FLOATTYPE *c,int n)
{
  int i;

  for (i=0;i<n;i++)
    c[i] = a[i]-b[i];
}

void vector_add(FLOATTYPE *a,FLOATTYPE *b,FLOATTYPE *c,FLOATTYPE fa,
                FLOATTYPE fb,int n)
{
  int i;

  for (i=0;i<n;i++)
    c[i] = fa*a[i]+fb*b[i];
}

void vector_multiply(FLOATTYPE **a,FLOATTYPE *b,FLOATTYPE *c,int n,int m)
{
  FLOATTYPE sum;
  int i,j;

  for (i=0;i<n;i++)
  {
    sum = 0.0;
    for (j=0;j<m;j++)
      sum += a[i][j]*b[j];
    c[i] = sum;
  }
}

void derivative_matrix(FLOATTYPE **D,int dy,int dx)
{
  int i,j,n=dx*dy;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      D[i][j]=0.0;
  for (i=0;i<n;i++)
  {
    D[i][i] = 1;
    D[i][(i+1)%n] = -1;
  }
}

void MGH_identity_matrix(FLOATTYPE **I,int n)
{
  int i,j;

  for (i=0;i<n;i++) for (j=0;j<n;j++) I[i][j]=(i==j)?1.0:0.0;
}

void regularization_matrix(FLOATTYPE **R,int n)
{
  int i,j;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      R[i][j]=0;
  for (i=0;i<n;i++)
  {
    R[i][(i==0)?n-1:i-1] = -1.0;
    R[i][i] = 2.1;
    R[i][(i==n-1)?0:i+1] = -1.0;
  }
}

void covariance_matrix(FLOATTYPE **R,int n)
{
  int i,j;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    {
      R[i][j]=1.0/(1.0+fabs(i-j));
      /*    printf("i=%d,j=%d,R[i][j]=%f\n",i,j,R[i][j]); */
    }
}


void mgh_ludcmp(FLOATTYPE **a,int n,int *indx,FLOATTYPE *d)
{
  int i,imax=0,j,k; /* imax = 0 to shut up gcc */
  FLOATTYPE big,dum,sum,temp;
  FLOATTYPE *vv,*vector();

  vv=MGH_vector(n);
  *d=1.0;
  for (i=0;i<n;i++)
  {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0)
      ErrorExit(ERROR_BADPARM, "LU decomposition: matrix in singular");
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++)
  {
    for (i=0;i<j;i++)
    {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++)
    {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big)
      {
        big=dum;
        imax=i;
      }
    }
    if (j != imax)
    {
      for (k=0;k<n;k++)
      {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n-1)
    {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}

void mgh_lubksb(FLOATTYPE **a,int n,int *indx,FLOATTYPE *b)
{
  int i,ii= -1,ip,j;
  FLOATTYPE sum;

  for (i=0;i<n;i++)
  {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>-1)
      for (j=ii;j<i;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n-1;i>=0;i--)
  {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

void inverse(FLOATTYPE **a,FLOATTYPE **y,int n)
{
  FLOATTYPE d,*col;
  int i,j,*indx;

  col = MGH_vector(n);
  indx = MGH_ivector(n);
  mgh_ludcmp(a,n,indx,&d);
  for (j=0;j<n;j++)
  {
    for (i=0;i<n;i++) col[i]=0.0;
    col[j]=1.0;
    mgh_lubksb(a,n,indx,col);
    for (i=0;i<n;i++) y[i][j]=col[i];
  }
}

FLOATTYPE determinant(FLOATTYPE **a,int n)
{
  FLOATTYPE d;
  int j,*indx;

  indx = MGH_ivector(n);
  mgh_ludcmp(a,n,indx,&d);
  print_matrix(a,n,n);
  for (j=0;j<n;j++) d *= a[j][j];
  return d;
}

FLOATTYPE MGH_svd(FLOATTYPE **A,FLOATTYPE **V,FLOATTYPE *z,int m,int n)
{
  int i,j,k,count;
  FLOATTYPE c,s,p,q,r,v,toll=0.1;

  MGH_identity_matrix(V,n);
  for (count=n*(n-1)/2;count>0;)
  {
    count=n*(n-1)/2;
    for (j=0;j<n-1;j++)
    {
      for (k=j+1;k<n;k++)
      {
        p=q=r=0;
        for (i=0;i<m;i++)
        {
          p += A[i][j]*A[i][k];
          q += A[i][j]*A[i][j];
          r += A[i][k]*A[i][k];
        }
        if ((q*r==0)||(p*p/(q*r)<toll)) count--;
        if (q<r)
        {
          c=0;
          s=1;
        }
        else
        {
          q = q-r;
          v = sqrt(4*p*p+q*q);
          c = sqrt((v+q)/(2*v));
          s = p/(v*c);
        }
        for (i=0;i<m;i++)
        {
          r = A[i][j];
          A[i][j] = r*c+A[i][k]*s;
          A[i][k] = -r*s+A[i][k]*c;
        }
        for (i=0;i<n;i++)
        {
          r = V[i][j];
          V[i][j] = r*c+V[i][k]*s;
          V[i][k] = -r*s+V[i][k]*c;
        }
      }
    }
    printf("count=%d\n",count);
    print_matrix(A,m,n);
  }
  for (j=0;j<n;j++)
  {
    q = 0;
    for (i=0;i<m;i++) q += A[i][j]*A[i][j];
    q = sqrt(q);
    z[j] = q;

    for (i=0;i<m;i++) A[i][j] /= q;

  }

  return (FLOATTYPE)(0.0); /* What is this *really* suppose to return? (RJW) */
}

static FLOATTYPE at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)):0.0))

#ifdef MAX
#undef MAX
#endif
static FLOATTYPE maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1):(maxarg2))

#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

void mgh_svdcmp(FLOATTYPE **a,FLOATTYPE *w,FLOATTYPE **v,int m,int n)
{
  int flag,i,its,j,jj,k,l=0,nm=0; /* Keep gcc happy */
  FLOATTYPE c,f,h,s,x,y,z;
  FLOATTYPE anorm=0.0,g=0.0,scale=0.0;
  FLOATTYPE *rv1;

  if (m<n)
    ErrorExit(ERROR_BADPARM,
              "Singular value decompostion: input matrix has too few rows");
  rv1=MGH_vector(n);
  for (i=0;i<n;i++)
  {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i<m)
    {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale)
      {
        for (k=i;k<m;k++)
        {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        if (i != n-1)
        {
          for (j=l;j<n;j++)
          {
            for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
            f=s/h;
            for (k=i;k<m;k++) a[k][j] += f*a[k][i];
          }
        }
        for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i < m && i != n-1)
    {
      for (k=l;k<n;k++) scale += fabs(a[i][k]);
      if (scale)
      {
        for (k=l;k<n;k++)
        {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
        if (i != m-1)
        {
          for (j=l;j<m;j++)
          {
            for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<n;k++) a[j][k] += s*rv1[k];
          }
        }
        for (k=l;k<n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }

  for (i=n-1;i>=0;i--)
  {
    if (i < n-1)
    {
      if (g)
      {
        for (j=l;j<n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<n;j++)
        {
          for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }

  for (i=n-1;i>=0;i--)
  {
    l=i+1;
    g=w[i];
    if (i < n-1)
      for (j=l;j<n;j++) a[i][j]=0.0;
    if (g)
    {
      g=1.0/g;
      if (i != n-1)
      {
        for (j=l;j<n;j++)
        {
          for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
          f=(s/a[i][i])*g;
          for (k=i;k<m;k++) a[k][j] += f*a[k][i];
        }
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    }
    else
    {
      for (j=i;j<m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }

  for (k=n-1;k>=0;k--)
  {
    for (its=1;its<=30;its++)
    {
      flag=1;
      for (l=k;l>=0;l--)
      {
        nm=l-1;
        if ((FLOATTYPE)(fabs(rv1[l])+anorm) == anorm)
        {
          flag=0;
          break;
        }
        if ((FLOATTYPE)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag)
      {
        c=0.0;
        s=1.0;
        for (i=l;i<=k;i++)
        {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if ((FLOATTYPE)(fabs(f)+anorm) == anorm) break;
          g=w[i];
          h=PYTHAG(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s=(-f*h);
          for (j=0;j<m;j++)
          {
            y=a[j][nm];
            z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l==k)
      {
        if (z < 0.0)
        {
          w[k] = -z;
          for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
        }
        break;
      }
      if (its==30)
        ErrorExit(ERROR_BADPARM,
                  "Singular value decomposition failed to converge");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;

      c=s=1.0;
      for (j=l;j<=nm;j++)
      {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=PYTHAG(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y=y*c;
        for (jj=0;jj<n;jj++)
        {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=PYTHAG(f,h);
        w[j]=z;
        if (z)
        {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=(c*g)+(s*y);
        x=(c*y)-(s*g);
        for (jj=0;jj<m;jj++)
        {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free(rv1);
}
