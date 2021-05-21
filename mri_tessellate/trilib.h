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


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef systemid
#include "systemid.h"
#endif

/*
   The Borland c compiler does not automatically convert from void*
   to float[3] *. Type casting in one statement does't succeed. So
   an array of type float[3] * has to be allocated in the following
   fashion: pnt=(PNT *) calloc.... To this end the following typedefs
   have been defined.
*/

typedef float PNT[3];
typedef int   DHK[3];


/* crosslin.c  */

int crossline(float *plane, int *itr, int *nlm, float (*rlm)[2],
              int npnt, float (*pnt)[3], int ndhk, int (*dhk)[3]);

/* lmoutr2d.c  */

int lmoutr2d(int *p1, int *p2, int *p3, int *r, float *lambda, float *mu);

/* lmoutr.c  */

int lmoutr(float *p1, float *p2, float *p3, float *r,
           float *lamba, float *mu, float *dif);

/* autoscal.c  */

void autoscal(float step, float *ss, char *f);

/* vecfun.c  */

float vecdot(float *r1, float *r2);
void veccross(float *r1, float *r2, float *r3);
float vecnorm(float *r);
int vecnormate(float *r);
float vecdist(float *r1, float *r2);
void vecmirror(float *r1, float *n, float *r2);
int linesect(float *p1, float *r1, float *p2, float *r2, float *lab, float *mu);

/* trinorm.c  */

void trinorm(float (*n)[3], int npnt, float (*pnt)[3],
             int ndhk, int (*dhk)[3], int do_normate);

/* checkasc.c  */

int check_asci(FILE *file);

/* dhkinp.c  */

int dhkinp(char *filename,int *npnt,int *ndhk,float (**pnt1)[3],int (**dhk1)[3]);

/* dhkout.c  */

int dhkout(char *filename,int npnt,int ndhk,float (*pnt1)[3],int (*dhk1)[3],
           int bin);

/* pntinp.c  */

int pntinp(char *filename,int *npnt, float (**pnt1)[3]);

/* nrutils.c  */

void nrerror (char *error_text);
float *vector(int nl, int nh);
void free_vector(float *v, int nl, int nh);
int *ivector(int nl, int nh);
void ifree_vector(int *v, int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void ifree_matrix(int **m, int nrl, int nrh, int ncl, int nch);

/* matrix3d.c  */

float ***matrix3d(int n3l, int n3h, int nrl, int nrh, int ncl, int nch);
void free_3dmatrix(float ***m, int n3l, int n3h, int nrl, int nrh, int ncl, int nch);

/* matrix.c  */

int putmatrix(char *name, float **m, int nrh, int nch, char *format);
int binputmatrix(char *name, float **m, int nrh, int nch);
float **getmatrix(char *name, int *nrh, int *nch);
int puttailmatrix(char *name, float **m, int nrh, int nch, char *format,
                  int nextra, char *extra);
int binputtailmatrix(char *name, float **m, int nrh, int nch,
                     int nextra, char *extra);
float **gettailmatrix(char *name, int *nrh, int *nch,
                      int *nextra, char **extra);
int ascimatrix(void);
int foreignmatrix(void);

/* rhoek.c   */

double rhoek(double r1[3], double r2[3], double r3[3], int *on_triangle);

/* solangle.c  */

int solangle(double r1[3], double r2[3], double r3[3],
             double *res1, double *res2, double *res3, double small);

/* rdds.c   */

void rdds(double r1[3], double r2[3], double r3[3],
          double *res1, double *res2, double *res3);

/* rxds.c   */

int rxds(double r1[3], double r2[3], double r3[3],
         double res1[3], double res2[3], double res3[3]);

/* ludec.c   */

void ludec(float **b, int npnt, float *cond);
void solvel(float **b, int npnt, float *x);
void solveu(float **b, int npnt, float *x);
void tsolvel(float **b, int npnt, float *x);
void tsolveu(float **b, int npnt, float *x);
void solvelu(float **b, int npnt, float *x);
int invertlu(float **b, int npnt, float **y);

/* marquard.c  */

int marquard(int nopt, int npar, float *v, float *phi, float *par,
             int (*f_bound)(int npar, float *par), float *res,
             void (*funeval) (int nopt, int npar, float *phi, float *par),
             void (*heseval) (int nopt, int npar, float **gtg, float *par,
                              float *err, float *gterr),
             int maxstp, float fstop, float *cond, int display,
             float *sdphi, float **cov);

/* svdcmp.c   */

void svdcmp(float **a, int m, int n, float *sigma, float **v);
void  sigdec(int m, int n, float **u, float *sigma, float **v);
void  svdsolve(float **u, int m, int n, float *sigma, float **v,
               float *b, float *x, float crit);
void  svdnsolve(float **u, int m, int n, float *sigma, float **v,
                float *b, float *x, int rang);

/* hfti.c    */

void hfti(float **a, int m, int n, float **b, int nb, float tau, int *krank);

/* exists.c   */

int exists(char *filename);
int overwrite(char *filename);

/* ask.c    */

int ask(char *prompt, char *line);

/* decouple.c   */

int decouple(int nopt, int nlpar, int nnpar, int ntim,
             float **v, float **phi, float *p, float **s,
             int (*f_bound)(int nnpar, float *p), float *res,
             void (*Aeval) (int nopt, int nlpar, int nnpar, float *p, float **A),
             void (*dAeval) (int nopt, int nlpar, int nnpar, float *p, float ***dA),
             int maxstp, float fstop, float *cond, int display, float *sdphi,
             float **cov, int truecov);

/* contour.c    */

typedef struct {
  int tr;
  float l;
  float m;
}
CONTOUR;

typedef struct {
  int tr;
  float l;
  float m;
  float l2;
  float m2;
}
DCONT;

void routlm(float *r, int tr, float l, float m,
            int npnt, float (*pnt)[3], int ndhk, int (*dhk)[3]);
int contour(float *p1, float *p2, int isp2, float *p3, int isp3, float *dir,
            int npnt, float (*pnt)[3], int ndhk, int (*dhk)[3], int *nlm,
            DCONT *c, int *p2tr);
int elecs(int nlm, DCONT *c, int i1, int i2,
          int npnt, float (*pnt)[3], int ndhk, int (*dhk)[3],
          int nel, float *pos, int isfrac, CONTOUR *el);

/* readtail.c    */

char *find_tail_item(int length, char tail[], char temp[]);
int tail_float(int length, char tail[], char temp[], float *value);
int tail_replace_item(int *length, char **tail, char temp[], char new_val[]);

/* trifun.c    */

void dtriag(float (*hp)[3], int (*ivl)[3], int *nvl, int ibcont[], int nb,
            int iocont[], int no);
void striag(float (*hp)[3], int (*ivl)[3], int *nvl, int ibcont[], int nb,
            int iocont[], int no);
void cdiag(float (*hp)[3], int (*ivl)[3], int nvl0, int nvl1, float *rl);
void pnt_cont2tri(float c[3], int ntrack, int nouter, int manual,
                  int log, int diag, int unreg, float t2,
                  int ncpnt, float (*cpnt)[3], int fpnt, float (*pnt)[3],
                  int fdhk, int (*dhk)[3], int *nnpnt, int *nndhk);

/* voisine.c    */

int voisine(int (*dhk)[3],int ndhk,int npnt,int ***nb1,int **nnb1,int **opn1,
            int ***nbtr1,int **nnbtr1,int nnbmax);

/* four1.c     */

void four1(float data[], unsigned long nn, int isign);
