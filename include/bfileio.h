/**
 * @brief Routines for handling bfile (bshort and bfloat) I/O
 *
 */
/*
 * Original Author: Douglas Greve
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

#ifndef BF_H_INC
#define BF_H_INC

#include <stdio.h>

#define BF_SHORT 1 /* bshort */
#define BF_FLOAT 2 /* bfloat */


/* for communicating error messages */
#define BFMSGLEN 2000
#ifdef BFILEIO_SRC
int  bferr = 0;
char bfmsg[BFMSGLEN] = {0};
#else
extern int  bferr;
extern char bfmsg[BFMSGLEN];
#endif

typedef struct
{
  int nrows;
  int ncols;
  int nslcs;
  int nfrms;
  float **slcdata; /* all data are represented as floats */
  int nrowcols;    /* rows*cols -- precompute for speed */
}
BF_DATA;

#ifdef BF_DEBUG
/* call function to get/set values */
#define BF_GETVAL bf_getval
#define BF_SETVAL bf_setval
#else
/* optimize getting/setting values */
#define BF_GETVAL(bfd,r,c,s,f)  \
   *(bfd->slcdata[s]+ c + r*bfd->ncols + f*bfd->nrowcols)
#define BF_SETVAL(val,bfd,r,c,s,f)  \
   *(bfd->slcdata[s]+ c + r * bfd->ncols + f * bfd->nrowcols) = val
#endif

/* --------- these are the most useful functions -----------------*/
BF_DATA * bf_ldvolume(char *stem);
int       bf_svvolume(BF_DATA *bfd, char *stem, int svendian, int svtype);
int       bf_iswritable(char *fname);

float bf_getval(BF_DATA *bfd, int r, int c, int s, int f);
int   bf_setval(float val, BF_DATA *bfd, int r, int c, int s, int f);
/* see also BF_SETVAL/BF_GETVAL for bf_setval/bf_getval  */

int bf_freebfd(BF_DATA **bfd);
int bf_rescale(BF_DATA *bfd, float min, float max);
int bf_get_minmax(BF_DATA *bfd, float *bfdmin, float *bfdmax);
/*----------------------------------------------------------------*/

BF_DATA * bf_allocbfd(int nrows, int ncols, int nslcs, int nfrms);

BF_DATA * bf_ldslice(char *stem, int slice);
int bf_svslice(BF_DATA *bfd, char *stem, int slice, int svendian, int svtype);

float * bf_ldbfile(char * bfname, int *nrows, int *ncols, int *nfrms);
int bf_svbfile(float *bfdata, char *bfname, int nrows, int ncols,
               int nfrms, int svendian);

int bf_getnslices(char *stem);

int bf_volume_exists(char *stem);
int bf_delete_volume(char *stem);

int bf_getarchendian(void);
int bf_readheader(char * hdrfile, int *nrows, int *ncols,
                  int *nfrms, int *endian);
int bf_writeheader(char * hdrfile, int nrows, int ncols,
                   int nfrms, int endian);
int bf_getbfiledim(char * bfname, int *nrows,
                   int *ncols, int *nfrms,
                   int *endian, int *type);
int bf_getvoltype(char *stem);
int bf_getvoldim(char *stem, int *nrows,  int *ncols, int *nslcs,
                 int *nfrms, int *endian, int *type);

BF_DATA * bf_preallocbfd(int nrows, int ncols, int nslcs, int nfrms);

int bf_gettypefromname(char *bfname);
char * bf_getstemfromname(char *bfname);

int bf_dumpvolinfo(FILE *fp, BF_DATA *bfd);

int bf_rcf2index(BF_DATA *bfd, int r, int c, int f);
int bf_index2rcf(BF_DATA *bfd, int index, int *r, int *c, int *f);

int byteswapbufdouble(void *buf, long int nbufbytes);
int byteswapbuffloat(void *buf, long int nbufbytes);
int byteswapbufshort(void *buf, long int nbufbytes);


#endif /*BF_H_INC*/
