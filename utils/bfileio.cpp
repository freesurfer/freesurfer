/**
 * @brief Routines for handling bfile (bshort and bfloat) I/O.
 *
 * Bfile names are assumed to have the following format:
 * stem_%03d.bext
 * where stem can include a path, and bext is either bshort or bfloat,
 * %03d indicates that the slice number is designated with a 3-digit,
 * zero-padded number. It is always assumed that the slices start with 0.
 *
 * Each Bfile has a corresponding header file with name stem_%03d.hdr.
 * This file has 4 numbers in it: nrows ncols nframes endianness.
 * Endianness = 1 for PCs and 0 for non-PCs.
 *
 * Bfile data within a slice are stored with column as fastest, row next,
 * and frame the slowest.
 *
 * Unless otherwise specified, functions return 0 if there where no
 *  errors and non-0 if there were errors.
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

#define BFILEIO_SRC

#include <string>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

#include "bfileio.h"

extern int errno;

/* ------------------------- */
int byteswapbufdouble(void *buf, long int nbufbytes)
{
  char *cbuf, c;
  long int n, nmax;

  nmax = nbufbytes;
  cbuf = (char *)buf;
  for (n = 0; n < nmax; n += 8) {
    c = *cbuf;
    *cbuf = *(cbuf + 7);
    *(cbuf + 7) = c;

    c = *(cbuf + 1);
    *(cbuf + 1) = *(cbuf + 6);
    *(cbuf + 6) = c;

    c = *(cbuf + 2);
    *(cbuf + 2) = *(cbuf + 5);
    *(cbuf + 5) = c;

    c = *(cbuf + 3);
    *(cbuf + 3) = *(cbuf + 4);
    *(cbuf + 4) = c;

    cbuf += 8;
  }
  return (0);
}

/* ------------------------- */
int byteswapbuffloat(void *buf, long int nbufbytes)
{
  char *cbuf, c;
  long int n, nmax;

  nmax = nbufbytes;
  cbuf = (char *)buf;
  for (n = 0; n < nmax; n += 4) {
    c = *cbuf;
    *cbuf = *(cbuf + 3);
    *(cbuf + 3) = c;

    c = *(cbuf + 1);
    *(cbuf + 1) = *(cbuf + 2);
    *(cbuf + 2) = c;

    cbuf += 4;
  }
  return (0);
}

/* ------------------------- */
int byteswapbufshort(void *buf, long int nbufbytes)
{
  char *cbuf, c;
  long int n, nmax;

  nmax = nbufbytes;
  cbuf = (char *)buf;
  for (n = 0; n < nmax; n += 2) {
    c = *cbuf;
    *cbuf = *(cbuf + 1);
    *(cbuf + 1) = c;
    cbuf += 2;
  }
  return (0);
}
/*---------------------------------------------------------
  Name: bf_getarchendian()
  Determines the endianness of the current computer
  architecture designated by either a 0 or 1, compatible with
  the designation in the .hdr file 0 = non-PC, 1 = PC.
  --------------------------------------------------------*/
int bf_getarchendian(void)
{
  int endian;
  short tmp = 1;
  char *ctmp;

  ctmp = (char *)(&tmp);
  if (*(ctmp + 1) == 1)
    endian = 0;
  else
    endian = 1;
  return (endian);
}
/*---------------------------------------------*/
char *bf_getstemfromname(char *bfname)
{
  extern int bferr;
  extern char bfmsg[BFMSGLEN];
  int len, stemlen;
  char *stem;
  bferr = 0;

  len = (int)strlen(bfname);
  if (len < 7) {
    sprintf(bfmsg, "bf_getstemfromname() %s, wrong format\n", bfname);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    return (NULL);
  }

  stemlen = len - 7;
  stem = (char *)calloc(stemlen + 1, sizeof(char));
  if (stem == NULL) {
    sprintf(bfmsg, "bf_getstemfromname() could not calloc\n");
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    return (NULL);
  }

  memmove(stem, bfname, stemlen);
  return (stem);
}
/*---------------------------------------------*/
int bf_gettypefromname(char *bfname)
{
  extern int bferr;
  extern char bfmsg[BFMSGLEN];
  int len;
  char c;
  bferr = 0;

  len = (int)strlen(bfname);
  if (len < 7) {
    sprintf(bfmsg, "bf_gettypefromname: %s, wrong format\n", bfname);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    return (0);
  }

  c = bfname[len - 5];
  if (c != 's' && c != 'f') {
    sprintf(bfmsg, "bf_gettypefromname: %s, wrong format\n", bfname);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    return (0);
  }

  if (c == 's') return (BF_SHORT);
  if (c == 'f') return (BF_FLOAT);

  return (0); /* should never get here */
}

/*---------------------------------------------*/
int bf_readheader(char *hdrfile, int *nrows, int *ncols, int *nfrms, int *endian)
{
  extern int bferr;
  extern char bfmsg[BFMSGLEN];
  FILE *fp;
  bferr = 0;

  /* open the header file */
  fp = fopen(hdrfile, "r");
  if (fp == NULL) {
    sprintf(bfmsg, "bf_readheader: could not open %s\n", hdrfile);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    perror("");
    return (1);
  }

  /* read the data */
  int scanned = fscanf(fp, "%d %d %d %d", nrows, ncols, nfrms, endian);
  if (scanned != 4) {
    fprintf(stderr, "Warning, scanned %d elements, expected 4\n", scanned);
  }
  fclose(fp);

  return (0);
}
/*---------------------------------------------*/
int bf_writeheader(char *hdrfile, int nrows, int ncols, int nfrms, int endian)
{
  extern int bferr;
  extern char bfmsg[BFMSGLEN];
  FILE *fp;
  bferr = 0;

  /* open the header file */
  fp = fopen(hdrfile, "w");
  if (fp == NULL) {
    sprintf(bfmsg, "bf_writeheader: could not open %s\n", hdrfile);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    perror("");
    return (1);
  }

  /* write the data */
  fprintf(fp, "%d %d %d %d\n", nrows, ncols, nfrms, endian);
  fclose(fp);

  return (0);
}
/*---------------------------------------------*/
int bf_getbfiledim(char *bfname, int *nrows, int *ncols, int *nfrms, int *endian, int *type)
{
  char *stem = NULL;
  char hdrfile[1000];
  int err;

  /* get the type (bshort/bfloat) from the name of the file */
  *type = bf_gettypefromname(bfname);
  if (*type == 0) return (1);

  /* construct the name of the header file from the bfile name*/
  stem = bf_getstemfromname(bfname);
  if (stem == NULL) return (1);
  sprintf(hdrfile, "%s.hdr", stem);
  free(stem);

  /* read the header */
  err = bf_readheader(hdrfile, nrows, ncols, nfrms, endian);
  if (err) return (1);

  return (0);
}
/*---------------------------------------------*/
float *bf_ldbfile(char *bfname, int *nrows, int *ncols, int *nfrms)
{
  extern int bferr;
  extern char bfmsg[BFMSGLEN];
  short *sdata = NULL;
  float *fdata = NULL;
  FILE *fp;
  int err, type, endian, archendian;
  size_t ntot, nread;
  bferr = 0;

  /* get endianness of current architecture */
  archendian = bf_getarchendian();

  /* get the dimensions of the bfile */
  err = bf_getbfiledim(bfname, nrows, ncols, nfrms, &endian, &type);
  if (err) return (NULL);

  /* total number of items in the file */
  ntot = (size_t)(*nrows) * (size_t)(*ncols) * (size_t)(*nfrms);

  /* open the data file */
  fp = fopen(bfname, "r");
  if (fp == NULL) {
    perror("bf_ldbfile");
    sprintf(bfmsg, "bf_ldbfile(): could not open %s\n", bfname);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    return (NULL);
  }

  /* create a buffer to hold the data */
  fdata = (float *)calloc(ntot, sizeof(float));
  if (fdata == NULL) {
    sprintf(bfmsg, "bf_ldbfile(): could not alloc float %d", (int)(ntot * sizeof(float)));
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    fclose(fp);
    return (NULL);
  }

  /* --------------------- bfloat ---------------------------*/
  if (type == BF_FLOAT) {
    /* read the data in */
    nread = fread(fdata, sizeof(float), ntot, fp);
    fclose(fp);
    if (nread != ntot) {
      perror("bf_ldbfile");
      sprintf(bfmsg, "bf_ldbfile(): error reading %s", bfname);
      bferr = 1;
      fprintf(stderr, "%s \n", bfmsg);
      fprintf(stderr, " ntoberead = %d, nread = %d\n", (int)ntot, (int)nread);
      free(fdata);
      return (NULL);
    }
    if (endian != archendian) /* swap bytes if necessary */
      byteswapbuffloat(fdata, (int)ntot * sizeof(float));
  }

  /* --------------------- bshort ---------------------------*/
  if (type == BF_SHORT) {
    /* create a temp short buf */
    sdata = (short *)calloc(ntot, sizeof(short));
    if (sdata == NULL) {
      sprintf(bfmsg, "bf_ldbfile(): could not alloc %d", (int)(ntot * sizeof(short)));
      bferr = 1;
      fprintf(stderr, "%s \n", bfmsg);
      free(fdata);
      return (NULL);
    }
    /* read the data in */
    nread = fread(sdata, sizeof(short), ntot, fp);
    fclose(fp);
    if (nread != ntot) {
      perror("bf_ldbfile");
      sprintf(bfmsg, "bf_ldbfile(): error reading %s", bfname);
      bferr = 1;
      fprintf(stderr, "%s \n", bfmsg);
      fprintf(stderr, " ntoberead = %d, nread = %d\n", (int)ntot, (int)nread);
      free(fdata);
      free(sdata);
      return (NULL);
    }
    if (endian != archendian) /* swap bytes if necessary */
      byteswapbufshort(sdata, ntot * sizeof(short));
    /* convert data from short to float */
    for (unsigned int n = 0; n < ntot; n++) fdata[n] = (float)sdata[n];
    free(sdata);
  }

  return (fdata);
}
/*---------------------------------------------*/
int bf_svbfile(float *bfdata, char *bfname, int nrows, int ncols, int nfrms, int svendian)
{
  extern int bferr;
  extern char bfmsg[BFMSGLEN];
  char *stem = NULL;
  short *sdata = NULL;
  float *fdata = NULL;
  char hdrfile[1000];
  FILE *fp;
  int type, archendian, err;
  int fdatadealloc = 0;
  size_t ntot, nwrote;
  bferr = 0;

  /* get endianness of current architecture */
  archendian = bf_getarchendian();

  /* total number of items in the file */
  ntot = (size_t)nrows * (size_t)ncols * (size_t)nfrms;

  /* construct the name of the header file */
  stem = bf_getstemfromname(bfname);
  if (stem == NULL) return (1);
  sprintf(hdrfile, "%s.hdr", stem);
  free(stem);

  err = bf_writeheader(hdrfile, nrows, ncols, nfrms, svendian);
  if (err) return (1);

  /* get the type (bshort or bfloat) from the file name */
  type = bf_gettypefromname(bfname);
  if (type == 0) return (1);

  /* open the output file */
  fp = fopen(bfname, "w");
  if (fp == NULL) {
    perror("");
    sprintf(bfmsg, "bf_svbfile(): could not open %s\n", bfname);
    bferr = 1;
    fprintf(stderr, "%s \n", bfmsg);
    return (1);
  }

  /*-------------------- bfloat ---------------------------*/
  if (type == BF_FLOAT) {
    if (svendian != archendian) {
      /* copy data to temp buf and swap bytes */
      fdata = (float *)calloc(ntot, sizeof(float));
      if (fdata == NULL) {
        sprintf(bfmsg, "bf_svbfile() could not alloc float %d\n", (int)ntot);
        bferr = 1;
        fprintf(stderr, "%s \n", bfmsg);
        return (1);
      }
      fdatadealloc = 1;
      memmove(fdata, bfdata, ntot * sizeof(float));
      byteswapbuffloat(fdata, ntot * sizeof(float));
    }
    else
      fdata = bfdata;

    /* write the data to the named file */
    nwrote = fwrite(fdata, sizeof(float), ntot, fp);
    if (fdatadealloc) free(fdata);

    /* check that the proper number of items were written */
    if (nwrote != ntot) {
      perror("bf_svbfile");
      sprintf(bfmsg, "bf_svbfile(): error writing to %s", bfname);
      bferr = 1;
      fprintf(stderr, "%s \n", bfmsg);
      fprintf(stderr, " ntobewritten = %d, nwritten = %d\n", (int)ntot, (int)nwrote);
      return (1);
    }
  }

  /*-------------------- bshort ---------------------------*/
  if (type == BF_SHORT) {
    /* copy data into a short buffer */
    sdata = (short *)calloc(ntot, sizeof(short));
    if (sdata == NULL) {
      sprintf(bfmsg, "bf_svbfile(): could not alloc short %d\n", (int)ntot);
      bferr = 1;
      fprintf(stderr, "%s \n", bfmsg);
      return (1);
    }
    for (unsigned int n = 0; n < ntot; n++) sdata[n] = (short)rint(bfdata[n]);

    /* swap bytes if necessary */
    if (svendian != archendian) byteswapbufshort(sdata, ntot * sizeof(short));

    /* write the data */
    nwrote = fwrite(sdata, sizeof(short), ntot, fp);
    free(sdata);

    /* check that the proper number of items were written */
    if (nwrote != ntot) {
      perror("bf_svbfile");
      sprintf(bfmsg, "bf_svbfile(): error writing to %s", bfname);
      bferr = 1;
      fprintf(stderr, "%s \n", bfmsg);
      fprintf(stderr, " ntobewritten = %d, nwritten = %d\n", (int)ntot, (int)nwrote);
      return (1);
    }
  }

  fclose(fp);

  return (0);
}
/*-------------------------------------------------
  Name: bf_getnslices(char *stem)
  Gets number of slices associatd with volume stem by
  counting the number of stem_%03d.hdr files, starting
  with 0.
  ---------------------------------------------------*/
int bf_getnslices(char *stem)
{
  FILE *fp;
  int nslices;
  char bfile[1000];

  memset(bfile, '\0', 1000);

  nslices = 0;
  while (1) {
    sprintf(bfile, "%s_%03d.hdr", stem, nslices);
    fp = fopen(bfile, "r");
    if (!fp) break;
    fclose(fp);
    nslices++;
  }
  /*printf("nslices = %d\n",nslices);*/
  return (nslices);
}

/*---------------------------------------------------------
  Name: bf_getvoltype(char *stem)
  Determines the type of the data (either short or float)
  based on the extension of stem_000. Returns 0 if the
  file does not exist.
  --------------------------------------------------------*/
int bf_getvoltype(char *stem)
{
  FILE *fp;
  char bfile[1000];

  memset(bfile, '\0', 1000);

  sprintf(bfile, "%s_000.%s", stem, "bshort");
  fp = fopen(bfile, "r");
  if (fp != NULL) {
    fclose(fp);
    return (BF_SHORT);
  }

  sprintf(bfile, "%s_000.%s", stem, "bfloat");
  fp = fopen(bfile, "r");
  if (fp != NULL) {
    fclose(fp);
    return (BF_FLOAT);
  }

  fprintf(stderr, "ERROR: bf_getvoltype: cannot find %s\n", stem);
  return (0);
}
/* -------------------------------------------------------- */
int bf_getvoldim(char *stem, int *nrows, int *ncols, int *nslcs, int *nfrms, int *endian, int *type)
{
  char hdrfile[1000];
  int err;

  memset(hdrfile, '\0', 1000);

  sprintf(hdrfile, "%s_000.hdr", stem);
  err = bf_readheader(hdrfile, nrows, ncols, nfrms, endian);
  if (err) return (1);

  *nslcs = bf_getnslices(stem);

  *type = bf_getvoltype(stem);
  if (*type == 0) return (1);

  return (0);
}

/*-----------------------------------------------*/
int bf_dumpvolinfo(FILE *fp, BF_DATA *bfd)
{
  fprintf(fp, "nrows  %d \n", bfd->nrows);
  fprintf(fp, "ncols  %d \n", bfd->ncols);
  fprintf(fp, "nslcs  %d \n", bfd->nslcs);
  fprintf(fp, "nfrms  %d \n", bfd->nfrms);
  return (1);
}

/* ------------------------------------------------- */
int bf_freebfd(BF_DATA **bfd)
{
  int n;

  /* Free the data */
  if ((*bfd)->slcdata != NULL) {
    for (n = 0; n < (*bfd)->nslcs; n++) {
      /* Free each slice of the data */
      if ((*bfd)->slcdata[n] != NULL) {
        free((*bfd)->slcdata[n]);
        (*bfd)->slcdata[n] = NULL;
      }
    }
    free((*bfd)->slcdata);
    (*bfd)->slcdata = NULL;
  }

  /* Finally free the BF_DATA */
  free(*bfd);
  *bfd = NULL;
  return (0);
}

/* -----------------------------------------------------------*/
BF_DATA *bf_preallocbfd(int nrows, int ncols, int nslcs, int nfrms)
{
  BF_DATA *bfd;

  if (nrows == 0 || ncols == 0 || nslcs == 0 || nfrms == 0) {
    fprintf(stderr, "bf_preallocbfd: cannot alloc zero-length volume\n");
    return (NULL);
  }

  bfd = (BF_DATA *)malloc(sizeof(BF_DATA));
  if (bfd == NULL) {
    fprintf(stderr, "bf_preallocbfd: could not alloc bfd\n");
    return (NULL);
  }

  bfd->nrows = nrows;
  bfd->ncols = ncols;
  bfd->nrowcols = bfd->nrows * bfd->ncols;
  bfd->nslcs = nslcs;
  bfd->nfrms = nfrms;

  /* allocate the array of pointers */
  bfd->slcdata = (float **)calloc(bfd->nslcs, sizeof(float *));
  if (bfd->slcdata == NULL) {
    fprintf(stderr, "bf_preallocbfd: could not alloc bfd->slcdata\n");
    free(bfd);
    return (NULL);
  }

  return (bfd);
}
/* -----------------------------------------------------------*/
int bf_iswritable(char *fname)
{
  extern int bferr;
  FILE *fp;
  char tmpstr[2000];
  bferr = 0;

  sprintf(tmpstr, "%s-doo-da-doo-da-day", fname);
  fp = fopen(tmpstr, "w");
  if (fp == NULL) return (0);
  fclose(fp);
  unlink(tmpstr);
  return (1);
}
/* -----------------------------------------------------------*/
int bf_volume_exists(char *stem)
{
  char fname[1000];
  FILE *fp;

  sprintf(fname, "%s_000.hdr", stem);
  fp = fopen(fname, "r");
  if (fp == NULL) return (0);

  fclose(fp);
  return (1);
}
/* -----------------------------------------------------------*/
int bf_delete_volume(char *stem)
{
  char fname[1000];
  int nrows, ncols, nslcs, nfrms, endian, type;
  int err, slice;
  std::string ext;

  if (!bf_volume_exists(stem)) {
    fprintf(stderr, "bf_delete_volume(): volume %s does not exist\n", stem);
    return (1);
  }

  err = bf_getvoldim(stem, &nrows, &ncols, &nslcs, &nfrms, &endian, &type);
  if (err) return (1);

  if (type == BF_SHORT) ext = "short";
  if (type == BF_FLOAT) ext = "float";

  for (slice = 0; slice < nslcs; slice++) {
    sprintf(fname, "%s_%03d.hdr", stem, slice);
    err = unlink(fname);
    sprintf(fname, "%s_%03d.%s", stem, slice, ext.c_str());
    err = unlink(fname);
  }
  return (0);
}
/* -----------------------------------------------------------*/
BF_DATA *bf_allocbfd(int nrows, int ncols, int nslcs, int nfrms)
{
  int slice, nperslice;
  BF_DATA *bfd;

  bfd = bf_preallocbfd(nrows, ncols, nslcs, nfrms);
  if (bfd == NULL) return (NULL);

  /* number of items that must be read for each slice */
  nperslice = bfd->nrows * bfd->ncols * bfd->nfrms;

  for (slice = 0; slice < bfd->nslcs; slice++) {
    /* alloc a slice's worth of float */
    bfd->slcdata[slice] = (float *)calloc(nperslice, sizeof(float));
    if (bfd->slcdata[slice] == NULL) {
      fprintf(stderr, "ERROR: bf_allocbfd: could not alloc %d\n", nperslice);
      bf_freebfd(&bfd);
      return (NULL);
    }
  }

  return (bfd);
}
/*---------------------------------------------------------*/
BF_DATA *bf_ldvolume(char *stem)
{
  BF_DATA *bfd;
  int nrows, ncols, nslcs, nfrms, endian, type;
  int err, slice;
  char bfname[1000];
  std::string ext;
  float *fdata;

  err = bf_getvoldim(stem, &nrows, &ncols, &nslcs, &nfrms, &endian, &type);
  if (err) return (NULL);

  if (type == BF_SHORT)
    ext = "bshort";
  else
    ext = "bfloat";

  bfd = bf_preallocbfd(nrows, ncols, nslcs, nfrms);
  if (bfd == NULL) return (NULL);

  for (slice = 0; slice < nslcs; slice++) {
    memset(bfname, '\0', 1000);
    sprintf(bfname, "%s_%03d.%s", stem, slice, ext.c_str());
    fdata = bf_ldbfile(bfname, &nrows, &ncols, &nfrms);
    if (fdata == NULL) {
      bf_freebfd(&bfd);
      return (NULL);
    }
    bfd->slcdata[slice] = fdata;
  }

  return (bfd);
}
/*---------------------------------------------------------*/
BF_DATA *bf_ldslice(char *stem, int slice)
{
  BF_DATA *bfd;
  int nrows, ncols, nfrms, endian, type;
  int err;
  char bfname[1000];
  std::string ext;
  float *fdata;

  memset(bfname, '\0', 1000);

  type = bf_getvoltype(stem);
  if (type == BF_SHORT)
    ext = "bshort";
  else
    ext = "bfloat";

  sprintf(bfname, "%s_%03d.%s", stem, slice, ext.c_str());

  err = bf_getbfiledim(bfname, &nrows, &ncols, &nfrms, &endian, &type);
  if (err) return (NULL);

  bfd = bf_preallocbfd(nrows, ncols, 1, nfrms);
  if (bfd == NULL) return (NULL);

  fdata = bf_ldbfile(bfname, &nrows, &ncols, &nfrms);
  if (fdata == NULL) {
    bf_freebfd(&bfd);
    return (NULL);
  }
  bfd->slcdata[0] = fdata;

  return (bfd);
}
/*---------------------------------------------------------*/
int bf_svvolume(BF_DATA *bfd, char *stem, int svendian, int svtype)
{
  int err, slice;
  char bfname[1000];
  std::string ext;
  float *fdata;

  memset(bfname, '\0', 1000);

  if (bf_volume_exists(stem)) {
    fprintf(stderr, "INFO: %s volume exists, deleting before saving new volume\n", stem);
    bf_delete_volume(stem);
  }

  if (svtype == BF_SHORT)
    ext = "bshort";
  else
    ext = "bfloat";

  for (slice = 0; slice < bfd->nslcs; slice++) {
    sprintf(bfname, "%s_%03d.%s", stem, slice, ext.c_str());
    fdata = bfd->slcdata[slice];
    err = bf_svbfile(fdata, bfname, bfd->nrows, bfd->ncols, bfd->nfrms, svendian);
    if (err) return (1);
  }

  return (0);
}
/*---------------------------------------------------------*/
int bf_svslice(BF_DATA *bfd, char *stem, int slice, int svendian, int svtype)
{
  int err;
  char bfname[1000];
  std::string ext;
  float *fdata;

  memset(bfname, '\0', 1000);

  if (svtype == BF_SHORT)
    ext = "bshort";
  else if (svtype == BF_FLOAT)
    ext = "bfloat";
  else {
    fprintf(stderr, "bf_svslice: unkown type %d\n", svtype);
    return (1);
  }

  fdata = bfd->slcdata[0];
  sprintf(bfname, "%s_%03d.%s", stem, slice, ext.c_str());
  err = bf_svbfile(fdata, bfname, bfd->nrows, bfd->ncols, bfd->nfrms, svendian);
  if (err) return (1);

  return (0);
}
/* ----------------------------------------------------------
  Name: bf_getval(*bfd, r, c, s, f)
  Returns value of data in row r, column c, slice s, and frame f.
  s refers to the slice in the volume (not with respect to the
  first slice loaded).
  Returns NaN if any of the subscripts are out of range.
  Note: for optimization, use the macro BF_GETVAL(). The macro
  does not do any error checking. If you get errors using the
  macro, define the macro BF_DEBUG; this will cause BF_GETVAL
  to point to bf_getval for easier debugging.
  -------------------------------------------------------------*/
float bf_getval(BF_DATA *bfd, int r, int c, int s, int f)
{
  float val;
  int i;

  if (s < 0 || s >= bfd->nslcs) {
    fprintf(stderr, "ERROR: bf_getval: slice %d out of bounds\n", s);
    return (-10000000000000.0);
  }

  if (r < 0 || r >= bfd->nrows) {
    fprintf(stderr, "ERROR: bf_getval: row %d out of bounds\n", r);
    return (-10000000000000.0);
  }

  if (c < 0 || c >= bfd->ncols) {
    fprintf(stderr, "ERROR: bf_getval: column %d out of bounds\n", r);
    return (-10000000000000.0);
  }

  if (f < 0 || f >= bfd->nfrms) {
    fprintf(stderr, "ERROR: bf_getval: frame %d out of bounds\n", r);
    return (-10000000000000.0);
  }

  i = bf_rcf2index(bfd, r, c, f);

  val = *(bfd->slcdata[s] + i);
  return (val);
}

/* ----------------------------------------------------------
  Name: bf_setval(val,*bfd, r, c, s, f)
  Sets the value of data in row r, column c, slice s, and frame f.
  s refers to the slice in the volume (not with respect to the
  first slice loaded).
  Returns 1 if any of the subscripts are out of range.
  Note: for optimization, use the macro BF_SETVAL(). The macro
  does not do any error checking. If you get errors using the
  macro, define the macro BF_DEBUG; this will cause BF_SETVAL
  to point to bf_setval for easier debugging.
  -------------------------------------------------------------*/
int bf_setval(float val, BF_DATA *bfd, int r, int c, int s, int f)
{
  int i;

  if (s < 0 || s >= bfd->nslcs) {
    fprintf(stderr, "ERROR: bf_setval: slice %d out of bounds\n", s);
    return (1);
  }

  if (r < 0 || r >= bfd->nrows) {
    fprintf(stderr, "ERROR: bf_setval: row %d out of bounds\n", r);
    return (1);
  }

  if (c < 0 || c >= bfd->ncols) {
    fprintf(stderr, "ERROR: bf_setval: column %d out of bounds\n", r);
    return (1);
  }

  if (f < 0 || f >= bfd->nfrms) {
    fprintf(stderr, "ERROR: bf_setval: frame %d out of bounds\n", r);
    return (1);
  }

  i = bf_rcf2index(bfd, r, c, f);

  *(bfd->slcdata[s] + i) = val;
  return (0);
}

/* ----------------------------------------------------------
  Name: bf_rcf2index(*bfd, r, c, f)
  Computes the index in a slice corresonding to row r, column c,
  and frame f.
  -------------------------------------------------------------*/
int bf_rcf2index(BF_DATA *bfd, int r, int c, int f)
{
  int index;

  index = c + r * bfd->ncols + f * bfd->nrowcols;
  return (index);
}
/* ----------------------------------------------------------
  Name: bf_index2rcf(*bfd, index, *r, *c, *f)
  Computes the row r, column c,  and frame f corresponding to
  index in a slice.
  -------------------------------------------------------------*/
int bf_index2rcf(BF_DATA *bfd, int index, int *r, int *c, int *f)

{
  int i = index;

  *f = (int)i / (bfd->nrowcols);
  i = i - *f * (bfd->nrowcols);

  *r = (int)i / (bfd->ncols);
  i = i - *r * bfd->ncols;

  *c = (int)i;

  return (0);
}
/* ----------------------------------------------------------
   bf_get_minmax() - finds the global minimum and maximum.
  ----------------------------------------------------------*/
int bf_get_minmax(BF_DATA *bfd, float *bfdmin, float *bfdmax)
{
  int r, c, s, f;
  float val;

  /* first, find the minimum and maximum */
  *bfdmin = BF_GETVAL(bfd, 0, 0, 0, 0);
  *bfdmax = BF_GETVAL(bfd, 0, 0, 0, 0);
  for (r = 0; r < bfd->nrows; r++) {
    for (c = 0; c < bfd->ncols; c++) {
      for (s = 0; s < bfd->nslcs; s++) {
        for (f = 0; f < bfd->nfrms; f++) {
          val = BF_GETVAL(bfd, r, c, s, f);
          if (*bfdmin > val) *bfdmin = val;
          if (*bfdmax < val) *bfdmax = val;
        }
      }
    }
  }
  return (0);
}
/* ----------------------------------------------------------
   bf_rescale() - rescales the data in a bfile data structure
   to the new min and max.
  ----------------------------------------------------------*/
int bf_rescale(BF_DATA *bfd, float min, float max)
{
  int r, c, s, f;
  float val, bfdmin, bfdmax, bfdrange, range;

  bf_get_minmax(bfd, &bfdmin, &bfdmax);
  bfdrange = bfdmax - bfdmin;
  range = max - min;

  /* now do the actual rescaling */
  for (r = 0; r < bfd->nrows; r++) {
    for (c = 0; c < bfd->ncols; c++) {
      for (s = 0; s < bfd->nslcs; s++) {
        for (f = 0; f < bfd->nfrms; f++) {
          val = (BF_GETVAL(bfd, r, c, s, f) - bfdmin) / bfdrange;
          val = range * val + min;
          BF_SETVAL(val, bfd, r, c, s, f);
        }
      }
    }
  }

  return (0);
}
/**************************************************************/
#if 0 /* exluded on 12/15/00 */
/* ----------------------------------------------------------
  Name: bf_read_register()
  Reads a registration file. While the regsitration file format
  is not part of the bfile format, the registration file often
  accompanies a bfile.

  subject -- name of subject as found in the data base
  ipr -- in-plane resolution
  bpr -- between-plane resolution
  intensity -- for the register program
  R - matrix to convert from xyz in COR space to xyz in BFile space,
      ie, xyzBF = R*xyzCOR
  -------------------------------------------------------------*/
int bf_read_register(char *regfile, char **subject, float *ipr,
                     float *bpr, float *intensity,  float R[4][4])
{
  FILE *fp;
  char tmp[1000];
  int r,c,n;

  fp = fopen(regfile,"r");
  if (fp==NULL)
  {
    perror("bf_read_register()");
    fprintf(stderr,"Could not open %s\n",regfile);
    return(1);
  }

  /* subject name */
  n = fscanf(fp,"%s",tmp);
  if (n != 1)
  {
    perror("bf_read_register()");
    fprintf(stderr,"Error reading subject from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  *subject = (char *) calloc(strlen(tmp)+2,sizeof(char));
  sprintf(*subject,"%s",tmp);

  /* in-plane resolution */
  n = fscanf(fp,"%f",ipr);
  if (n != 1)
  {
    perror("bf_read_register()");
    fprintf(stderr,"Error reading ipr from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* between-plane resolution */
  n = fscanf(fp,"%f",bpr);
  if (n != 1)
  {
    perror("bf_read_register()");
    fprintf(stderr,"Error reading bpr from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* intensity*/
  n = fscanf(fp,"%f",intensity);
  if (n != 1)
  {
    perror("bf_read_register()");
    fprintf(stderr,"Error reading intensity from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* regisration matrix */
  for (r=0;r<4;r++)
  {
    for (c=0;c<4;c++)
    {
      n = fscanf(fp,"%f",&R[r][c]);
      /*printf("%d %d %f\n",r,c,R[r][c]);*/
      if (n != 1)
      {
        perror("bf_read_register()");
        fprintf(stderr,"Error reading R[%d][%d] from %s\n",r,c,regfile);
        fclose(fp);
        return(1);
      }
    }
  }

  fclose(fp);

  return(0);
}
/* ----------------------------------------------------
  bf_qmtx() -- Computes the quantization matrix Q.
  xyz = Q*rcs
  --------------------------------------------------------*/
int bf_qmtx(int nrows, int ncols, int nslcs,
            float rowres, float colres, float slcres,
            float Q[4][4])
{
  int r,c;

  /* first, zero the matrix */
  for (r=0; r<4; r++) for (c=0; c<4; c++) Q[r][c] = 0.0;

  Q[0][0] = -colres;
  Q[0][3] = colres*(ncols-1)/2;

  Q[1][2] = slcres;
  Q[1][3] = -slcres*(nslcs-1)/2;

  Q[2][1] = -rowres;
  Q[2][3] = rowres*(nrows-1)/2;

  Q[3][3] = 1.0;

  return(0);
}


int bf_print_4x4(FILE *fp, char *fmt, float M[4][4])
{
  int r,c;
  char fmt2[100];

  sprintf(fmt2,"%s ",fmt);

  for (r=0; r<4; r++)
  {
    for (c=0; c<4; c++) fprintf(fp,fmt2,M[r][c]);
    fprintf(fp,"\n");
  }

  return(0);
}
#endif
