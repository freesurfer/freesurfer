/**
 * @brief utilities for manipulating statistical volumes
 *
 */
/*
 * Original Authors: Bruce Fischl and Doug Greve
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

#define _STATS_SRC

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "const.h"
#include "diag.h"
#include "error.h"
#include "fsgdf.h"
#include "machine.h"
#include "matrix.h"
#include "mghendian.h"
#include "mri.h"
#include "mri_identify.h"
#include "mrinorm.h"
#include "proto.h"
#include "registerio.h"
#include "stats.h"
#include "transform.h"
#include "utils.h"

extern const char *Progname;

#define REG_ROWS 4
#define REG_COLS 4
#define STRUCT_DIM 256

MATRIX *StatLoadTalairachXFM(const char *subjid, const char *xfmfile);

/*--------------------------------------------------------------------*/
// Load output of asegstats2table or aparcstats2table.
STAT_TABLE *LoadStatTable(const char *statfile)
{
  STAT_TABLE *st;
  FILE *fp;
  char tmpstr[100000];
  int r, c, n;

  fp = fopen(statfile, "r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n", statfile);
    return (NULL);
  }

  st = (STAT_TABLE *)calloc(sizeof(STAT_TABLE), 1);
  st->filename = strcpyalloc(statfile);

  // Read in the first line
  fgets(tmpstr, 100000, fp);
  st->ncols = gdfCountItemsInString(tmpstr) - 1;
  if (st->ncols < 1) {
    printf("ERROR: format:  %s\n", statfile);
    return (NULL);
  }
  printf("Found %d data colums\n", st->ncols);

  // Count the number of rows
  st->nrows = 0;
  while (fgets(tmpstr, 100000, fp) != NULL) st->nrows++;
  printf("Found %d data rows\n", st->nrows);
  fclose(fp);

  st->colnames = (char **)calloc(st->ncols, sizeof(char *));
  st->rownames = (char **)calloc(st->nrows, sizeof(char *));

  // OK, now read everything in
  fp = fopen(statfile, "r");

  // Read the measure
  fscanf(fp, "%s", tmpstr);
  st->measure = strcpyalloc(tmpstr);

  // Read the column headers
  for (c = 0; c < st->ncols; c++) {
    fscanf(fp, "%s", tmpstr);
    st->colnames[c] = strcpyalloc(tmpstr);
  }

  // Alloc the data
  st->data = (double **)calloc(st->nrows, sizeof(double *));
  for (r = 0; r < st->nrows; r++) st->data[r] = (double *)calloc(st->ncols, sizeof(double));

  // Read each row
  for (r = 0; r < st->nrows; r++) {
    fscanf(fp, "%s", tmpstr);
    st->rownames[r] = strcpyalloc(tmpstr);
    for (c = 0; c < st->ncols; c++) {
      n = fscanf(fp, "%lf", &(st->data[r][c]));
      if (n != 1) {
        printf("ERROR: format: %s at row %d, col %d\n", statfile, r, c);
        return (NULL);
      }
    }
    // printf("%s %lf\n",st->rownames[r],st->data[r][st->ncols-1]);
  }
  fclose(fp);

  st->mri = MRIallocSequence(st->ncols, 1, 1, MRI_FLOAT, st->nrows);
  st->mri->xsize = 1;
  st->mri->ysize = 1;
  st->mri->zsize = 1;
  st->mri->x_r = 1;
  st->mri->x_a = 0;
  st->mri->x_s = 0;
  st->mri->y_r = 0;
  st->mri->y_a = 1;
  st->mri->y_s = 0;
  st->mri->z_r = 0;
  st->mri->z_a = 0;
  st->mri->z_s = 1;

  for (r = 0; r < st->nrows; r++)
    for (c = 0; c < st->ncols; c++) MRIsetVoxVal(st->mri, c, 0, 0, r, st->data[r][c]);

  return (st);
}

STAT_TABLE *AllocStatTable(int nrows, int ncols)
{
  STAT_TABLE *st;
  int r;

  st = (STAT_TABLE *)calloc(sizeof(STAT_TABLE), 1);
  st->nrows = nrows;
  st->ncols = ncols;
  st->colnames = (char **)calloc(st->ncols, sizeof(char *));
  st->rownames = (char **)calloc(st->nrows, sizeof(char *));

  st->data = (double **)calloc(st->nrows, sizeof(double *));
  for (r = 0; r < st->nrows; r++) st->data[r] = (double *)calloc(st->ncols, sizeof(double));

  return (st);
}

// Write output equivalant of asegstats2table or aparcstats2table.
int WriteStatTable(const char *fname, STAT_TABLE *st)
{
  FILE *fp;
  int err;

  fp = fopen(fname, "w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", fname);
    exit(1);
  }
  err = PrintStatTable(fp, st);
  return (err);
}

int PrintStatTable(FILE *fp, STAT_TABLE *st)
{
  int r, c;

  fprintf(fp, "%-33s ", st->measure);
  for (c = 0; c < st->ncols; c++) fprintf(fp, "%s ", st->colnames[c]);
  fprintf(fp, "\n");
  for (r = 0; r < st->nrows; r++) {
    fprintf(fp, "%-33s ", st->rownames[r]);
    for (c = 0; c < st->ncols; c++) fprintf(fp, "%7.3lf ", st->data[r][c]);
    fprintf(fp, "\n");
  }
  return (0);
}

STAT_TABLE *InitStatTableFromMRI(MRI *mri_in, const char *tablefile)
// sets data from mri
// also (if passed) reads in cols, rows and measure from tablefile (but not the data)
{
  int r, c, ncols, nrows;
  FILE *fp;
  char tmpstr[100000];

  STAT_TABLE *st = AllocStatTable(mri_in->nframes, mri_in->width);
  st->mri = MRIcopy(mri_in, NULL);

  for (r = 0; r < st->nrows; r++)
    for (c = 0; c < st->ncols; c++) st->data[r][c] = MRIgetVoxVal(st->mri, c, 0, 0, r);

  if (tablefile == NULL || strcmp(tablefile, "") == 0) return st;

  // Process template table file:
  fp = fopen(tablefile, "r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n", tablefile);
    return (NULL);
  }
  // Read in the first line
  fgets(tmpstr, 100000, fp);
  ncols = gdfCountItemsInString(tmpstr) - 1;
  printf("Found %d data colums\n", ncols);
  if (ncols != st->ncols) {
    printf("ERROR: Col numbers do not agree in MRI and:  %s\n", tablefile);
    return (NULL);
  }

  // Count the number of rows
  nrows = 0;
  while (fgets(tmpstr, 100000, fp) != NULL) nrows++;
  printf("Found %d data rows\n", nrows);
  if (nrows < st->nrows) {
    printf("ERROR: Not enough row headers for MRI in:  %s\n", tablefile);
    return (NULL);
  }
  if (nrows > st->nrows) {
    printf("WARNING: Too many row headers for MRI in:  %s, will crop ...\n", tablefile);
  }
  fclose(fp);

  // OK, now read everything in
  fp = fopen(tablefile, "r");

  // Read the measure
  fscanf(fp, "%s", tmpstr);
  st->measure = strcpyalloc(tmpstr);

  // Read the column headers
  for (c = 0; c < st->ncols; c++) {
    fscanf(fp, "%s", tmpstr);
    st->colnames[c] = strcpyalloc(tmpstr);
  }

  // Read each row header
  for (r = 0; r < st->nrows; r++) {
    fscanf(fp, "%s", tmpstr);
    st->rownames[r] = strcpyalloc(tmpstr);
    fgets(tmpstr, 100000, fp);  // skip till end of line
  }
  fclose(fp);

  return st;
}

// Stuff below here is not used anymore (really?)

/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
fMRI_REG *StatReadRegistration(const char *fname)
{
  int float2int, err;
  fMRI_REG *reg;
  char *subject;

  reg = (fMRI_REG *)calloc(1, sizeof(fMRI_REG));
  err = regio_read_register(
      fname, &subject, &reg->in_plane_res, &reg->slice_thickness, &reg->brightness_scale, &reg->mri2fmri, &float2int);
  if (err) return (NULL);
  int req = snprintf(reg->name, 100, "%s", subject);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  free(subject);
  reg->fmri2mri = MatrixInverse(reg->mri2fmri, NULL);
  return (reg);
}

/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
int StatFreeRegistration(fMRI_REG **preg)
{
  fMRI_REG *reg;

  reg = *preg;
  MatrixFree(&reg->fmri2mri);
  MatrixFree(&reg->mri2fmri);
  free(reg);
  reg = NULL;  // yes, this leaks a small amount of memory
  return (NO_ERROR);
}

/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
SV *StatReadVolume(const char *prefix)
{
  char path[STRLEN], fname[STRLEN], line[MAX_LINE_LEN], *cp;
  STAT_VOLUME *sv;
  FILE *fp;
  unsigned int nitems;
  int dof_mean, dof_sigma, event_number, slice_number, which_alloc, width, height, nframes, nslices, t, event,
      x, y, z;
  float *buf, fval;
  int DatVersion, DOF;
  float TER;
  char *regfile = NULL;

  FileNamePath(prefix, path);
  sv = (SV *)calloc(1, sizeof(SV));
  if (!sv) ErrorExit(ERROR_NOMEMORY, "StatReadVolume(%s): could not allocate sv", prefix);

  /* read in register.dat */
  if (regfile != NULL) {
    int req = snprintf(fname, STRLEN, "%s", regfile);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  } else {
    int req = snprintf(fname, STRLEN, "%s/register.dat", path);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }

  sv->reg = StatReadRegistration(fname);
  if (!sv->reg) return (NULL);

  /* read the selavg/selxavg dat file, if it exists */
  int req = snprintf(fname, STRLEN, "%s.dat", prefix);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fp = fopen(fname, "r");
  which_alloc = ALLOC_MEANS;
  if (fp) /* means there are time points and means and sigmas */
  {
    fgetl(line, MAX_LINE_LEN - 1, fp);
    sscanf(line, "%*s %f", &sv->tr);
    fgetl(line, MAX_LINE_LEN - 1, fp);
    sscanf(line, "%*s %f", &sv->timewindow);
    fgetl(line, MAX_LINE_LEN - 1, fp);
    sscanf(line, "%*s %f", &sv->prestim);
    fgetl(line, MAX_LINE_LEN - 1, fp);
    sscanf(line, "%*s %d", &sv->nevents);
    fgetl(line, MAX_LINE_LEN - 1, fp);
    sscanf(line, "%*s %d", &sv->time_per_event);
    if (fscanf(fp, "%*s %d", &DatVersion) != EOF) {
      fscanf(fp, "%*s %f", &TER);
      fscanf(fp, "%*s %d", &DOF);
      printf("SelXAvg Format TER = %g, DOF =  %d\n", TER, DOF);
      sv->voltype = 2;
      fprintf(stderr, "INFO: detected volume %s as type selxavg\n", prefix);
    }
    else {
      sv->voltype = 1;
      fprintf(stderr, "INFO: detected volume %s as type selavg\n", prefix);
    }
    fclose(fp);
    which_alloc |= ALLOC_STDS;
  }
  else {
    /*fprintf(stderr,"WARNING: %s: StatReadVolume():\n",Progname);
      fprintf(stderr,"%s does not exist\n",fname);*/
    fprintf(stderr, "INFO: detected volume %s as type raw\n", prefix);
    sv->nevents = 1;
    sv->time_per_event = 0; /* will be filled in later by .hdr file */
    sv->voltype = 0;
  }

  if (sv->nevents > MAX_EVENTS) {
    fprintf(stderr, "ERROR: %s, StatReadVolume():\n", Progname);
    fprintf(stderr, "Number of events (%d) exceeds maximum allowed (%d)\n", sv->nevents, MAX_EVENTS);
    exit(1);
  }

  if (sv->voltype == 1) {
    /* read in the dof file */
    sprintf(fname, "%s_000.dof", prefix);
    fp = fopen(fname, "r");
    if (fp) {
      while ((cp = fgetl(line, MAX_LINE_LEN - 1, fp)) != NULL) {
        sscanf(cp, "%d %d %d", &event_number, &dof_mean, &dof_sigma);
        sv->mean_dofs[event_number] = (float)dof_mean;
        sv->std_dofs[event_number] = (float)dof_sigma;
      }
      fclose(fp);
    }
    else {
      fprintf(stderr, "WARNING: %s: StatReadVolume():\n", Progname);
      fprintf(stderr, "%s does not exist\n", fname);
    }
  }
  else {
    if (sv->voltype == 0) DOF = 1; /* for raw type */
    for (event_number = 0; event_number < sv->nevents; event_number++) {
      sv->mean_dofs[event_number] = (float)DOF + 1;
      sv->std_dofs[event_number] = (float)DOF;
    }
  }

  /* count # of slices */
  slice_number = 0;
  do {
    sprintf(fname, "%s_%3.3d.hdr", prefix, slice_number);
    fp = fopen(fname, "r");
    if (fp) /* this is a valid slice */
    {
      fclose(fp);
      slice_number++;
    }
  } while (fp);

  sv->nslices = nslices = slice_number;

  /* first read .hdr file to get slice dimensions */
  sprintf(fname, "%s_%3.3d.hdr", prefix, 0);
  fp = fopen(fname, "r");
  if (!fp) ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s", fname);
  fscanf(fp, "%d %d %d", &width, &height, &nframes);
  if (!sv->time_per_event) sv->time_per_event = nframes; /* no global .dat file */
  fclose(fp);

  nframes = sv->time_per_event * sv->nevents;
  if (which_alloc & ALLOC_STDS) nframes *= 2.0f;
  if (Gdiag & DIAG_SHOW) fprintf(stderr, "reading %d slices, %d images per slice\n", nslices, nframes);

  buf = (float *)calloc(width * height, sizeof(float));
  if (!buf) ErrorExit(ERROR_NO_MEMORY, "StatReadVolume: could not allocate buffer");

  nframes = sv->time_per_event; /* one frame per time point */

  StatAllocVolume(sv, sv->nevents, width, height, nslices, sv->time_per_event, which_alloc);

  /* read it after nevents */
  if (stricmp(sv->reg->name, "talairach") && stricmp(sv->reg->name, "spherical")) 
    StatReadTransform(sv, sv->reg->name);

  /* read in the actual data */
  nitems = width * height;
  for (z = 0; z < nslices; z++) {
    /* first read .hdr file to get slice dimensions */
    sprintf(fname, "%s_%3.3d.hdr", prefix, z);
    fp = fopen(fname, "r");
    if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "StatReadVolume: could not open %s", fname));
    fscanf(fp, "%d %d %d", &width, &height, &nframes);
    sv->slice_width = width;
    sv->slice_height = height;
    fclose(fp);

    /* now read .bfloat file to parse actual data */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z);
    fp = fopen(fname, "r");
    if (!fp) {
      sprintf(fname, "%s_%3.3d.bshort", prefix, z);
      if (!fp) {
        ErrorReturn(NULL, (ERROR_NOFILE, "StatReadVolume: could not open %s", fname));
      }
      else
        fprintf(stderr, "ERROR: %s does not support bshort volumes\n", Progname);
    }

    for (event = 0; event < sv->nevents; event++) {
      /* read slice of means for each time point */
      for (t = 0; t < sv->time_per_event; t++) {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(NULL,
                      (ERROR_BADFILE,
                       "StatReadVolume: could not read "
                       "%dth slice",
                       z));
        for (y = 0; y < height; y++) {
          for (x = 0; x < width; x++) {
            fval = buf[y * width + x];
#if (BYTE_ORDER == LITTLE_ENDIAN)
            fval = swapFloat(fval);
#endif
            MRIFseq_vox(sv->mri_avgs[event], x, y, z, t) = fval;
          }
        }
      }

      /* read slice of standard deviations for each time point */
      if (sv->mri_stds[event])
        for (t = 0; t < sv->time_per_event; t++) {
          if (fread(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(NULL,
                        (ERROR_BADFILE,
                         "StatReadVolume: could not read "
                         "%dth slice",
                         z));
          for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
              fval = buf[y * width + x];
#if (BYTE_ORDER == LITTLE_ENDIAN)
              fval = swapFloat(fval);
#endif
              MRIFseq_vox(sv->mri_stds[event], x, y, z, t) = fval;
            }
          }
        }
    }

    fclose(fp);

    if (!sv->mri_avg_dofs[event]) continue;

    /* if dof files exist for each slice, read them in */
    sprintf(fname, "%s_dof_%3.3d.bfloat", prefix, z);
    fp = fopen(fname, "r");
    if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "StatReadVolume: could not open %s", fname));

    for (event = 0; event < sv->nevents; event++) {
      /* read slice of mean dofs for each time point */
      for (t = 0; t < sv->time_per_event; t++) {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(NULL,
                      (ERROR_BADFILE,
                       "StatReadVolume: could not read "
                       "%dth slice",
                       z));
        for (y = 0; y < height; y++) {
          for (x = 0; x < width; x++) {
            fval = buf[y * width + x];
#if (BYTE_ORDER == LITTLE_ENDIAN)
            fval = swapFloat(fval);
#endif
            MRIFseq_vox(sv->mri_avg_dofs[event], x, y, z, t) = fval;
          }
        }
      }

      /* read slice of standard deviation dofs for each time point */
      if (sv->mri_stds[event])
        for (t = 0; t < sv->time_per_event; t++) {
          if (fread(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(NULL,
                        (ERROR_BADFILE,
                         "StatReadVolume: could not read "
                         "%dth slice",
                         z));
          for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
              fval = buf[y * width + x];
#if (BYTE_ORDER == LITTLE_ENDIAN)
              fval = swapFloat(fval);
#endif
              MRIFseq_vox(sv->mri_std_dofs[event], x, y, z, t) = fval;
            }
          }
        }
    }

    fclose(fp);
  }

  free(buf);
  return (sv);
}

/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
SV *StatReadVolume2(const char *prefix)
{
  char path[STRLEN], fname[STRLEN], line[MAX_LINE_LEN];
  STAT_VOLUME *sv;
  FILE *fp;
  int event_number, which_alloc, nframes, t;
  int event, x, y, z, f, DatVersion, DOF;
  float fval, TER;
  char *regfile = NULL, *hfile = NULL;
  MRI *h;

  FileNamePath(prefix, path);
  sv = (SV *)calloc(1, sizeof(SV));
  if (!sv) ErrorExit(ERROR_NOMEMORY, "StatReadVolume(%s): could not allocate sv", prefix);

  /* read in register.dat */
  if (regfile != NULL) {
    int req = snprintf(fname, STRLEN, "%s", regfile);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  } else {
    int req = snprintf(fname, STRLEN, "%s/register.dat", path);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }

  sv->reg = StatReadRegistration(fname);
  if (!sv->reg) return (NULL);

  /* read the selavg/selxavg dat file, if it exists */
  sprintf(fname, "%s.dat", prefix);
  fp = fopen(fname, "r");
  which_alloc = ALLOC_MEANS | ALLOC_STDS;
  if (!fp) {
    printf("ERROR: could not open %s\n", fname);
    return (NULL);
  }
  fgetl(line, MAX_LINE_LEN - 1, fp);
  sscanf(line, "%*s %f", &sv->tr);
  fgetl(line, MAX_LINE_LEN - 1, fp);
  sscanf(line, "%*s %f", &sv->timewindow);
  fgetl(line, MAX_LINE_LEN - 1, fp);
  sscanf(line, "%*s %f", &sv->prestim);
  fgetl(line, MAX_LINE_LEN - 1, fp);
  sscanf(line, "%*s %d", &sv->nevents);
  fgetl(line, MAX_LINE_LEN - 1, fp);
  sscanf(line, "%*s %d", &sv->time_per_event);
  if (fscanf(fp, "%*s %d", &DatVersion) != EOF) {
    fscanf(fp, "%*s %f", &TER);
    fscanf(fp, "%*s %d", &DOF);
    printf("SelXAvg Format TER = %g, DOF =  %d\n", TER, DOF);
    sv->voltype = 2;
    fprintf(stderr, "INFO: detected volume %s as type selxavg\n", prefix);
  }
  fclose(fp);

  if (sv->nevents > MAX_EVENTS) {
    fprintf(stderr, "ERROR: %s, StatReadVolume():\n", Progname);
    fprintf(stderr, "Number of events (%d) exceeds maximum allowed (%d)\n", sv->nevents, MAX_EVENTS);
    exit(1);
  }

  if (sv->voltype == 0) DOF = 1; /* for raw type */
  for (event_number = 0; event_number < sv->nevents; event_number++) {
    sv->mean_dofs[event_number] = (float)DOF + 1;
    sv->std_dofs[event_number] = (float)DOF;
  }

  hfile = IDnameFromStem(prefix);
  if (hfile == NULL) return (NULL);
  h = MRIread(hfile);
  if (h == NULL) return (NULL);

  sv->nslices = h->depth;

  nframes = sv->time_per_event * sv->nevents;
  if (which_alloc & ALLOC_STDS) nframes *= 2;

  StatAllocVolume(sv, sv->nevents, h->width, h->height, h->depth, sv->time_per_event, which_alloc);

  /* read it after nevents */
  if (stricmp(sv->reg->name, "talairach") && stricmp(sv->reg->name, "spherical")) {
    StatReadTransform(sv, sv->reg->name);
  }

  f = 0;
  for (event = 0; event < sv->nevents; event++) {
    /* read slice of means for each time point */
    for (t = 0; t < sv->time_per_event; t++) {
      for (z = 0; z < h->depth; z++) {
        for (y = 0; y < h->height; y++) {
          for (x = 0; x < h->width; x++) {
            fval = MRIgetVoxVal(h, x, y, z, f);
            MRIFseq_vox(sv->mri_avgs[event], x, y, z, t) = fval;
          }  // x
        }    // y
      }      // z
      printf("Avgs event %d, t %d, f = %2d\n", event, t, f);
      f++;
    }  // t
    /* read slice of std for each time point */
    for (t = 0; t < sv->time_per_event; t++) {
      for (z = 0; z < h->depth; z++) {
        for (y = 0; y < h->height; y++) {
          for (x = 0; x < h->width; x++) {
            fval = MRIgetVoxVal(h, x, y, z, f);
            MRIFseq_vox(sv->mri_stds[event], x, y, z, t) = fval;
          }  // x
        }    // y
      }      // z
      printf("Stds event %d, t %d, f = %2d\n", event, t, f);
      f++;
    }  // t
  }    // event

  for (event = 0; event < sv->nevents; event++) {
    MRIsetResolution(sv->mri_avgs[event], sv->reg->in_plane_res, sv->reg->in_plane_res, sv->reg->slice_thickness);
    if (sv->mri_stds[event])
      MRIsetResolution(sv->mri_stds[event], sv->reg->in_plane_res, sv->reg->in_plane_res, sv->reg->slice_thickness);
  }

  return (sv);
}

/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
int StatFree(SV **psv)
{
  SV *sv;
  int event, width, height, nslices;

  sv = *psv;
  *psv = NULL;

  width = sv->slice_width;
  height = sv->slice_height;
  nslices = sv->nslices;
  for (event = 0; event < sv->nevents; event++) {
    MRIfree(&sv->mri_avgs[event]);
    if (sv->voltype != 0) MRIfree(&sv->mri_stds[event]);
    if (sv->mri_avg_dofs[event]) MRIfree(&sv->mri_avg_dofs[event]);
    if (sv->mri_std_dofs[event]) MRIfree(&sv->mri_std_dofs[event]);
  }

  delete_general_transform(&sv->transform);
  StatFreeRegistration(&sv->reg);
  free(sv);

  return (NO_ERROR);
}
/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
#if 0
MRI *
StatVolumeToTalairach(SV *sv, MRI *mri, int resolution)
{
  if (!mri)
  {
    int dim ;

    dim = nint((float)STRUCT_DIM / (float)resolution) ;
    mri = MRIallocSequence(dim, dim, dim, MRI_FLOAT, ) ;
  }

  return(mri) ;
}
#endif
/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
STAT_VOLUME *StatAllocVolume(SV *sv, int nevents, int width, int height, int nslices, int time_points, int which_alloc)
{
  int event;

  if (!sv) {
    sv = (SV *)calloc(1, sizeof(SV));
    if (!sv) ErrorExit(ERROR_NOMEMORY, "StatAllocVolume: could not allocate sv");

    sv->reg = (fMRI_REG *)calloc(1, sizeof(fMRI_REG));
    sv->reg->in_plane_res = sv->reg->slice_thickness = 1.0f;
    strcpy(sv->reg->name, "none");
    sv->reg->fmri2mri = MatrixIdentity(4, NULL);
    sv->reg->mri2fmri = MatrixIdentity(4, NULL);
    sv->time_per_event = time_points;
    sv->nevents = nevents;
  }

  for (event = 0; event < sv->nevents; event++) {
    sv->mri_avgs[event] = MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points);
    if (!sv->mri_avgs[event]) ErrorExit(ERROR_NO_MEMORY, "StatAllocVolume: could not allocate volume");
    sv->mri_avgs[event]->x_r = -1.0;
    sv->mri_avgs[event]->y_s = -1.0;
    sv->mri_avgs[event]->z_a = +1.0;
    sv->mri_avgs[event]->c_r = 0.0;
    sv->mri_avgs[event]->c_a = 0.0;
    sv->mri_avgs[event]->c_s = 0.0;

    if (which_alloc & ALLOC_STDS) {
      sv->mri_stds[event] = MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points);
      if (!sv->mri_stds[event])
        ErrorExit(ERROR_NO_MEMORY,
                  "StatAllocVolume: could not allocate "
                  "volume");
    }
    if (which_alloc & ALLOC_DOFS) {
      sv->mri_avg_dofs[event] = MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points);
      if (!sv->mri_avg_dofs[event]) ErrorExit(ERROR_NO_MEMORY, "StatAllocVolume: could not allocate volume");
      sv->mri_std_dofs[event] = MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points);
      if (!sv->mri_std_dofs[event]) ErrorExit(ERROR_NO_MEMORY, "StatAllocVolume: could not allocate volume");
    }
  }
  return (sv);
}
/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
SV *StatAllocStructuralVolume(SV *sv, float fov, float resolution, const char *name)
{
  SV *sv_tal;
  int width, height, depth, event;

  width = height = depth = nint(fov / resolution);
  sv_tal = StatAllocVolume(
      NULL, sv->nevents, width, height, depth, sv->time_per_event, ALLOC_MEANS | ALLOC_STDS | ALLOC_DOFS);

  sv_tal->voltype = sv->voltype;
  strcpy(sv_tal->reg->name, name);
  sv_tal->nslices = depth;
  sv_tal->slice_width = width;
  sv_tal->slice_height = height;
  sv_tal->reg->slice_thickness = sv_tal->reg->in_plane_res = resolution;
  sv_tal->reg->brightness_scale = 2.0f;

  for (event = 0; event < sv->nevents; event++) {
    MRIsetResolution(sv_tal->mri_avgs[event], resolution, resolution, resolution);
    MRIsetResolution(sv_tal->mri_stds[event], resolution, resolution, resolution);
    MRIsetResolution(sv_tal->mri_avg_dofs[event], resolution, resolution, resolution);
    MRIsetResolution(sv_tal->mri_std_dofs[event], resolution, resolution, resolution);
    sv_tal->mri_avgs[event]->x_r = -1.0;
    sv_tal->mri_avgs[event]->y_s = -1.0;
    sv_tal->mri_avgs[event]->z_a = +1.0;
    sv_tal->mri_avgs[event]->c_r = 0.0;
    sv_tal->mri_avgs[event]->c_a = 0.0;
    sv_tal->mri_avgs[event]->c_s = 0.0;
  }

  sv_tal->timewindow = sv->timewindow;
  sv_tal->prestim = sv->prestim;
  sv_tal->tr = sv->tr;
  sv_tal->timewindow = sv->timewindow;
  return (sv_tal);
}
/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
int StatAccumulateSurfaceVolume(SV *sv_surf, SV *sv, MRI_SURFACE *mris)
{
  int x, y, z, width, height, depth, event, t, xv, yv, zv, swidth, sheight, sdepth, vno;
  double xf, yf, zf, xr, yr, zr;
  float mean, surf_mean, std, surf_std, surf_dof, dof, xoff, yoff, zoff, sxoff, syoff, szoff, xs, ys, zs;
  VECTOR *v_struct, *v_func;
  MRI *mri_avg, *mri_std, *mri_ctrl;
  VERTEX *vertex;

  v_func = VectorAlloc(4, MATRIX_REAL);
  v_struct = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v_func, 4) = VECTOR_ELT(v_struct, 4) = 1.0f;

  width = sv_surf->mri_avgs[0]->width;
  height = sv_surf->mri_avgs[0]->height;
  depth = sv_surf->mri_avgs[0]->depth;

  swidth = sv->mri_avgs[0]->width;
  sheight = sv->mri_avgs[0]->height;
  sdepth = sv->mri_avgs[0]->depth;

  xoff = (float)(width - 1) / 2.0f;
  yoff = (float)(height - 1) / 2.0f;
  zoff = (float)(depth - 1) / 2.0f;

  sxoff = (float)(sv->slice_width - 1) / 2.0f;
  syoff = (float)(sv->slice_height - 1) / 2.0f;
  szoff = (float)(sv->nslices - 1) / 2.0f;

  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR);
  mri_avg = MRIallocSequence(width, height, depth, MRI_FLOAT, 2);
  MRIcopyHeader(sv_surf->mri_avgs[0], mri_avg);
  mri_std = MRIallocSequence(width, height, depth, MRI_FLOAT, 2);
  MRIcopyHeader(sv_surf->mri_stds[0], mri_std);

  if (sv_surf->nevents != sv->nevents)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "StatAccumulateTalairachVolume: inconsistent # of events "
                 "(%d vs %d)",
                 sv_surf->nevents,
                 sv->nevents));

  for (event = 0; event < sv_surf->nevents; event++) {
    sv_surf->mean_dofs[event] += sv->mean_dofs[event];
    sv_surf->std_dofs[event] += sv->std_dofs[event];
    if (Gdiag & DIAG_SHOW) fprintf(stderr, "\rprocessing event %d of %d....", event + 1, sv_surf->nevents);

    for (t = 0; t < sv_surf->time_per_event; t++) {
      MRIclear(mri_avg);
      MRIclear(mri_std);
      MRIclear(mri_ctrl);
      /*
        sample from the surface, through the strucural space,
        into the functional
        one.
      */

      /*
        first build average volume for this subjects.
        This is to avoid sampling
        from the volume onto the surface which would take forever since it
        would require searching the entire surface for the nearest point for
        every point in the volume.
      */
      for (vno = 0; vno < mris->nvertices; vno++) {
        vertex = &mris->vertices[vno];

        /*
          read the functional data from the coordinates in 'folded' space,
          then write them into the structural volume using the
          canonical coordinates.
        */
        xs = vertex->x;
        ys = vertex->y;
        zs = vertex->z;

        /* transform it into functional space */
        VECTOR3_LOAD(v_struct, xs, ys, zs);
        MatrixMultiply(sv->reg->mri2fmri, v_struct, v_func);

        /* transform it into a (functional) voxel coordinate */
        MRIworldToVoxel(
            sv->mri_avgs[event], (double)V3_X(v_func), (double)V3_Y(v_func), (double)V3_Z(v_func), &xf, &yf, &zf);
        xv = nint(xf);
        yv = nint(yf);
        zv = nint(zf);

        if (xv >= 0 && xv < swidth && yv >= 0 && yv < sheight && zv >= 0 && zv < sdepth) {
          /* convert from canonical surface coordinate to voxel coordinate */
          xs = vertex->cx;
          ys = vertex->cy;
          zs = vertex->cz;
          MRIworldToVoxel(sv_surf->mri_avgs[event], xs, ys, zs, &xr, &yr, &zr);
          x = nint(xr);
          y = nint(yr);
          z = nint(zr);
          if ((vno == 161) || (x == 15 && y == 20 && z == 3 && t == 0)) DiagBreak();

          if ((x < 0) || (y < 0) || (z < 0) || (z >= mri_avg->depth) || (y >= mri_avg->height) || (x >= mri_avg->width))
            continue;

          /*
            at this point (xv,yv,zv) is a functional coordinate,
            and (x,y,z) is the corresponding structural coordinate.
          */
          /* update means */
          surf_mean = MRIFseq_vox(mri_avg, x, y, z, 0);
          surf_dof = MRIFseq_vox(mri_avg, x, y, z, 1);
          mean = MRIFseq_vox(sv->mri_avgs[event], xv, yv, zv, t);
          surf_mean = (surf_mean * surf_dof + mean) / (surf_dof + 1);
          MRIFseq_vox(mri_avg, x, y, z, 0) = surf_mean;
          MRIFseq_vox(mri_avg, x, y, z, 1) = ++surf_dof;
          MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED;

#if 0
          if (x == 15 && y == 20 && z == 3 && t == 0)
            fprintf(stderr, "adding mean %2.3f, avg = %2.3f\n",mean,surf_mean);
#endif

          /* update stds */
          surf_std = MRIFseq_vox(mri_std, x, y, z, 0);
          surf_dof = MRIFseq_vox(mri_std, x, y, z, 1);
          std = MRIFseq_vox(sv->mri_stds[event], xv, yv, zv, t);

          /* work with variances so things are linear */
          surf_std *= surf_std;
          std *= std;
          surf_std = sqrt((surf_std * surf_dof + std) / (surf_dof + 1));
          MRIFseq_vox(mri_std, x, y, z, 0) = surf_std;
          MRIFseq_vox(mri_std, x, y, z, 1) = ++surf_dof;
        }
      }
      if (getenv("SOAP_STATS")) {
        int niter = atoi(getenv("SOAP_STATS"));

        fprintf(stderr, "performing soap bubble for %d iterations...\n", niter);
        /*        MRIsoapBubbleExpand(mri_avg, mri_ctrl, mri_avg, 1) ;*/
        MRIbuildVoronoiDiagram(mri_avg, mri_ctrl, mri_avg);
        MRIsoapBubble(mri_avg, mri_ctrl, mri_avg, niter, -1);
        MRIbuildVoronoiDiagram(mri_avg, mri_ctrl, mri_avg);
        /*        MRIsoapBubbleExpand(mri_std, mri_ctrl, mri_std, 1) ;*/
        MRIsoapBubble(mri_std, mri_ctrl, mri_std, niter, -1);
      }

      /*
        now use the intermediate structural volumes to update the
        cross-subject average volume.
      */
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          for (x = 0; x < width; x++) {
            if (x == 15 && y == 20 && z == 3 && t == 0) DiagBreak();

            /* update means */
            surf_mean = MRIFseq_vox(sv_surf->mri_avgs[event], x, y, z, t);
            surf_dof = MRIFseq_vox(sv_surf->mri_avg_dofs[event], x, y, z, t);
            mean = MRIFvox(mri_avg, x, y, z);
            dof = sv->mean_dofs[event];
            surf_mean = (surf_mean * surf_dof + mean * dof) / (surf_dof + dof);
            surf_dof += dof;
            MRIFseq_vox(sv_surf->mri_avg_dofs[event], x, y, z, t) = surf_dof;
            MRIFseq_vox(sv_surf->mri_avgs[event], x, y, z, t) = surf_mean;

            /* update stds */
            surf_std = MRIFseq_vox(sv_surf->mri_stds[event], x, y, z, t);
            surf_dof = MRIFseq_vox(sv_surf->mri_std_dofs[event], x, y, z, t);
            std = MRIFvox(mri_std, x, y, z);
            dof = sv->std_dofs[event];

            /* work with variances so things are linear */
            surf_std *= surf_std;
            std *= std;
            surf_std = sqrt((surf_std * surf_dof + std * dof) / (surf_dof + dof));
            surf_dof += dof;
            MRIFseq_vox(sv_surf->mri_std_dofs[event], x, y, z, t) = surf_dof;
            MRIFseq_vox(sv_surf->mri_stds[event], x, y, z, t) = surf_std;
          }
        }
      }
    }
  }

  MRIfree(&mri_std);
  MRIfree(&mri_avg);

  if (Gdiag & DIAG_SHOW) fprintf(stderr, "done.\n");
  VectorFree(&v_func);
  VectorFree(&v_struct);
  return (NO_ERROR);
}
/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
int StatAccumulateTalairachVolume(SV *sv_tal, SV *sv)
{
  int x, y, z, width, height, depth, event, t, xv, yv, zv, swidth, sheight, sdepth;
  // double   xf, yf, zf ;
  float mean, tal_mean, std, tal_std, tal_dof, dof, xoff, yoff, zoff, sxoff, syoff, szoff;
  VECTOR *v_struct, *v_func;
  MRI *mri_avg, *mri_std;
  MATRIX *Tfunc, *Ttal;
  MATRIX *Mcor2tal, *Mtal2cor, *Vtal2func;
  MATRIX *Vtal, *Vfunc;
  float xf2, yf2, zf2;
  extern int stats_fixxfm, statnorm_float2int;

  if (!sv) {
    fprintf(stderr, "ERROR: %s: StatAccumulateTalairachVolume():\n", Progname);
    fprintf(stderr, "Input stat volume is null\n");
    fprintf(stderr, "%s: %d\n", __FILE__, __LINE__);
    exit(1);
  }

  printf("INFO: statnorm_float2int = %d\n", statnorm_float2int);

  v_func = VectorAlloc(4, MATRIX_REAL);
  v_struct = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v_func, 4) = VECTOR_ELT(v_struct, 4) = 1.0f;

  width = sv_tal->mri_avgs[0]->width;
  height = sv_tal->mri_avgs[0]->height;
  depth = sv_tal->mri_avgs[0]->depth;

  swidth = sv->mri_avgs[0]->width;
  sheight = sv->mri_avgs[0]->height;
  sdepth = sv->mri_avgs[0]->depth;

  xoff = (float)(sv_tal->mri_avgs[0]->width - 1) / 2.0f;
  yoff = (float)(sv_tal->mri_avgs[0]->height - 1) / 2.0f;
  zoff = (float)(sv_tal->mri_avgs[0]->depth - 1) / 2.0f;

  sxoff = (float)(sv->slice_width - 1) / 2.0f;
  syoff = (float)(sv->slice_height - 1) / 2.0f;
  szoff = (float)(sv->nslices - 1) / 2.0f;

  if (sv_tal->nevents != sv->nevents)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "StatAccumulateTalairachVolume: inconsistent # of events "
                 "(%d vs %d)",
                 sv_tal->nevents,
                 sv->nevents));

  /*---------------------------------------------------------------*/
  /* This section was added to circumvent MRItalairachVoxelToWorld
     and other functions. Instead, it allows loading of the file
     pointed to by stats_talxfm.*/
  Tfunc = MRIxfmCRS2XYZtkreg(sv->mri_avgs[0]);
  Ttal = MRIxfmCRS2XYZ(sv_tal->mri_avgs[0], 0);
  Mcor2tal = StatLoadTalairachXFM(sv->reg->name, stats_talxfm);
  if (stats_fixxfm) {
    printf("INFO: devolving talairach.xfm\n");
    DevolveXFM(sv->reg->name, Mcor2tal, stats_talxfm);
  }
  Mtal2cor = MatrixInverse(Mcor2tal, NULL);
  Vtal2func = MatrixInverse(Tfunc, NULL);
  MatrixMultiply(Vtal2func, sv->reg->mri2fmri, Vtal2func);
  MatrixMultiply(Vtal2func, Mtal2cor, Vtal2func);
  MatrixMultiply(Vtal2func, Ttal, Vtal2func);

  printf("Tfunc %g %g %g ----------------------\n",
         sv->mri_avgs[0]->xsize,
         sv->mri_avgs[0]->ysize,
         sv->mri_avgs[0]->zsize);
  MatrixPrint(stdout, Tfunc);
  printf("Ttal ---------------------------------------------\n");
  MatrixPrint(stdout, Ttal);
  printf("TalXFM ---------------------------------------------\n");
  MatrixPrint(stdout, Mcor2tal);
  printf("TalVox2FuncVox ---------------------------------------------\n");
  MatrixPrint(stdout, Vtal2func);
  printf("---------------------------------------------\n");

  Vtal = MatrixAlloc(4, 1, MATRIX_REAL);
  Vtal->rptr[4][1] = 1;
  Vfunc = MatrixAlloc(4, 1, MATRIX_REAL);
  Vfunc->rptr[4][1] = 1;
  /*---------------------------------------------------------------*/

  printf("Resampling - nevents = %d, nperevent = %d\n", sv_tal->nevents, sv_tal->time_per_event);
  for (event = 0; event < sv_tal->nevents; event++) {
    printf("Event: %2d\n", event);
    fflush(stdout);
    sv_tal->mean_dofs[event] += sv->mean_dofs[event];
    sv_tal->std_dofs[event] += sv->std_dofs[event];
    mri_avg = sv_tal->mri_avgs[event];
    mri_std = sv_tal->mri_stds[event];
    mri_avg->linear_transform = sv->mri_avgs[event]->linear_transform;
    mri_avg->inverse_linear_transform = sv->mri_avgs[event]->inverse_linear_transform;

    /* Go through each col, row, and slice in the tal volume */
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          Vtal->rptr[1][1] = x;
          Vtal->rptr[2][1] = y;
          Vtal->rptr[3][1] = z;
          MatrixMultiply(Vtal2func, Vtal, Vfunc);
          xf2 = Vfunc->rptr[1][1];
          yf2 = Vfunc->rptr[2][1];
          zf2 = Vfunc->rptr[3][1];

          /* convert from float 2 int */
          switch (statnorm_float2int) {
            case FLT2INT_ROUND:
              xv = nint(xf2);
              yv = nint(yf2);
              zv = nint(zf2);
              break;
            case FLT2INT_TKREG:
              xv = (int)(xf2);
              yv = (int)ceil(yf2);
              zv = (int)(zf2);
              break;
            case FLT2INT_FLOOR:
              xv = (int)(xf2);
              yv = (int)(yf2);
              zv = (int)(zf2);
              break;
            default:
              printf("ERROR: float2int code %d unrecognized\n", statnorm_float2int);
              exit(1);
              break;
          }

          if (xv >= 0 && xv < swidth && yv >= 0 && yv < sheight && zv >= 0 && zv < sdepth) {
            for (t = 0; t < sv_tal->time_per_event; t++) {
              /* update means */
              tal_mean = MRIFseq_vox(sv_tal->mri_avgs[event], x, y, z, t);
              mean = MRIFseq_vox(sv->mri_avgs[event], xv, yv, zv, t);
              dof = sv->mean_dofs[event];
              tal_dof = MRIFseq_vox(sv_tal->mri_avg_dofs[event], x, y, z, t);
              tal_mean = (tal_mean * tal_dof + mean * dof) / (tal_dof + dof);
              tal_dof += dof;
              MRIFseq_vox(sv_tal->mri_avg_dofs[event], x, y, z, t) = tal_dof;
              MRIFseq_vox(sv_tal->mri_avgs[event], x, y, z, t) = tal_mean;

              if (sv->voltype != 0) {
                tal_std = MRIFseq_vox(sv_tal->mri_stds[event], x, y, z, t);
                std = MRIFseq_vox(sv->mri_stds[event], xv, yv, zv, t);
                dof = sv->std_dofs[event];
                tal_dof = MRIFseq_vox(sv_tal->mri_std_dofs[event], x, y, z, t);

                /* work with variances so things are linear */
                tal_std *= tal_std;
                std *= std;
                tal_std = sqrt((tal_std * tal_dof + std * dof) / (tal_dof + dof));
                tal_dof += dof;
                MRIFseq_vox(sv_tal->mri_std_dofs[event], x, y, z, t) = tal_dof;
                MRIFseq_vox(sv_tal->mri_stds[event], x, y, z, t) = tal_std;
              }
            }
          }
        }
      }
    }
  }
  printf("\n");

  if (Gdiag & DIAG_SHOW) fprintf(stderr, "done.\n");
  VectorFree(&v_func);
  VectorFree(&v_struct);
  return (NO_ERROR);
}
/*--------------------------------------------------------------
  ----------------------------------------------------------------*/
int StatWriteVolume(SV *sv, const char *prefix)
{
  char path[STRLEN], fname[STRLEN];
  FILE *fp;
  unsigned int nitems;
  int event_number, width, height, nslices, t, event, x, y, z, nframes;
  float *buf, fval;

  width = sv->slice_width;
  height = sv->slice_height;
  nslices = sv->nslices;
  FileNamePath(prefix, path);
  int req = snprintf(fname, STRLEN, "%s/register.dat", path);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  StatWriteRegistration(sv->reg, fname);

  if (sv->voltype != 0) /* not a raw stats file (sel averaged) */
  {
    /* write the global header file */
    int req = snprintf(fname, STRLEN, "%s.dat", prefix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fp = fopen(fname, "w");
    if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "StatWriteVolume: could not open dat file %s", fname));
    fprintf(fp, "tr %f\n", sv->tr);
    fprintf(fp, "timewindow %f\n", sv->timewindow);
    fprintf(fp, "prestim %f\n", sv->prestim);
    fprintf(fp, "nbins %d\n", sv->nevents);
    fprintf(fp, "perevent %d\n", sv->time_per_event);
    fclose(fp);
  }

  buf = (float *)calloc(sv->slice_width * sv->slice_height, sizeof(float));
  if (!buf) ErrorExit(ERROR_NO_MEMORY, "StatWriteVolume: could not allocate buffer");

  if (sv->voltype == 0)
    nframes = sv->time_per_event * sv->nevents;
  else
    nframes = 2 * sv->time_per_event * sv->nevents;

  /* write out the dof files, the .hdr files and the .bfloat data files */
  nitems = width * height;
  for (z = 0; z < sv->nslices; z++) {
    if (sv->voltype == 1) {
      /* write out .dof file */
      sprintf(fname, "%s_%3.3d.dof", prefix, z);
      fp = fopen(fname, "w");
      if (!fp)
        ErrorReturn(ERROR_NOFILE,
                    (ERROR_NOFILE,
                     "StatWriteVolume: could not open "
                     "dof file %s",
                     fname));
      for (event_number = 0; event_number < sv->nevents; event_number++)
        fprintf(fp, "%d %2.0f %2.0f\n", event_number, sv->mean_dofs[event_number], sv->std_dofs[event_number]);
      fclose(fp);
    }

    /* write out .hdr file */
    sprintf(fname, "%s_%3.3d.hdr", prefix, z);
    fp = fopen(fname, "w");
    if (!fp) ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s", fname);
    fprintf(fp, "%d %d %d 0\n", sv->slice_width, sv->slice_height, nframes);
    fclose(fp);

    /* now write out actual data into .bfloat files */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z);
    fp = fopen(fname, "w");
    if (!fp) ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s", fname);

    for (event = 0; event < sv->nevents; event++) {
      /* write slice of means for each time point */
      for (t = 0; t < sv->time_per_event; t++) {
        for (y = 0; y < sv->slice_height; y++) {
          for (x = 0; x < sv->slice_width; x++) {
            fval = MRIFseq_vox(sv->mri_avgs[event], x, y, z, t);
#if (BYTE_ORDER == LITTLE_ENDIAN)
            fval = swapFloat(fval);
#endif
            buf[y * width + x] = fval;
          }
        }
        if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(ERROR_BADFILE,
                      (ERROR_BADFILE,
                       "StatWriteVolume: could not write "
                       "%dth slice",
                       z));
      }

      if (sv->voltype != 0) {
        /* write slice of standard deviations for each time point */
        for (t = 0; t < sv->time_per_event; t++) {
          for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
              fval = MRIFseq_vox(sv->mri_stds[event], x, y, z, t);
#if (BYTE_ORDER == LITTLE_ENDIAN)
              fval = swapFloat(fval);
#endif
              buf[y * width + x] = fval;
            }
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE,
                        (ERROR_BADFILE,
                         "StatWriteVolume: could not write "
                         "%dth slice",
                         z));
        }
      }
    }

    fclose(fp);

    if (sv->mri_avg_dofs[0]) /* keeping track of dofs on a per-voxel basis */
    {
#if 0 /* don't write out _dof_ file */
      /* now write out dofs on a per-voxel basis */
      sprintf(fname, "%s_dof_%3.3d.bfloat", prefix, z) ;
      fp = fopen(fname, "w") ;
      if (!fp)
        ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s",fname) ;

      for (event = 0 ; event < sv->nevents ; event++)
      {
        /* write slice of means for each time point */
        for (t = 0 ; t < sv->time_per_event ; t++)
        {
          for (y = 0 ; y < sv->slice_height ; y++)
          {
            for (x = 0 ; x < sv->slice_width ; x++)
            {
              fval = MRIFseq_vox(sv->mri_avg_dofs[event],x,y,z,t) ;
#if (BYTE_ORDER == LITTLE_ENDIAN)
              fval = swapFloat(fval) ;
#endif
              buf[y*width+x] = fval ;
            }
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE,
                        (ERROR_BADFILE, "StatWriteVolume: could not write "
                         "%dth slice",z)) ;
        }

        /* write slice of standard deviations for each time point */
        for (t = 0 ; t < sv->time_per_event ; t++)
        {
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
            {
              fval = MRIFseq_vox(sv->mri_std_dofs[event],x,y,z,t) ;
#if (BYTE_ORDER == LITTLE_ENDIAN)
              fval = swapFloat(fval) ;
#endif
              buf[y*width+x] = fval ;
            }
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE,
                        (ERROR_BADFILE, "StatWriteVolume: could not write "
                         "%dth slice", z)) ;
        }
      }

      fclose(fp) ;
#endif
    }
  }

  /* now remove old files which may have had more slices */
  for (; z < sv->nslices * 4; z++) {
    /* write out .dof file */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z);
    unlink(fname);
  }
  free(buf);
  return (NO_ERROR);
}
/*------------------------------------------------------------------
  -------------------------------------------------------------------*/
int StatWriteRegistration(fMRI_REG *reg, const char *fname)
{
  FILE *fp;
  int row, col;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "StatWriteRegistration: could not open %s", fname));
  fprintf(fp, "%s\n", reg->name);
  fprintf(fp, "%f\n", reg->in_plane_res);
  fprintf(fp, "%f\n", reg->slice_thickness);
  fprintf(fp, "%f\n", reg->brightness_scale);
  for (row = 1; row <= REG_ROWS; row++) {
    for (col = 1; col <= REG_COLS; col++) fprintf(fp, "%f  ", reg->fmri2mri->rptr[row][col]);

    fprintf(fp, "\n");
  }
  fclose(fp);
  return (NO_ERROR);
}
/*-------------------------------------------------------------------
  -------------------------------------------------------------------*/
int StatReadTransform(STAT_VOLUME *sv, const char *name)
{
  char *sd, subjects[STRLEN], fname[STRLEN];
  int event;

  /* read in the Talairach transform file */
  sd = getenv("SUBJECTS_DIR");
  if (sd == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined\n");
    exit(1);
  }
  strcpy(subjects, sd);
  int req = snprintf(fname, STRLEN, "%s/%s/mri/transforms/talairach.xfm", subjects, name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  if (input_transform_file(fname, &sv->transform) != OK)
    ErrorPrintf(ERROR_NO_FILE, "%s: could not read xform file '%s'\n", Progname, fname);
  sv->linear_transform = get_linear_transform_ptr(&sv->transform);
  sv->inverse_linear_transform = get_inverse_linear_transform_ptr(&sv->transform);

  for (event = 0; event < sv->nevents; event++) {
    MRIsetTransform(sv->mri_avgs[event], &sv->transform);
    if (sv->mri_stds[event]) MRIsetTransform(sv->mri_stds[event], &sv->transform);
  }
  return (NO_ERROR);
}
/*------------------------------------------------------------------------
  ------------------------------------------------------------------------*/
int StatVolumeExists(const char *prefix)
{
  char fname[STRLEN];
  FILE *fp;

  int req = snprintf(fname, STRLEN, "%s_%3.3d.bfloat", prefix, 0);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fp = fopen(fname, "r");
  if (!fp) return (0);
  fclose(fp);
  return (1);
}
/*------------------------------------------------------------------------
  StatLoadTalairachXFM() - loads the matrix from the given xfm file.
  If the transform is Vox2Vox, converts it to RAS2RAS. The file is
  assumed to exist in SUBJECTS_DIR/subjid/mri/transforms/xfmfile.
  ------------------------------------------------------------------------*/
MATRIX *StatLoadTalairachXFM(const char *subjid, const char *xfmfile)
{
  char subjects[STRLEN], fname[STRLEN];
  MATRIX *Mcor2tal;
  LTA *lta;
  char *cp;
  FILE *fp;

  cp = getenv("SUBJECTS_DIR");
  if (cp) {
    strcpy(subjects, cp);
  } else {
    strcpy(subjects, "~inverse/subjects");
  }
  int req = snprintf(fname, STRLEN, "%s/%s/mri/transforms/%s", subjects, subjid, xfmfile);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: could not open %s for reading \n", fname);
    exit(1);
  }

  lta = LTAread(fname);
  if (lta->type == LINEAR_VOX_TO_VOX) {
    printf("INFO: converting LTA to RAS\n");
    LTAvoxelTransformToCoronalRasTransform(lta);
  }
  Mcor2tal = lta->xforms[0].m_L;

  return (Mcor2tal);
}
FS_STATS *FSstatsRead(char *fname)
{
  FS_STATS *stats;
  char line[MAX_LINE_LEN], *cp, name[STRLEN];
  FILE *fp;
  int n;

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "FSstatsRead(%s): could not open file", fname));

  stats = (FS_STATS *)calloc(1, sizeof(FS_STATS));

  while ((cp = fgetl(line, MAX_LINE_LEN, fp)) != NULL) stats->nlabels++;
  rewind(fp);
  stats->labels = (FS_STAT *)calloc(stats->nlabels, sizeof(FS_STAT));
  for (n = 0; n < stats->nlabels; n++) {
    cp = fgetl(line, MAX_LINE_LEN, fp);
    sscanf(cp,
           "%*d %d %d %lf %s %lf %lf %lf %lf",
           &stats->labels[n].label,
           &stats->labels[n].nvoxels,
           &stats->labels[n].volume,
           name,
           &stats->labels[n].int_mean,
           &stats->labels[n].int_std,
           &stats->labels[n].int_min,
           &stats->labels[n].int_max);
    stats->labels[n].name = (char *)calloc(strlen(name) + 1, sizeof(char));
    strcpy(stats->labels[n].name, name);
  }

  fclose(fp);
  return (stats);
}

int PrintSegStat(FILE *fp, SEGSTAT *segstat)
{
  int n, c;
  char tmpstr[1000];

  fprintf(fp, "# TableCol  1 ColHeader Index \n");
  fprintf(fp, "# TableCol  1 FieldName Index \n");
  fprintf(fp, "# TableCol  1 Units     NA \n");
  segstat->ColHeaders[0] = strcpyalloc("Index");

  fprintf(fp, "# TableCol  2 ColHeader SegId \n");
  fprintf(fp, "# TableCol  2 FieldName Segmentation Id\n");
  fprintf(fp, "# TableCol  2 Units     NA\n");
  segstat->ColHeaders[1] = strcpyalloc("SegId");
  if (!segstat->IsSurf) {
    fprintf(fp, "# TableCol  3 ColHeader NVoxels \n");
    fprintf(fp, "# TableCol  3 FieldName Number of Voxels\n");
    fprintf(fp, "# TableCol  3 Units     unitless\n");
    segstat->ColHeaders[2] = strcpyalloc("NVoxels");
    fprintf(fp, "# TableCol  4 ColHeader Volume_mm3\n");
    fprintf(fp, "# TableCol  4 FieldName Volume\n");
    fprintf(fp, "# TableCol  4 Units     mm^3\n");
    segstat->ColHeaders[3] = strcpyalloc("Volume_mm3");
  }
  else {
    fprintf(fp, "# TableCol  3 ColHeader NVertices \n");
    fprintf(fp, "# TableCol  3 FieldName Number of Vertices\n");
    fprintf(fp, "# TableCol  3 Units     unitless\n");
    segstat->ColHeaders[2] = strcpyalloc("NVertices");
    fprintf(fp, "# TableCol  4 ColHeader Area_mm2\n");
    fprintf(fp, "# TableCol  4 FieldName Area\n");
    fprintf(fp, "# TableCol  4 Units     mm^2\n");
    segstat->ColHeaders[3] = strcpyalloc("Area_mm2");
  }
  c = 5;
  fprintf(fp, "# TableCol %2d ColHeader StructName\n", c);
  fprintf(fp, "# TableCol %2d FieldName Structure Name\n", c);
  fprintf(fp, "# TableCol %2d Units     NA\n", c);
  segstat->ColHeaders[c - 1] = strcpyalloc("StructName");
  c++;

  if (segstat->DoIntensity) {
    fprintf(fp, "# TableCol %2d ColHeader %sMean \n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d FieldName Intensity %sMean\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d Units     %s\n", c, segstat->InIntensityUnits);
    sprintf(tmpstr, "%sMean", segstat->InIntensityName);
    segstat->ColHeaders[c - 1] = strcpyalloc(tmpstr);
    c++;

    fprintf(fp, "# TableCol %2d ColHeader %sStdDev\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d FieldName Itensity %sStdDev\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d Units     %s\n", c, segstat->InIntensityUnits);
    sprintf(tmpstr, "%sStdDev", segstat->InIntensityName);
    segstat->ColHeaders[c - 1] = strcpyalloc(tmpstr);
    c++;

    fprintf(fp, "# TableCol %2d ColHeader %sMin\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d FieldName Intensity %sMin\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d Units     %s\n", c, segstat->InIntensityUnits);
    sprintf(tmpstr, "%sMin", segstat->InIntensityName);
    segstat->ColHeaders[c - 1] = strcpyalloc(tmpstr);
    c++;

    fprintf(fp, "# TableCol %2d ColHeader %sMax\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d FieldName Intensity %sMax\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d Units     %s\n", c, segstat->InIntensityUnits);
    sprintf(tmpstr, "%sMax", segstat->InIntensityName);
    segstat->ColHeaders[c - 1] = strcpyalloc(tmpstr);
    c++;

    fprintf(fp, "# TableCol %2d ColHeader %sRange\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d FieldName Intensity %sRange\n", c, segstat->InIntensityName);
    fprintf(fp, "# TableCol %2d Units     %s\n", c, segstat->InIntensityUnits);
    sprintf(tmpstr, "%sRange", segstat->InIntensityName);
    segstat->ColHeaders[c - 1] = strcpyalloc(tmpstr);
    c++;
  }
  if (segstat->DoSNR) {
    fprintf(fp, "# TableCol %2d ColHeader SNR\n", c);
    fprintf(fp, "# TableCol %2d FieldName Intensity SNR\n", c);
    fprintf(fp, "# TableCol %2d Units     none\n", c);
    sprintf(tmpstr, "%sSNR", segstat->InIntensityName);
    segstat->ColHeaders[c - 1] = strcpyalloc(tmpstr);
    c++;
  }

  fprintf(fp, "# NRows %d \n", segstat->nentries);
  fprintf(fp, "# NTableCols %d \n", c - 1);

  fprintf(fp, "# ColHeaders  ");
  for (n = 0; n < c - 1; n++) fprintf(fp, "%s ", segstat->ColHeaders[n]);
  fprintf(fp, " \n");
  for (n = 0; n < segstat->nentries; n++) {
    fprintf(fp, "%3d %3d  %8d %10.1f  ", n + 1, segstat->entry[n].id, segstat->entry[n].nhits, segstat->entry[n].vol);

    if (segstat->UseName)
      fprintf(fp, "%-30s ", segstat->entry[n].name);
    else
      fprintf(fp, "Seg%04d ", segstat->entry[n].id);

    if (segstat->DoIntensity) {
      fprintf(fp,
              "%10.4f %10.4f %10.4f %10.4f %10.4f ",
              segstat->entry[n].mean,
              segstat->entry[n].std,
              segstat->entry[n].min,
              segstat->entry[n].max,
              segstat->entry[n].range);
      if (segstat->DoSNR) fprintf(fp, "%10.4f ", segstat->entry[n].snr);
    }
    fprintf(fp, "\n");
  }

  return (0);
}
