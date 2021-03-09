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

#include "diag.h"
#include "error.h"
#include "fio.h"
#include "macros.h"
#include "vlabels.h"

extern const char *Progname;

VOXEL_LABELS_IMAGE *VLalloc(int width, int height, int depth, float resolution)
{
  VOXEL_LABELS_IMAGE *vli;
  VOXEL_LABELS ***vl;
  int x, y, z;
  VL *v;

  vli = (VLI *)calloc(1, sizeof(VLI));
  if (!vli) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate VLI *", Progname);
  vli->width = width;
  vli->height = height;
  vli->depth = depth;
  vli->resolution = resolution;
  vl = vli->vl = (VL ***)calloc(width, sizeof(VL **));
  if (!vl) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate VL ***", Progname);
  for (x = 0; x < width; x++) {
    vl[x] = (VL **)calloc(height, sizeof(VL *));
    if (!vl[x]) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate  VL ** %d", Progname, x);
    for (y = 0; y < height; y++) {
      vl[x][y] = (VL *)calloc(depth, sizeof(VL));
      if (!vl[x][y]) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate  VL * %d,%d", Progname, x, y);
    }
  }

  // initialize pointers to be null
  for (x = 0; x < vli->width; x++)
    for (y = 0; y < vli->height; y++)
      for (z = 0; z < vli->depth; z++) {
        v = &vli->vl[x][y][z];
        v->labels = 0;
        v->counts = 0;
      }

  return (vli);
}

int VLfree(VLI **pvli)
{
  VLI *vli;
  int x, y, z;
  VL ***vl;
  VL *v;

  vli = *pvli;
  *pvli = NULL;
  vl = vli->vl;

  // freeup labels and counts first
  for (x = 0; x < vli->width; x++)
    for (y = 0; y < vli->height; y++)
      for (z = 0; z < vli->depth; z++) {
        v = &vli->vl[x][y][z];
        if (v->labels) free(v->labels);
        if (v->counts) free(v->counts);
      }

  // then free up arrays
  for (x = 0; x < vli->width; x++) {
    for (y = 0; y < vli->height; y++) {
      free(vl[x][y]);
    }
    free(vl[x]);
  }
  free(vl);
  free(vli);
  return (NO_ERROR);
}

int VLwrite(VLI *vli, char *fname)
{
  FILE *fp;
  int x, y, z, n;
  long here, there;
  VL *vl;

  fp = fopen(fname, "wb");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "VLwrite: could not open %s", fname));
  fwriteInt(VL_MAGIC, fp);
  fwriteInt(vli->width, fp);
  fwriteInt(vli->height, fp);
  fwriteInt(vli->depth, fp);
  fwriteFloat(vli->resolution, fp);

  /* leave room for lookup table*/
  there = ftell(fp); /* first spot to write offset into */
  fseek(fp, vli->width * sizeof(long), SEEK_CUR);

  for (x = 0; x < vli->width; x++) {
    here = ftell(fp);
    fseek(fp, there, SEEK_SET);
    fwriteInt((int)here, fp);
    fseek(fp, here, SEEK_SET);
    for (y = 0; y < vli->height; y++) {
      for (z = 0; z < vli->depth; z++) {
        vl = &vli->vl[x][y][z];
        fwriteShort(vl->nlabels, fp);
        for (n = 0; n < vl->nlabels; n++) {
          fputc(vl->labels[n], fp);
          fwriteShort(vl->counts[n], fp);
        }
      }
    }
  }

  fclose(fp);
  return (NO_ERROR);
}
VLI *VLread(char *fname)
{
  FILE *fp;
  int x, y, z, n, width, height, depth, magic;
  float resolution;
  VL *vl;
  VLI *vli;

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "VLwrite: could not open %s", fname));
  magic = freadInt(fp);
  if (magic != VL_MAGIC)
    ErrorReturn(NULL, (ERROR_BADFILE, "VLread(%s): not a VL file (magic %x != %x)", fname, magic, VL_MAGIC));

  width = freadInt(fp);
  height = freadInt(fp);
  depth = freadInt(fp);
  resolution = freadFloat(fp);
  vli = VLalloc(width, height, depth, resolution);
  if (!vli) ErrorReturn(NULL, (Gerror, "VLread(%s) failed", fname));

  /* leave room for lookup table*/
  fseek(fp, vli->width * sizeof(long), SEEK_CUR);

  for (x = 0; x < vli->width; x++) {
    for (y = 0; y < vli->height; y++) {
      for (z = 0; z < vli->depth; z++) {
        vl = &vli->vl[x][y][z];
        vl->nlabels = freadShort(fp);
        vl->labels = (unsigned char *)calloc(vl->nlabels, sizeof(unsigned char));
        vl->counts = (unsigned short *)calloc(vl->nlabels, sizeof(unsigned short));
        if (!vl->labels || !vl->counts)
          ErrorExit(ERROR_NOMEMORY, "VLread(%s): could not allocate voxel %d,%d,%d", fname, x, y, z);

        for (n = 0; n < vl->nlabels; n++) {
          vl->labels[n] = fgetc(fp);
          vl->counts[n] = freadShort(fp);
        }
      }
    }
  }

  fclose(fp);
  return (vli);
}

VL *VLreadVoxel(char *fname, int x, int y, int z, VL *vl) { return (vl); }

int VLnormalize(VLI *vli)
{
  int x, y, z, n;
  VL *vl;
  float pct, total;

  for (x = 0; x < vli->width; x++) {
    for (y = 0; y < vli->height; y++) {
      for (z = 0; z < vli->depth; z++) {
        vl = &vli->vl[x][y][z];
        for (total = 0.0, n = 0; n < vl->nlabels; n++) total += (float)vl->counts[n];

        for (n = 0; n < vl->nlabels; n++) {
          pct = 100.0f * (float)vl->counts[n] / total;
          vl->counts[n] = nint(pct);
        }
      }
    }
  }
  return (NO_ERROR);
}
