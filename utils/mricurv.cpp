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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "MC.h"
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "mricurv.h"

#define SMOOTH_MMT 0.25
#define NITERSMOOTH 2
#define MAX_CURV 1
#define AVERAGE_CURV 2
#define MODE MAX_CURV

static float computeLocalCurvature(int *ref_tab, int connectivity);

float Nbhcurvature(Nbh *nbh, int connectivity)
{
  int a, b, c, ref;
  int reference_table[8];
  float curv;

  // building reference table for the 8 cubes
  memset(reference_table, 0, 8 * sizeof(int));
  for (ref = 0, c = 0; c < 2; c++)
    for (b = 0; b < 2; b++)
      for (a = 0; a < 2; a++, ref++) {
        if ((*nbh)[a][b][c] == 1) reference_table[ref] += 1;
        if ((*nbh)[a + 1][b][c] == 1) reference_table[ref] += 2;
        if ((*nbh)[a][b + 1][c] == 1) reference_table[ref] += 4;
        if ((*nbh)[a + 1][b + 1][c] == 1) reference_table[ref] += 8;
        if ((*nbh)[a][b][c + 1] == 1) reference_table[ref] += 16;
        if ((*nbh)[a + 1][b][c + 1] == 1) reference_table[ref] += 32;
        if ((*nbh)[a][b + 1][c + 1] == 1) reference_table[ref] += 64;
        if ((*nbh)[a + 1][b + 1][c + 1] == 1) reference_table[ref] += 128;
      }

  curv = computeLocalCurvature(reference_table, connectivity);
  return curv;
}

float MRIcurvature(MRI *mri, int i, int j, int k, int label, int connectivity)
{
  int a, b, c, ref;
  int reference_table[8];
  float curv;

  // building reference table for the 8 cubes
  memset(reference_table, 0, 8 * sizeof(int));
  for (ref = 0, c = -1; c < 1; c++)
    for (b = -1; b < 1; b++)
      for (a = -1; a < 1; a++, ref++) {
        if (MRIvox(mri, i + a, j + b, k + c) == label) reference_table[ref] += 1;
        if (MRIvox(mri, i + a + 1, j + b, k + c) == label) reference_table[ref] += 2;
        if (MRIvox(mri, i + a, j + b + 1, k + c) == label) reference_table[ref] += 4;
        if (MRIvox(mri, i + a + 1, j + b + 1, k + c) == label) reference_table[ref] += 8;
        if (MRIvox(mri, i + a, j + b, k + c + 1) == label) reference_table[ref] += 16;
        if (MRIvox(mri, i + a + 1, j + b, k + c + 1) == label) reference_table[ref] += 32;
        if (MRIvox(mri, i + a, j + b + 1, k + c + 1) == label) reference_table[ref] += 64;
        if (MRIvox(mri, i + a + 1, j + b + 1, k + c + 1) == label) reference_table[ref] += 128;
      }

  curv = computeLocalCurvature(reference_table, connectivity);
  return curv;
}

typedef struct VerTex
{
  float x, y, z;
  int f[20];
  int fnum;
  int v[20];
  int vnum;
  int marked;
} VerTex;
typedef struct FaCe
{
  int v[3];
  float nx, ny, nz, area;
} FaCe;

static void computeFaceProperties(VerTex *vertex, int nvt, FaCe *face, int nfc);
static void computeCurvature(VerTex *vertex, int nvt, FaCe *face, int nfc, int *ref_tab, int nb, float *curv);

static float fx[12] = {0.5, 0, 0.5, 1, 0, 1, 0, 1, 0.5, 0, 0.5, 1};
static float fy[12] = {0, 0.5, 1, 0.5, 0, 0, 1, 1, 0, 0.5, 1, 0.5};
static float fz[12] = {0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1};

#define MAX_VERTICES 100
#define MAX_FACES 40
static float computeLocalCurvature(int *ref_tab, int connectivity)
{
  int number_of_vertices;
  VerTex vertex[MAX_VERTICES];
  int number_of_faces, reference, refdst;
  FaCe face[MAX_FACES];
  int kept_vertex_table[6], kvt_nbr, save_vertex_indice[30];
  int *Case, vt[12], vind[12], fnbr, n, m, l;
  float x, y, z, xtmp[6], ytmp[6], ztmp[6], curvature[6], maxcurv;

  number_of_vertices = 0;
  number_of_faces = 0;
  memset(save_vertex_indice, -1, sizeof(int) * 30);

  kvt_nbr = 0;
  // first deal with the cube #1
  reference = ref_tab[0];
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = y = z = 0;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;
    // allocate all vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    if (vt[7]) kept_vertex_table[kvt_nbr++] = vind[7];
    if (vt[10]) kept_vertex_table[kvt_nbr++] = vind[10];
    if (vt[11]) kept_vertex_table[kvt_nbr++] = vind[11];
    // save vertices
    if (vt[2]) save_vertex_indice[0] = vind[2];
    if (vt[3]) save_vertex_indice[10] = vind[3];
    if (vt[5]) save_vertex_indice[20] = vind[5];
    if (vt[6]) save_vertex_indice[22] = vind[6];
    if (vt[7]) save_vertex_indice[24] = vind[7];
    if (vt[8]) save_vertex_indice[2] = vind[8];
    if (vt[9]) save_vertex_indice[12] = vind[9];
    if (vt[10]) save_vertex_indice[4] = vind[10];
    if (vt[11]) save_vertex_indice[14] = vind[11];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #2
  reference = (ref_tab[1]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 1;
    y = z = 0;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 1) {
          vind[n] = save_vertex_indice[10];
          continue;
        }
        if (n == 4) {
          vind[n] = save_vertex_indice[20];
          continue;
        }
        if (n == 6) {
          vind[n] = save_vertex_indice[24];
          continue;
        }
        if (n == 9) {
          vind[n] = save_vertex_indice[14];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    if (vt[10]) kept_vertex_table[kvt_nbr++] = vind[10];
    // save vertices
    if (vt[2]) save_vertex_indice[1] = vind[2];
    if (vt[7]) save_vertex_indice[26] = vind[7];
    if (vt[8]) save_vertex_indice[3] = vind[8];
    if (vt[10]) save_vertex_indice[5] = vind[10];
    if (vt[11]) save_vertex_indice[16] = vind[11];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #3
  reference = (ref_tab[2]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 0;
    y = 1;
    z = 0;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 0) {
          vind[n] = save_vertex_indice[0];
          continue;
        }
        if (n == 4) {
          vind[n] = save_vertex_indice[22];
          continue;
        }
        if (n == 5) {
          vind[n] = save_vertex_indice[24];
          continue;
        }
        if (n == 8) {
          vind[n] = save_vertex_indice[4];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    if (vt[11]) kept_vertex_table[kvt_nbr++] = vind[11];
    // save vertices
    if (vt[3]) save_vertex_indice[11] = vind[3];
    if (vt[7]) save_vertex_indice[28] = vind[7];
    if (vt[9]) save_vertex_indice[13] = vind[9];
    if (vt[10]) save_vertex_indice[6] = vind[10];
    if (vt[11]) save_vertex_indice[15] = vind[11];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #4
  reference = (ref_tab[3]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 1;
    y = 1;
    z = 0;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 0) {
          vind[n] = save_vertex_indice[1];
          continue;
        }
        if (n == 1) {
          vind[n] = save_vertex_indice[11];
          continue;
        }
        if (n == 4) {
          vind[n] = save_vertex_indice[24];
          continue;
        }
        if (n == 5) {
          vind[n] = save_vertex_indice[26];
          continue;
        }
        if (n == 6) {
          vind[n] = save_vertex_indice[28];
          continue;
        }
        if (n == 8) {
          vind[n] = save_vertex_indice[5];
          continue;
        }
        if (n == 9) {
          vind[n] = save_vertex_indice[15];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    // none
    // save vertices
    if (vt[10]) save_vertex_indice[7] = vind[10];
    if (vt[11]) save_vertex_indice[17] = vind[11];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #5
  reference = (ref_tab[4]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 0;
    y = 0;
    z = 1;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 0) {
          vind[n] = save_vertex_indice[2];
          continue;
        }
        if (n == 1) {
          vind[n] = save_vertex_indice[12];
          continue;
        }
        if (n == 2) {
          vind[n] = save_vertex_indice[4];
          continue;
        }
        if (n == 3) {
          vind[n] = save_vertex_indice[14];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    if (vt[7]) kept_vertex_table[kvt_nbr++] = vind[7];
    // save vertices
    if (vt[5]) save_vertex_indice[21] = vind[5];
    if (vt[6]) save_vertex_indice[23] = vind[6];
    if (vt[7]) save_vertex_indice[25] = vind[7];
    if (vt[10]) save_vertex_indice[8] = vind[10];
    if (vt[11]) save_vertex_indice[18] = vind[11];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #6
  reference = (ref_tab[5]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 1;
    y = 0;
    z = 1;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 0) {
          vind[n] = save_vertex_indice[3];
          continue;
        }
        if (n == 1) {
          vind[n] = save_vertex_indice[14];
          continue;
        }
        if (n == 2) {
          vind[n] = save_vertex_indice[5];
          continue;
        }
        if (n == 3) {
          vind[n] = save_vertex_indice[16];
          continue;
        }
        if (n == 4) {
          vind[n] = save_vertex_indice[21];
          continue;
        }
        if (n == 6) {
          vind[n] = save_vertex_indice[25];
          continue;
        }
        if (n == 9) {
          vind[n] = save_vertex_indice[18];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    // none
    // save vertices
    if (vt[7]) save_vertex_indice[27] = vind[7];
    if (vt[10]) save_vertex_indice[9] = vind[10];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #7
  reference = (ref_tab[6]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 0;
    y = 1;
    z = 1;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;
    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 0) {
          vind[n] = save_vertex_indice[4];
          continue;
        }
        if (n == 1) {
          vind[n] = save_vertex_indice[13];
          continue;
        }
        if (n == 2) {
          vind[n] = save_vertex_indice[6];
          continue;
        }
        if (n == 3) {
          vind[n] = save_vertex_indice[15];
          continue;
        }
        if (n == 4) {
          vind[n] = save_vertex_indice[23];
          continue;
        }
        if (n == 5) {
          vind[n] = save_vertex_indice[25];
          continue;
        }
        if (n == 8) {
          vind[n] = save_vertex_indice[8];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    // none
    // save vertices
    if (vt[7]) save_vertex_indice[29] = vind[7];
    if (vt[11]) save_vertex_indice[19] = vind[11];
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }

  // cube #8
  reference = (ref_tab[7]);
  if (reference && reference != 255) {
    // position in the cube 3*3*3
    x = 1;
    y = 1;
    z = 1;

    switch (connectivity) {
      case 1:
        Case = MC6p[reference];
        break;
      case 2:
        Case = MC18[reference];
        break;
      case 3:
        Case = MC6[reference];
        break;
      case 4:
        Case = MC26[reference];
        break;
      default:
        Case = MC6p[reference];
        break;
    }
    // number of faces
    fnbr = 0;
    while (Case[3 * fnbr] >= 0) fnbr++;

    memset(vt, 0, 12 * sizeof(int));
    memset(vind, -1, sizeof(int));
    for (n = 0; n < 3 * fnbr; n++) vt[Case[n]]++;

    // find and allocate vertex
    for (n = 0; n < 12; n++)
      if (vt[n]) {
        if (n == 0) {
          vind[n] = save_vertex_indice[5];
          continue;
        }
        if (n == 1) {
          vind[n] = save_vertex_indice[15];
          continue;
        }
        if (n == 2) {
          vind[n] = save_vertex_indice[7];
          continue;
        }
        if (n == 3) {
          vind[n] = save_vertex_indice[17];
          continue;
        }
        if (n == 4) {
          vind[n] = save_vertex_indice[25];
          continue;
        }
        if (n == 5) {
          vind[n] = save_vertex_indice[27];
          continue;
        }
        if (n == 6) {
          vind[n] = save_vertex_indice[29];
          continue;
        }
        if (n == 8) {
          vind[n] = save_vertex_indice[9];
          continue;
        }
        if (n == 9) {
          vind[n] = save_vertex_indice[19];
          continue;
        }

        vertex[number_of_vertices].x = x + fx[n];
        vertex[number_of_vertices].y = y + fy[n];
        vertex[number_of_vertices].z = z + fz[n];
        vertex[number_of_vertices].fnum = 0;
        vind[n] = number_of_vertices++;
      }
    // save kept vertices
    // none
    // save vertices
    // none
    // allocate faces
    for (n = 0; n < fnbr; n++) {
      face[number_of_faces].v[0] = vind[Case[3 * n]];
      face[number_of_faces].v[1] = vind[Case[3 * n + 1]];
      face[number_of_faces++].v[2] = vind[Case[3 * n + 2]];
    }
  }
  if (kvt_nbr == 0) return 0;

  // now allocate face neighbors
  for (n = 0; n < number_of_faces; n++) {
    for (m = 0; m < 3; m++) vertex[face[n].v[m]].f[vertex[face[n].v[m]].fnum++] = n;
  }

  // now allocate vertex neighbors for the kept vertices only
  for (m = 0; m < number_of_vertices; m++) vertex[m].vnum = 0;
  for (n = 0; n < kvt_nbr; n++) {
    for (m = 0; m < number_of_vertices; m++) vertex[m].marked = 0;

    reference = kept_vertex_table[n];
    vertex[reference].marked = 1;
    for (m = 0; m < vertex[reference].fnum; m++)
      for (l = 0; l < 3; l++) {
        refdst = face[vertex[reference].f[m]].v[l];
        if (vertex[refdst].marked == 0) {
          vertex[reference].v[vertex[reference].vnum++] = refdst;
          vertex[refdst].marked = 1;
        }
      }
  }

#if 1
  // now smooth the surface (only the kept vertices)
  l = NITERSMOOTH;
  while (l) {
    for (n = 0; n < kvt_nbr; n++) {
      reference = kept_vertex_table[n];
      xtmp[n] = ytmp[n] = ztmp[n] = 0;
      for (m = 0; m < vertex[reference].vnum; m++) {
        xtmp[n] += vertex[vertex[reference].v[m]].x;
        ytmp[n] += vertex[vertex[reference].v[m]].y;
        ztmp[n] += vertex[vertex[reference].v[m]].z;
      }
      xtmp[n] /= (float)vertex[reference].vnum;
      ytmp[n] /= (float)vertex[reference].vnum;
      ztmp[n] /= (float)vertex[reference].vnum;
    }
    // update
    for (n = 0; n < kvt_nbr; n++) {
      reference = kept_vertex_table[n];
      vertex[reference].x = vertex[reference].x + SMOOTH_MMT * (xtmp[n] - vertex[reference].x);
      vertex[reference].y = vertex[reference].y + SMOOTH_MMT * (ytmp[n] - vertex[reference].y);
      vertex[reference].z = vertex[reference].z + SMOOTH_MMT * (ztmp[n] - vertex[reference].z);
    }
    l--;
  }
#endif

#if 0
  fprintf(stderr,"\nWE HAVE:");
  fprintf(stderr," %d vertices and %d faces ",number_of_vertices,number_of_faces);
  for (n=0;n<number_of_faces;n++)
  {
    fprintf(stderr,"\nface %d (%d,%d,%d):",n,face[n].v[0],
            face[n].v[1],face[n].v[2]);
  }

  for (n=0;n<number_of_vertices;n++)
  {
    fprintf(stderr,"\nvertex %d (%3.2f,%3.2f,%3.2f) -> %d,%d: ",n,vertex[n].x,
            vertex[n].y,vertex[n].z,vertex[n].fnum,vertex[n].vnum);
    for (m=0;m<vertex[n].fnum;m++)
      fprintf(stderr,"f#%d ",vertex[n].f[m]);
    for (m=0;m<vertex[n].vnum;m++)
      fprintf(stderr,"\n     ->v%d (%3.2f,%3.2f,%3.2f)",vertex[n].v[m],vertex[vertex[n].v[m]].x
              ,vertex[vertex[n].v[m]].y,vertex[vertex[n].v[m]].z);
  }
  fprintf(stderr,"\n %d: ",kvt_nbr);
  for (n=0;n<kvt_nbr;n++)
    fprintf(stderr," %d ",kept_vertex_table[n]);
  fprintf(stderr,"\n");
#endif

  // now compute the n(<6) curvatures
  computeFaceProperties(vertex, number_of_vertices, face, number_of_faces);
#if 0
  for (n=0;n<number_of_faces;n++)
  {
    fprintf(stderr,"\nface #%d, %f %f, %f, %f ",n,face[n].area,face[n].nx,
            face[n].ny,face[n].nz);
    if (isnan(face[n].nx)||isnan(face[n].ny)||isnan(face[n].nz)||isnan(face[n].area))
      pause();
  }
#endif

  computeCurvature(vertex, number_of_vertices, face, number_of_faces, kept_vertex_table, kvt_nbr, curvature);

  if (MODE == AVERAGE_CURV) {
    maxcurv = 0;
    for (n = 0; n < kvt_nbr; n++) maxcurv += curvature[n];
    maxcurv /= (float)kvt_nbr;
  }
  else {
    maxcurv = 0;
    m = 0;
    for (n = 0; n < kvt_nbr; n++)
      if (fabs(curvature[n]) > maxcurv) {
        maxcurv = fabs(curvature[n]);
        m = n;
      }
  }
  return curvature[m];
}

static void computeCurvature(VerTex *vertex, int nvt, FaCe *face, int nfc, int *ref_tab, int nb, float *curv)
{
  int n, m, reference;
  VECTOR *v_n, *v_e1, *v_e2, *v;
  float nx, ny, nz, area, dx, dy, dz, y, r2, u1, u2, YR2, R4;

  v_n = VectorAlloc(3, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v = VectorAlloc(3, MATRIX_REAL);
  for (n = 0; n < nb; n++) {
    reference = ref_tab[n];
    // first need to compute normal
    nx = ny = nz = area = 0;
    for (m = 0; m < vertex[reference].fnum; m++) {
      nx += face[vertex[reference].f[m]].nx * face[vertex[reference].f[m]].area;
      ny += face[vertex[reference].f[m]].ny * face[vertex[reference].f[m]].area;
      nz += face[vertex[reference].f[m]].nz * face[vertex[reference].f[m]].area;
      area += face[vertex[reference].f[m]].area;
    }
    nx /= area;
    ny /= area;
    nz /= area;
    VECTOR_LOAD(v_n, nx, ny, nz);
    // now need to compute the tangent plane!
    VECTOR_LOAD(v, ny, nz, nx);
    V3_CROSS_PRODUCT(v_n, v, v_e1);
    if ((V3_LEN_IS_ZERO(v_e1))) {
      if (nz != 0)
        VECTOR_LOAD(v, ny, -nz, nx)
      else if (ny != 0)
        VECTOR_LOAD(v, -ny, nz, nx)
      else
        VECTOR_LOAD(v, ny, nz, -nx);
      V3_CROSS_PRODUCT(v_n, v, v_e1);
    }
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2);
    V3_NORMALIZE(v_e1, v_e1);
    V3_NORMALIZE(v_e2, v_e2);
    // finally compute curvature by fitting a 1-d quadratic r->a*r*r: curv=2*a

    for (YR2 = 0, R4 = 0, m = 0; m < vertex[reference].vnum; m++) {
      dx = vertex[vertex[reference].v[m]].x - vertex[reference].x;
      dy = vertex[vertex[reference].v[m]].y - vertex[reference].y;
      dz = vertex[vertex[reference].v[m]].z - vertex[reference].z;
      VECTOR_LOAD(v, dx, dy, dz);

      y = V3_DOT(v, v_n);
      u1 = V3_DOT(v_e1, v);
      u2 = V3_DOT(v_e2, v);
      r2 = u1 * u1 + u2 * u2;
      YR2 += y * r2;
      R4 += r2 * r2;
    }
    curv[n] = 2 * YR2 / R4;
  }

  VectorFree(&v);
  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
}

static void computeFaceProperties(VerTex *vertex, int nvt, FaCe *face, int nfc)
{
  int n;
  float nx, ny, nz, area, dx1, dx2, dy1, dy2, dz1, dz2;

  for (n = 0; n < nfc; n++) {
    dx1 = vertex[face[n].v[1]].x - vertex[face[n].v[0]].x;
    dy1 = vertex[face[n].v[1]].y - vertex[face[n].v[0]].y;
    dz1 = vertex[face[n].v[1]].z - vertex[face[n].v[0]].z;

    dx2 = vertex[face[n].v[2]].x - vertex[face[n].v[0]].x;
    dy2 = vertex[face[n].v[2]].y - vertex[face[n].v[0]].y;
    dz2 = vertex[face[n].v[2]].z - vertex[face[n].v[0]].z;

    nx = dy1 * dz2 - dy2 * dz1;
    ny = dx2 * dz1 - dx1 * dz2;
    nz = dx1 * dy2 - dx2 * dy1;

    area = sqrt(nx * nx + ny * ny + nz * nz);

    face[n].area = area / 2;
    face[n].nx = nx / area;
    face[n].ny = ny / area;
    face[n].nz = nz / area;
  }
}
