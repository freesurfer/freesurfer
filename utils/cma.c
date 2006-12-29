/**
 * @file  cma.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:31 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "error.h"

/* see ch notebook 2 */

#include "cma.h"

#ifndef TRUE
#define TRUE (1==1)
#endif

#ifndef FALSE
#define FALSE (0==1)
#endif

extern int errno;

int CMAfreeOutlineField(CMAoutlineField **of)
{

  CMAoutlineField *ofp;
  int i;

  ofp = *of;

  if (ofp->claim_field != NULL)
  {
    for (i = 0;i < ofp->height;i++)
    {
      if (ofp->claim_field[i] != NULL)
        free(ofp->claim_field[i]);
    }
    free(ofp->claim_field);
  }

  if (ofp->fill_field != NULL)
  {
    for (i = 0;i < ofp->height;i++)
    {
      if (ofp->fill_field[i] != NULL)
        free(ofp->fill_field[i]);
    }
    free(ofp->fill_field);
  }

  if (ofp->outline_points_field != NULL)
  {
    for (i = 0;i < ofp->height;i++)
    {
      if (ofp->outline_points_field[i] != NULL)
        free(ofp->outline_points_field[i]);
    }
    free(ofp->outline_points_field);
  }

  *of = NULL;

  return(NO_ERROR);

} /* end CMAfreeOutlineField() */

CMAoutlineField *CMAoutlineFieldAlloc(int width, int height)
{

  CMAoutlineField *of;
  int i;

  of = (CMAoutlineField *)malloc(sizeof(CMAoutlineField));
  if (of == NULL)
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating structure"));

  of->claim_field = NULL;
  of->fill_field = NULL;
  of->outline_points_field = NULL;
  of->width = width;
  of->height = height;

  of->claim_field = (CMAoutlineClaim **)malloc(height * sizeof(CMAoutlineClaim *));
  if (of->claim_field == NULL)
  {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating claim field"));
  }
  memset(of->claim_field, 0x00, height * sizeof(CMAoutlineClaim *));

  of->fill_field = (unsigned char **)malloc(height * sizeof(unsigned char *));
  if (of->fill_field == NULL)
  {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating fill field"));
  }
  memset(of->fill_field, 0x00, height * sizeof(unsigned char *));

  of->outline_points_field = (unsigned char **)malloc(height * sizeof(unsigned char *));
  if (of->outline_points_field == NULL)
  {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating outline points field"));
  }
  memset(of->outline_points_field, 0x00, height * sizeof(unsigned char *));

  for (i = 0;i < height;i++)
  {

    of->claim_field[i] = (CMAoutlineClaim *)malloc(width * sizeof(CMAoutlineClaim));
    if (of->claim_field[i] == NULL)
    {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating claim field"));
    }
    memset(of->claim_field[i], 0x00, width * sizeof(CMAoutlineClaim));

    of->fill_field[i] = (unsigned char *)malloc(width * sizeof(unsigned char));
    if (of->fill_field[i] == NULL)
    {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating fill field"));
    }
    memset(of->fill_field[i], 0x00, width * sizeof(unsigned char));

    of->outline_points_field[i] = (unsigned char *)malloc(width * sizeof(unsigned char));
    if (of->outline_points_field[i] == NULL)
    {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating outline points field"));
    }
    memset(of->outline_points_field[i], 0x00, width * sizeof(unsigned char));

  }

  return(of);

} /* end CMAoutlineFieldAlloc() */

int CMAclearFillField(CMAoutlineField *field)
{

  int i;

  for (i = 0;i < field->height;i++)
    memset(field->fill_field[i], 0x00, field->width * sizeof(unsigned char));

  return(NO_ERROR);

} /* end CMAclearFillField() */

/* fills with CMA_FILL_INTERIOR to non-zero values */
int CMAfill(CMAoutlineField *field, short seed_x, short seed_y)
{

  if (seed_x < 0 || seed_x >= field->width)
    return(NO_ERROR);

  if (seed_y < 0 || seed_y >= field->height)
    return(NO_ERROR);

  if (field->fill_field[seed_y][seed_x] != 0)
    return(NO_ERROR);

  field->fill_field[seed_y][seed_x] = CMA_FILL_INTERIOR;

  CMAfill(field, seed_x + 1, seed_y    );
  CMAfill(field, seed_x - 1, seed_y    );
  CMAfill(field, seed_x    , seed_y + 1);
  CMAfill(field, seed_x    , seed_y - 1);

  return(NO_ERROR);

} /* end CMAfill() */

int CMAclaimPoints(CMAoutlineField *field, short label, short *points, int n_points, short seed_x, short seed_y)
{

  int i, j;
  short x, y;

  if (label < 0 || label > MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): label out of range (label is %d, MAX_CMA_LABEL is %d)", label, MAX_CMA_LABEL));

  if (seed_x < 0 || seed_x >= field->width)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): seed point out of range (seed_x = %d, field width = %d)", seed_x, field->width));
  if (seed_y < 0 || seed_y >= field->height)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): seed point out of range (seed_y = %d, field height = %d)", seed_y, field->height));

  CMAclearFillField(field);

  for (i = 0;i < n_points;i++)
  {
    x = points[2*i];
    y = points[2*i+1];
    if (x < 0 || x >= field->width)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): outline point out of range (x)"));
    if (y < 0 || y >= field->height)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): outline point out of range (y)"));
    field->fill_field[y][x] = CMA_FILL_OUTLINE;
    field->outline_points_field[y][x] = 1;
  }

  CMAfill(field, seed_x, seed_y);

  for (i = 0;i < field->width;i++)
  {
    for (j = 0;j < field->height;j++)
    {

      if (field->fill_field[j][i] == CMA_FILL_INTERIOR)
      {
        field->claim_field[j][i].n_claims = 1;
        field->claim_field[j][i].interior_claim_flag = TRUE;
        field->claim_field[j][i].claim_labels[0] = label;
      }

      if (field->fill_field[j][i] == CMA_FILL_OUTLINE)
      {
        if (field->claim_field[j][i].n_claims < MAX_OUTLINE_CLAIMS)
        {
          field->claim_field[j][i].claim_labels[field->claim_field[j][i].n_claims] = label;
          field->claim_field[j][i].n_claims++;
          field->claim_field[j][i].interior_claim_flag = FALSE;
        }
      }

    }
  }

  return(NO_ERROR);

} /* end CMAclaimPoints() */

int CMAvalueClaims(CMAoutlineClaim *claim)
{

  int i;

  for (i = 0;i < MAX_OUTLINE_CLAIMS;i++)
    claim->claim_values[i] = 0.0;

  if (claim->n_claims == 0)
  {
    claim->no_label_claim = 0.5;
  }
  else if (claim->n_claims == 1)
  {
    if (claim->interior_claim_flag == 1)
      claim->claim_values[0] = 1.0;
    else
    {
      claim->claim_values[0] = 0.5;
      claim->no_label_claim = 0.5;
    }
  }
  else
  {
    float ct = 1.0 / (float)claim->n_claims;
    for (i = 0;i < claim->n_claims;i++)
      claim->claim_values[i] = ct;
  }

  return(NO_ERROR);

} /* end CMAvalueClaims() */

int CMAvalueAllClaims(CMAoutlineField *field)
{

  int i, j;

  for (i = 0;i < field->width;i++)
    for (j = 0;j < field->height;j++)
      CMAvalueClaims(&(field->claim_field[j][i]));

  return(NO_ERROR);

} /* end CMAvalueAllClaims() */

int CMAaddWeightedTotals(CMAoutlineClaim *claim, float weight, float *claim_totals)
{

  int i;

  /* we've checked the label range in CMAclaimPoints() */

  for (i = 0;i < claim->n_claims;i++)
    claim_totals[claim->claim_labels[i]] += weight * claim->claim_values[i];

  claim_totals[0] += weight * claim->no_label_claim;

  return(NO_ERROR);

} /* end CMAaddWeightedTotals() */

/* returns the label with the greatest claim to a given pixel */
short CMAtotalClaims(CMAoutlineField *field, int x, int y)
{

  float claim_totals[MAX_CMA_LABEL+1];
  float best_claim;
  short best_index;
  int i;

  if (x < 0 || x >= field->width)
    ErrorReturn(-1, (ERROR_BADPARM, "CMAtotalClaims(): x out of range (x = %d, field->width = %d)", x, field->width));

  if (y < 0 || y >= field->height)
    ErrorReturn(-1, (ERROR_BADPARM, "CMAtotalClaims(): y out of range (y = %d, field->height = %d)", y, field->height));

  for (i = 0;i < MAX_CMA_LABEL;i++)
    claim_totals[i] = 0.0;

  /* center point (x, y checks done above) */
  CMAaddWeightedTotals(&(field->claim_field[y][x]), OTL_CLAIM_WEIGHT_CENTER, claim_totals);

  /* immediately adjoining points */
  if (x-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y  ][x-1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x+1 < field->width)
    CMAaddWeightedTotals(&(field->claim_field[y  ][x+1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (y-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y-1][x  ]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (y+1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y+1][x  ]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);

  /* diagonally adjoining points */
  if (x-1 >= 0 && y-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y-1][x-1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x-1 >= 0 && y+1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y+1][x-1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x+1 < field->width && y-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y-1][x+1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x+1 < field->width && y+1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y+1][x+1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);

  /* find the highest claim and its index */
  best_claim = claim_totals[0];
  best_index = 0;
  for (i = 1;i <= MAX_CMA_LABEL;i++)
  {
    if (claim_totals[i] > best_claim)
    {
      best_claim = claim_totals[i];
      best_index = i;
    }
  }

  return(best_index);

} /* end CMAtotalClaims() */

int CMAassignLabels(CMAoutlineField *field)
{

  int i, j;

  CMAclearFillField(field);

  CMAvalueAllClaims(field);

  for (i = 0;i < field->width;i++)
    for (j = 0;j < field->height;j++)
      field->fill_field[j][i] = CMAtotalClaims(field, i, j);

  return(NO_ERROR);

} /* end CMAassignLabels() */

int CMAzeroOutlines(CMAoutlineField *field)
{

  int i, j;

  for (i = 0;i < field->width;i++)
    for (j = 0;j < field->height;j++)
      if (field->outline_points_field[j][i])
        field->fill_field[j][i] = 0;

  return(NO_ERROR);

} /* end CMAzeroOutlines() */

/* eof */
