/**
 * @file  AFNI.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:58 $
 *    $Revision: 1.2 $
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


//
// AFNI.h
//
// define afni header structure

typedef struct
{
  // mandatory attributes
  int dataset_rank[2];
  int dataset_dimensions[3];
  char typestring[16];
  int scene_data[3];
  int orient_specific[3];
  float origin[3];
  float delta[3];
  // almost mandatory attributes
  char *idcode_string;
  int numchars;
  char byteorder_string[10];
  float *brick_stats;
  int numstats;
  int *brick_types;
  int numtypes;
  float *brick_float_facs;
  int numfacs;
}
AFNI_HEADER, AF;

MRI *afniRead(char *fname, int read_volume);
int afniWrite(MRI *mri, char *fname);
int readAFNIHeader(FILE *fp, AF *paf);
void AFinit(AF *pAF);
void AFclean(AF *pAF);
void printAFNIHeader(AF *pAF);
