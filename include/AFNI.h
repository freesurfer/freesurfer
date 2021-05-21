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

MRI *afniRead(const char *fname, int read_volume);
int afniWrite(MRI *mri,const char *fname);
int readAFNIHeader(FILE *fp, AF *paf);
void AFinit(AF *pAF);
void AFclean(AF *pAF);
void printAFNIHeader(AF *pAF);
