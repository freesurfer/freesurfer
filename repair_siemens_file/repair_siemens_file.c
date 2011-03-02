/**
 * @file  repair_siemens_file.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.11 $
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


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "const.h"
#include "error.h"
#include "version.h"
#include "mghendian.h"

#define OLD_APPEND  ".orig"
#define HEADER_LENGTH  6144

char *Progname;

#ifndef Darwin
#ifndef SunOS
#ifndef Windows_NT
extern void swab(const void *from, void *to, size_t n);
#endif
#endif
#endif

int repair_file(char *fname);

void usage(void) {

  fprintf(stderr, "usage: %s <siemens file> ...\n", Progname);
  exit(1);

} /* end usage() */

int main(int argc, char *argv[]) {

  int i;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: repair_siemens_file.c,v 1.11 2011/03/02 00:04:35 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  /* ----- get the basename of the executable ----- */
  Progname = strrchr(argv[0], '/');
  Progname = (Progname == NULL ? argv[0] : Progname + 1);

  if (argc < 2) {
    usage();
  }

  for (i = 1;i < argc;i++) {
    repair_file(argv[i]);
  }

  exit(0);

} /* end main() */

int repair_file(char *fname) {

  FILE *fp;
  short rows, cols;
  short bits_per_voxel, bytes_per_voxel;
  int file_length;
  char header[HEADER_LENGTH];
  char *data;
  int data_bytes;
  char new_fname[STRLEN];

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "can't open file %s for reading\n", fname);
    return(ERROR_BADFILE);
  }

  fseek(fp, 0, SEEK_END);
  file_length = ftell(fp);

  fseek(fp, 4994, SEEK_SET);
  fread(&rows, sizeof(short), 1, fp);
  fseek(fp, 4996, SEEK_SET);
  fread(&cols, sizeof(short), 1, fp);
  fseek(fp, 5024, SEEK_SET);
  fread(&bits_per_voxel, sizeof(short), 1, fp);

#if (BYTE_ORDER == LITTLE_ENDIAN)
#if defined(SunOS)
  swab((const char *)&rows, (char *)&rows, 2);
  swab((const char *)&cols, (char *)&cols, 2);
  swab((const char *)&bits_per_voxel, (char *)&bits_per_voxel, 2);
#else
  swab(&rows, &rows, 2);
  swab(&cols, &cols, 2);
  swab(&bits_per_voxel, &bits_per_voxel, 2);
#endif
#endif

  bytes_per_voxel = bits_per_voxel / 8;

  data_bytes = rows * cols * bytes_per_voxel;

  if (file_length == HEADER_LENGTH + data_bytes) {
    fprintf(stderr, "file %s is the correct length, skipping\n", fname);
    fclose(fp);
    return(ERROR_BADPARM);
  }

  if (file_length < HEADER_LENGTH + data_bytes) {
    fprintf(stderr, "file %s is too short; can't repair\n", fname);
    fclose(fp);
    return(ERROR_BADPARM);
  }

  fseek(fp, 0, SEEK_SET);
  if (fread(&header, 1, HEADER_LENGTH, fp) != HEADER_LENGTH) {
    fprintf(stderr, "error reading header from file %s\n", fname);
    fclose(fp);
    return(ERROR_BADPARM);
  }

  data = (char *)malloc(data_bytes);
  if (data == NULL) {
    fprintf(stderr, "error allocating data memory for file %s\n", fname);
    fclose(fp);
    return(ERROR_NOMEMORY);
  }

  fseek(fp, -data_bytes, SEEK_END);
  if (fread(data, 1, data_bytes, fp) != data_bytes) {
    fprintf(stderr, "error reading data from file %s\n", fname);
    free(data);
    fclose(fp);
    return(ERROR_BADPARM);
  }

  fclose(fp);

  strcpy(new_fname, fname);
  strcat(new_fname, OLD_APPEND);
  if (rename(fname, new_fname) == -1) {
    fprintf(stderr, "error moving file %s to %s\n", fname, new_fname);
    free(data);
    return(ERROR_BADPARM);
  }

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "error opening file %s for writing\n", new_fname);
    free(data);
    return(ERROR_BADPARM);
  }

  if (fwrite(&header, 1, HEADER_LENGTH, fp) != HEADER_LENGTH) {
    fprintf(stderr, "error writing header to file %s\n", fname);
    free(data);
    fclose(fp);
    return(ERROR_BADPARM);
  }

  if (fwrite(data, 1, data_bytes, fp) != data_bytes) {
    fprintf(stderr, "error writing data to file %s\n", fname);
    free(data);
    fclose(fp);
    return(ERROR_BADPARM);
  }

  fclose(fp);

  free(data);

  return(NO_ERROR);

} /* end repair_file() */

/* EOF */
