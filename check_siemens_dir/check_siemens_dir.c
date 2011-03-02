/**
 * @file  check_siemens_dir.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.10 $
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
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include "const.h"
#include "version.h"
#include "mghendian.h"

#define MAX_FILES  10000

struct file {
  char fname[STRLEN];
  int exam, series, image, length;
  short rows, cols, bytes_per_voxel;
};

char *Progname;

#ifndef Darwin
#ifndef SunOS
#ifndef Windows_NT
extern void swab(const void *from, void *to, size_t n);
#endif
#endif
#endif

void check_directory(DIR *dp, char *dir_name);

void usage(void) {

  fprintf(stderr, "usage: %s <siemens direcotry> ...\n", Progname);
  exit(1);

} /* end usage() */

int main(int argc, char *argv[]) {

  int i;
  struct stat stat_buf;
  DIR *dp;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: check_siemens_dir.c,v 1.10 2011/03/02 00:04:01 nicks Exp $", 
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
    if (stat(argv[i], &stat_buf) == -1) {
      printf("can't stat %s\n", argv[i]);
    } else if (!S_ISDIR(stat_buf.st_mode)) {
      printf("%s isn't a directory\n", argv[i]);
    } else {
      if ((dp = opendir(argv[i])) == NULL) {
        printf("error opening directory %s\n", argv[i]);
      } else {
        printf("%s\n", argv[i]);
        check_directory(dp, argv[i]);
        closedir(dp);
      }
    }
  }

  exit(0);

} /* end main() */

void check_directory(DIR *dp, char *dir_name) {

  struct dirent *de;
  int exam, series, image;
  int fname_length;
  struct file files[MAX_FILES];
  int i, j;
  char full_fname[STRLEN];
  int last_series;
  int first_image;
  int last_image;
  FILE *fp;
  short rows, cols, bits_per_voxel;

  for (i = 0;i < MAX_FILES;i++) {
    files[i].image = -1;
  }

  i = 0;
  while ((de = readdir(dp)) != NULL) {
    if (i == MAX_FILES) {
      printf("maximum number of files exceeded (MAX_FILES = %d)\n", MAX_FILES);
      return;
    }
    fname_length = strlen(de->d_name);
    if (fname_length > 4) {
      if (strcmp(&(de->d_name)[fname_length-4], ".ima") == 0) {
        if (sscanf(de->d_name, "%d-%d-%d.ima", &exam, &series, &image) == 3) {
          if (image >= MAX_FILES || image <= 0) {
            printf("image number out of range (file is %s, MAX_FILES is %d)\n",
       de->d_name, MAX_FILES);
          } else {
            sprintf(full_fname, "%s/%s", dir_name, de->d_name);
            strcpy(files[image].fname, full_fname);
            fp = fopen(full_fname, "r");
            if (fp == NULL) {
              printf("error opening file %s; can't check file length\n", 
         full_fname);
              files[image].length = -1;
            } else {
              fseek(fp, 0, SEEK_END);
              files[image].length = ftell(fp);
              fseek(fp, 4994, SEEK_SET);
              fread(&rows, sizeof(short), 1, fp);
              fseek(fp, 4996, SEEK_SET);
              fread(&cols, sizeof(short), 1, fp);
              fseek(fp, 5024, SEEK_SET);
              fread(&bits_per_voxel, sizeof(short), 1, fp);
              fclose(fp);
              // #ifdef Linux
#if (BYTE_ORDER == LITTLE_ENDIAN)
#ifdef SunOS
              swab((const char *)&rows, (char *)&rows, 2);
              swab((const char *)&cols, (char *)&cols, 2);
              swab((const char *)&bits_per_voxel, (char *)&bits_per_voxel, 2);
#else
              swab(&rows, &rows, 2);
              swab(&cols, &cols, 2);
              swab(&bits_per_voxel, &bits_per_voxel, 2);
#endif
#endif
              files[image].rows = rows;
              files[image].cols = cols;
              files[image].bytes_per_voxel = bits_per_voxel / 8;
            }
            files[image].exam = exam;
            files[image].series = series;
            files[image].image = image;
            i++;
          }
        }
      }
    }
  }

  if (i == 0) {
    printf("directory contains no .ima files\n");
    return;
  }

  exam = -1;
  last_series = -1;
  for (i = 0;i < MAX_FILES;i++) {
    if (files[i].image != -1) {
      if (exam == -1)
        exam = files[i].exam;
      else {
        if (files[i].exam != exam) {
          printf("multiple exams in directory\n");
          return;
        }
      }
      if (files[i].series > last_series)
        last_series = files[i].series;
    }
  }

  i = 1;
  if (files[1].image == -1) {
    printf("exam does not start with image 1\n");
    for (;files[i].image != -1;i++);
  }

  first_image = i;

  /* ----- check for missing series ----- */

  last_series = files[first_image].series;

  if (last_series != 1)
    printf("exam does not start with series 1\n");

  for (i = first_image;i < MAX_FILES;i++) {
    if (files[i].image != -1) {
      /* --- ignore jumps back in series -- 
   we'll get these in the out of place image checks --- */
      if (files[i].series == last_series + 1)
        last_series = files[i].series;
      else if (files[i].series > last_series) {
        printf("missing series:");
        for (j = last_series+1;j < files[i].series;j++)
          printf(" %d", j);
        printf("\n");
        last_series = files[i].series;
      }
    }
  }

  /* ----- check for missing images ----- */
  last_image = first_image;
  for (i = first_image;i < MAX_FILES;i++) {
    if (files[i].image != -1) {
      if (files[i].image == last_image + 1)
        last_image = files[i].image;
      else if (files[i].image > last_image + 2) {
        printf("missing images: %d to %d\n", last_image+1, files[i].image - 1);
        last_image = files[i].image;
      } else if (files[i].image > last_image) {
        printf("missing image: %d\n", last_image+1);
        last_image = files[i].image;
      }
    }
  }


  /* ----- check for out of place images ----- */

  last_series = files[first_image].series;

  for (i = first_image;i < MAX_FILES;i++) {
    if (files[i].image != -1) {
      if (files[i].series < last_series) {
        printf("image out of place: %d-%d-%d.ima\n", 
         files[i].exam, files[i].series, files[i].image);
      } else
        last_series = files[i].series;
    }
  }


  /* ----- check file sizes ----- */
  for (i = first_image;i < MAX_FILES;i++) {
    if (files[i].image != -1 && files[image].length != -1) {
      if (files[i].length != 
    files[i].rows * files[i].cols * files[i].bytes_per_voxel + 6144) {
        printf("bad file size for file %d-%d-%d.ima:\n", 
         files[i].exam, files[i].series, files[i].image);
        printf("  6144 (hdr) + %d rows * %d cols * "
         "%d bpv != %d (file length)\n", 
         files[i].rows, 
         files[i].cols, 
         files[i].bytes_per_voxel, 
         files[i].length);
      }
    }
  }

  return;

} /* end check_directory() */

/* EOF */
