/**
 * @brief Routine to check .license file
 *
 */
/*
 * Original Author: Bruce Fischl
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

#ifndef Darwin
#include <gnu/libc-version.h>
#endif

#include <chklc.h>
#include <const.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include "diag.h"

extern char *crypt(const char *, const char *);

static const char *errmsg =
    "--------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer environment FREESURFER_HOME is not defined.\n"
    "  If you are outside the NMR-Martinos Center, please set this\n"
    "  variable to the location where you installed FreeSurfer.\n"
    "  If you are inside the NMR-Martinos Center, please source\n"
    "  the standard environment. If you need to install FreeSurfer,\n"
    "  go to: http://surfer.nmr.mgh.harvard.edu\n"
    "--------------------------------------------------------------------------\n";

static const char *licmsg =
    "--------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer license file %s not found.\n"
    "  If you are outside the NMR-Martinos Center,\n"
    "  go to http://surfer.nmr.mgh.harvard.edu/registration.html to \n"
    "  get a valid license file (it's free).\n"
    "  If you are inside the NMR-Martinos Center,\n"
    "  make sure to source the standard environment.\n"
    "  A path to an alternative license file can also be\n"
    "  specified with the FS_LICENSE environmental variable.\n"
    "--------------------------------------------------------------------------\n";

static const char *licmsg2 =
    "--------------------------------------------------------------------------\n"
    "ERROR: Invalid FreeSurfer license key found in license file %s\n"
    "  If you are outside the NMR-Martinos Center,\n"
    "  go to http://surfer.nmr.mgh.harvard.edu/registration.html to \n"
    "  get a valid license file (it's free).\n"
    "  If you are inside the NMR-Martinos Center,\n"
    "  make sure to source the standard environment.\n"
    "--------------------------------------------------------------------------\n";

static const char *permission_msg =
    "---------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer license file %s exists "
    "but you do not have read permission.\n"
    "Try running:\n\n"
    "  chmod a+r %s\n"
    "---------------------------------------------------------------------------\n";

static const char *isdir_msg =
    "---------------------------------------------------------------------------\n"
    "ERROR: FS_LICENSE environment variable points to a folder not a file\n"
    "---------------------------------------------------------------------------\n";


void chklc(void)
{
  char dirname[STRLEN], *cp, *alt;
  FILE *lfile = NULL;
  char *email;
  char *magic;
  char *key;
  char *key2;
  char *gkey;
  char *lfilename;
  char str[STRLEN];
  char *crypt_gkey;
  static int first_time = 1;

  sprintf(str, "S%sER%sRONT%sOR", "URF", "_F", "DO");
  if (getenv(str) != NULL) return;

  cp = getenv("FREESURFER_HOME");
  if (cp == NULL) {
    fprintf(stderr, "%s", errmsg);
#ifdef Darwin
    fprintf(stderr, "\n");
    fprintf(stderr, "%s", "Attempting to use /Applications/freesurfer directory.\n");
    strncpy(dirname, "/Applications/freesurfer", STRLEN);
#else
    exit(-1);
#endif
  }
  else {
    strncpy(dirname, cp, STRLEN-1);
  }

  lfilename = (char *)calloc(1, 512);
  email = (char *)calloc(1, 512);
  magic = (char *)calloc(1, 512);
  key = (char *)calloc(1, 512);
  key2 = (char *)calloc(1, 512);
  gkey = (char *)calloc(1, 1024);

  // check if alternative license path is provided:
  alt = getenv("FS_LICENSE");
  if (alt != NULL) {
    strncpy(lfilename, alt, 511);	// leave a nul on the end
    if (Gdiag_no > 0 && first_time) printf("Trying license file %s\n", lfilename);
    lfile = fopen(lfilename, "r");
    if (lfile == NULL) {
      if (errno == EACCES) {
        printf(permission_msg, lfilename, lfilename);
        exit(-1);
      }
      fprintf(stderr, licmsg, lfilename);
      exit(-1);
    }
    // make sure that the path is not a directory
    struct stat path_stat;
    stat(alt, &path_stat);
    if S_ISDIR(path_stat.st_mode) {;
      puts(isdir_msg);
      exit(-1);
    }
  }

  // check for license in FREESURFER_HOME:
  if (lfile == NULL) {
    auto cx = snprintf(lfilename, 511, "%s/lic%s", dirname, "ense.txt");
    if( (cx<0) || (cx>511) ) {
      std::cerr << __FUNCTION__
		<< ": snprintf returned error on line "
		<< __LINE__ << std::endl;
    }
    if (Gdiag_no > 0 && first_time) printf("Trying license file %s\n", lfilename);
    lfile = fopen(lfilename, "r");
  }
  if (lfile == NULL) {
    if (errno == EACCES) {
      printf(permission_msg, lfilename, lfilename);
      exit(-1);
    }
    auto cx = snprintf(lfilename, 511, "%s/.lic%s", dirname, "ense");
    if( (cx<0) || (cx>511) ) {
      std::cerr << __FUNCTION__
		<< ": snprintf returned error on line "
		<< __LINE__ << std::endl;
    }
    if (Gdiag_no > 0 && first_time) printf("Now trying license file %s\n", lfilename);
    lfile = fopen(lfilename, "r");
  }
  if (lfile == NULL) {
    if (errno == EACCES) {
      printf(permission_msg, lfilename, lfilename);
      exit(-1);
    }
    fprintf(stderr, licmsg, lfilename);
    exit(-1);
  }

  if (fscanf(lfile, "%s\n%s\n%s\n%s\n", email, magic, key, key2) != 4) {
    fprintf(stderr, "error parsing license file, expected 4 values\n");
  }

  sprintf(gkey, "%s.%s", email, magic);

  if (Gdiag_no > 0 && first_time) {
    printf("email %s\n", email);
    printf("magic %s\n", magic);
    printf("key   %s\n", key);
    printf("key2  %s\n", key2);
    printf("gkey  %s\n", gkey);
  }

  // This code is meant to provide backwards compatibility
  // of freesurfer license checking. Unfortunately previous
  // freesurfer license keys used an improper salt value of
  // "*C" which is not valid because it should be purely
  // alpha-numeric. New license files have a 4th line with
  // a key generated using a proper salt.

  if (strcmp(key2, "") != 0) {
    // We have a 4 line license file.
    if (Gdiag_no > 0 && first_time) printf("4 line license file\n");
    strcpy(key, key2);
    crypt_gkey = crypt(gkey, "FS");
    if (crypt_gkey == NULL) {
      printf("ERROR: crypt() returned null with 4-line file\n");
      exit(1);
    }
  }
  else {
    // We have a 3 line license file.
    if (Gdiag_no > 0 && first_time) printf("3 line license file\n");
#ifdef Darwin
    // On Darwin systems the key produced with a salt of '*C'
    // is different than that produced on Linux. So to be backwards
    // compatible we must avoid the call to crypt.
    crypt_gkey = key;
#else
    cmp_glib_version();
    crypt_gkey = crypt(gkey, "*C");
#endif
  }

  if (Gdiag_no > 0 && first_time) printf("crypt_gkey %s\n", crypt_gkey);

  if (strcmp(key, crypt_gkey) != 0) {
    fprintf(stderr, licmsg2, lfilename);
    exit(-1);
  }

  free(email);
  free(magic);
  free(key);
  free(key2);
  free(gkey);
  free(lfilename);
  fclose(lfile);

  if (Gdiag_no > 0 && first_time) printf("chklc() done\n");
  first_time = 0;
  return;
}

//  Unfortunately we need separate, but nearly identicle, license checking code
//  for freeview. This is because the license checking above will exit(-1)
//  if the license check fails. But freeview is a clickable application on mac
//  which might not have a terminal window open, so we want a message to appear
//  if the license check fails. The only way to do this is to either modify
//  all existing calls to chklc throughout the code, or have separate license
//  checking code for freeview. The first scares me, so Im going with the
//  latter.
//
//  Also there are difference with the way the FREESURFER_HOME environment variable
//  is handled. It doesnt need to be defined for freeview to operate. So perhaps
//  this somewhat redundant code is a tad more justified.
//
//  return value:   0 - failed
//                  1 - passed
//
//  if failed, error msg will be returned in msg. make sure msg is pre-allocated with enough space
int chklc2(char *msg)
{
  char dirname[STRLEN], *cp, *alt;
  FILE *lfile = NULL;
  char *email;
  char *magic;
  char *key;
  char *key2;
  char *gkey;
  char *lfilename;
  char str[STRLEN];
  char *crypt_gkey;

  sprintf(str, "S%sER%sRONT%sOR", "URF", "_F", "DO");
  if (getenv(str) != NULL) return 1;

  cp = getenv("FREESURFER_HOME");
  if (cp == NULL) {
    fprintf(stderr, "%s", errmsg);
#ifdef Darwin
    fprintf(stderr, "\n");
    fprintf(stderr, "%s", "Attempting to use the /Applications/freesurfer directory.\n");
    strncpy(dirname, "/Applications/freesurfer", STRLEN);
#else
    exit(-1);
#endif
  }
  else {
    strncpy(dirname, cp, STRLEN-1);
  }

  lfilename = (char *)calloc(1, 512);
  email = (char *)calloc(1, 512);
  magic = (char *)calloc(1, 512);
  key = (char *)calloc(1, 512);
  key2 = (char *)calloc(1, 512);
  gkey = (char *)calloc(1, 1024);

  // check if alternative license path is provided:
  alt = getenv("FS_LICENSE");
  if (alt != NULL) {
    strncpy(lfilename, alt, 511);	// leave a nul on the end
    lfile = fopen(lfilename, "r");
    if (lfile == NULL) {
      if (errno == EACCES) {
        fprintf(stderr, permission_msg, lfilename, lfilename);
        if (msg) sprintf(msg, permission_msg, lfilename, lfilename);
        return 0;
      }
      fprintf(stderr, licmsg, lfilename);
      if (msg) sprintf(msg, licmsg, lfilename);
      return 0;
    }
  }

  // check for license in FREESURFER_HOME:
  if (lfile == NULL) {
    auto cx = snprintf(lfilename, 511, "%s/.lic%s", dirname, "ense");
    if( (cx<0) || (cx>511) ) {
      std::cerr << __FUNCTION__
		<< ": snprintf returned error on line "
		<< __LINE__ << std::endl;
    }
    lfile = fopen(lfilename, "r");
  }
  if (lfile == NULL) {
    if (errno == EACCES) {
      fprintf(stderr, permission_msg, lfilename, lfilename);
      if (msg) sprintf(msg, permission_msg, lfilename, lfilename);
      return 0;
    }
    auto cx = snprintf(lfilename, 511, "%s/lic%s", dirname, "ense.txt");
    if( (cx<0) || (cx>511) ) {
      std::cerr << __FUNCTION__
		<< ": snprintf returned error on line "
		<< __LINE__ << std::endl;
    }
    lfile = fopen(lfilename, "r");
  }
  if (lfile == NULL) {
    if (errno == EACCES) {
      fprintf(stderr, permission_msg, lfilename, lfilename);
      if (msg) sprintf(msg, permission_msg, lfilename, lfilename);
      return 0;
    }
    fprintf(stderr, licmsg, lfilename);
    if (msg) sprintf(msg, licmsg, lfilename);
    return 0;
  }

  if (fscanf(lfile, "%s\n%s\n%s\n%s\n", email, magic, key, key2) != 4) {
    fprintf(stderr, "error parsing license file, expected 4 values\n");
    return 0;
  }

  sprintf(gkey, "%s.%s", email, magic);

  // This code is meant to provide backwards compatibility
  // of freesurfer license checking. Unfortunately previous
  // freesurfer license keys used an improper salt value of
  // "*C" which is not valid because it should be purely
  // alpha-numeric. New license files have a 4th line with
  // a key generated using a proper salt.

  if (strcmp(key2, "") != 0) {
    // We have a 4 line license file.
    strcpy(key, key2);
    crypt_gkey = crypt(gkey, "FS");
  }
  else {
// We have a 3 line license file.
#ifdef Darwin
    // On Darwin systems the key produced with a salt of '*C'
    // is different than that produced on Linux. So to be backwards
    // compatible we must avoid the call to crypt.
    crypt_gkey = key;
#else
    cmp_glib_version();
    crypt_gkey = crypt(gkey, "*C");
#endif
  }

  if (strcmp(key, crypt_gkey) != 0) {
    fprintf(stderr, licmsg2, lfilename);
    if (msg) sprintf(msg, licmsg2, lfilename);
    return 0;
  }

  free(email);
  free(magic);
  free(key);
  free(key2);
  free(gkey);
  free(lfilename);
  fclose(lfile);

  return 1;
}

#ifndef Darwin
void cmp_glib_version(void)
{
  int i;
  const char *GNU_LIBC_VERSION_MAX = "2.15";
  int glibc_max[2], glibc_current[2];
  static const char *new_license_msg =
      "--------------------------------------------------------------------------\n"
      "GNU libc version: %s\n"
      "ERROR: Systems running GNU glibc version greater than %s\n"
      "  require a newly formatted license file (it's free). Please\n"
      "  download a new one from the following page:\n"
      "  http://surfer.nmr.mgh.harvard.edu/registration.html\n"
      "--------------------------------------------------------------------------\n";

  sscanf(GNU_LIBC_VERSION_MAX, "%d.%d", &glibc_max[0], &glibc_max[1]);
  sscanf(gnu_get_libc_version(), "%d.%d", &glibc_current[0], &glibc_current[1]);

  for (i = 0; i < 2; i++) {
    if (glibc_current[i] > glibc_max[i]) {
      printf(new_license_msg, gnu_get_libc_version(), GNU_LIBC_VERSION_MAX);
      exit(-1);
    }
  }
}
#endif
