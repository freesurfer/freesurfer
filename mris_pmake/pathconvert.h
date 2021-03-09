/**
 * @brief convert directory path specifications in *nix.
 *
 * 'pathconvert.h' declares two functions, 'abs2rel' and 'rel2abs'
 * that convert directory path specifications in *nix.
 *
 * These functions were originally defined by Shigio Yamaguchi are
 * copyright 1999 Tama Communications Corporation and are released
 * under the LGPL.
 */
/*
 * Original Authors: Rudolph Pienaar / Christian Haselgrove
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


#ifndef __PATHCONVERT_H__
#define __PATHCONVERT_H__

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>

char    *abs2rel __P((const char *, const char *, char *, size_t));
char    *rel2abs __P((const char *, const char *, char *, size_t));

#endif // __PATHCONVERT_H__
