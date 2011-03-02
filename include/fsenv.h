/**
 * @file  fsenv.h
 * @brief Load, set freesurfer environment variables
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.5 $
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

#ifndef FSENV_INC
#define FSENV_INC

#include <sys/utsname.h>
#include "colortab.h"

typedef struct
{
  char *FREESURFER_HOME;
  char *SUBJECTS_DIR;
  char *user;      // current user
  char *date;      // current date and time
  char *cwd;       // current working directory
  char *hostname;  // eg, icebox (same as nodename)
  char *sysname;   // eg, Linux
  char *machine;   // eg, i686
  COLOR_TABLE *ctab; // FREESURFER_HOME/FreeSurferColorLUT.txt
}
FSENV;

const char *FSENVsrcVersion(void);
FSENV *FSENVgetenv(void);
int FSENVprintenv(FILE *fp,FSENV *env);
int FSENVsetSUBJECTS_DIR(char *SUBJECTS_DIR);
char *FSENVgetSUBJECTS_DIR(void);
int FSENVfree(FSENV **ppenv);

#endif //#ifndef FSENV_INC
