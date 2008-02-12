/**
 * @file  fsenv.h
 * @brief Load, set freesurfer environment variables
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/02/12 00:41:45 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2006-2008,
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
