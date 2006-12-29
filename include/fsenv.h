/**
 * @file  fsenv.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
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


// $Id: fsenv.h,v 1.3 2006/12/29 02:08:59 nicks Exp $
// Load, set freesurfer environment variables

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
int FSENVfree(FSENV **ppenv);

#endif //#ifndef FSENV_INC
