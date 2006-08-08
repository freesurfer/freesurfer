// $Id: fsenv.h,v 1.2 2006/08/08 20:12:41 greve Exp $
// Load, set freesurfer environment variables

#ifndef FSENV_INC
#define FSENV_INC

#include <sys/utsname.h>
#include "colortab.h"

typedef struct {
  char *FREESURFER_HOME;
  char *SUBJECTS_DIR;
  char *user;      // current user
  char *date;      // current date and time
  char *cwd;       // current working directory
  char *hostname;  // eg, icebox (same as nodename)
  char *sysname;   // eg, Linux
  char *machine;   // eg, i686
  COLOR_TABLE *ctab; // FREESURFER_HOME/FreeSurferColorLUT.txt
} FSENV;

const char *FSENVsrcVersion(void);
FSENV *FSENVgetenv(void);
int FSENVprintenv(FILE *fp,FSENV *env);
int FSENVsetSUBJECTS_DIR(char *SUBJECTS_DIR);
int FSENVfree(FSENV **ppenv);

#endif //#ifndef FSENV_INC
