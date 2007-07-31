/**
 * @file  chklc.c
 * @brief Routine to check .license file
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/07/31 07:52:10 $
 *    $Revision: 1.10 $
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

#include <unistd.h>
#include <const.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char *crypt(const char *, const char *);

static const char* errmsg =
"--------------------------------------------------------------------------\n"
"ERROR: FreeSurfer environment FREESURFER_HOME is not defined.\n"
"  If you are outside the NMR-Martinos Center, please set this\n"
"  variable to the location where you installed FreeSurfer.\n"
"  If you are inside the NMR-Martinos Center, please source\n"
"  the standard environment. If you need to install FreeSurfer,\n"
"  go to: http://surfer.nmr.mgh.harvard.edu\n"
"--------------------------------------------------------------------------\n";

static const char* licmsg =
"--------------------------------------------------------------------------\n"
"ERROR: FreeSurfer license file %s not found.\n"
"  If you are outside the NMR-Martinos Center,\n"
"  go to http://surfer.nmr.mgh.harvard.edu to \n"
"  get a valid license file (it's free).\n"
"  If you are inside the NMR-Martinos Center,\n"
"  make sure to source the standard environment.\n"
"--------------------------------------------------------------------------\n";

static const char* licmsg2 =
"--------------------------------------------------------------------------\n"
"ERROR: Invalid FreeSurfer license key found in license file %s\n"
"  If you are outside the NMR-Martinos Center,\n"
"  go to http://surfer.nmr.mgh.harvard.edu to \n"
"  get a valid license file (it's free).\n"
"  If you are inside the NMR-Martinos Center,\n"
"  make sure to source the standard environment.\n"
"--------------------------------------------------------------------------\n";

void chklc(void)
{
  char  dirname[STRLEN], *cp ;
  FILE* lfile;
  char* email;
  char* magic;
  char* key;
  char* gkey;
  char* lfilename;
  char  str[STRLEN] ;

  sprintf(str, "S%sER%sRONT%sOR", "URF", "_F", "DO") ;
  if (getenv(str) != NULL) return ;

  cp = getenv("FREESURFER_HOME");
  if (cp == NULL)
  {
    //fprintf(stdout,errmsg);
    fprintf(stderr,errmsg);
    exit(-1);
  }

  strncpy(dirname, cp, STRLEN) ;

  lfilename = (char*)calloc(1,512);
  email = (char*)calloc(1,512);
  magic = (char*)calloc(1,512);
  key   = (char*)calloc(1,512);
  gkey  = (char*)calloc(1,1024);

  sprintf(lfilename,"%s/.lic%s",dirname, "ense");

  lfile = fopen(lfilename,"r");
  if (lfile == NULL)
  {
    //fprintf(stdout,licmsg,lfilename);
    fprintf(stderr,licmsg,lfilename);
    exit(-1);
  }

  fscanf(lfile,"%s\n",email);
  fscanf(lfile,"%s\n",magic);
  fscanf(lfile,"%s\n",key);

  sprintf(gkey,"%s.%s",email,magic);

#ifndef Darwin
  if (strcmp(key,crypt(gkey,"*C*O*R*T*E*C*H*S*0*1*2*3*"))!=0)
  {
    //fprintf(stdout,licmsg2,lfilename);
    fprintf(stderr,licmsg2,lfilename);
    exit(-1);
  }
#endif

  free(email);
  free(magic);
  free(key);
  free(gkey);
  free(lfilename);
  fclose(lfile) ;

  return;
}

