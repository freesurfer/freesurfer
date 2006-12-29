/**
 * @file  chklc.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:30 $
 *    $Revision: 1.9 $
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

extern char     *crypt(const char *, const char *);

void chklc(void)
{
  /*#ifndef Darwin*/
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
    printf("-----------------------------------------------------\n");
    printf("FreeSurfer Error: environment FREESURFER_HOME is not defined.\n"
           "  If you are outside the NMR-Martinos Center, please set this\n"
           "  variable to the location where you installed FreeSurfer.\n"
           "  If you are inside the NMR-Martinos Center, please source\n"
           "  the standard environment.\n");
    printf(
      "  If you need to install FreeSurfer, go to surfer.nmr.mgh.harvard.edu\n");
    fprintf(stderr,"-----------------------------------------------------\n");
    fprintf(stderr,
            "FreeSurfer Error: environment FREESURFER_HOME is not defined.\n"
            "  If you are outside the NMR-Martinos Center, please set this\n"
            "  variable to the location where you installed FreeSurfer.\n"
            "  If you are inside the NMR-Martinos Center, please source\n"
            "  the standard environment.\n");
    fprintf(stderr,
            "  If you need to install FreeSurfer, go to surfer.nmr.mgh.harvard.edu\n");
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
    printf("-----------------------------------------------------\n");
    printf("FreeSurfer license file %s not found.\n",lfilename);
    printf("  If you are outside the NMR-Martinos center,\n"
           "    go to http://surfer.nmr.mgh.harvard.edu to \n"
           "    get a valid license file (it's free).\n");
    printf("  If you are inside the NMR-Martinos center,\n"
           "    make sure to source the standard environment.\n");
    fprintf(stderr,"-----------------------------------------------------\n");
    fprintf(stderr,"FreeSurfer license file %s not found.\n",lfilename);
    fprintf(stderr,
            "  If you are outside the NMR-Martinos center,\n"
            "    go to http://surfer.nmr.mgh.harvard.edu to \n"
            "    get a valid license file (it's free).\n");
    fprintf(stderr,
            "  If you are inside the NMR-Martinos center,\n"
            "    make sure to source the standard environment.\n");
    exit(-1);
  }

  fscanf(lfile,"%s\n",email);
  fscanf(lfile,"%s\n",magic);
  fscanf(lfile,"%s\n",key);

  sprintf(gkey,"%s.%s",email,magic);

#ifndef Darwin
  if (strcmp(key,crypt(gkey,"*C*O*R*T*E*C*H*S*0*1*2*3*"))!=0)
  {
    printf("-----------------------------------------------------\n");
    printf("FreeSurfer: Invalid license key found in license file %s \n",
           lfilename);
    printf("  If you are outside the NMR-Martinos center,\n"
           "  go to http://surfer.nmr.mgh.harvard.edu to \n"
           "  to get a valid license file (it's free).\n");
    printf("  If you are inside the NMR-Martinos center,\n"
           "  make sure to source the standard environment.\n");
    fprintf(stderr,"-----------------------------------------------------\n");
    fprintf(stderr,
            "FreeSurfer: Invalid license key found in license file %s \n",
            lfilename);
    fprintf(stderr,
            "  If you are outside the NMR-Martinos center,\n"
            "  go to http://surfer.nmr.mgh.harvard.edu to \n"
            "  to get a valid license file (it's free).\n");
    fprintf(stderr,
            "  If you are inside the NMR-Martinos center,\n"
            "  make sure to source the standard environment.\n");
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
  /*#endif*/
}

