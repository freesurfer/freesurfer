/**
 * @file  fsenv.c
 * @brief Utils to get, set, and print freesurfer environment vars
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:44 $
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pwd.h>
#include <time.h>
#include "fsenv.h"
#include "utils.h"
#include "version.h"

/* --------------------------------------------- */
// Return the CVS version of this file.
const char *FSENVsrcVersion(void)
{
  return("$Id: fsenv.c,v 1.5 2011/03/02 00:04:44 nicks Exp $");
}

FSENV *FSENVgetenv(void)
{
  FSENV *fsenv;
  char *pc, tmpstr[2000];
  struct utsname uts;

  fsenv = (FSENV *) calloc(sizeof(FSENV),1);

  pc = getenv("FREESURFER_HOME");
  if (pc == NULL)
  {
    printf("FREESURFER_HOME not defined\n");
    return(NULL);
  }
  fsenv->FREESURFER_HOME = strcpyalloc(pc);

  pc = getenv("SUBJECTS_DIR");
  if (pc == NULL)
  {
    printf("SUBJECTS_DIR not defined\n");
    return(NULL);
  }
  fsenv->SUBJECTS_DIR = strcpyalloc(pc);
  fsenv->user = strcpyalloc(VERuser());

  // Current working directory
  getcwd(tmpstr,2000);
  fsenv->cwd = strcpyalloc(tmpstr);

  // Kernel information
  uname(&uts);
  fsenv->sysname  = strcpyalloc(uts.sysname);
  fsenv->machine  = strcpyalloc(uts.machine);
  fsenv->hostname = strcpyalloc(uts.nodename);

  // Load the default color table
  sprintf(tmpstr,"%s/FreeSurferColorLUT.txt",fsenv->FREESURFER_HOME);
  fsenv->ctab = CTABreadASCII(tmpstr);
  if (fsenv->ctab == NULL)
  {
    printf("ERROR: reading %s\n",tmpstr);
    return(NULL);
  }

  // Get time and date at the time this function was called
  fsenv->date = VERcurTimeStamp();

  return(fsenv);
}
/*-----------------------------------------------*/
int FSENVfree(FSENV **ppenv)
{
  FSENV *env = *ppenv;
  free(env->FREESURFER_HOME);
  free(env->SUBJECTS_DIR);
  free(env->user);
  free(env->date);
  free(env->cwd);
  free(env->hostname);
  free(env->sysname);
  free(env->machine);
  CTABfree(&env->ctab);
  free(*ppenv);
  *ppenv = NULL;
  return(0);
}


/*-----------------------------------------------*/
int FSENVprintenv(FILE *fp,FSENV *env)
{
  fprintf(fp,"FREESURFER_HOME %s\n",env->FREESURFER_HOME);
  fprintf(fp,"SUBJECTS_DIR %s\n",env->SUBJECTS_DIR);
  fprintf(fp,"user %s\n",env->user);
  fprintf(fp,"date %s\n",env->date);
  fprintf(fp,"cwd %s\n",env->cwd);
  fprintf(fp,"hostname %s\n",env->hostname);
  fprintf(fp,"sysname %s\n",env->sysname);
  fprintf(fp,"machine %s\n",env->machine);
  return(0);
}


/*-----------------------------------------------*/
char *FSENVgetSUBJECTS_DIR(void)
{
  char *pc = getenv("SUBJECTS_DIR");
  if (pc == NULL)
  {
    printf("FSENVgetSUBJECTS_DIR: SUBJECTS_DIR not defined\n");
    return(NULL);
  }
  return strcpyalloc(pc);
}


/*-----------------------------------------------*/
int FSENVsetSUBJECTS_DIR(char *SUBJECTS_DIR)
{
  setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
  return(0);
}






