// $Id: fsenv.c,v 1.1 2006/08/08 20:04:17 greve Exp $


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
const char *FSENVsrcVersion(void) { 
  return("$Id: fsenv.c,v 1.1 2006/08/08 20:04:17 greve Exp $");
}

FSENV *FSENVgetenv(void)
{
  FSENV *fsenv;
  char *pc, tmpstr[2000];
  struct utsname uts;

  fsenv = (FSENV *) calloc(sizeof(FSENV),1); 

  pc = getenv("FREESURFER_HOME");
  if(pc == NULL){
    printf("FREESURFER_HOME not defined\n");
    return(NULL);
  }
  fsenv->FREESURFER_HOME = strcpyalloc(pc);

  pc = getenv("SUBJECTS_DIR");
  if(pc == NULL){
    printf("SUBJECTS_DIR not defined\n");
    return(NULL);
  }
  fsenv->SUBJECTS_DIR = strcpyalloc(pc);
  fsenv->user = strcpyalloc(VERuser());

  getcwd(tmpstr,2000);
  fsenv->cwd = strcpyalloc(tmpstr);

  uname(&uts);
  fsenv->sysname  = strcpyalloc(uts.sysname);
  fsenv->machine  = strcpyalloc(uts.machine);
  fsenv->hostname = strcpyalloc(uts.nodename);

  sprintf(tmpstr,"%s/FreeSurferColorLUT.txt",fsenv->FREESURFER_HOME);
  fsenv->ctab = CTABreadASCII(tmpstr);
  if(fsenv->ctab == NULL){
    printf("ERROR: reading %s\n",tmpstr);
    return(NULL);
  }

  fsenv->date = VERcurTimeStamp();

  return(fsenv);
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
int FSENVsetSUBJECTS_DIR(char *SUBJECTS_DIR)
{
  setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
  return(0);
}






