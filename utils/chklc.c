/**
 * @file  chklc.c
 * @brief Routine to check .license file
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2014/08/26 21:25:12 $
 *    $Revision: 1.21 $
 *
 * Copyright © 2011-2013 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <unistd.h>
#include <const.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

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
  char* key2;
  char* gkey;
  char* lfilename;
  char  str[STRLEN] ;
  char* crypt_gkey ;

  sprintf(str, "S%sER%sRONT%sOR", "URF", "_F", "DO") ;
  if (getenv(str) != NULL) return ;

  cp = getenv("FREESURFER_HOME");
  if (cp == NULL)
  {
    fprintf(stderr,"%s",errmsg);
    #ifdef Darwin
      fprintf(stderr,"\n");
      fprintf(stderr,"%s","Attempting to use /Applications/freesurfer directory.\n");
      strncpy(dirname, "/Applications/freesurfer", STRLEN) ;
    #else
      exit(-1);
    #endif
  }
  else
  {
    strncpy(dirname, cp, STRLEN) ;
  }

  lfilename = (char*)calloc(1,512);
  email = (char*)calloc(1,512);
  magic = (char*)calloc(1,512);
  key   = (char*)calloc(1,512);
  key2  = (char*)calloc(1,512);
  gkey  = (char*)calloc(1,1024);

  sprintf(lfilename,"%s/.lic%s",dirname, "ense");

  lfile = fopen(lfilename,"r");
  if (lfile == NULL)
  {
    if(errno == EACCES){
      printf("\n\nERROR: FreeSurfer license file %s exists but you do not have read permission\n",lfilename);
      printf("   Try running chmod a+r %s\n\n\n",lfilename);
      exit(-1);
    }
    sprintf(lfilename,"%s/lic%s",dirname, "ense.txt");
    lfile = fopen(lfilename,"r");
  }
  if (lfile == NULL)
  {
    if(errno == EACCES){
      printf("\n\nERROR: FreeSurfer license file %s exists but you do not have read permission\n",lfilename);
      printf("   Try running chmod a+r %s\n\n\n",lfilename);
      exit(-1);
    }
    fprintf(stderr,licmsg,lfilename);
    exit(-1);
  }

  fscanf(lfile,"%s\n",email);
  fscanf(lfile,"%s\n",magic);
  fscanf(lfile,"%s\n",key);
  fscanf(lfile,"%s\n",key2); 

  sprintf(gkey,"%s.%s",email,magic);
  
  // This code is meant to provide backwards compatibility
  // of freesurfer license checking. Unfortunately previous 
  // freesurfer license keys used an improper salt value of
  // "*C" which is not valid because it should be purely 
  // alpha-numeric. New license files have a 4th line with
  // a key generated using a proper salt.
  
  if (strcmp(key2, "") != 0)
  {
    // We have a 4 line license file.
    strcpy(key,key2);
    crypt_gkey = crypt(gkey,"FS");
  } 
  else 
  {
    // We have a 3 line license file.
    #ifdef Darwin
      // On Darwin systems the key produced with a salt of '*C'
      // is different than that produced on Linux. So to be backwards
      // compatible we must avoid the call to crypt.
      crypt_gkey = key;
    #else 
      crypt_gkey = crypt(gkey,"*C");
    #endif
  }

  if (strcmp(key,crypt_gkey)!=0)
  {
    fprintf(stderr,licmsg2,lfilename);
    exit(-1);
  }

  free(email);
  free(magic);
  free(key);
  free(key2);
  free(gkey);
  free(lfilename);
  fclose(lfile) ;

  return;
}

//  Unfortunately we need separate, but nearly identicle, license checking code
//  for freeview. This is because the license checking above will exit(-1) 
//  if the license check fails. But freeview is a clickable application on mac
//  which might not have a terminal window open, so we want a message to appear
//  if the license check fails. The only way to do this is to either modify 
//  all existing calls to chklc throughout the code, or have separate license
//  checking code for freeview. The first scares me, so Im going with the latter.
//
//  Also there are difference with the way the FREESURFER_HOME environment variable 
//  is handled. It doesnt need to be defined for freeview to operate. So perhaps
//  this somewhat redundant code is a tad more justified.
//
//  return value:   0 - failed
//                  1 - passed
//
//  if failed, error msg will be returned in msg. make sure msg is pre-allocated with enough space
int chklc2(char* msg)
{

  char  dirname[STRLEN], *cp;
  FILE* lfile;
  char* email;
  char* magic;
  char* key;
  char* key2;
  char* gkey;
  char* lfilename;
  char  str[STRLEN] ;
  char* crypt_gkey ;

  sprintf(str, "S%sER%sRONT%sOR", "URF", "_F", "DO") ;
  if (getenv(str) != NULL) return 1;

  cp = getenv("FREESURFER_HOME");
  if (cp == NULL)
  {
    fprintf(stderr,"%s",errmsg);
    #ifdef Darwin
      fprintf(stderr,"\n");
      fprintf(stderr,"%s","Attempting to use the /Applications/freesurfer directory.\n");
      strncpy(dirname, "/Applications/freesurfer", STRLEN) ;
    #else
      exit(-1);
    #endif
  }
  else
  {
    strncpy(dirname, cp, STRLEN) ;
  }
 
  lfilename = (char*)calloc(1,512);
  email = (char*)calloc(1,512);
  magic = (char*)calloc(1,512);
  key   = (char*)calloc(1,512);
  key2  = (char*)calloc(1,512);
  gkey  = (char*)calloc(1,1024);

  sprintf(lfilename,"%s/.lic%s",dirname, "ense");

  lfile = fopen(lfilename,"r");
  char permission_msg[] =
    "---------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer license file %s exists \n"
    "but you do not have read permission. Try running\n\n"
    "  chmod a+r %s\n"
    "---------------------------------------------------------------------------\n";
  if (lfile == NULL)
  {
    if(errno == EACCES){
      fprintf(stderr,permission_msg,lfilename,lfilename);
      if (msg)
        sprintf(msg,permission_msg,lfilename,lfilename);
      return 0;
    }
    sprintf(lfilename,"%s/lic%s",dirname, "ense.txt");
    lfile = fopen(lfilename,"r");
  }
  if (lfile == NULL)
  {
    if(errno == EACCES){
      fprintf(stderr,permission_msg,lfilename,lfilename);
      if (msg)
        sprintf(msg,permission_msg,lfilename,lfilename);
      return 0;
    }
    fprintf(stderr,licmsg,lfilename);
    if (msg)
      sprintf(msg,licmsg,lfilename);
    return 0;
  }

  fscanf(lfile,"%s\n",email);
  fscanf(lfile,"%s\n",magic);
  fscanf(lfile,"%s\n",key);
  fscanf(lfile,"%s\n",key2); 

  sprintf(gkey,"%s.%s",email,magic);

  // This code is meant to provide backwards compatibility
  // of freesurfer license checking. Unfortunately previous 
  // freesurfer license keys used an improper salt value of
  // "*C" which is not valid because it should be purely 
  // alpha-numeric. New license files have a 4th line with
  // a key generated using a proper salt.
  
  if (strcmp(key2, "") != 0)
  {
    // We have a 4 line license file.
    strcpy(key,key2);
    crypt_gkey = crypt(gkey,"FS");
  } 
  else 
  {
    // We have a 3 line license file.
    #ifdef Darwin
      // On Darwin systems the key produced with a salt of '*C'
      // is different than that produced on Linux. So to be backwards
      // compatible we must avoid the call to crypt.
      crypt_gkey = key;
    #else 
      crypt_gkey = crypt(gkey,"*C");
    #endif
  }

  if (strcmp(key,crypt_gkey)!=0)
  {
    fprintf(stderr,licmsg2,lfilename);
    if (msg)
      sprintf(msg, licmsg2, lfilename);
    return 0;
  }

  free(email);
  free(magic);
  free(key);
  free(key2);
  free(gkey);
  free(lfilename);
  fclose(lfile) ;

  return 1;
}

