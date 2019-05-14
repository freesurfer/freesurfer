/**
 * @file  tridec.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.7 $
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "version.h"

const char* Progname = "tridec";

int
main(int argc,char *argv[]) {

  FILE *fp1,*fp2,*fp3;
  char fname[STRLEN],fpref[STRLEN],pname[STRLEN];
  int vnum1, vnum2, cnum2;
  char line[STRLEN];
  char *data_dir,*mri_dir;
  int i;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: tridec.c,v 1.7 2011/03/02 00:04:41 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (argc!=5) {
    printf("Usage: %s subject_name fine_file ico_file out_file\n",argv[0]);
    exit(10);
  }

  data_dir = getenv("SUBJECTS_DIR");
  if (data_dir==NULL) {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(0);
  }
  mri_dir = getenv("FREESURFER_HOME");
  if (mri_dir==NULL) {
    printf("environment variable FREESURFER_HOME undefined (use setenv)\n");
    exit(0);
  }

  sprintf(pname,"%s",argv[1]);
  sprintf(fpref,"%s/%s",data_dir,pname);

  sprintf(fname,"%s/bem/%s",fpref,argv[2]);
  fp1 = fopen(fname,"r");
  if (fp1==NULL) {
    printf("File %s not found.\n",fname);
    exit(0);
  }
  sprintf (fname,"%s/lib/bem/%s",mri_dir,argv[3]);
  fp2 = fopen(fname,"r");
  if (fp2==NULL) {
    printf("File %s not found.\n",fname);
    exit(0);
  }
  sprintf(fname,"%s/bem/%s",fpref,argv[4]);
  fp3 = fopen(fname,"w");
  if (fp3==NULL) {
    printf("can't write to file %s\n",fname);
    exit(0);
  }

  fgets(line,STRLEN,fp2);
  sscanf(line,"%d",&vnum2);
  fgets(line,STRLEN,fp1);
  sscanf(line,"%d",&vnum1);
  fprintf(fp3,"%d\n",vnum2);

  for (i=0;i<vnum2;i++) {
    fgets(line,STRLEN,fp1);
    fprintf(fp3,"%s",line);
  }

  for (i=0;i<vnum2;i++) {
    fgets(line,STRLEN,fp2);
  }

  fgets(line,STRLEN,fp2);
  sscanf(line,"%d",&cnum2);
  fprintf(fp3,"%d\n",cnum2);

  for (i=0;i<cnum2;i++) {
    fgets(line,STRLEN,fp2);
    fprintf(fp3,"%s",line);
  }

  printf("written triangle file %s\n",fname);

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  return(0) ;
}
