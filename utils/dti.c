// $Id: dti.c,v 1.1 2006/09/08 22:21:26 greve Exp $

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
#include "fio.h"
#include "mri.h"
#include "DICOMRead.h"

/* --------------------------------------------- */
// Return the CVS version of this file.
const char *DTIsrcVersion(void) { 
  return("$Id: dti.c,v 1.1 2006/09/08 22:21:26 greve Exp $");
}


/*-----------------------------------------------------------------
  DTIparamsFromSiemensAscii() - reads in diffusion parameters from the
  Siemens ASCII header stored in fname. fname may be a siemens dicom
  file or an infodump file as produced by mri_probedicom run on a
  siemens dicom file. 
  -----------------------------------------------------------------*/
int DTIparamsFromSiemensAscii(char *fname, float *bValue, 
			      int *nAcq, int *nDir, int *DiffMode)
{
  char *tag, *pc;

  if(!fio_FileExistsReadable(fname)){
    printf("ERROR: cannot read %s\n",fname);
    return(1);
  }

  tag = "sDiffusion.alBValue[1]";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(pc, "%f", bValue);
  free(pc);

  tag = "sWiPMemBlock.alFree[8]";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(pc, "%d", nAcq);
  free(pc);

  tag = "sDiffusion.lDiffDirections";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(fname, "%d", nDir);
  free(pc);

  tag = "sWiPMemBlock.alFree[1]";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(pc, "%d", DiffMode);
  free(pc);

  return(0);
}
