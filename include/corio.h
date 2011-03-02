/**
 * @file  corio.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.3 $
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




#ifndef CORIO_H_INC
#define CORIO_H_INC

#define CORVAL(ppCOR,row,col,slc) *(ppCOR[slc]+col+row*256)

int               free_cor(unsigned char ***pppCOR);
unsigned char ** alloc_cor(void);
unsigned char **    ld_cor(char *cordir);
int                 sv_cor(unsigned char **COR, char *cordir);
int setcorval(unsigned char val, unsigned char ** COR,
              int row, int col, int slc);
unsigned char getcorval(unsigned char ** COR,
                        int row, int col, int slc);
int cordir_iswritable(char *cordir);


#endif
