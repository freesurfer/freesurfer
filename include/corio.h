/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
