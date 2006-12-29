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
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.2 $
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
