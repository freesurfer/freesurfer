/**
 * @file  MRIio_old.h
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


/* header file for MRI I/O */
/* 2/1/91 - AD */

#ifndef MRIIO_OLD_H
#define MRIIO_OLD_H

char *lmalloc(unsigned long size) ;
char *lcalloc(size_t nmemb,size_t size) ;
void buffer_to_image(unsigned char *buf,unsigned char**im,int ysize,int xsize);
void image_to_buffer(unsigned char **im,unsigned char*buf,int ysize,int xsize);
void file_name(char *fpref, char *fname, int num, char *form) ;


#endif

