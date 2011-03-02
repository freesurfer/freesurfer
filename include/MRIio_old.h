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

