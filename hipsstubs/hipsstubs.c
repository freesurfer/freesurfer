/**
 * @file  hipsstubs.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.11 $
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

#include "error.h"

#if 0
/* wrote replacements for these (in hipsrepl.c) */
int init_header(void) ;
int init_header(void) {
  return(-1);
}
int h_copy(void) ;
int h_copy(void) {
  return(-1);
}
int free_header(void) ;
int free_header(void) {
  return(-1);
}
#endif

#if 1
int perr(void) ;
int perr(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_canny(void) ;
int h_canny(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_enlarge(void) ;
int h_enlarge(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_reduce(void) ;
int h_reduce(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_invert(void) ;
int h_invert(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
#include <stdlib.h>
int h_tob(void) ;
int h_tob(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_minmax(void) ;
int h_minmax(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_linscale(void) ;
int h_linscale(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int update_header(void) ;
int update_header(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int free_hdrcon(void) ;
int free_hdrcon(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int *memalloc(size_t, size_t)  ;
int *memalloc(size_t a, size_t b) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(calloc(a,b));
}
int *halloc(size_t a, size_t b) ;
int *halloc(size_t a, size_t b) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(calloc(a,b)) ;
}
int fwrite_header(void) ;
int fwrite_header(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int fwrite_image(void) ;
int fwrite_image(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int fread_header(void) ;
int fread_header(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int fread_image(void) ;
int fread_image(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_tof(void) ;
int h_tof(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_toc(void) ;
int h_toc(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_fourtr(void) ;
int h_fourtr(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_invfourtr(void) ;
int h_invfourtr(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_divscale(void) ;
int h_divscale(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_mul(void) ;
int h_mul(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_flipquad(void) ;
int h_flipquad(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_softthresh(void) ;
int h_softthresh(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_tod(void) ;
int h_tod(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_todc(void) ;
int h_todc(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_toi(void) ;
int h_toi(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_mulscale(void) ;
int h_mulscale(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_div(void) ;
int h_div(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_add(void) ;
int h_add(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_diff(void) ;
int h_diff(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int alloc_histo(void) ;
int alloc_histo(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_clearhisto(void) ;
int h_clearhisto(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_histo(void) ;
int h_histo(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_histoeq(void) ;
int h_histoeq(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_pixmap(void) ;
int h_pixmap(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_entropycnt(void) ;
int h_entropycnt(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int h_entropy(void) ;
int h_entropy(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
#if 0
unsigned char *hmalloc(unsigned long i) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(0) ;
}
int WinShowImage(void) ;
int WinShowImage(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int WinSetName(void) ;
int WinSetName(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int WinCreate(void) ;
int WinCreate(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
#endif
int h_median(void) ;
int h_median(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}
int hformatname(void) ;
int hformatname(void) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(-1);
}

#if 1
int hipserrno = 0 ;

int hips_hclip = 0 ;
int hips_lclip = 0 ;
int hips_oldhdr = 0 ;
int hips_cplxtor = 0 ;
int hips_rtocplx = 0 ;
#endif

void MachThreadInit(void) ;
void MachThreadInit(void) {}
void MachThreadGetTid(void) ;
void MachThreadGetTid(void) {}
void MachThreadKill(void) ;
void MachThreadKill(void) {}
void MachThreadStart(void) ;
void MachThreadStart(void) {}
void MachThreadSuspend(void) ;
void MachThreadSuspend(void) {}
void MachThreadResume(void) ;
void MachThreadResume(void) {}
void MachThreadSleep(void) ;
void MachThreadSleep(void) {}
void MachThreadYield(void) ;
void MachThreadYield(void) {}
void MachThreadExit(void) ;
void MachThreadExit(void) {}


#endif
/* getparam(hd,name,fmt,count,valuepointer) */

int getparam(void *hd, ...) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(0) ;
}
int setparam(void *hd, ...) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(0) ;
}
void *findparam(void *hd, char *name) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(NULL) ;
}
int clearparam(void *hd, char *name) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(0) ;
}


int h_applylut(void *hdi,void *hdo,int count,unsigned char *lut) {
  ErrorExit(ERROR_UNSUPPORTED, "HIPS unsupported") ;
  return(0) ;
}




