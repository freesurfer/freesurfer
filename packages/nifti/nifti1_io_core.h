/** \file nifti1_io_core.h
    \brief Data structures for using nifti1_io API.
           - Written by Bob Cox, SSCC NIMH
           - Revisions by Rick Reynolds, SSCC NIMH
 */
#ifndef _NIFTI_IO_CORE_HEADER_
#define _NIFTI_IO_CORE_HEADER_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "nifti1.h"                  /*** NIFTI-1 header specification ***/

/*=================*/
#ifdef  __cplusplus
extern "C" {
#endif
/*=================*/

/*****===================================================================*****/
/*****         File nifti1_io.h == Declarations for nifti1_io.c          *****/
/*****...................................................................*****/
/*****            This code is released to the public domain.            *****/
/*****...................................................................*****/
/*****  Author: Robert W Cox, SSCC/DIRP/NIMH/NIH/DHHS/USA/EARTH          *****/
/*****  Date:   August 2003                                              *****/
/*****...................................................................*****/
/*****  Neither the National Institutes of Health (NIH), nor any of its  *****/
/*****  employees imply any warranty of usefulness of this software for  *****/
/*****  any purpose, and do not assume any liability for damages,        *****/
/*****  incidental or otherwise, caused by any use of this document.     *****/
/*****===================================================================*****/

/* 
   Modified by: Mark Jenkinson (FMRIB Centre, University of Oxford, UK)
   Date: July/August 2004 

      Mainly adding low-level IO and changing things to allow gzipped files
      to be read and written
      Full backwards compatability should have been maintained

   Modified by: Rick Reynolds (SSCC/DIRP/NIMH, National Institutes of Health)
   Date: December 2004

      Modified and added many routines for I/O.
*/

/********************** Some sample data structures **************************/

typedef struct {                   /** 4x4 matrix struct **/
  float m[4][4] ;
} mat44 ;

typedef struct {                   /** 3x3 matrix struct **/
  float m[3][3] ;
} mat33 ;



mat44 nifti_mat44_inverse( mat44 R ) ;

mat33 nifti_mat33_inverse( mat33 R ) ;
mat33 nifti_mat33_polar  ( mat33 A ) ;
float nifti_mat33_rownorm( mat33 A ) ;
float nifti_mat33_colnorm( mat33 A ) ;
float nifti_mat33_determ ( mat33 R ) ;
mat33 nifti_mat33_mul    ( mat33 A , mat33 B ) ;

void  nifti_swap_2bytes ( size_t n , void *ar ) ;
void  nifti_swap_4bytes ( size_t n , void *ar ) ;
void  nifti_swap_8bytes ( size_t n , void *ar ) ;
void  nifti_swap_16bytes( size_t n , void *ar ) ;
void  nifti_swap_Nbytes ( size_t n , int siz , void *ar ) ;

void  swap_nifti_header ( struct nifti_1_header *h , int is_nifti ) ;


void nifti_mat44_to_quatern( mat44 R ,
                             float *qb, float *qc, float *qd,
                             float *qx, float *qy, float *qz,
                             float *dx, float *dy, float *dz, float *qfac ) ;

mat44 nifti_quatern_to_mat44( float qb, float qc, float qd,
                              float qx, float qy, float qz,
                              float dx, float dy, float dz, float qfac );

mat44 nifti_make_orthog_mat44( float r11, float r12, float r13 ,
                               float r21, float r22, float r23 ,
                               float r31, float r32, float r33  ) ;



/* Orientation codes that might be returned from nifti_mat44_to_orientation().*/

#define NIFTI_L2R  1    /* Left to Right         */
#define NIFTI_R2L  2    /* Right to Left         */
#define NIFTI_P2A  3    /* Posterior to Anterior */
#define NIFTI_A2P  4    /* Anterior to Posterior */
#define NIFTI_I2S  5    /* Inferior to Superior  */
#define NIFTI_S2I  6    /* Superior to Inferior  */

void nifti_mat44_to_orientation( mat44 R , int *icod, int *jcod, int *kcod ) ;

/*=================*/
#ifdef  __cplusplus
}
#endif
/*=================*/

#endif /* _NIFTI_IO_CORE_HEADER_ */
