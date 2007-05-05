/**
 * @file  librms.h
 *
 */
/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/05/05 00:00:07 $
 *    $Revision: 1.2 $
 *
 * Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
 * Washington University, Mallinckrodt Institute of Radiology.
 * All Rights Reserved.
 *
 * This software may not be reproduced, copied, or distributed without 
 * written permission of Washington University. For further information 
 * contact A. Z. Snyder.
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/************/
/* imgpad.f */
/************/
int	npad_	(int *n, int *m);
void	imgpad_ (float *imag, int *nx, int *ny, int *nz, float *imagp, int *nxp, int *nyp, int *nzp);
void	imgdap_ (float *imag, int *nx, int *ny, int *nz, float *imagp, int *nxp, int *nyp, int *nzp);

/*************/
/* gauss3d.f */
/*************/
void	gauss3d_ (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);

/***************/
/* img2lmask.c */
/***************/
void	img2lmask (int *pnx, int *pny, int *pnz, float *imag, int *mask, float *mmppix, float *pfhalf, float *pcrit);

/***************/
/* param6opr.f */
/***************/
void	param2warp_ (int *mode, float *param, float *a);
void	warp2param_ (int *mode, float *a, float *param);
void	img2vrt_ (float *mmppix, float *center, float *vox2ras);
void	vrt2img_ (float *mmppix, float *center, float *ras2vox);

/************/
/* matopr.f */
/************/
void	matmul_	(float *a, float *b, float *c, int *n);

/************************/
/* fftsol.f or fftsun.f */
/************************/
void fft_   (float *a, float *b, int *nseg, int *n, int *nspn, int *isn);
void realt_ (float *a, float *b, int *nseg, int *n, int *nspn, int *isn);
