/*$Header: /space/repo/1/dev/dev/talairach_avi/librms.h,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: librms.h,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/**/

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
