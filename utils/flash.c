#include <math.h>
#include "macros.h"
#include "flash.h"



double
dFlash_dT1(double T1, double PD, double TR, double flip_angle, double TE)
{
	double e1, numer, denom ;

	e1 = exp(TR/T1) ;
	numer = e1 * PD * TR * (cos(flip_angle)-1)*sin(flip_angle) ;
	denom = T1 * (e1 - cos(flip_angle)) ; denom *= denom ;
	if (DZERO(denom))
		denom = 0.0001 ;
	return(numer/denom) ;
}

double
dFlash_dPD(double T1, double PD, double TR, double flip_angle, double TE)
{
	double e1, numer, denom ;

	e1 = exp(TR/T1) ;
	numer = (e1-1) * sin(flip_angle) ;
	denom = e1 - cos(flip_angle) ;
	if (DZERO(denom))
		denom = 0.0001 ;
	return(numer/denom) ;
}

double
FLASHforwardModel(double T1, double PD, double TR, double flip_angle, double TE)
{
  double FLASH, E1 ;
  double  CFA, SFA ;


  CFA = cos(flip_angle) ; SFA = sin(flip_angle) ;
  E1 = exp(-TR/T1) ;
      
  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}
MRI *
MRIparameterMapsToFlash(MRI *mri_src, MRI *mri_dst, double *TRs, double *TEs, double *FAs, int nflash)
{
	int    x, y, z, n ;
	double T1, PD ;
	Real   val ;

	if (!mri_dst)
		mri_dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth,
																	 mri_src->type, nflash) ;

	for (x = 0 ; x < mri_src->width ; x++)
	{
		for (y = 0 ; y < mri_src->height ; y++)
		{
			for (z = 0 ; z < mri_src->depth ; z++)
			{
				MRIsampleVolumeFrame(mri_src, x, y, z, 0, &val); T1 = val ;
				MRIsampleVolumeFrame(mri_src, x, y, z, 1, &val); PD = val ;
				
				for (n = 0 ; n < nflash ; n++)
				{
					val = FLASHforwardModel(T1, PD, TRs[n], FAs[n], TEs[n]) ;
					MRISseq_vox(mri_dst, x, y, z, n) = val ;
				}
			}
		}
	}


	return(mri_dst) ;
}

