#include "tukey.h"

double
tukey_biweight(double residual, double C)
{
	double p ;

	if (abs(residual) > C)
		return(C*C/2) ;
	else
	{
		p = residual/C ; p *= p ;
		p = 1-p ; p = p*p*p ;
		return(((C*C)/2)*(1-p)) ;
	}
}
