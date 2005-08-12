#include <stdio.h>

#include "fio.h"
#include "error.h"
#include "tags.h"

int
TAGreadStart(FILE *fp, double *plen)
{
	int  tag ;
	
	tag = freadInt(fp) ;
	*plen = freadDouble(fp) ;
	
	return(tag) ;
}

int
TAGwriteStart(FILE *fp, int tag, long *phere)
{
	fwriteInt(tag, fp) ;
	*phere = ftell(fp) ;
	fwriteDouble(0.0, fp) ;
	
	return(NO_ERROR) ;
}
int
TAGwriteEnd(FILE *fp, long there)
{
	long here ;

	here = ftell(fp) ;
	fseek(fp, there, SEEK_SET) ;
	fwriteDouble((double)(here-(there+sizeof(double))), fp) ;
	fseek(fp, here, SEEK_SET) ;

	return(NO_ERROR) ;
}

