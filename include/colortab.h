#ifndef COLORTAB_H
#define COLORTAB_H

#include "const.h"

typedef struct
{
	int   index ;
	char  name[STRLEN] ;
	int   r ;
	int   g ;
	int   b ;
	int   flag ;
} COLOR_TABLE_ENTRY, CTE ;

typedef struct
{
	CTE   *bins ;
	int   nbins ;
	char  fname[STRLEN] ;  /* file it was read from */
} COLOR_TABLE, CT ;

COLOR_TABLE *CTABread(char *fname) ;
int         CTABfree(COLOR_TABLE **pct) ;
int         CTABwriteInto(FILE *fp, COLOR_TABLE *ct) ;
COLOR_TABLE *CTABreadFrom(FILE *fp) ;



#endif
