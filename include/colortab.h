#ifndef COLORTAB_H
#define COLORTAB_H

#include "const.h"
#include "stdio.h"

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

int         CTABcolorToIndex(COLOR_TABLE *pct, int r, int g, int b, int*index);
int         CTABindexToColor(COLOR_TABLE *pct, int index, int*r, int*g, int*b);
int         CTABindexToAnnotation(COLOR_TABLE *pct, int index, int *pannot);
int         CTABcopyName(COLOR_TABLE *pct, int index, char *name);
int         CTABannotationToIndex(COLOR_TABLE *ctab, int annotation) ;
int         CTABnameToIndex(COLOR_TABLE *ctab, char *name) ;
int         CTABnameToAnnotation(COLOR_TABLE *ctab, char *name) ;
int         CTABindexToItemNo(COLOR_TABLE *ctab, int index);
int         CTABwriteTxt(char *fname, COLOR_TABLE *ctab);
int         CTABprint(FILE *fp, COLOR_TABLE *ctab);

#endif
