#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "error.h"
#include "colortab.h"
#include "fio.h"
#include "proto.h"

COLOR_TABLE *
CTABread(char *fname)
{
	COLOR_TABLE *ct ;
	COLOR_TABLE_ENTRY *cte ;
	char        line[STRLEN], *cp ;
	int         nbins ;
	FILE        *fp ;

	fp = fopen(fname, "r") ;
	if (fp == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABread(%s): could not open file", fname)) ;

	nbins = 0 ;
	while ((cp = fgetl(line, STRLEN, fp)) != NULL)
		nbins++ ;

	if (nbins <= 0)
	{
		fclose(fp) ;
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABread(%s): badly formed file", fname)) ;
	}

	ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE)) ;
	if (ct == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABread(%s): could not allocate table", fname)) ;
	ct->bins = (COLOR_TABLE_ENTRY *)calloc(nbins, sizeof(COLOR_TABLE_ENTRY)) ;
	if (ct->bins == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABread(%s): could not allocate %d bin table", fname, nbins)) ;
	strcpy(ct->fname, fname) ;
	ct->nbins = nbins ;  

	rewind(fp) ; nbins = 0 ;
	while ((cp = fgetl(line, STRLEN, fp)) != NULL)
	{
		cte = &ct->bins[nbins] ;
		if (sscanf(line, "%d %s %d %d %d %d", &cte->index,cte->name,&cte->r,&cte->g,&cte->b,&cte->trans) != 6)
		{
			printf("CTABread(%s): could not scan 6 parms from line %s (%d)\n",fname,cp, nbins) ;
			continue ;
		}
		nbins++ ;
	}


	fclose(fp) ;
	return(ct) ;
}

int
CTABfree(COLOR_TABLE **pct)
{
	COLOR_TABLE *ct = *pct ;
	*pct = NULL ;
	free(ct->bins) ;
	free(ct) ;
	return(NO_ERROR) ;
}

int
CTABwriteInto(FILE *fp, COLOR_TABLE *ct)
{
	int                i ;
	COLOR_TABLE_ENTRY  *cte ;

	fwriteInt(ct->nbins, fp) ;
	fwriteInt(strlen(ct->fname)+1, fp) ;
	fwrite(ct->fname, sizeof(char), strlen(ct->fname)+1, fp) ;
	for (i = 0 ; i < ct->nbins ; i++)
	{
		cte = &ct->bins[i] ;
		fwriteInt(strlen(cte->name)+1, fp) ;
		fwrite(cte->name, sizeof(char), strlen(cte->name)+1, fp) ;
		fwriteInt(cte->r, fp) ;
		fwriteInt(cte->g, fp) ;
		fwriteInt(cte->b, fp) ;
		fwriteInt(cte->trans, fp) ;
	}
	return(NO_ERROR) ;
}

COLOR_TABLE *
CTABreadFrom(FILE *fp)
{
	COLOR_TABLE        *ct ;
	int                i, nbins, len ;
	COLOR_TABLE_ENTRY  *cte ;

	nbins = freadInt(fp) ;
	ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE)) ;
	if (ct == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABreadFrom: could not allocate table")) ;
	ct->bins = (COLOR_TABLE_ENTRY *)calloc(nbins, sizeof(COLOR_TABLE_ENTRY)) ;
	if (ct->bins == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABreadFrom: could not allocate %d bin table", nbins)) ;
	ct->nbins = nbins ;

	len = freadInt(fp) ;
	fread(ct->fname, sizeof(char), len, fp) ;
	for (i = 0 ; i < ct->nbins ; i++)
	{
		cte = &ct->bins[i] ;
		len = freadInt(fp) ;
		fread(cte->name, sizeof(char), len, fp) ;
		cte->r = freadInt(fp) ;
		cte->g = freadInt(fp) ;
		cte->b = freadInt(fp) ;
		cte->trans   = freadInt(fp) ;
		cte->index = i ;
	}
	return(ct) ;
}

int
CTABcolorToIndex(COLOR_TABLE *pct, int r, int g, int b, int*index)
{
  CTE *bin;
  int nbin;

  if (NULL == pct)
    return(ERROR_BAD_PARM);

  nbin = 0;
  while (nbin < pct->nbins)
    {
      bin = &(pct->bins[nbin]);
      if (bin->r == r &&
	  bin->g == g &&
	  bin->b == b)
	{
	  *index = nbin;
	  return(NO_ERROR);
	}
      
      nbin++;
    }
  
  *index = -1;
  return(ERROR_BAD_PARM);
}

int
CTABindexToColor(COLOR_TABLE *pct, int index, int*r, int*g, int*b)
{
  CTE *bin;

  if (NULL == pct)
    return(ERROR_BAD_PARM);

  if (index < 0 || index >= pct->nbins)
    return(ERROR_BAD_PARM);
    
  bin = &(pct->bins[index]);
  *r = bin->r;
  *g = bin->g;
  *b = bin->b;

  return(NO_ERROR);
}

int
CTABindexToAnnotation(COLOR_TABLE *ct, int index, int *pannot)
{
  CTE *bin;
	int  annotation, i ;

  if (NULL == ct)
    return(ERROR_BAD_PARM);

	for (i = 0 ; i < ct->nbins ; i++)
	{
		if (ct->bins[i].index == index)
			break ;
	}
	if (i >= ct->nbins)
		ErrorReturn(ERROR_BAD_PARM, (ERROR_BAD_PARM, "CTABindexToAnnotation(%s, %d): could not find index in table", ct->fname, index)) ;

  bin = &(ct->bins[i]);
	annotation = (bin->b << 16) + (bin->g << 8) + bin->r ;
	*pannot = annotation ;

  return(NO_ERROR);
}

int
CTABannotationToIndex(COLOR_TABLE *ctab, int annotation)
{
  int   r, g, b, index ;

	if (ctab == NULL)
		return(-1) ;
	r = annotation & 0x0000ff ;
	g = (annotation >> 8) & 0x0000ff ;
	b = (annotation >> 16) & 0x0000ff ;
	CTABcolorToIndex(ctab, r, g, b, &index) ;
	return(index) ;
}

int
CTABalphaLevel(COLOR_TABLE *pct, int index, int *alpha)
{
  CTE *bin;

  if (NULL == pct)
    return(ERROR_BAD_PARM);

  if (index < 0 || index >= pct->nbins)
    return(ERROR_BAD_PARM);
    
  /* We return 255 - the transparency value here to make it an alpha
     level. */
  bin = &(pct->bins[index]);
  *alpha = 255 - bin->trans;

  return(NO_ERROR);
}

int
CTABcopyName(COLOR_TABLE *pct, int index, char *name)
{
  CTE *bin;

  if (NULL == pct)
    return(ERROR_BAD_PARM);

  if (index < 0 || index >= pct->nbins)
    return(ERROR_BAD_PARM);
    
  bin = &(pct->bins[index]);
  strcpy (name, bin->name);

  return(NO_ERROR);
}

int
CTABnameToIndex(COLOR_TABLE *ctab, char *name)
{
  CTE *bin;
  int nbin;

  if (NULL == ctab)
    return(-1);

  nbin = 0;
  while (nbin < ctab->nbins)
	{
		bin = &(ctab->bins[nbin]);
		if (!stricmp(name, bin->name))
		{
			return(nbin) ;
		}
      
		nbin++;
	}
  
  return(-1);
}

int
CTABnameToAnnotation(COLOR_TABLE *ctab, char *name)
{
  CTE *bin;
  int nbin, annotation ;

  if (NULL == ctab)
    return(-1);

  nbin = 0;
  while (nbin < ctab->nbins)
	{
		bin = &(ctab->bins[nbin]);
		if (!stricmp(name, bin->name))
		{
			annotation = (bin->b << 16) + (bin->g << 8) + bin->r ;
			return(annotation) ;
		}
      
		nbin++;
	}
  
  return(-1);
}

/*--------------------------------------------------------
  CTABindexToItemNo() - returns the 0-based item number in the list
  corresponding to the index. Unforutnately, the term "index" would be
  a better name for the item number and "id" would be a better name
  for index, but I thought it would be best to stick with the notation
  used throught out the file. 
  --------------------------------------------------------*/
int CTABindexToItemNo(COLOR_TABLE *ctab, int index)
{
  int n;
  
  for(n=0; n < ctab->nbins; n++)
    if(ctab->bins[n].index == index) return(n);

  return(-1);
}
/*----------------------------------------------------------
  CTABprint() - print color table to the terminal
  ----------------------------------------------------------*/
int CTABprint(FILE *fp, COLOR_TABLE *ctab)
{
  int n;
  
  for(n=0; n < ctab->nbins; n++)
    fprintf(fp,"%3d  %-30s  %3d %3d %3d   0\n",n,ctab->bins[n].name,
	    ctab->bins[n].r,ctab->bins[n].g,ctab->bins[n].b);
  return(0);
}

/*----------------------------------------------------------
  CTABwriteTxt() - print color table to a file.
  ----------------------------------------------------------*/
int CTABwriteTxt(char *fname, COLOR_TABLE *ctab)
{
  FILE *fp;
  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: could not open %s for writing\n",fname);
    return(1);
  }
  CTABprint(fp, ctab);
  fclose(fp);
  return(0);
}

COLOR_TABLE *
CTABalloc(int nbins)
{
	COLOR_TABLE        *ct ;
	int                i ;
	COLOR_TABLE_ENTRY  *cte ;

	ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE)) ;
	if (ct == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABalloc(%d): could not allocate table", nbins)) ;
	ct->bins = (COLOR_TABLE_ENTRY *)calloc(nbins, sizeof(COLOR_TABLE_ENTRY)) ;
	if (ct->bins == NULL)
		ErrorReturn(NULL, (ERROR_NOFILE, "CTABalloc: could not allocate %d bin table", nbins)) ;
	strcpy(ct->fname, "none") ;
	ct->nbins = nbins ;  

	for (i = 0 ; i < nbins ; i++)
	{
		cte = &ct->bins[i] ;
		cte->r = nint(randomNumber(0, 255)) ;
		cte->g = nint(randomNumber(0, 255)) ;
		cte->b = nint(randomNumber(0, 255)) ;
		cte->index = i ;
		sprintf(cte->name, "cluster%d", i) ;
	}

	return(ct) ;
}


