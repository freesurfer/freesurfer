#ifndef LPAFILE_H
#define LPAFILE_H


#define NPOINTS   4

typedef struct
{
  int  xp[NPOINTS] ;    /* x coordinates of corners */
  int  yp[NPOINTS] ;    /* y coordinates of corners */
  int  xc ;             /* x coordinate of centroid */
  int  yc ;             /* y coordinate of centroid */
  long fpos ;           /* position in answer file */
} LP_BOX ;

typedef struct
{
  char    fname[100] ;     /* name of the file */
  FILE    *fp ;
  char    **filelist ;     /* array of pointers */
  int     nfiles ;         /* size of filelist */
  LP_BOX  *coords ;
  int     last_written ;   /* index of last written lp_box */
  int     current ;
} LP_ANSWER_FILE, LPAF ;

LPAF *LPAFcreate(char *out_fname, int argc, char *argv[]) ;
int  LPAFwrite(LPAF *lpaf, int current) ;
int  LPAFread(LPAF *lpaf, int current) ;
int  LPAFset(LPAF *lpaf, int current, int *xp, int *yp, int xc, int yc) ;

#endif
