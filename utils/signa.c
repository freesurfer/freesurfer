#include <stdio.h>
#include <stdlib.h>

#include "machine.h"
#include "mri.h"
#include "utils.h"
#include "error.h"
#include "signa.h"
#include "mghendian.h"
/*#include "idbm_hdr_def.h"*/

static int orderShortBuffer(short *sbuf, int nbytes) ;


#if 0
MRI *
read_signa(char *fname, char *h, float scale)
{
  int i,j;
  unsigned long k;
  FILE *fp;
  int min=10000,max=0,n;
  unsigned short *buf ;
  MRI *mri ;
  HINFO  header ;

  mri = MRIalloc(256,256,256, MRI_SHORT) ;
  fp = fopen(fname,"rb");
  if (fp==NULL) 
    ErrorReturn(NULL, (ERROR_NOFILE, "read_signa(%s): file not found.\n",
                       fname)) ;

  fread(h,sizeof(char),HLENGTH,fp);
  
  fread(buf,sizeof(short),IMGSIZE*IMGSIZE,fp);
  k=0;
  for (i=0;i<IMGSIZE;i++)
    for (j=0;j<IMGSIZE;j++)
    {
      if (buf[k]<min) min=buf[k];
      if (buf[k]>max) max=buf[k];
      n = (int)(buf[k]/scale+0.5);
      im[i][j] = (n<0)?0:(n>255)?255:n;
      k++;
    }

  fclose(fp);
  printf("File %s read\n",fname);
  printf("min=%d,max=%d\n",min,max);
  return(mri) ;
}
#endif
/*------------------------------------------------------------*/
/* Get Header (Character) */
char *ghc( char *destin, char *header, int offset, int byte_length )
{
  strncpy( destin, header+(2*offset), byte_length );
  return(destin);
}
/*------------------------------------------------------------*/
/* Get Header (integer) */
int ghi( char *header, int offset )
{
  int num;
  int i;
  int byte_length = 2;

  num = 0 ;
  for ( i=0; i < byte_length ; i++ )
  {
   num = ( num << 8) | *(header + (2*offset) + i )  ;
  }

  return(num);
}
/*------------------------------------------------------------*/
/* Get Header (float) */
float ghf( char *header, int offset )
{
#define sign_bit 020000000000
#define dmantissa 077777777
#define dexponent 0177
#define dmantlen 24
#define smantissa 037777777
#define sexponent 0377
#define smantlen 23

  union { float ff; long ii; char buf[4]; } thing ;
  long k;
  long dg_exp, dg_sign, dg_mantissa;
  long sun_exp;

/* read in to char string */
  thing.buf[0] = *(header + (2*offset) + 0 );
  thing.buf[1] = *(header + (2*offset) + 1 );
  thing.buf[2] = *(header + (2*offset) + 2 );
  thing.buf[3] = *(header + (2*offset) + 3 );
  thing.ff = orderFloatBytes(thing.ff) ;

/* convert mv floating point numbers to sun float - from GE's mvtsunf */
  dg_exp = ( thing.ii >> 24 ) & dexponent;
  dg_sign = thing.ii  & sign_bit;
  dg_mantissa = ( thing.ii  & dmantissa ) << 8;

  sun_exp = 4 * ( dg_exp - 64 );
  while ( (k = dg_mantissa & sign_bit ) == 0 && dg_mantissa != 0) {
      sun_exp--;
      dg_mantissa = dg_mantissa << 1;
  }

  sun_exp += 126;
  if (sun_exp < 0) sun_exp = 0;
  if (sun_exp > 255) sun_exp = 255;

  dg_mantissa = dg_mantissa<<1;

  thing.ii = 0;
  thing.ii = ( dg_sign ) |
                   ( sun_exp << smantlen ) |
                   ( (dg_mantissa >> 9) & smantissa );

  return(thing.ff/1000);
}
/*------------------------------------------------------------*/

/* Get header info */
int
get_signa_header_info(char *h,HINFO *hinfo)
{
  hinfo->x = ghi( h, IHDR_START+IHDR_X ) ;
  hinfo->y = ghi( h, IHDR_START+IHDR_Y ) ;
  hinfo->tr = ghf( h, IHDR_START+IHDR_TR ) ;
  hinfo->ti = ghf( h, IHDR_START+IHDR_TI ) ;
  hinfo->te = ghf( h, IHDR_START+IHDR_TE ) ;
  hinfo->strtx = ghf( h, IHDR_START+IHDR_STRTX ) ;
  hinfo->endx = ghf( h, IHDR_START+IHDR_ENDX ) ;
  hinfo->strty = ghf( h, IHDR_START+IHDR_STRTY ) ;
  hinfo->endy = ghf( h, IHDR_START+IHDR_ENDY ) ;
  hinfo->strtz = ghf( h, IHDR_START+IHDR_STRTZ ) ;
  hinfo->endz = ghf( h, IHDR_START+IHDR_ENDZ ) ;
  hinfo->locatn = ghf( h, IHDR_START+IHDR_LOCATN ) ;  
  hinfo->fov = ghf( h, SEHDR_START+SEHDR_FOV ) ;
  hinfo->psiz = ghf( h, IHDR_START+IHDR_PIXSIZ ) ;
  hinfo->thick = ghf( h, IHDR_START+IHDR_THICK ) ;
  hinfo->ptype = ghi( h, SEHDR_START+SEHDR_PTYPE ) ;
  hinfo->imnr0 = 1 ;
  hinfo->imnr1 = ghi( h, SEHDR_START+SEHDR_IALLOC ) ;
  return(NO_ERROR) ;
}
int
is_signa(char *fname)
{
  HINFO  header ;
  char   h[HLENGTH+2] ;
  FILE   *fp ;
  int    ret ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    return(0) ;

  if ((ret = fread(h,sizeof(char),HLENGTH,fp)) != HLENGTH)
    return(0) ;
  get_signa_header_info((char *)&h, &header) ;
  fclose(fp) ;

  return(header.x > 2 && header.x < 1024*10 &&
         header.y > 2 && header.y < 1024*10 &&
         header.tr > 0 && header.tr < 1000000 &&
         header.fov > .0002 && header.fov < 300 &&
         (fabs(header.strtx) + fabs(header.endx) > .000001) &&
         (fabs(header.strtz) + fabs(header.endy) > .000001) &&
         (fabs(header.strty) + fabs(header.endz) > .000001) &&
         fabs(header.ti) < 100000 && fabs(header.te) < 10000) ;
}

MRI *
signaRead(char *fname, int read_volume_flag)
{
  HINFO  header ;
  char   h[HLENGTH+2] ;
  FILE   *fp ;
  int    ret, i ;
  MRI    *mri ;
  char   path[STRLEN] ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(0, (ERROR_NOFILE, "is_signa(%s): could not open file",fname)) ;

  if ((ret = fread(h,sizeof(char),HLENGTH,fp)) != HLENGTH)
    ErrorReturn(0, (ERROR_BADFILE, 
                    "is_signa(%s): could not read %d byte header",
                    fname, HLENGTH)) ;
  get_signa_header_info((char *)&h, &header) ;

  if (header.imnr1 < 0)
  {
    fprintf(stderr, "WARNING: # of slices=%d in header - assuming 124...\n",
            header.imnr1) ;
    header.imnr1 = 124 ;
  }
  mri = MRIalloc(header.x, header.y, header.imnr1-header.imnr0+1, MRI_SHORT) ;
  if (!mri)
    ErrorReturn(NULL, (ERROR_NOMEMORY, 
                       "signaRead(%s): could not read %dx%dx%d volume",
                       fname, header.x,header.y, header.imnr1-header.imnr0+1));

  mri->xstart = 1000*header.strtx ; mri->xend = 1000*header.endx ;
  mri->ystart = 1000*header.strty ; mri->yend = 1000*header.endy ;
  mri->zstart = 1000*header.strtz ; mri->zend = 1000*header.endz ;
  mri->tr = header.tr ; mri->ti = header.ti ; mri->te = header.te ;
  mri->location = 1000*header.locatn ;
  mri->fov = 1000*header.fov ;
  mri->xsize = mri->ysize = mri->ps = 1000*header.psiz;
  mri->zsize = mri->thick = 1000*header.thick;
  // no orientation info and thus sets to coronal
  setDirectionCosine(mri, MRI_CORONAL); 

  fclose(fp) ;
  if (!read_volume_flag)
    return(mri) ;

  FileNamePath(fname, path) ;
  for (i = mri->imnr0 ; i <= mri->imnr1 ; i++)
  {
    sprintf(fname, "%s/I.%03d", path, i) ;

    fp = fopen(fname, "rb") ;
    if (!fp)
      ErrorReturn(0, (ERROR_NOFILE, "read_signa(%s): could not open file",
                      fname)) ;

    fseek(fp, HLENGTH, SEEK_SET) ;
    fread(&MRISvox(mri, 0, 0, i-mri->imnr0),
          sizeof(short),mri->width*mri->height,fp);
    orderShortBuffer(&MRISvox(mri,0,0,i-mri->imnr0),mri->width*mri->height) ;
    fclose(fp) ;
  }

  /* now figure out what the image #s start and end at */
  return(mri) ;
}

static int
orderShortBuffer(short *sbuf, int nbytes)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  int   i ;

  for (i = 0 ; i < nbytes ; i++, sbuf++)
    *sbuf = orderShortBytes(*sbuf) ;
#endif
  return(NO_ERROR) ;
}

