#ifndef MATFILE_H
#define MATFILE_H
       
#include "matrix.h"

typedef struct
{
    long  type ;
    long  mrows ;
    long  ncols ;
    long  imagf ;
    long  namlen ;
} MATHD ;

typedef struct
{
    long  type ;
    long  mrows ;
    long  ncols ;
    long  imagf ;
    long  namlen ;
    char  *data ;
    char  *idata ;
} MATFILE ;

char    *MatReadHeader(FILE *fp, MATFILE *mf) ;
MATFILE *MatFileRead(const char *fname, int type) ;
MATRIX  *MatlabRead(const char *fname) ;
int     MatlabWrite(MATRIX *mat, const char *fname, char *name) ;
int     MatFileWrite(const char *fname,
                    float *data, int rows, int cols, char *name) ;
int Matlab_Install_printf( int (*new_printf)(const char *szFormat, ...) );

#define MAT_BYTE     0
#define MAT_DOUBLE   1
#define MAT_INT      2
#define MAT_SHORT    3
#define MAT_FLOAT    4

#define MATFILE_PC      0000
#define MATFILE_SPARC   1000

#define MATFILE_DOUBLE  00
#define MATFILE_FLOAT   10
#define MATFILE_LONG    20
#define MATFILE_SHORT   30
#define MATFILE_USHORT  40
#define MATFILE_BYTE    50

#define MATFILE_FULL    0
#define MATFILE_TEXT    1
#define MATFILE_SPARSE  2


#endif

