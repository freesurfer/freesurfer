#ifndef MATFILE_H
#define MATFILE_H
       
#ifdef Darwin
typedef unsigned char  Byte;
#endif

#include "znzlib.h"  
#include "matrix.h"
#include "machine.h"

typedef struct
{
    long32  type ;
    long32  mrows ;
    long32  ncols ;
    long32  imagf ;
    long32  namlen ;
} MATHD ;

typedef struct
{
    /*common*/
    long32  type;
    long32  mrows ;
    long32  ncols ;
    long32  imagf;
    long32  namlen ;
    int     version ;
    char  *data ;
    char  *idata ;
    /*only 5*/  
    char endian;
} MATFILE ;

typedef struct
{
  char *mfile;
  int nvars;
  char *varname[1000];
  MATRIX *varmtx[1000];
} MATFILECONTENTS, MLFC;

char *MatReadHeader0(FILE *fp, MATFILE *mf);
char    *MatReadHeader(FILE *fp, MATFILE *mf, long32 *compressed) ;
char    *znzMatReadHeader(FILE *fp, MATFILE *mf, char **data) ;
MATFILE *MatFileRead(const char *fname, int type) ;
MATRIX  *MatlabRead(const char *fname) ;
MATRIX  *MatlabRead2(const char *fname) ;
int     MatlabWrite(MATRIX *mat, const char *fname, char *name) ;
int     MatFileWrite(const char *fname,
                    float *data, int rows, int cols, char *name) ;
int Matlab_Install_printf( int (*new_printf)(const char *szFormat, ...) );
MLFC *ReadMatlabFileContents(const char *fname);
int   MLFCprint(FILE *fp, MLFC *mlfc);
int MLFCfree(MLFC **ppmlfc);
MATRIX *ReadMatlabFileVariable(char *fname, char *varname);


#define MAT_BYTE     0
#define MAT_DOUBLE   1
#define MAT_INT      2
#define MAT_SHORT    3
#define MAT_FLOAT    4

#define MATFILE_PC      0000
#define MATFILE_SPARC   1000

#define MATFILE_PC5     'I' 
#define MATFILE_SPARC5  'M' 

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

