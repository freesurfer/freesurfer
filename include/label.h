#ifndef LABEL_H
#define LABEL_H

typedef struct
{
  int    n_points ;           /* # of points in area */
  char   name[100] ;          /* name of label file */
  char   subject_name[100] ;  /* name of subject */
  int    *vno ;       /* indices into surface tessellation */
  float  *x ;
  float  *y ;
  float  *z ;
} LABEL ;


int LabelFree(LABEL **parea) ;
int LabelDump(FILE *fp, LABEL *area) ;
LABEL *LabelRead(char *subject_name, char *label_name) ;
int LabelWrite(LABEL *area, char *fname) ;
int LabelToCanonical(LABEL *area, MRI_SURFACE *mris) ;
int LabelFromCanonical(LABEL *area, MRI_SURFACE *mris) ;
int LabelToFlat(LABEL *area, MRI_SURFACE *mris) ;

#endif
