#ifndef BOX_H
#define BOX_H

typedef struct
{
  int  x0 ;
  int  y0 ;
  int  z0 ;
  int  x1 ;
  int  y1 ;
  int  z1 ;
  int  width ;
  int  height ;
  int  depth ;
} BOX ;


int   BoxPrint(BOX *box, FILE *fp) ;
int   BoxExpand(BOX *box_src, BOX *box_dst, int dx, int dy, int dz) ;

#endif
