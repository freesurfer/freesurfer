/* @(#)pr_line.h 1.15 88/02/08 SMI */

/*
 * Copyright 1986 by Sun Microsystems, Inc.
 */


#ifndef pr_line_h_DEFINED
#define pr_line_h_DEFINED

#define POLY_CLOSE ((u_char *) 1)
#define POLY_DONTCLOSE ((u_char *) 0)

extern short pr_tex_dotted[];
extern short pr_tex_dashed[];
extern short pr_tex_dashdot[];
extern short pr_tex_dashdotdotted[];
extern short pr_tex_longdashed[];


typedef
struct pr_texture
{
  short *pattern;
  short offset;
  struct pr_texture_options
  {
unsigned startpoint :
    1,
endpoint :
    1,
balanced :
    1,
givenpattern :
    1,
res_fat :
    1,
res_poly :
    1,
res_mvlist :
    1,
res_right :
    1,
res_close :
    1,
res_cliprt :
    1;
  }
  options;
  short res_polyoff;
  short res_oldpatln;
  short res_fatoff;
  short *res_patfat;
  short res_numsegs;
}
Pr_texture;


typedef
struct pr_brush
{
  int width;
}
Pr_brush;

#endif /*pr_line_h_DEFINED */



