/**
 * @file  TexFont.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */



/* Copyright (c) Mark J. Kilgard, 1997. */

/* This program is freely distributable without licensing fees  and is
   provided without guarantee or warrantee expressed or  implied. This
   program is -not- in the public domain. */

#ifndef __TEXFONT_H__
#define __TEXFONT_H__

#include <GL/gl.h>

#define TXF_FORMAT_BYTE   0
#define TXF_FORMAT_BITMAP 1

typedef struct
{
  unsigned short c;       /* Potentially support 16-bit glyphs. */
  unsigned char width;
  unsigned char height;
  signed char xoffset;
  signed char yoffset;
  signed char advance;
  char dummy;           /* Space holder for alignment reasons. */
  short x;
  short y;
}
TexGlyphInfo;

typedef struct
{
  GLfloat t0[2];
  GLshort v0[2];
  GLfloat t1[2];
  GLshort v1[2];
  GLfloat t2[2];
  GLshort v2[2];
  GLfloat t3[2];
  GLshort v3[2];
  GLfloat advance;
}
TexGlyphVertexInfo;

typedef struct
{
  GLuint texobj;
  int tex_width;
  int tex_height;
  int max_ascent;
  int max_descent;
  int num_glyphs;
  int min_glyph;
  int range;
  unsigned char *teximage;
  TexGlyphInfo *tgi;
  TexGlyphVertexInfo *tgvi;
  TexGlyphVertexInfo **lut;
}
TexFont;

extern char *txfErrorString(void);

extern TexFont *txfLoadFont(
    char *filename);

extern void txfUnloadFont(
    TexFont * txf);

extern GLuint txfEstablishTexture(
    TexFont * txf,
    GLuint texobj,
    GLboolean setupMipmaps);

extern void txfBindFontTexture(
    TexFont * txf);

extern void txfGetStringMetrics(
    TexFont * txf,
    char *string,
    int len,
    int *width,
    int *max_ascent,
    int *max_descent);

extern void txfRenderGlyph(
    TexFont * txf,
    int c);

extern void txfRenderString(
    TexFont * txf,
    char *string,
    int len);

extern void txfRenderFancyString(
    TexFont * txf,
    char *string,
    int len);

extern int txfInFont(TexFont * txf, int c) ;

#endif /* __TEXFONT_H__ */
