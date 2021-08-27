/**
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef xTypes_H
#define xTypes_H

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef NULL
#define NULL 0
#endif

typedef unsigned char tBoolean;
typedef long          tSignature;

typedef struct
{
  int mnX, mnY;
}
xPoint2n, *xPoint2nRef;

typedef struct
{
  float mfX, mfY;
}
xPoint2f, *xPoint2fRef;

typedef struct
{
  int mnRed, mnGreen, mnBlue;
}
xColor3n, *xColor3nRef;

typedef struct
{
  float mfRed, mfGreen, mfBlue;
}
xColor3f, *xColor3fRef;

typedef struct
{
  float mfX, mfY, mfZ;
}
xVector3f, *xVector3fRef;

typedef enum
{
  tAxis_X = 0,
  tAxis_Y,
  tAxis_Z,
  knNumAxes
} tAxis;

typedef enum
{
  xColr_tComponent_None = 0,
  xColr_tComponent_Red,
  xColr_tComponent_Green,
  xColr_tComponent_Blue,
  xColr_knNumComponents
} xColr_tComponent;

void  xColr_SetFloat            ( xColor3fRef iColor,
                                  float ifRed, float ifGreen, float ifBlue );
void  xColr_SetFloatComponent   ( xColor3fRef iColor,
                                  xColr_tComponent iComponent,
                                  float ifValue );
float xColr_GetFloatComponent   ( xColor3fRef iColor,
                                  xColr_tComponent iComponent );
void  xColr_SetInt              ( xColor3nRef iColor,
                                  int ifRed, int ifGreen, int ifBlue );
void  xColr_SetIntComponent     ( xColor3nRef iColor,
                                  xColr_tComponent iComponent,
                                  int inValue );
int   xColr_GetIntComponent     ( xColor3nRef iColor,
                                  xColr_tComponent iComponent );
void xColr_PackFloatArray  ( xColor3fRef iColor,
                             float*      iafColor );

void xColr_HilightComponent ( xColor3fRef      iColor,
                              xColr_tComponent iComponent );

#define xColr_ExpandFloat(iColor) (iColor)->mfRed,(iColor)->mfGreen,(iColor)->mfBlue

#endif
