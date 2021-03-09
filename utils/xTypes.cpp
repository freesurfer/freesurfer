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

#include "xTypes.h"

#define xColr_kfHilightAmt 0.3
#define xColr_kfMinHighlightDistance 0.3

void xColr_SetFloat(xColor3fRef iColor, float ifRed, float ifGreen, float ifBlue)
{
  iColor->mfRed = ifRed;
  iColor->mfGreen = ifGreen;
  iColor->mfBlue = ifBlue;
}

void xColr_SetFloatComponent(xColor3fRef iColor, xColr_tComponent iComponent, float ifValue)
{
  switch (iComponent) {
    case xColr_tComponent_Red:
      iColor->mfRed = ifValue;
      break;
    case xColr_tComponent_Green:
      iColor->mfGreen = ifValue;
      break;
    case xColr_tComponent_Blue:
      iColor->mfBlue = ifValue;
      break;
    default:
      return;
  }
}

float xColr_GetFloatComponent(xColor3fRef iColor, xColr_tComponent iComponent)
{
  switch (iComponent) {
    case xColr_tComponent_Red:
      return iColor->mfRed;
      break;
    case xColr_tComponent_Green:
      return iColor->mfGreen;
      break;
    case xColr_tComponent_Blue:
      return iColor->mfBlue;
      break;
    default:
      return 0;
  }
}

void xColr_SetInt(xColor3nRef iColor, int inRed, int inGreen, int inBlue)
{
  iColor->mnRed = inRed;
  iColor->mnGreen = inGreen;
  iColor->mnBlue = inBlue;
}

void xColr_SetIntComponent(xColor3nRef iColor, xColr_tComponent iComponent, int inValue)
{
  switch (iComponent) {
    case xColr_tComponent_Red:
      iColor->mnRed = inValue;
      break;
    case xColr_tComponent_Green:
      iColor->mnGreen = inValue;
      break;
    case xColr_tComponent_Blue:
      iColor->mnBlue = inValue;
      break;
    default:
      return;
  }
}

int xColr_GetIntComponent(xColor3nRef iColor, xColr_tComponent iComponent)
{
  switch (iComponent) {
    case xColr_tComponent_Red:
      return iColor->mnRed;
      break;
    case xColr_tComponent_Green:
      return iColor->mnGreen;
      break;
    case xColr_tComponent_Blue:
      return iColor->mnBlue;
      break;
    default:
      return 0;
  }
}

void xColr_PackFloatArray(xColor3fRef iColor, float *iafColor)
{
  iafColor[0] = iColor->mfRed;
  iafColor[1] = iColor->mfGreen;
  iafColor[2] = iColor->mfBlue;
}

void xColr_HilightComponent(xColor3fRef iColor, xColr_tComponent iComponent)
{
  if (NULL == iColor) {
    return;
  }

  switch (iComponent) {
    case xColr_tComponent_Red:
      iColor->mfRed += xColr_kfHilightAmt;
      if (iColor->mfRed > 1.0) {
        iColor->mfRed = 1.0;
      }
      if (iColor->mfRed - iColor->mfGreen < xColr_kfMinHighlightDistance) {
        iColor->mfGreen = iColor->mfRed - xColr_kfMinHighlightDistance;
      }
      if (iColor->mfRed - iColor->mfBlue < xColr_kfMinHighlightDistance) {
        iColor->mfBlue = iColor->mfRed - xColr_kfMinHighlightDistance;
      }
      break;
    case xColr_tComponent_Green:
      iColor->mfGreen += xColr_kfHilightAmt;
      if (iColor->mfGreen > 1.0) {
        iColor->mfGreen = 1.0;
      }
      if (iColor->mfGreen - iColor->mfRed < xColr_kfMinHighlightDistance) {
        iColor->mfRed = iColor->mfGreen - xColr_kfMinHighlightDistance;
      }
      if (iColor->mfGreen - iColor->mfBlue < xColr_kfMinHighlightDistance) {
        iColor->mfBlue = iColor->mfGreen - xColr_kfMinHighlightDistance;
      }
      break;
    case xColr_tComponent_Blue:
      iColor->mfBlue += xColr_kfHilightAmt;
      if (iColor->mfBlue > 1.0) {
        iColor->mfBlue = 1.0;
      }
      if (iColor->mfBlue - iColor->mfGreen < xColr_kfMinHighlightDistance) {
        iColor->mfGreen = iColor->mfBlue - xColr_kfMinHighlightDistance;
      }
      if (iColor->mfBlue - iColor->mfRed < xColr_kfMinHighlightDistance) {
        iColor->mfRed = iColor->mfBlue - xColr_kfMinHighlightDistance;
      }
      break;
    default:
      return;
  }
}
