#include "xTypes.h"

#define xColr_kfHilightAmt 0.3
#define xColr_kfMinHighlightDistance 0.3

void xColr_Set ( xColor3fRef iColor, 
      float ifRed, float ifGreen, float ifBlue ) {
  iColor->mfRed   = ifRed;
  iColor->mfGreen = ifGreen;
  iColor->mfBlue  = ifBlue;
}

void xColr_PackFloatArray ( xColor3fRef iColor,
          float*      iafColor ) {
  iafColor[0] = iColor->mfRed;
  iafColor[1] = iColor->mfGreen;
  iafColor[2] = iColor->mfBlue;
}

void xColr_HilightComponent ( xColor3fRef      iColor,
            xColr_tComponent iComponent ) {

  if( NULL == iColor )
    return;

  switch( iComponent ) {
  case xColr_tComponent_Red:
    iColor->mfRed += xColr_kfHilightAmt;
    if( iColor->mfRed > 1.0 )
      iColor->mfRed = 1.0;
    if( iColor->mfRed - iColor->mfGreen < xColr_kfMinHighlightDistance )
      iColor->mfGreen = iColor->mfRed - xColr_kfMinHighlightDistance;
    if( iColor->mfRed - iColor->mfBlue < xColr_kfMinHighlightDistance )
      iColor->mfBlue = iColor->mfRed - xColr_kfMinHighlightDistance;
    break;
  case xColr_tComponent_Green:
    iColor->mfGreen += xColr_kfHilightAmt;
    if( iColor->mfGreen > 1.0 )
      iColor->mfGreen = 1.0;
    if( iColor->mfGreen - iColor->mfRed < xColr_kfMinHighlightDistance )
      iColor->mfRed = iColor->mfGreen - xColr_kfMinHighlightDistance;
    if( iColor->mfGreen - iColor->mfBlue < xColr_kfMinHighlightDistance )
      iColor->mfBlue = iColor->mfGreen - xColr_kfMinHighlightDistance;
    break;
  case xColr_tComponent_Blue:
    iColor->mfBlue += xColr_kfHilightAmt;
    if( iColor->mfBlue > 1.0 )
      iColor->mfBlue = 1.0;
    if( iColor->mfBlue - iColor->mfGreen < xColr_kfMinHighlightDistance )
      iColor->mfGreen = iColor->mfBlue - xColr_kfMinHighlightDistance;
    if( iColor->mfBlue - iColor->mfRed < xColr_kfMinHighlightDistance )
      iColor->mfRed = iColor->mfBlue - xColr_kfMinHighlightDistance;
    break;
  default:
    return;    
  }
}
