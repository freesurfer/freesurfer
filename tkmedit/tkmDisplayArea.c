/**
 * @file  tkmDisplayArea.c
 * @brief Graphics and UI interaction for data display.
 *
 * Manages a single pane of display in the window. Has slots for
 * viewing different kinds of data. Handles UI events in the panel.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/01 01:41:22 $
 *    $Revision: 1.144 $
 *
 * Copyright (C) 2002-2011, CorTechs Labs, Inc. (La Jolla, CA) and
 * The General Hospital Corporation (Boston, MA).
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer/CorTechs Software License Agreement' contained
 * in the file 'license.cortechs.txt' found in the FreeSurfer distribution,
 * and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferCorTechsLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "tkmDisplayArea.h"
#include "tkmMeditWindow.h"
#include "tkmFunctionalVolume.h"
#include "xUtilities.h"
#include "cma.h"
#include "utils.h"
#include "error.h"
#include "proto.h" //  nint

/* i'm not sure what to do about these y flips. it seems that whenever we're
   using a point that's going to go into the buffer to be drawn to the screen,
   we should flip when not in horizontal view. when using regular gl drawing
   commands to draw to the screen, only do it in horizontal orientation. */
#define BUFFER_Y_FLIP(y) ( this->mOrientation != mri_tOrientation_Horizontal? \
                          (this->mnVolumeSizeY - (y)) : (y) )
#define GLDRAW_Y_FLIP(y) ( this->mOrientation == mri_tOrientation_Horizontal? \
                          (this->mnVolumeSizeY - (y)) : (y) )
#define GLDRAW_Y_FLIP_FLOAT(y) \
                         ( this->mOrientation == mri_tOrientation_Horizontal? \
                          ((float)this->mnVolumeSizeY - (y)) : (y) )

//#define BUFFER_Y_FLIP(y) y
//#define GLDRAW_Y_FLIP(y) y

#define Y_FLIP(y)        (this->mnVolumeSizeY - (y))

/* these describe the direction a volume voxel goes in the buffer in the y
   direction. to handle stuff like:

   bufferPt.mnY = Y_FLIP(bufferPt.mnY);
   for( nY = bufferPt.mnY; nY < bufferPt.mnY + this->mnZoomLevel; nY++ ) ...

   use

   bufferPt.mnY = Y_FLIP(bufferPt.mnY);
   for( nY = bufferPt.mnY;
   BUFFER_Y_LT(nY,BUFFER_Y_INC(bufferPt.mnY,this->mnZoomLevel));
   nY = BUFFER_Y_INC(bufferPt.mnY,1) )...
*/

#define BUFFER_Y_LT(y,i)  (this->mOrientation != mri_tOrientation_Horizontal? \
                          ( (y) > (i) ) : ( (y) < (i) ))
#define BUFFER_Y_INC(y,i) (this->mOrientation != mri_tOrientation_Horizontal? \
                           ( (y) - (i) ) : ( (y) + (i) ) )


/* tool and brush info is static */
static DspA_tTool                sTool  = DspA_tTool_SelectVoxels;
static DspA_tBrushSettings       sBrush;
static DspA_tSegBrushSettings   sSegBrush;
static DspA_tFloodSelectSettings sFloodSelectSettings;

/* cursor info too */
static xColor3f     sCursorColor = {
                                     1, 0, 0
                                   };
static DspA_tMarker sCursorShape = DspA_tMarker_Crosshair;

/* static focused display. */
static tkmDisplayAreaRef sFocusedDisplay = NULL;

char *DspA_ksaErrorStrings [DspA_knNumErrorCodes] = {
      "No error",
      "Allocation failed.",
      "Invalid ptr to object, was probably null.",
      "Invalid signature, object was not valid.",
      "Invalid paramter, probably out of bounds.",
      "Invalid cursor.",
      "Invalid display flag.",
      "Invalid orientation.",
      "Invalid volume voxel.",
      "Invalid buffer point.",
      "Invalid screen point.",
      "Error accessing surface.",
      "Out of memory.",
      "Error accessing control points.",
      "Error accessing selection.",
      "Error accessing parent window.",
      "Error accessing functional volume.",
      "Error accessing surface cache list.",
      "Error accessing head point list.",
      "Error accessing list.",
      "Couldn't find a closest voxel.",
      "Error opening file.",
      "Invalid error code."
    };

char *DspA_ksaOrientation [mri_knNumOrientations] = {
      "Coronal",
      "Horizontal",
      "Sagittal"
    };

char *DspA_ksaSurface [Surf_knNumVertexSets] = {
      "Main",
      "Original",
      "Pial"
    };

char *DspA_ksaDisplaySet [DspA_knNumDisplaySets] = {
      "cursor",
      "mouseover"
    };

DspA_tErr DspA_New ( tkmDisplayAreaRef* oppWindow,
                     tkmMeditWindowRef  ipWindow ) {

  DspA_tErr         eResult      = DspA_tErr_NoErr;
  tkmDisplayAreaRef this         = NULL;
  int               nFlag        = 0;
  int               nVolume      = 0;
  int               nSegVolume   = 0;
  int               nSurface     = 0;
  int               nVertexSet   = 0;
  xColor3f          color;

  /* allocate us. */
  this = (tkmDisplayAreaRef) malloc( sizeof(tkmDisplayArea) );
  if ( NULL == this ) {
    eResult = DspA_tErr_AllocationFailed;
    goto error;
  }

  /* set the signature */
  this->mSignature = DspA_kSignature;

  /* set the parent window. */
  this->mpWindow             = ipWindow;
  this->mID                  = -1;
  this->mLocationInSuper.mnX = 0;
  this->mLocationInSuper.mnY = 0;

  /* set default window data */
  this->mnWidth     = 512;
  this->mnHeight    = 512;

  /* frame buffer starts out null until we get a volume */
  this->mpFrameBuffer       = NULL;
  this->mfFrameBufferScaleX = 1.0;
  this->mfFrameBufferScaleY = 1.0;

  /* allocate our voxels */
  xVoxl_New( &this->mpLastCursor );
  xVoxl_New( &this->mpCursor );
  xVoxl_New( &this->mpMouseLocationAnaIdx );
  xVoxl_New( &this->mpZoomCenter );
  xVoxl_New( &this->mpOriginalZoomCenter );

  /* stuff in default values for display states. */
  this->mOrientation           = mri_tOrientation_Coronal;
  this->mnZoomLevel            = 1;
  this->mnHilitedVertexIndex   = -1;
  this->mHilitedSurface        = Surf_tVertexSet_None;
  this->mbSliceChanged         = TRUE;
  this->mnVolumeSizeX           = 0;
  this->mnVolumeSizeY           = 0;
  this->mnVolumeSizeZ           = 0;
  this->mpSelectedHeadPoint     = NULL;
  this->mnSegmentationVolumeIndex        = -1;

  /* Set all line widths to 1. */
  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ )
    for ( nVertexSet = 0; nVertexSet < Surf_knNumVertexSets; nVertexSet++ )
      DspA_SetSurfaceLineWidth( this, nSurface, nVertexSet, 1 );

  /* Set default colors for surface lines. */
  xColr_SetFloat( &color, 1, 1, 0 );
  DspA_SetSurfaceLineColor( this, tkm_tSurfaceType_Main,
                            Surf_tVertexSet_Main, &color );
  DspA_SetSurfaceLineColor( this, tkm_tSurfaceType_Aux,
                            Surf_tVertexSet_Main, &color );
  xColr_SetFloat( &color, 0, 1, 0 );
  DspA_SetSurfaceLineColor( this, tkm_tSurfaceType_Main,
                            Surf_tVertexSet_Original, &color );
  DspA_SetSurfaceLineColor( this, tkm_tSurfaceType_Aux,
                            Surf_tVertexSet_Original, &color );
  xColr_SetFloat( &color, 1, 0, 0 );
  DspA_SetSurfaceLineColor( this, tkm_tSurfaceType_Main,
                            Surf_tVertexSet_Pial, &color );
  DspA_SetSurfaceLineColor( this, tkm_tSurfaceType_Aux,
                            Surf_tVertexSet_Pial, &color );

  this->mfSegmentationAlpha = 1.0;
  this->mfDTIAlpha = 1.0;
  this->mfFuncOverlayAlpha = 1.0;


  /* all our display flags start out false. */
  for ( nFlag = 0; nFlag < DspA_knNumDisplayFlags; nFlag++ )
    this->mabDisplayFlags[nFlag] = FALSE;

  /* null ptrs for display data. */
  for ( nVolume = 0; nVolume < tkm_knNumVolumeTypes; nVolume++ ) {
    this->mpVolume[nVolume] = NULL;
  }
  for ( nSegVolume = 0; nSegVolume < tkm_knNumSegTypes; nSegVolume++ ) {
    this->mSegmentationVolume[nSegVolume] = NULL;
    this->mSegmentationColorTable[nSegVolume] = NULL;
  }
  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {
    this->mpSurface[nSurface]       = NULL;
    this->maSurfaceLists[nSurface]  = NULL;
  }
  this->mpFunctionalVolume      = NULL;
  this->mpControlPoints         = NULL;
  this->mpSelection             = NULL;
  this->mHeadPoints             = NULL;
  this->mGCAVolume              = NULL;
  this->mGCATransform           = NULL;
  this->mVLI1                   = NULL;
  this->mVLI2                   = NULL;
  this->mpDTIVolume             = NULL;

  /* No line vertices yet. */
  this->mLineVertex1.mnX = this->mLineVertex1.mnY = -1;
  this->mLineVertex2.mnX = this->mLineVertex2.mnY = -1;
  this->mNumLineVoxels = 0;
  this->mLineDistance = 0;

  /* set default brush info */
  sBrush.mnRadius    = 1;
  sBrush.mShape      = DspA_tBrushShape_Square;
  sBrush.mb3D        = FALSE;
  sBrush.mb3DFill    = FALSE;
  sBrush.mFuzzy     = 0;
  sBrush.mDistance  = 0;
  DspA_SetBrushInfoToDefault( this, DspA_tBrush_EditOne );
  DspA_SetBrushInfoToDefault( this, DspA_tBrush_EditTwo );

  /* default seg brush info */
  sSegBrush.mNewValue    = 0;
  sSegBrush.mb3D         = FALSE;
  sSegBrush.mSrc         = tkm_tVolumeTarget_MainAna;
  sSegBrush.mFuzzy      = 0;
  sSegBrush.mDistance   = 0;
  sSegBrush.mnPaintValue = 0;
  sSegBrush.mnEraseValue = 0;

  /* default flood select info */
  sFloodSelectSettings.mb3D       = FALSE;
  sFloodSelectSettings.mSrc       = tkm_tVolumeTarget_MainAna;
  sFloodSelectSettings.mFuzzy    = 0;
  sFloodSelectSettings.mDistance = 0;

  /* set default cursor color */
  color.mfRed   = 1.0;
  color.mfGreen = 0.0;
  color.mfBlue  = 0.0;
  DspA_SetCursorColor( this, &color );

  /* return window. */
  *oppWindow = this;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_New: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_Delete ( tkmDisplayAreaRef* ioppWindow ) {

  DspA_tErr         eResult      = DspA_tErr_NoErr;
  tkmDisplayAreaRef this         = NULL;

  /* get us */
  this = *ioppWindow;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* delete frame buffer */
  if ( NULL != this->mpFrameBuffer )
    free( this->mpFrameBuffer );

  /* delete our voxels */
  xVoxl_Delete( &this->mpLastCursor );
  xVoxl_Delete( &this->mpCursor );
  xVoxl_Delete( &this->mpMouseLocationAnaIdx );
  xVoxl_Delete( &this->mpZoomCenter );
  xVoxl_Delete( &this->mpOriginalZoomCenter );

  /* trash the signature */
  this->mSignature = 0x1;

  /* delete us */
  free( this );

  /* return null */
  *ioppWindow = NULL;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_Delete: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetPosition ( tkmDisplayAreaRef this,
                             xPoint2n          iLocation,
                             int               inWidth,
                             int               inHeight ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set our location */
  this->mLocationInSuper = iLocation;

  /* set our size */
  this->mnWidth  = inWidth;
  this->mnHeight = inHeight;

  /* set our scale */
  if ( this->mnVolumeSizeX > 0 ) {

    this->mfFrameBufferScaleX =
      (float)this->mnWidth  / (float)this->mnVolumeSizeX;
    this->mfFrameBufferScaleY =
      (float)this->mnHeight / (float)this->mnVolumeSizeY;
  }

  /* rebuild our current frame and redraw. */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetPosition: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetID ( tkmDisplayAreaRef this, int inID ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the id */
  this->mID = inID;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetPosition: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_UpdateWindowTitle ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult           = DspA_tErr_NoErr;
  char      sTitle[STRLEN]       = "";
  char      sSubjectName[STRLEN] = "";
  char      sVolumeName[STRLEN]  = "";
  char      sAuxVolumeName[STRLEN]  = "";

  DebugEnterFunction( ("DspA_UpdateWindowTitle( this=%p )", this) );

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* get the title information */
  DebugNote( ("Getting main subject name") );
  Volm_CopySubjectName( this->mpVolume[tkm_tVolumeType_Main],
                        sSubjectName, sizeof(sSubjectName) );
  DebugNote( ("Getting main volume name") );
  Volm_CopyVolumeName( this->mpVolume[tkm_tVolumeType_Main],
                       sVolumeName, sizeof(sVolumeName) );
  if ( NULL != this->mpVolume[tkm_tVolumeType_Aux] ) {
    DebugNote( ("Getting aux volume name") );
    Volm_CopyVolumeName( this->mpVolume[tkm_tVolumeType_Aux],
                         sAuxVolumeName, sizeof(sAuxVolumeName) );
  }

  /* if we don't have an aux volume */
  if ( NULL == this->mpVolume[tkm_tVolumeType_Aux] ) {

    /* just use the subject and volume name */
    sprintf( sTitle, "%s: %s", sSubjectName, sVolumeName );

  } else {

    /* else see which one is displayed. use that name first and then
       the other name in parens. */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      sprintf( sTitle, "%s: %s (** %s **)",
               sSubjectName, sVolumeName, sAuxVolumeName );
    } else {
      sprintf( sTitle, "%s: ** %s ** (%s)",
               sSubjectName, sVolumeName, sAuxVolumeName );
    }
  }

  /* set the window name. */
  DebugNote( ("Setting composed window name") );
  MWin_SetWindowTitle( this->mpWindow, sTitle );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_UpdateWindowTitle: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_SetVolume ( tkmDisplayAreaRef this,
                           mriVolumeRef      ipVolume,
                           int               inSizeX,
                           int               inSizeY,
                           int               inSizeZ ) {

  DspA_tErr eResult                        = DspA_tErr_NoErr;
  xVoxelRef pCenter                        = NULL;
  char      sTclArguments[tkm_knTclCmdLen] = "";
  char      sVolumeName[tkm_knNameLen]     = "";
  int         nSize ;

  DebugEnterFunction( ("DspA_SetVolume( this=%p, ipVolume=%p, "
                       "inSize=%d,%d,%d )",
                       this, ipVolume, inSizeX, inSizeY, inSizeZ) );

  xVoxl_New( &pCenter );

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the main volume */
  this->mpVolume[tkm_tVolumeType_Main] = ipVolume;

  /* save the volume size */
  this->mnVolumeSizeX = inSizeX;
  this->mnVolumeSizeY = inSizeY;
  this->mnVolumeSizeZ = inSizeZ;

  /* Get the real volume dimensions. inSize{X,Y,Z} are actually the
     screen or window space dimensions. These real dimensions are used
     for bounds checking.  */
  Volm_GetDimensions( this->mpVolume[tkm_tVolumeType_Main],
                      &this->mnVolumeDimensionX, &this->mnVolumeDimensionY,
                      &this->mnVolumeDimensionZ );


  /* Calculate the min and max volume index in screen coords. This is
     the screen space that corresponds to volume voxels. Coordinates
     outside of thie range are invalid. Used in VerifyVolumeVoxel. Use
     min/max to get a minimum range of 0-256, this is because some
     volumes have weird spacing issues or something. I don't know. */
  this->mnMinVolumeIndexX = MIN( 0,
                                 floor((float)inSizeX/2.0) -
                                 floor((float)this->mnVolumeDimensionX/2.0) );
  this->mnMinVolumeIndexY = MIN( 0,
                                 floor((float)inSizeY/2.0) -
                                 floor((float)this->mnVolumeDimensionY/2.0) );
  this->mnMinVolumeIndexZ = floor((float)inSizeZ/2.0) -
                            floor((float)this->mnVolumeDimensionZ/2.0);

  this->mnMaxVolumeIndexX = MAX( inSizeX,
                                 ceil((float)inSizeX/2.0) +
                                 ceil((float)this->mnVolumeDimensionX/2.0) );
  this->mnMaxVolumeIndexY = MAX( inSizeY,
                                 ceil((float)inSizeY/2.0) +
                                 ceil((float)this->mnVolumeDimensionY/2.0) );
  this->mnMaxVolumeIndexZ = ceil((float)inSizeZ/2.0) +
                            ceil((float)this->mnVolumeDimensionZ/2.0);

  /* This ugly little hack is for single slice images that show up in
     slice 127 but should probably go into slice 128 as the calculated
     bounds suggest. Ah well. If it isn't one of these images, calc
     the min/max as usual. */
  if ( 128 == this->mnMinVolumeIndexZ &&
       129 == this->mnMaxVolumeIndexZ ) {
    this->mnMinVolumeIndexZ = 127;
  } else {
    this->mnMinVolumeIndexZ = MIN( 0, this->mnMinVolumeIndexZ );
    this->mnMaxVolumeIndexZ = MAX( inSizeZ, this->mnMaxVolumeIndexZ );
  }


  /* if we alreayd have a frame buffer, delete it */
  if ( NULL == this->mpFrameBuffer ) {
    free( this->mpFrameBuffer );
    this->mpFrameBuffer = NULL;
  }

  /* Allocate a new one. NOTES: Normally mnVolumeSize{X,Y,Z} will all
     be equal to each other, so this shouldn't be an issue. */
  nSize = MAX( MAX(this->mnVolumeSizeX, this->mnVolumeSizeY),
               this->mnVolumeSizeZ) ;
  this->mpFrameBuffer = 
    (GLubyte*) malloc( nSize * nSize * DspA_knNumBytesPerPixel );
  if ( NULL == this->mpFrameBuffer ) {
    eResult = DspA_tErr_AllocationFailed;
    goto error;
  }

  /* set inital values for our buffer scale */
  this->mfFrameBufferScaleX =
    (float)this->mnWidth  / (float)this->mnVolumeSizeX;
  this->mfFrameBufferScaleY =
    (float)this->mnHeight  / (float)this->mnVolumeSizeY;

  /* initialize surface point lists. */
  eResult = DspA_InitSurfaceLists_( this, nSize * Surf_knNumVertexSets *
                                    mri_knNumOrientations );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* update window title */
  DspA_UpdateWindowTitle( this );

  /* send volume name */
  Volm_CopyVolumeName( this->mpVolume[tkm_tVolumeType_Main],
                       sVolumeName, sizeof(sVolumeName) );
  sprintf( sTclArguments, "\"%s value\"", sVolumeName );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );

  /* get the center of the volume */
  xVoxl_Set( pCenter, floor(this->mnVolumeDimensionX/2),
             floor(this->mnVolumeDimensionY/2),
             floor(this->mnVolumeDimensionZ/2) );
  xVoxl_Set( pCenter, floor(this->mnVolumeSizeX/2),
             floor(this->mnVolumeSizeY/2),
             floor(this->mnVolumeSizeZ/2) );

  /* set cursor and zoom center to middle of volume if not already set */
  if (xVoxl_GetX(this->mpCursor)==0 &&
      xVoxl_GetY(this->mpCursor)==0 &&
      xVoxl_GetZ(this->mpCursor)==0) {
    eResult = DspA_SetCursor( this, pCenter );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    eResult = DspA_SetZoomCenter( this, pCenter );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  } else {
    eResult = DspA_SetCursor( this, this->mpCursor );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }


  /* show cursor. */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Cursor, TRUE );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* show volume. */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Anatomical, TRUE );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set dirty flag and redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetVolume: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pCenter );

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_SetAuxVolume ( tkmDisplayAreaRef this,
                              mriVolumeRef      ipVolume,
                              int               inSizeX,
                              int               inSizeY,
                              int               inSizeZ ) {

  DspA_tErr eResult                        = DspA_tErr_NoErr;
  char      sVolumeName[tkm_knNameLen]     = "";
  char      sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction( ("DspA_SetAuxVolume( this=%p, ipVolume=%p, "
                       "inSizeX=%d, inSizeY=%d, inSizeZ=%d)",
                       this, ipVolume, inSizeX, inSizeY, inSizeZ) );

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make sure this size is the same as the main volume size. */
  if ( inSizeX != this->mnVolumeSizeX ||
       inSizeY != this->mnVolumeSizeY ||
       inSizeZ != this->mnVolumeSizeZ ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* save the aux volume */
  this->mpVolume[tkm_tVolumeType_Aux] = ipVolume;

  /* if we got a volume... */
  if ( NULL != this->mpVolume[tkm_tVolumeType_Aux] ) {

    /* send volume name to tk window */
    Volm_CopyVolumeName( this->mpVolume[tkm_tVolumeType_Aux],
                         sVolumeName, sizeof(sVolumeName) );
    xUtil_snprintf( sTclArguments, sizeof(sTclArguments),
                    "\"%s value\"", sVolumeName );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeName, sTclArguments );

    /* show the volume value */
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxValue, "1" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeOptions, "1" );

    /* if we're currently showing the aux volume, the slice has changed */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      this->mbSliceChanged = TRUE;
    }

    /* if we're focused, send the new information for the cursor */
    if ( sFocusedDisplay == this ) {
      DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
                                       this->mpCursor );
    }
  } else {

    /* hide the volume value */
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxValue, "0" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeOptions, "0" );

    /* don't show the aux volume */
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_AuxVolume, FALSE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* update window title */
  DspA_UpdateWindowTitle( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetAuxVolume: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_SetSegmentationVolume ( tkmDisplayAreaRef this,
                                       tkm_tSegType      iType,
                                       mriVolumeRef      iVolume ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  tBoolean  bHaveVolume         = FALSE;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  DebugAssertThrowX( (iType >= 0 && iType < tkm_knNumSegTypes),
                     eResult, DspA_tErr_InvalidParameter );

  /* save the group */
  this->mSegmentationVolume[iType] = iVolume;

  /* turn stuff on or off based on if we have one. */
  if ( this->mSegmentationVolume[iType] != NULL ) {
    bHaveVolume = TRUE;
  }

  /* show the appropriate seg label and highlight the appropriate
     volume. */
  sprintf( sTclArguments, "%d", (int)bHaveVolume );
  switch ( iType ) {
  case tkm_tSegType_Main:
    tkm_SendTclCommand( tkm_tTclCommand_ShowSegLabel, sTclArguments );
    if ( bHaveVolume ) {
      DspA_SetDisplayFlag( this,
                           DspA_tDisplayFlag_SegmentationVolumeOverlay, TRUE );
      DspA_SetDisplayFlag( this,
                           DspA_tDisplayFlag_AuxSegmentationVolume, FALSE );
    }
    break;
  case tkm_tSegType_Aux:
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxSegLabel, sTclArguments );
    if ( bHaveVolume ) {
      DspA_SetDisplayFlag( this,
                           DspA_tDisplayFlag_SegmentationVolumeOverlay, TRUE );
      DspA_SetDisplayFlag( this, DspA_tDisplayFlag_AuxSegmentationVolume,
                           TRUE );
    }
    break;
  default:
    break;
  }
  tkm_SendTclCommand( tkm_tTclCommand_ShowSegmentationOptions,
                      sTclArguments );

  /* if we're focused, send the new information for the cursor */
  if ( sFocusedDisplay == this ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
                                     this->mpCursor );
  }

  /* redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSegmentationVolume: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSegmentationColorTable ( tkmDisplayAreaRef this,
    tkm_tSegType      iType,
    COLOR_TABLE*      iCTAB ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  DebugAssertThrowX( (iType >= 0 && iType < tkm_knNumSegTypes),
                     eResult, DspA_tErr_InvalidParameter );

  /* save the table */
  this->mSegmentationColorTable[iType] = iCTAB;

  /* if we have a segmentation displayed, redraw */
  if ( NULL != this->mSegmentationVolume[iType] &&
       this->mabDisplayFlags[DspA_tDisplayFlag_SegmentationVolumeOverlay] ) {
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSegmentationColorTable: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSurface ( tkmDisplayAreaRef this,
                            tkm_tSurfaceType  iType,
                            mriSurfaceRef     ipSurface ) {

  DspA_tErr        eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the surface */
  this->mpSurface[iType] = ipSurface;

  /* purge all the surface lists. */
  DspA_PurgeSurfaceLists_( this );

  /* set slice dirty flag and redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* set surface to null */
  this->mpSurface[iType] = NULL;

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSurface: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetOverlayVolume ( tkmDisplayAreaRef      this,
                                  tkmFunctionalVolumeRef ipVolume) {

  DspA_tErr eResult        = DspA_tErr_NoErr;
  FunV_tErr eFunctional    = FunV_tErr_NoError;
  tBoolean  bOverlayLoaded = FALSE;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the volume */
  this->mpFunctionalVolume = ipVolume;

  /* see if there's data... */
  eFunctional = FunV_IsOverlayPresent( this->mpFunctionalVolume,
                                       &bOverlayLoaded );
  if ( FunV_tErr_NoError != eFunctional ) {
    DspA_Signal( "DspA_SetOverlayVolume", __LINE__,
                 DspA_tErr_ErrorAccessingFunctionalVolume );
  }

  /* turn functional data and color scale bar on if there is */
  if ( bOverlayLoaded ) {
    eResult =
      DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay, TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    eResult =
      DspA_SetDisplayFlag(this,DspA_tDisplayFlag_FunctionalColorScaleBar,TRUE);
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* if we're focused, send the new information for the cursor */
  if ( sFocusedDisplay == this ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
                                     this->mpCursor );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetOverlayVolume: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetControlPointsSpace ( tkmDisplayAreaRef this,
                                       x3DListRef        ipVoxels ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the control points space. */
  this->mpControlPoints = ipVoxels;

  /* turn control points on. */
  if ( NULL != this->mpControlPoints ) {

    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_ControlPoints,
                                   TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetControlPointsSpace: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSelectionSpace( tkmDisplayAreaRef this,
                                  mriVolumeRef      ipVolume ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the selected voxels */
  this->mpSelection = ipVolume;

  /* turn selection display on. */
  if ( NULL != this->mpSelection ) {

    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Selection, TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSelectionSpace: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetHeadPointList ( tkmDisplayAreaRef   this,
                                  mriHeadPointListRef iList ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  tBoolean  bHaveList          = FALSE;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the selected voxels */
  this->mHeadPoints = iList;

  /* turn selection display on. */
  if ( NULL != this->mHeadPoints ) {

    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_HeadPoints, TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    bHaveList = TRUE;
  }

  sprintf( sTclArguments, "%d", (int)bHaveList );
  tkm_SendTclCommand( tkm_tTclCommand_ShowHeadPointLabel, sTclArguments );

  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetHeadPointList: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetVLIs  ( tkmDisplayAreaRef this,
                          VLI*              iVLI1,
                          VLI*              iVLI2,
                          char*             isVLI1_name,
                          char*             isVLI2_name) {

  DspA_tErr eResult            = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the VLIs and TRANSFORM */
  this->mVLI1    = iVLI1;
  this->mVLI2    = iVLI2;
  strcpy(this->isVLI1_name, isVLI1_name) ;
  strcpy(this->isVLI2_name, isVLI2_name) ;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetVLIs: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_SetDTIVolume  ( tkmDisplayAreaRef this,
                               mriVolumeRef        iVolume ) {


  DspA_tErr eResult             = DspA_tErr_NoErr;
  char       sTclArguments[tkm_knTclCmdLen] = "";

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Save the volume. */
  this->mpDTIVolume = iVolume;

  /* Show the DTI options if we got a volume.  */
  sprintf( sTclArguments, "%d", (int)(NULL != iVolume) );
  tkm_SendTclCommand( tkm_tTclCommand_ShowDTIOptions, sTclArguments );

  /* turn DTI display on. */
  if ( NULL != this->mpDTIVolume ) {
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_DTIOverlay, TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetDTIVolume: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_SetGCA ( tkmDisplayAreaRef   this,
                        GCA*                iVolume,
                        TRANSFORM*          iTransform ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the GCA and LTA */
  this->mGCAVolume    = iVolume;
  this->mGCATransform = iTransform;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetGCA: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetCursor ( tkmDisplayAreaRef this,
                           xVoxelRef          ipCursor ) {

  DspA_tErr             eResult     = DspA_tErr_NoErr;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  HPtL_tHeadPointRef    pHeadPoint  = NULL;
  int                   nSlice      = 0;
  xPoint2f              planePt;
  xVoxel                MRIIdx;

  planePt.mfX=0;
  planePt.mfY=0;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the cursor */
  eResult = DspA_VerifyVolumeVoxel_( this, ipCursor );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* allow the functional display to respond. */
  if ( NULL != this->mpFunctionalVolume ) {
    Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                             ipCursor, &MRIIdx );
    eFunctional = FunV_MRIIdxClicked( this->mpFunctionalVolume, &MRIIdx );
  }

  /* Allow a point to be selected in the GDF volume. */
  tkm_SelectGDFMRIIdx( &MRIIdx );

  /* get our current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* Copy the current cursor into the last cursor. */
  xVoxl_Copy( this->mpLastCursor, this->mpCursor );

  /* set the cursor */
  xVoxl_Copy( this->mpCursor, ipCursor );

  /* if the new slice number is diffrent, set our dirty slice flag. */
  if ( DspA_GetCurrentSliceNumber_( this ) != nSlice ) {
    this->mbSliceChanged = TRUE;
  }

  // WRONG ===========================================================
  // the convention is that the pixel center is (*.0, *.0) and thus
  // no need to add.
  /* if cursor is .0 .0, change to .5 .5 so that it will draw in the
     center of a voxel on screen */
  DspA_ConvertVolumeToPlane_( this, this->mpCursor, this->mOrientation,
                              &planePt, &nSlice );
  if ( planePt.mfX == (float)(int)planePt.mfX &&
       planePt.mfY == (float)(int)planePt.mfY ) {
    //  planePt.mfX += 0.5;
    //  planePt.mfY += 0.5;
    DspA_ConvertPlaneToVolume_( this, &planePt, nSlice,
                                this->mOrientation, this->mpCursor );
  }

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* notify the window that the cursor has changed. */
    MWin_CursorChanged( this->mpWindow, this, this->mpCursor );

    /* send the information for this point */
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
                                     this->mpCursor );

  }

  /* if we have head point data... */
  if (  NULL != this->mHeadPoints ) {

    /* find a head point to select. look in flattened space if we're
       looking at the max int proj.*/
    tkm_GetHeadPoint( this->mpCursor, this->mOrientation,
                      this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj],
                      &pHeadPoint );

    this->mpSelectedHeadPoint = pHeadPoint;
  }

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetCursor(%d,%d,%d): %s\n",
                 eResult, xVoxl_ExpandInt(ipCursor),
                 DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertAndSetCursor ( tkmDisplayAreaRef this,
                                     mri_tCoordSpace   iFromSpace,
                                     xVoxelRef         ipCoord ) {

  DspA_tErr  eResult     = DspA_tErr_NoErr;
  Volm_tErr  eVolume     = Volm_tErr_NoErr;
  xVoxel     MRIIdx;
  xVoxel     anaIdx;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* convert the coord to the right space. */
  switch ( iFromSpace ) {
  case mri_tCoordSpace_VolumeIdx:
    // src may not be (256,256,256) so that we convert into "normalized" coords
    eVolume =
      Volm_ConvertMRIIdxToIdx(this->mpVolume[tkm_tVolumeType_Main],
                              ipCoord, &anaIdx);
    break;
  case mri_tCoordSpace_SurfaceRAS:
    if ( tkm_UseRealRAS() ) {
      eVolume = Volm_ConvertRASToIdx( this->mpVolume[tkm_tVolumeType_Main],
                                      ipCoord, &anaIdx );
    } else {
      eVolume =
        Volm_ConvertSurfaceRASToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                                        ipCoord, &MRIIdx );
      eVolume =
        Volm_ConvertMRIIdxToIdx( this->mpVolume[tkm_tVolumeType_Main],
                                 &MRIIdx, &anaIdx );
    }
    break;
  case mri_tCoordSpace_RAS:
    eVolume = Volm_ConvertRASToIdx( this->mpVolume[tkm_tVolumeType_Main],
                                    ipCoord, &anaIdx );
    break;
  case mri_tCoordSpace_Talairach:
    eVolume = Volm_ConvertTalToIdx( this->mpVolume[tkm_tVolumeType_Main],
                                    ipCoord, &anaIdx );
    break;
  default:
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  if ( Volm_tErr_NoErr != eVolume )
    goto error;
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the cursor */
  eResult = DspA_SetCursor( this, &anaIdx );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_ConvertAndSetCursor(%d, %d,%d,%d): %s\n",
                 eResult, (int)iFromSpace, xVoxl_ExpandInt(ipCoord),
                 DspA_GetErrorString(eResult) ) );
  }
  if ( Volm_tErr_NoErr != eVolume ) {
    DebugPrint( ("Error %d in DspA_ConvertAndSetCursor(%d, %d,%d,%d): %s\n",
                 eVolume, (int)iFromSpace, xVoxl_ExpandInt(ipCoord),
                 Volm_GetErrorString(eVolume) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSlice ( tkmDisplayAreaRef this,
                          int               inSlice ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xVoxelRef  pCursor = NULL;

  xVoxl_New( &pCursor );

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* ignore if this is already our slice number */
  if ( inSlice == DspA_GetCurrentSliceNumber_( this ) )
    goto cleanup;

  /* copy the cursor. */
  xVoxl_Copy( pCursor, this->mpCursor );

  /* change the slice */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_SetZ( pCursor, inSlice );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_SetY( pCursor, inSlice );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetX( pCursor, inSlice );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
  }

  /* set the cursor. */
  eResult = DspA_SetCursor( this, pCursor );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* update mouse over information if we're focused */
  if ( sFocusedDisplay == this ) {
    eResult = DspA_SendMouseInfoToTcl( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSlice: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pCursor );

  return eResult;
}

DspA_tErr DspA_SetOrientation ( tkmDisplayAreaRef this,
                                mri_tOrientation  iOrientation ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  int       nSlice             = 0;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the orientation */
  if ( iOrientation <= mri_tOrientation_None
       || iOrientation >= mri_knNumOrientations ) {
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
  }

  /* if the orientations are different, our slice will have changed. */
  if ( this->mOrientation != iOrientation ) {
    this->mbSliceChanged = TRUE;
  }

  /* set the orientation */
  this->mOrientation = iOrientation;

  /* when zoomed, set the zoom center to the cursor so it appears to
     reorient around the cursor (zoom centers are only 2d and not valid
     when switching planes. */
  if ( this->mnZoomLevel != DspA_knMinZoomLevel ) {
    eResult = DspA_SetZoomCenterToCursor( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* notify of change */
    MWin_OrientationChanged( this->mpWindow, this, iOrientation );

    /* send the orientation */
    sprintf( sTclArguments, "%d", (int)this->mOrientation );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateOrientation, sTclArguments );

    /* send the new slice number */
    nSlice = DspA_GetCurrentSliceNumber_( this );
    sprintf( sTclArguments, "%d", nSlice );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeSlice, sTclArguments );
  }

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetOrientation: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetZoomLevel ( tkmDisplayAreaRef this,
                              int               inLevel ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  int       nNewLevel          = 0;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the zoom level. */
  nNewLevel = this->mnZoomLevel;
  if ( inLevel >= DspA_knMinZoomLevel
       && inLevel <= DspA_knMaxZoomLevel ) {
    nNewLevel = inLevel;
  }

  DspA_SetZoomCenterToCursor( this );

  /* set our zoom level. */
  this->mnZoomLevel = nNewLevel;

  /* this requires a rebuild of the frame buffer. */
  this->mbSliceChanged = TRUE;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* notify the window that the zoom level has changed. */
    MWin_ZoomLevelChanged( this->mpWindow, this, this->mnZoomLevel );

    /* send zoom level update. */
    sprintf( sTclArguments, "%d", (int)this->mnZoomLevel );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateZoomLevel, sTclArguments );
  }

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetZoomLevel: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetZoomCenter ( tkmDisplayAreaRef this,
                               xVoxelRef          ipCenter ) {

  DspA_tErr eResult  = DspA_tErr_NoErr;
  int       nX       = 0;
  int       nY       = 0;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the center */
  eResult = DspA_VerifyVolumeVoxel_( this, ipCenter );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    nX = xVoxl_GetX( ipCenter );
    nY = xVoxl_GetY( ipCenter );
    break;
  case mri_tOrientation_Horizontal:
    nX = xVoxl_GetX( ipCenter );
    nY = xVoxl_GetZ( ipCenter );
    break;
  case mri_tOrientation_Sagittal:
    nX = xVoxl_GetZ( ipCenter );
    nY = xVoxl_GetY( ipCenter );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* set the center */
  xVoxl_SetX( this->mpZoomCenter, nX );
  xVoxl_SetY( this->mpZoomCenter, nY );

  /* should rebuild frame */
  this->mbSliceChanged = TRUE;

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetZoomCenter: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetZoomCenterToCursor ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the center to our cursor */
  DspA_SetZoomCenter( this, this->mpCursor );

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetZoomCenterToCursor: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_HiliteSurfaceVertex ( tkmDisplayAreaRef this,
                                     Surf_tVertexSet  inSurface,
                                     int               inVertex ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the values */
  this->mHilitedSurface        = inSurface;
  this->mnHilitedVertexIndex   = inVertex;

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HiliteSurfaceVertex: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_SetDisplayFlag ( tkmDisplayAreaRef this,
                                DspA_tDisplayFlag iWhichFlag,
                                tBoolean          ibNewValue ) {

  DspA_tErr    eResult               = DspA_tErr_NoErr;
  char         sTclArguments[STRLEN] = "";
  tBoolean     bNewValue             = FALSE;
  tkm_tSegType segType               = tkm_tSegType_Main;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* check the flag */
  if ( iWhichFlag <= DspA_tDisplayFlag_None
       || iWhichFlag >= DspA_knNumDisplayFlags ) {
    eResult = DspA_tErr_InvalidDisplayFlag;
    goto error;
  }

  /* get the new value. */
  bNewValue = ibNewValue;

  /* verify the flag. */
  switch ( iWhichFlag ) {

  case DspA_tDisplayFlag_AuxVolume:

    /* if no aux volume, set to false. */
    if ( NULL == this->mpVolume[tkm_tVolumeType_Aux] )
      bNewValue = FALSE;

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    /* change the window title. */
    DspA_UpdateWindowTitle( this );

    break;

  case DspA_tDisplayFlag_Anatomical:

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_AuxSegmentationVolume:

    /* if no aux seg volume, set to false. */
    if ( NULL == this->mSegmentationVolume[tkm_tSegType_Aux] )
      bNewValue = FALSE;

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_SegmentationVolumeOverlay:
  case DspA_tDisplayFlag_SegLabelVolumeCount:

    if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume] ) {
      segType = tkm_tSegType_Aux;
    } else {
      segType = tkm_tSegType_Main;
    }

    /* if no segmentation, set to false. */
    if ( NULL == this->mSegmentationVolume[segType] )
      bNewValue = FALSE;

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_MainSurface:
  case DspA_tDisplayFlag_OriginalSurface:
  case DspA_tDisplayFlag_PialSurface:
  case DspA_tDisplayFlag_DisplaySurfaceVertices:

    /* if no surface, set to false. */
    if ( NULL == this->mpSurface[tkm_tSurfaceType_Main] )
      bNewValue = FALSE;

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_InterpolateSurfaceVertices:

    /* if no surface, set to false. */
    if ( NULL == this->mpSurface[tkm_tSurfaceType_Main] )
      bNewValue = FALSE;

    /* if the flag is different, purge lists and set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue ) {
      this->mbSliceChanged = TRUE;
      eResult = DspA_PurgeSurfaceLists_( this );
      if ( DspA_tErr_NoErr != eResult )
        goto error;
    }

    break;

  case DspA_tDisplayFlag_FunctionalOverlay:
  case DspA_tDisplayFlag_FunctionalColorScaleBar:
  case DspA_tDisplayFlag_MaskToFunctionalOverlay:

    /* if no func data, set to false. */
    if ( NULL == this->mpFunctionalVolume ) {
      bNewValue = FALSE;
    }

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_MaskFunctionalOverlayToAux:

    /* if no func data or aux data, set to false. */
    if ( NULL == this->mpFunctionalVolume ||
         NULL == this->mpVolume[tkm_tVolumeType_Aux] ) {
      bNewValue = FALSE;
    }

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;
  case DspA_tDisplayFlag_HistogramPercentChange:

    /* if no VLI data, set to false. */
    if ( NULL == this->mVLI1 ) {
      bNewValue = FALSE;
    }

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_Selection:
  case DspA_tDisplayFlag_MaxIntProj:
  case DspA_tDisplayFlag_UndoVolume:

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_HeadPoints:

    /* check to see if we have a list */
    if ( NULL == this->mHeadPoints )
      bNewValue = FALSE;

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_DTIOverlay:

    /* if no DTI data, set to false. */
    if ( NULL == this->mpDTIVolume ) {
      bNewValue = FALSE;
    }

    /* if the flag is different, set dirty flag */
    if ( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  case DspA_tDisplayFlag_VerboseGCADump:
    /* Just set it, don't need to redraw as it only affects text
       output. */
    break;

  default:
    break;
  }

  /* set the value */
  this->mabDisplayFlags[iWhichFlag] = bNewValue;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* notify the window that the flag has changed. */
    MWin_DisplayFlagChanged( this->mpWindow, this, iWhichFlag,
                             this->mabDisplayFlags[iWhichFlag] );

    /* send the tcl update. */
    sprintf( sTclArguments, "%d %d", (int)iWhichFlag, (int)bNewValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateDisplayFlag, sTclArguments );
  }

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetDisplayFlag %d %d: %s\n",
                 eResult, iWhichFlag, ibNewValue,
                 DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ToggleDisplayFlag ( tkmDisplayAreaRef this,
                                   DspA_tDisplayFlag iWhichFlag ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* check the flag */
  if ( iWhichFlag <= DspA_tDisplayFlag_None
       || iWhichFlag >= DspA_knNumDisplayFlags ) {
    eResult = DspA_tErr_InvalidDisplayFlag;
    goto error;
  }

  /* set the flag to the opposite of what it is now */
  eResult = DspA_SetDisplayFlag( this, iWhichFlag,
                                 !(this->mabDisplayFlags[iWhichFlag]) );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_ToggleDisplayFlag: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetTool ( tkmDisplayAreaRef this,
                         DspA_tTool        iTool ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Set the tool */
  sTool = iTool;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d", (int)sTool );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateTool, sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetTool: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetBrushTarget ( tkmDisplayAreaRef this,
                                DspA_tBrushTarget iTarget ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* check the input before setting the info */
  if ( iTarget > DspA_tBrushTarget_None
       && iTarget <= DspA_knNumBrushTargets ) {
    sBrush.mTarget = iTarget;
  }

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d", (int)sBrush.mTarget );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushTarget, sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetBrushTarget: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetBrushShape ( tkmDisplayAreaRef this,
                               int               inRadius,
                               DspA_tBrushShape  iShape,
                               tBoolean          ib3D ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* check the input before setting the info */
  if ( inRadius >= DspA_knMinBrushRadius
       && inRadius <= DspA_knMaxBrushRadius ) {
    sBrush.mnRadius = inRadius;
  }

  if ( iShape >= 0
       && iShape < DspA_knNumBrushShapes ) {
    sBrush.mShape = iShape;
  }

  /* Set the brush info */
  sBrush.mb3D = ib3D;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d",
              (int)sBrush.mnRadius, (int)sBrush.mShape, (int)sBrush.mb3D );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushShape, sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetBrushShape: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_SetAnatomicalFillInfo ( tkmDisplayAreaRef this,
                                       tBoolean          ib3DFill,
                                       float             iFuzzy,
                                       float             iDistance ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Set the brush info */
  sBrush.mb3DFill   = ib3DFill;
  sBrush.mFuzzy     = iFuzzy;
  sBrush.mDistance  = iDistance;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %f %f",
              (int)sBrush.mb3DFill, sBrush.mFuzzy,
              sBrush.mDistance );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAnatomicalFillInfo,
                        sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetBrushShape: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetBrushInfo ( tkmDisplayAreaRef this,
                              DspA_tBrush       iBrush,
                              DspA_tBrushInfoRef iInfo ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make sure the brush is in bounds */
  if ( iBrush < 0 ||
       iBrush >= DspA_knNumBrushes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set the brush theshold info */
  sBrush.mInfo[ iBrush ] = *iInfo;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %f %f %f %d %d",
              (int)iBrush,
              sBrush.mInfo[iBrush].mLow,
              sBrush.mInfo[iBrush].mHigh,
              sBrush.mInfo[iBrush].mNewValue,
              sBrush.mInfo[iBrush].mMode,
              sBrush.mInfo[iBrush].mCloneSource );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushInfo, sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetBrushInfo: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetCursorColor ( tkmDisplayAreaRef this,
                                xColor3fRef       iColor ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  if ( NULL == iColor ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set cursor color and notify of change */
  sCursorColor = *iColor;
  DspA_Redraw_( this );

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%f %f %f",
              sCursorColor.mfRed, sCursorColor.mfGreen, sCursorColor.mfBlue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateCursorColor, sTclArguments );
  }


  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetCursorColor: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSurfaceLineWidth ( tkmDisplayAreaRef this,
                                     tkm_tSurfaceType  iSurface,
                                     Surf_tVertexSet   iVertexSet,
                                     int               inWidth ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  if ( iSurface < 0 ||
       iSurface >= tkm_knNumSurfaceTypes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  if ( iVertexSet <= Surf_tVertexSet_None ||
       iVertexSet >= Surf_knNumVertexSets ||
       inWidth < 0 ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set surface width and notify of change */
  this->manSurfaceLineWidth[iSurface][iVertexSet] = inWidth;
  DspA_Redraw_( this );

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d", (int)iSurface, (int)iVertexSet,
              this->manSurfaceLineWidth[iSurface][iVertexSet] );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceLineWidth,
                        sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSurfaceLineWidth: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSurfaceLineColor ( tkmDisplayAreaRef this,
                                     tkm_tSurfaceType  iSurface,
                                     Surf_tVertexSet   iVertexSet,
                                     xColor3fRef       iColor ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  if ( iSurface < 0 ||
       iSurface >= tkm_knNumSurfaceTypes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  if ( iVertexSet <= Surf_tVertexSet_None ||
       iVertexSet >= Surf_knNumVertexSets ||
       NULL == iColor ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set surface color and notify of change */
  this->maSurfaceLineColor[iSurface][iVertexSet] = *iColor;
  DspA_Redraw_( this );

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf 
      ( sTclArguments, "%d %d %f %f %f", (int)iSurface, (int)iVertexSet,
        xColr_ExpandFloat( &(this->maSurfaceLineColor[iSurface][iVertexSet])));
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceLineColor,
                        sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSurfaceLineColor: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetFloodSelectParams ( tkmDisplayAreaRef          this,
                                      DspA_tFloodSelectSettings* iSettings ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the brush theshold info */
  sFloodSelectSettings = *iSettings;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %f %f",
              (int)sFloodSelectSettings.mb3D,
              (int)sFloodSelectSettings.mSrc,
              sFloodSelectSettings.mFuzzy,
              sFloodSelectSettings.mDistance );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateFloodSelectParams,
                        sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetFloodSelectParams: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetCursorShape ( tkmDisplayAreaRef this,
                                DspA_tMarker      iShape ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  if ( iShape <= DspA_tMarker_None ||
       iShape >= DspA_knNumMarkers ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set cursor shape and notify of change */
  sCursorShape = iShape;
  DspA_Redraw_( this );

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d", sCursorShape );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateCursorShape, sTclArguments );
  }


  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetCursorShape: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetSegBrushInfo ( tkmDisplayAreaRef        this,
                                 DspA_tSegBrushSettings* iSettings ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  if ( iSettings->mSrc < tkm_tVolumeType_Main ||
       iSettings->mSrc >= tkm_knNumVolumeTargets ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set brush data */
  sSegBrush = *iSettings;

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf( sTclArguments, "%d %d %d %f %f",
             sSegBrush.mnPaintValue, sSegBrush.mb3D,
             sSegBrush.mSrc, sSegBrush.mFuzzy, sSegBrush.mDistance );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSegBrushInfo, sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSegBrushInfo: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ChangeSliceBy ( tkmDisplayAreaRef this,
                               int               inDelta ) {

  DspA_tErr  eResult = DspA_tErr_NoErr;
  xVoxelRef  pCursor = NULL;
  float      nSlice  = 0;

  xVoxl_New( &pCursor );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* apply the increment to the proper part of the cursor. wrap around
     if necessary. */
  nSlice = DspA_GetCurrentSliceNumber_( this );
  nSlice += inDelta;
  while ( nSlice < 0 )
    nSlice += this->mnVolumeSizeZ;
  while ( nSlice >= this->mnVolumeSizeZ )
    nSlice -= this->mnVolumeSizeZ;

  /* set the slice */
  eResult = DspA_SetSlice( this, nSlice );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* notify window of change if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {
    MWin_SliceChanged( this->mpWindow, this, inDelta );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_ChangeSliceBy: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pCursor );

  return eResult;
}

DspA_tErr DspA_SetSegmentationAlpha ( tkmDisplayAreaRef this,
                                      float             ifAlpha ) {

  DspA_tErr  eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Set the alpha. */
  if ( this->mfSegmentationAlpha != ifAlpha ) {

    this->mfSegmentationAlpha = ifAlpha;

    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSegmentationAlpha: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetDTIAlpha ( tkmDisplayAreaRef this,
                             float             ifAlpha ) {

  DspA_tErr  eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Set the alpha. */
  if ( this->mfDTIAlpha != ifAlpha ) {

    this->mfDTIAlpha = ifAlpha;

    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetDTIAlpha: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetDTIAxisForComponent ( tkmDisplayAreaRef this,
                                        tkm_tAxis         iAxis,
                                        xColr_tComponent  iComponent ) {

  DspA_tErr  eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Set the axis. */
  if ( this->maDTIAxisForComponent[iComponent] != iAxis ) {

    this->maDTIAxisForComponent[iComponent] = iAxis;

    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetDTIAxisForComponent: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetFuncOverlayAlpha ( tkmDisplayAreaRef this,
                                     float             ifAlpha ) {

  DspA_tErr  eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Set the alpha. */
  if ( this->mfFuncOverlayAlpha != ifAlpha ) {

    this->mfFuncOverlayAlpha = ifAlpha;

    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetFuncOverlayAlpha: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_Focus ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* if we're already focused, return. */
  if ( this == sFocusedDisplay )
    goto cleanup;

  /* blur the last focused display if there was one. */
  if ( sFocusedDisplay ) {
    eResult = DspA_Blur ( sFocusedDisplay );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* focus on us. */
  sFocusedDisplay = this;

  /* send our info to the tk window. */
  DspA_SendViewStateToTcl_( this );

  /* update our window title. */
  DspA_UpdateWindowTitle( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_Focus: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_Blur ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_Blur: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_GetPosition ( tkmDisplayAreaRef this,
                             xPoint2nRef       opLocation,
                             int*              onWidth,
                             int*              onHeight ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* return the info */
  *opLocation = this->mLocationInSuper;
  *onWidth     = this->mnWidth;
  *onHeight    = this->mnHeight;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetPosition: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;

}

DspA_tErr DspA_GetSlice ( tkmDisplayAreaRef this,
                          int*              onSlice ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xVoxelRef  pCursor = NULL;

  xVoxl_New( &pCursor );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* return current slice */
  *onSlice = DspA_GetCurrentSliceNumber_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetSlice: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pCursor );

  return eResult;
}

DspA_tErr DspA_SetCursorToCenterOfSelectionVolume ( tkmDisplayAreaRef this ) {

  DspA_tErr  eResult     = DspA_tErr_NoErr;
  tBoolean   bGotLabel   = FALSE;
  int        nDimensionX = 0;
  int        nDimensionY = 0;
  int        nDimensionZ = 0;
  float      value       = 0;
  xVoxel     MRIIdx;
  xVoxel     anaIdx;
  int        nMinX       = 0;
  int        nMaxX       = 0;
  int        nMinY       = 0;
  int        nMaxY       = 0;
  int        nMinZ       = 0;
  int        nMaxZ       = 0;

  DebugEnterFunction( ("DspA_SetCursorToCenterOfSelectionVolume( this=%p )",
                       this) );

  DebugNote( ("Verifying self") );
  eResult = DspA_Verify( this );
  DebugAssertThrow( (eResult == DspA_tErr_NoErr) );

  /* Go through the selection volume and get the bounds of the
     selection. */
  bGotLabel = FALSE;
  nMinX = nMinY = nMinZ = 999999;
  nMaxX = nMaxY = nMaxZ = 0;
  Volm_GetDimensions( this->mpSelection, &nDimensionX,
                      &nDimensionY, &nDimensionZ );
  xVoxl_Set( &MRIIdx, 0, 0, 0 );
  while ( xVoxl_IncrementUntilLimits( &MRIIdx, nDimensionX-1,
                                      nDimensionY-1, nDimensionZ-1 )) {

    Volm_GetValueAtMRIIdx_( this->mpSelection, &MRIIdx, &value );
    if ( 0 != value ) {

      bGotLabel = TRUE;
      if ( xVoxl_GetZ(&MRIIdx) > nMaxZ ) nMaxZ = xVoxl_GetZ(&MRIIdx);
      if ( xVoxl_GetY(&MRIIdx) > nMaxY ) nMaxY = xVoxl_GetY(&MRIIdx);
      if ( xVoxl_GetX(&MRIIdx) > nMaxX ) nMaxX = xVoxl_GetX(&MRIIdx);
      if ( xVoxl_GetZ(&MRIIdx) < nMinZ ) nMinZ = xVoxl_GetZ(&MRIIdx);
      if ( xVoxl_GetY(&MRIIdx) < nMinY ) nMinY = xVoxl_GetY(&MRIIdx);
      if ( xVoxl_GetX(&MRIIdx) < nMinX ) nMinX = xVoxl_GetX(&MRIIdx);
    }
  }

  if ( bGotLabel ) {

    /* Set the cursor to the center of those bounds. */
    xVoxl_Set( &MRIIdx,
               nMinX + ((nMaxX - nMinX) / 2),
               nMinY + ((nMaxY - nMinY) / 2),
               nMinZ + ((nMaxZ - nMinZ) / 2) );
    Volm_ConvertMRIIdxToIdx( this->mpSelection, &MRIIdx, &anaIdx );
    DspA_SetCursor( this, &anaIdx );
  }

  DebugCatch;
  DebugCatchError( eResult, DspA_tErr_NoErr, DspA_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_HandleEvent ( tkmDisplayAreaRef this,
                             xGWin_tEventRef   ipEvent ) {

  DspA_tErr   eResult     = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  switch ( ipEvent->mType ) {

  case xGWin_tEventType_Draw:

    eResult = DspA_HandleDraw_ ( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    break;

  case xGWin_tEventType_MouseUp:

    eResult = DspA_HandleMouseUp_( this, ipEvent );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    break;

  case xGWin_tEventType_MouseDown:

    eResult = DspA_HandleMouseDown_( this, ipEvent );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    break;

  case xGWin_tEventType_MouseMoved:

    /* This function will return an error when the mouse point is out
       of bounds, but handles it correctly, so we don't need to pass
       it up. Clear the error code if so. */
    eResult = DspA_HandleMouseMoved_( this, ipEvent );
    if ( DspA_tErr_NoErr != eResult &&
         DspA_tErr_InvalidVolumeVoxel != eResult ) {
      goto error;
    } else {
      eResult = DspA_tErr_NoErr;
    }

    break;

  case xGWin_tEventType_KeyDown:

    eResult = DspA_HandleKeyDown_( this, ipEvent );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    break;

  default:
    break;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleEvent: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_HandleMouseUp_ ( tkmDisplayAreaRef this,
                                xGWin_tEventRef   ipEvent ) {

  DspA_tErr       eResult       = DspA_tErr_NoErr;
  xPoint2n        bufferPt      = {0,0};
  xVoxelRef       pVolumeVox    = NULL;
  xVoxel          MRIIdx;
  int             nVolValue     = 0;
  tkm_tVolumeType volumeType = tkm_tVolumeType_Main;
  int             nSegIndex     = 0;
  tkm_tSegType    segType       = tkm_tSegType_Main;
  DspA_tSegBrushSettings segBrush;
  tBoolean        bSelect       = FALSE;
  xVoxel          lineVox;
  xVoxel          lineVox1;
  xVoxel          lineVox2;
  xVoxel          lineIdx1;
  xVoxel          lineIdx2;
  xVoxel          lineRAS1;
  xVoxel          lineRAS2;

  lineVox.mfX=0;
  lineVox.mfY=0;

  xVoxl_New( &pVolumeVox );

  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  switch ( sTool ) {

  case DspA_tTool_Navigate:

    if ( !ipEvent->mbCtrlKey ) {

      /* if there was little delta... */
      if ( this->mTotalDelta.mfX > -1.0 && this->mTotalDelta.mfX < 1.0 &&
           this->mTotalDelta.mfY > -1.0 && this->mTotalDelta.mfY < 1.0 ) {

        /* do a single slice up or down or zoom in or out depending on what
           vertical side of the screen we clicked on */
        switch ( ipEvent->mButton ) {
        case 2:
          if ( GLDRAW_Y_FLIP(bufferPt.mnY) > 128 ) {
            DspA_SetSlice( this,this->mnOriginalSlice - 1 );
          } else {
            DspA_SetSlice( this,this->mnOriginalSlice + 1 );
          }
          break;
        case 3:
          if ( GLDRAW_Y_FLIP(bufferPt.mnY) > 128 ) {
            DspA_SetZoomLevel( this,this->mnOriginalZoomLevel - 1 );
          } else {
            DspA_SetZoomLevel( this,this->mnOriginalZoomLevel + 1 );
          }
          break;
        }
      }
    }
    break;


  case DspA_tTool_SelectVoxels:

    /* If select tool and shift key, do a flood select. If button 2,
       select, and if button 3, unselect. */
    if ( ipEvent->mbShiftKey &&
         (2 == ipEvent->mButton || 3 == ipEvent->mButton) ) {
      if ( 2 == ipEvent->mButton ) {
        bSelect = TRUE;
      } else if ( 3 == ipEvent->mButton ) {
        bSelect = FALSE;
      }

      tkm_FloodSelect( pVolumeVox,
                       sFloodSelectSettings.mb3D,
                       sFloodSelectSettings.mSrc,
                       sFloodSelectSettings.mFuzzy,
                       sFloodSelectSettings.mDistance,
                       bSelect );
      this->mbSliceChanged = TRUE;
    }
    break;


  case DspA_tTool_EditVoxels:

#if 0
    /* If Shift-edit, restore-undo around the clicked point. Don't do
       this now as it was really confusing people.*/
    if ( ipEvent->mbShiftKey ) {
      tkm_RestoreUndoVolumeAroundMRIIdx( pVolumeVox );
    }
#endif

    /* Shift-2 or Shift-3 does a fill. */
    if ( ipEvent->mbShiftKey &&
         !ipEvent->mbCtrlKey &&
         !ipEvent->mbAltKey ) {

      switch ( ipEvent->mButton ) {
      case 2:
      case 3:
        /* button determines whether we're using paint or erase value */
        if ( 2 == ipEvent->mButton ) {
          nVolValue = sBrush.mInfo[DspA_tBrush_EditOne].mNewValue;
        } else if ( 3 == ipEvent->mButton ) {
          nVolValue = sBrush.mInfo[DspA_tBrush_EditTwo].mNewValue;
        }

        if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
          volumeType = tkm_tVolumeType_Aux;
        } else {
          volumeType = tkm_tVolumeType_Main;
        }

        tkm_FloodFillAnatomicalVolume( volumeType,
                                       pVolumeVox,
                                       nVolValue,
                                       sBrush.mb3DFill,
                                       sBrush.mFuzzy,
                                       sBrush.mDistance );
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
        break;
      }
    }



    break;


  case DspA_tTool_EditSegmentation:

    if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume] ) {
      segType = tkm_tSegType_Aux;
    } else {
      segType = tkm_tSegType_Main;
    }

    /* No modifications, its a normal edit. */
    if ( !ipEvent->mbShiftKey
         && !ipEvent->mbCtrlKey
         && !ipEvent->mbAltKey ) {

      switch ( ipEvent->mButton ) {
      case 2:
      case 3:
        /* button determines whether we're using paint or erase value */
        if ( 2 == ipEvent->mButton ) {
          sSegBrush.mNewValue = sSegBrush.mnPaintValue;
        } else if ( 3 == ipEvent->mButton ) {
          sSegBrush.mNewValue = sSegBrush.mnEraseValue;
        }

        /* edit the seg volume */
        sSegBrush.mDest = segType;
        eResult = DspA_BrushVoxels_( this, pVolumeVox,
                                     NULL, DspA_EditSegmentationVoxels_ );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* editing requires us to rebuild buffer. */
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
      }
    }

    /* shift-ctrl-2 sucks the color */
    if ( ipEvent->mbShiftKey &&
         ipEvent->mbCtrlKey &&
         !ipEvent->mbAltKey &&
         ipEvent->mButton == 2 ) {

      /* get the color and set our brush info with the same settings
         except for the new color */
      tkm_GetSegLabel( segType, pVolumeVox, &nSegIndex, NULL );
      segBrush = sSegBrush;
      segBrush.mnPaintValue = nSegIndex;
      segBrush.mDest = segType;
      DspA_SetSegBrushInfo( this, &segBrush );
    }

    /* Shift-2 or Shift-3 does a fill. */
    if ( ipEvent->mbShiftKey &&
         !ipEvent->mbCtrlKey &&
         !ipEvent->mbAltKey ) {

      switch ( ipEvent->mButton ) {
      case 2:
      case 3:
        /* button determines whether we're using paint or erase value */
        if ( 2 == ipEvent->mButton ) {
          sSegBrush.mNewValue = sSegBrush.mnPaintValue;
        } else if ( 3 == ipEvent->mButton ) {
          sSegBrush.mNewValue = sSegBrush.mnEraseValue;
        }

        tkm_FloodFillSegmentation( segType, pVolumeVox,
                                   sSegBrush.mNewValue,
                                   sSegBrush.mb3D, sSegBrush.mSrc,
                                   sSegBrush.mFuzzy, sSegBrush.mDistance );
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
        break;
      }
    }

    /* send the cursor info again to update the new seg label volume
       after editing. */
    eResult = DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
              pVolumeVox );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;


  case DspA_tTool_EditCtrlPts:

    if ( !ipEvent->mbCtrlKey &&
         !ipEvent->mbAltKey ) {

      Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                               pVolumeVox, &MRIIdx );

      /* if button 2, make a new point here. */
      if ( 2 == ipEvent->mButton ) {

        tkm_MakeControlPoint( &MRIIdx );

        /* if button 3, delete the nearest point. */
      } else if ( 3 == ipEvent->mButton ) {

        tkm_RemoveControlPointWithinDist( &MRIIdx, this->mOrientation, 3 );
      }
    }
    break;


  case DspA_tTool_Line:

    /* If button 2, set one of the vertices. If button three, clear
       both vertices. */
    if ( !ipEvent->mbCtrlKey &&
         !ipEvent->mbAltKey ) {

      DspA_NormalizeVoxel_( pVolumeVox, this->mOrientation, &lineVox );

      if ( 2 == ipEvent->mButton ) {

        if ( this->mLineVertex1.mnX < 0 || this->mLineVertex1.mnY < 0 ) {
          this->mLineVertex1.mnX = xVoxl_GetX( &lineVox );
          this->mLineVertex1.mnY = xVoxl_GetY( &lineVox );

        } else {

          if ( this->mLineVertex2.mnX >= 0 || this->mLineVertex2.mnY >= 0 ) {
            this->mLineVertex1 = this->mLineVertex2;
          }
          this->mLineVertex2.mnX = xVoxl_GetX( &lineVox );
          this->mLineVertex2.mnY = xVoxl_GetY( &lineVox );
          ;
        }


        /* Calculate the distance. */
        xVoxl_Set( &lineVox1, this->mLineVertex1.mnX,
                   this->mLineVertex1.mnY, 0 );
        DspA_UnnormalizeVoxel_( &lineVox1, this->mOrientation, &lineIdx1 );
        Volm_ConvertIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                              &lineIdx1, &lineRAS1 );

        xVoxl_Set( &lineVox2, this->mLineVertex2.mnX,
                   this->mLineVertex2.mnY, 0 );
        DspA_UnnormalizeVoxel_( &lineVox2, this->mOrientation, &lineIdx2 );
        Volm_ConvertIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                              &lineIdx2, &lineRAS2 );

        this->mLineDistance =
          sqrt( (xVoxl_GetFloatX(&lineRAS2) - xVoxl_GetFloatX(&lineRAS1)) *
                (xVoxl_GetFloatX(&lineRAS2) - xVoxl_GetFloatX(&lineRAS1))  +
                (xVoxl_GetFloatY(&lineRAS2) - xVoxl_GetFloatY(&lineRAS1)) *
                (xVoxl_GetFloatY(&lineRAS2) - xVoxl_GetFloatY(&lineRAS1))  +
                (xVoxl_GetFloatZ(&lineRAS2) - xVoxl_GetFloatZ(&lineRAS1)) *
                (xVoxl_GetFloatZ(&lineRAS2) - xVoxl_GetFloatZ(&lineRAS1)) );

      } else if ( 3 == ipEvent->mButton ) {

        this->mLineVertex1.mnX = -1;
        this->mLineVertex1.mnY = -1;
        this->mLineVertex2 = this->mLineVertex1;
        this->mLineDistance = 0;
      }
      DspA_BuildLineToolVoxelList_( this );
      this->mbSliceChanged = TRUE;
    }
    break;


  default:
    break;
  }

  /* If this isn't the nav tool, or is, but the cntrl key is down, set
     the cursor. */
  if ( !ipEvent->mbShiftKey &&
       (!(DspA_tTool_Navigate == sTool) ||
        (DspA_tTool_Navigate == sTool && ipEvent->mbCtrlKey)) ) {

    /* set the cursor. */
    eResult = DspA_SetCursor( this, pVolumeVox );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* if ctrl (only) was down... */
  if ( ipEvent->mbCtrlKey &&
       !ipEvent->mbShiftKey ) {

    /* set zoom center to cursor */
    eResult = DspA_SetZoomCenterToCursor( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;

    /* also zoom in or out. */
    if ( 1 == ipEvent->mButton ) {

      /* zoom in */
      eResult = DspA_SetZoomLevel( this, this->mnZoomLevel*2 );
      if ( DspA_tErr_NoErr != eResult )
        goto error;

    } else if ( 3 == ipEvent->mButton ) {

      /* zoom out */
      eResult = DspA_SetZoomLevel( this, this->mnZoomLevel/2 );
      if ( DspA_tErr_NoErr != eResult )
        goto error;
    }
  }



  /* Most things want a mouse up after we're done, so schedule one and
     let the draw function figure out if we don't need it. */
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleMouseUp_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleMouseDown_ ( tkmDisplayAreaRef this,
                                  xGWin_tEventRef   ipEvent ) {

  DspA_tErr          eResult      = DspA_tErr_NoErr;
  xPoint2n           bufferPt     = {0,0};
  xVoxelRef          pVolumeVox   = NULL;
  DspA_tBrush        brushAction  = DspA_tBrush_None;
  DspA_tSelectAction selectAction = DspA_tSelectAction_None;

  xVoxl_New( &pVolumeVox );

  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* If shift button-1 was down, we're going to do an interactive
     brightness/contrast modification. */
  if ( 1 == ipEvent->mButton &&
       ipEvent->mbShiftKey &&
       !ipEvent->mbCtrlKey &&
       !ipEvent->mbAltKey ) {

    /* Record the location of the click down. */
    this->mLastClick = ipEvent->mWhere;
    this->mTotalDelta.mfX = 0;
    this->mTotalDelta.mfY = 0;

    /* Get the original brightness and contrast. */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      Volm_GetBrightnessAndContrast( this->mpVolume[tkm_tVolumeType_Aux],
                                     &this->mfOriginalBrightness,
                                     &this->mfOriginalContrast );
    } else {
      Volm_GetBrightnessAndContrast( this->mpVolume[tkm_tVolumeType_Main],
                                     &this->mfOriginalBrightness,
                                     &this->mfOriginalContrast );
    }

    goto cleanup;
  }

  switch ( sTool ) {

  case DspA_tTool_Navigate:

    this->mLastClick = ipEvent->mWhere;
    this->mTotalDelta.mfX = 0;
    this->mTotalDelta.mfY = 0;
    xVoxl_Copy( this->mpOriginalZoomCenter, this->mpZoomCenter );
    this->mnOriginalSlice = DspA_GetCurrentSliceNumber_( this );
    this->mnOriginalZoomLevel = this->mnZoomLevel;
    break;


  case DspA_tTool_SelectVoxels:

    switch ( ipEvent->mButton ) {
      /* button 2 or 3 select tool with no modifiers: */
    case 2:
    case 3:
      if ( !ipEvent->mbShiftKey
           && !ipEvent->mbCtrlKey
           && !ipEvent->mbAltKey ) {

        /* button determines the select action */
        if ( 2 == ipEvent->mButton ) {
          selectAction = DspA_tSelectAction_Select;
        } else if ( 3 == ipEvent->mButton ) {
          selectAction = DspA_tSelectAction_Deselect;
        }

        /* brush the voxels */
        eResult = DspA_BrushVoxels_( this, pVolumeVox, (void*)&selectAction,
                                     DspA_SelectVoxels_ );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* selecting requires us to rebuild buffer. */
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
      }
      break;
    }
    break;


  case DspA_tTool_EditVoxels:

    switch ( ipEvent->mButton ) {
      /* button 2 or 3 edit tool with no modifiers: */
    case 2:
    case 3:
      if ( !ipEvent->mbShiftKey
           && !ipEvent->mbCtrlKey
           && !ipEvent->mbAltKey ) {

        /* clear the undo list. */
        tkm_ClearUndoList();

        /* button determines the brush action */
        if ( 2 == ipEvent->mButton ) {
          brushAction = DspA_tBrush_EditOne;
        } else if ( 3 == ipEvent->mButton ) {
          brushAction = DspA_tBrush_EditTwo;
        }

        /* brush the voxels */
        eResult = DspA_BrushVoxels_( this, pVolumeVox, (void*)&brushAction,
                                     DspA_BrushVoxelsInThreshold_ );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* editing requires us to rebuild buffer. */
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
      }
      break;
    }
    break;


  case DspA_tTool_EditSegmentation:

    switch ( ipEvent->mButton ) {
      /* button 2 or 3 edit seg tool with no modifiers: */
    case 2:
    case 3:
      if ( !ipEvent->mbShiftKey
           && !ipEvent->mbCtrlKey
           && !ipEvent->mbAltKey ) {

        tkm_ClearUndoList();
      }
      break;
    }
    break;

  default:
    break;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleMouseDown_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleMouseMoved_ ( tkmDisplayAreaRef this,
                                   xGWin_tEventRef   ipEvent ) {
  DspA_tErr          eResult      = DspA_tErr_NoErr;
  xPoint2n           bufferPt     = {0,0};
  xVoxel             anaIdx;
  DspA_tBrush        brushAction  = DspA_tBrush_None;
  DspA_tSelectAction selectAction = DspA_tSelectAction_None;
  xPoint2f           delta;
  xPoint2f           newCenterPt;
  xVoxel             newCenterIdx;
  int                nNewSlice    = 0;
  int                nNewZoomLevel= 0;
  tkm_tSegType       segType      = tkm_tSegType_Main;
  float              newBrightness= 0;
  float              newContrast  = 0;

  DebugEnterFunction( ("DspA_HandleMouseMoved_( this=%p, ipEvent=%p )",
                       this, ipEvent) );

  /* For some reason we get MouseMoved events for points at the same
     value as the window width or height. So check that here and if
     that's the case, skip. */
  if ( ipEvent->mWhere.mnX == this->mnWidth ||
       ipEvent->mWhere.mnY == this->mnHeight )
    goto cleanup;

  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, &anaIdx );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Check if this is a valid voxel. If so, send the information about
     this point to Tcl for the mouseover info pane. */
  eResult = DspA_VerifyVolumeVoxel_( this, &anaIdx );
  if ( DspA_tErr_NoErr == eResult ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Mouseover,
                                     &anaIdx );
    xVoxl_Copy( this->mpMouseLocationAnaIdx, &anaIdx );
  } else {
    goto cleanup;
  }

  /* save this mouse location */
  this->mMouseLocation = ipEvent->mWhere;

  /* if a button isn't down, skip this */
  if ( ipEvent->mButton == 0 )
    goto cleanup;

  /* If shift button-1 was down, we're going to do an interactive
     brightness/contrast modification. */
  if ( 1 == ipEvent->mButton &&
       ipEvent->mbShiftKey &&
       !ipEvent->mbCtrlKey &&
       !ipEvent->mbAltKey ) {

    /* Get the delta. */
    delta.mfX = (float)(this->mLastClick.mnX - ipEvent->mWhere.mnX) / 
      (float)this->mnZoomLevel / 2.0;
    delta.mfY = (float)(this->mLastClick.mnY - ipEvent->mWhere.mnY) / 
      (float)this->mnZoomLevel / 2.0;

    /* flip y if horizontal cuz our freaking screen is upside down */
    if ( this->mOrientation == mri_tOrientation_Horizontal )
      delta.mfY = -delta.mfY;

    /* add to the total delta */
    this->mTotalDelta.mfX += delta.mfX;
    this->mTotalDelta.mfY += delta.mfY;

    /* save this mouse position */
    this->mLastClick = ipEvent->mWhere;

    /* Delta brightness and contrast is a factor of the delta x and y,
       respectively. */
    newBrightness = this->mfOriginalBrightness +
                    (-this->mTotalDelta.mfX / 512.0);
    newContrast = this->mfOriginalContrast +
                  (this->mTotalDelta.mfY / 512.0 * 30.0);


    /* Apply the changes and redraw. */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      tkm_SetVolumeBrightnessContrast( tkm_tVolumeType_Aux,
                                       newBrightness, newContrast );
    } else {
      tkm_SetVolumeBrightnessContrast( tkm_tVolumeType_Main,
                                       newBrightness, newContrast );
    }

    /* Redraw the buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  switch ( sTool ) {

  case DspA_tTool_Navigate:

    if ( !ipEvent->mbCtrlKey ) {

      /* figure out how much we've moved */
      delta.mfX = (float)(this->mLastClick.mnX - ipEvent->mWhere.mnX) / 
        (float)this->mnZoomLevel / 2.0;
      delta.mfY = (float)(this->mLastClick.mnY - ipEvent->mWhere.mnY) / 
        (float)this->mnZoomLevel / 2.0;

      /* flip y if horizontal cuz our freaking screen is upside down */
      if ( this->mOrientation == mri_tOrientation_Horizontal )
        delta.mfY = -delta.mfY;

      /* add to the total delta */
      this->mTotalDelta.mfX += delta.mfX;
      this->mTotalDelta.mfY += delta.mfY;

      /* save this mouse position */
      this->mLastClick = ipEvent->mWhere;

      switch ( ipEvent->mButton ) {

        /* button one is 'drag' the volume view. add the rounded delta
           to the original center and set the center */
      case 1:
        newCenterPt.mfX = this->mTotalDelta.mfX +
                          xVoxl_GetX( this->mpOriginalZoomCenter );
        newCenterPt.mfY = this->mTotalDelta.mfY +
                          xVoxl_GetY( this->mpOriginalZoomCenter );
        DspA_ConvertPlaneToVolume_( this, &newCenterPt,
                                    DspA_GetCurrentSliceNumber_( this ),
                                    this->mOrientation, &newCenterIdx );
        DspA_SetZoomCenter( this, &newCenterIdx );
        break;

        /* button 2 is move the slice. add the rounded y delta to the
           slice */
      case 2:
        nNewSlice = this->mnOriginalSlice + rint(this->mTotalDelta.mfY);
        if ( nNewSlice >= 0 && nNewSlice < this->mnVolumeSizeZ ) {
          DspA_SetSlice( this, nNewSlice );
        }
        break;

        /* button 3 is zoom. add the rounded delta to the zoom level */
      case 3:
        nNewZoomLevel = this->mnOriginalZoomLevel +rint(this->mTotalDelta.mfY);
        if ( nNewZoomLevel >= DspA_knMinZoomLevel
             && nNewZoomLevel <= DspA_knMaxZoomLevel ) {
          DspA_SetZoomLevel( this, nNewZoomLevel );
        }
        break;
      }
    }
    break;


  case DspA_tTool_SelectVoxels:

    switch ( ipEvent->mButton ) {
      /* button 2 or 3 select tool with no modifiers: */
    case 2:
    case 3:
      if ( !ipEvent->mbShiftKey
           && !ipEvent->mbCtrlKey
           && !ipEvent->mbAltKey ) {

        /* button determines the select action */
        if ( 2 == ipEvent->mButton ) {
          selectAction = DspA_tSelectAction_Select;
        } else if ( 3 == ipEvent->mButton ) {
          selectAction = DspA_tSelectAction_Deselect;
        }

        /* brush the voxels */
        eResult = DspA_BrushVoxels_( this, &anaIdx, (void*)&selectAction,
                                     DspA_SelectVoxels_ );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* selecting requires us to rebuild buffer. */
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
      }
      break;
    }
    break;


  case DspA_tTool_EditVoxels:

    switch ( ipEvent->mButton ) {
      /* button 2 or 3 edit tool with no modifiers: */
    case 2:
    case 3:
      if ( !ipEvent->mbShiftKey
           && !ipEvent->mbCtrlKey
           && !ipEvent->mbAltKey ) {

        /* button determines the brush action */
        if ( 2 == ipEvent->mButton ) {
          brushAction = DspA_tBrush_EditOne;
        } else if ( 3 == ipEvent->mButton ) {
          brushAction = DspA_tBrush_EditTwo;
        }

        /* brush the voxels */
        eResult = DspA_BrushVoxels_( this, &anaIdx, (void*)&brushAction,
                                     DspA_BrushVoxelsInThreshold_ );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* editing requires us to rebuild buffer. */
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
      }
      break;
    }
    break;


  case DspA_tTool_EditSegmentation:

    switch ( ipEvent->mButton ) {
      /* button 2 or 3 edit seg tool with no modifiers: */
    case 2:
    case 3:
      if ( !ipEvent->mbShiftKey
           && !ipEvent->mbCtrlKey
           && !ipEvent->mbAltKey ) {

        if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume] ) {
          segType = tkm_tSegType_Aux;
        } else {
          segType = tkm_tSegType_Main;
        }

        /* button determines whether we're using paint or erase value */
        if ( 2 == ipEvent->mButton ) {
          sSegBrush.mNewValue = sSegBrush.mnPaintValue;
        } else if ( 3 == ipEvent->mButton ) {
          sSegBrush.mNewValue = sSegBrush.mnEraseValue;
        }

        /* edit the seg volume */
        sSegBrush.mDest = segType;
        eResult = DspA_BrushVoxels_( this, &anaIdx,
                                     NULL, DspA_EditSegmentationVoxels_ );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* editing requires us to rebuild buffer. */
        this->mbSliceChanged = TRUE;
        DspA_Redraw_( this );
      }
      break;
    }
    break;
  default:
    break;
  }


  /* Clear the error flag. */
  eResult = DspA_tErr_NoErr;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleMouseMoved_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_HandleKeyDown_ ( tkmDisplayAreaRef this,
                                xGWin_tEventRef   ipEvent ) {

  DspA_tErr eResult        = DspA_tErr_NoErr;
  MWin_tErr eWindowResult  = MWin_tErr_NoErr;
  FunV_tErr eFunctional    = FunV_tErr_NoError;
  tBoolean  bOverlayLoaded = FALSE;

  /* Ctrl key combos */
  if ( ipEvent->mbCtrlKey && !ipEvent->mbAltKey ) {
    switch ( ipEvent->mKey ) {
    case '1':
      eResult =
        DspA_SetDisplayFlag( this, DspA_tDisplayFlag_AuxVolume, FALSE );
      break;
    case '2':
      eResult =
        DspA_SetDisplayFlag( this, DspA_tDisplayFlag_AuxVolume, TRUE );
      break;
    case 'a':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_Anatomical );
      break;
    case 'c':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_Cursor );
      break;
    case 'd':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_DTIOverlay );
      break;
    case 'f':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay );
      break;
    case 'g':
      eResult =
        DspA_ToggleDisplayFlag( this, 
                                DspA_tDisplayFlag_SegmentationVolumeOverlay );
      break;
    case 'i':
      eResult =
        DspA_ToggleDisplayFlag( this, 
                                DspA_tDisplayFlag_InterpolateSurfaceVertices );
      break;
    case 'k':
      eWindowResult = MWin_ToggleLinkedCursorFlag( this->mpWindow );
      if ( MWin_tErr_NoErr != eWindowResult ) {
        eResult = DspA_tErr_ErrorAccessingWindow;
      }
      break;
    case 'm':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_MainSurface );
      break;
    case 'o':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_OriginalSurface );
      break;
    case 'p':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_PialSurface );
      break;
    case 's':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_Selection );
      break;
    case 't':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_ControlPoints );
      break;
    case 'v':
      eResult =
        DspA_ToggleDisplayFlag( this, 
                                DspA_tDisplayFlag_DisplaySurfaceVertices );
      break;
    case 'z':
      tkm_RestoreUndoList();
      break;
    case '+':
      if ( this->mnZoomLevel < DspA_knMaxZoomLevel ) {
        eResult = DspA_SetZoomLevel( this, this->mnZoomLevel * 2 );
      }
      break;
    case '-':
      if ( this->mnZoomLevel > 1 ) {
        eResult = DspA_SetZoomLevel( this, this->mnZoomLevel / 2 );
      }
      break;
    }
  }

  /* Plain keys */
  if ( !ipEvent->mbCtrlKey && !ipEvent->mbAltKey ) {
    switch ( ipEvent->mKey ) {
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      eResult = DspA_SetBrushShape( this, (int)(ipEvent->mKey - '0'),
                                    sBrush.mShape, sBrush.mb3D );
      break;
    case 'a':
      eResult = DspA_SetTool( this, DspA_tTool_EditVoxels );
      break;
    case 'g':
      eResult = DspA_SetTool( this, DspA_tTool_EditSegmentation );
      break;
    case 'l':
      eResult = DspA_SetTool( this, DspA_tTool_Line );
      break;
    case 'n':
      eResult = DspA_SetTool( this, DspA_tTool_Navigate );
      break;
    case 's':
      eResult = DspA_SetTool( this, DspA_tTool_SelectVoxels );
      break;
    case 't':
      eResult = DspA_SetTool( this, DspA_tTool_EditCtrlPts );
      break;
    case 'x':
      eResult = DspA_SetOrientation( this, mri_tOrientation_Sagittal );
      break;
    case 'y':
      eResult = DspA_SetOrientation( this, mri_tOrientation_Horizontal );
      break;
    case 'z':
      eResult = DspA_SetOrientation( this, mri_tOrientation_Coronal );
      break;
    case xGWin_tKey_UpArrow:
    case xGWin_tKey_RightArrow:
      eResult = DspA_ChangeSliceBy( this, 1 );
      DspA_Redraw_( this );
      break;
    case xGWin_tKey_DownArrow:
    case xGWin_tKey_LeftArrow:
      eResult = DspA_ChangeSliceBy( this, -1 );
      DspA_Redraw_( this );
      break;
    case '+':
      if ( NULL != this->mpFunctionalVolume ) {
        eFunctional = FunV_ChangeTimePointBy( this->mpFunctionalVolume, 1 );
        if ( FunV_tErr_NoError != eFunctional ) {
          eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
        }
      }
      break;
    case '-':
      if ( NULL != this->mpFunctionalVolume ) {
        eFunctional = FunV_ChangeTimePointBy( this->mpFunctionalVolume, -1 );
        if ( FunV_tErr_NoError != eFunctional ) {
          eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
        }
      }
      break;
    }
  }

  /* Alt-key combos */
  if ( !ipEvent->mbCtrlKey && ipEvent->mbAltKey ) {
    switch ( ipEvent->mKey ) {
    case 'a':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_AuxVolume );
      break;
    case 'c':
      eResult =
        DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_AuxVolume );
      break;
    case 'f':
      if ( NULL != this->mpFunctionalVolume ) {
        FunV_IsOverlayPresent( this->mpFunctionalVolume, &bOverlayLoaded );
        if ( bOverlayLoaded ) {
          if ( 1 == this->mabDisplayFlags[DspA_tDisplayFlag_Anatomical] ) {
            DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Anatomical, 0 );
            DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay,1 );
          } else {
            DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Anatomical, 1 );
            DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay,0 );
          }
        }
      }
      break;
    case 'g':
      eResult =
        DspA_ToggleDisplayFlag( this, 
                                DspA_tDisplayFlag_AuxSegmentationVolume );
      break;
    case 'm':
      tkm_ShowNearestSurfaceVertex( Surf_tVertexSet_Main );
      break;
    case 'o':
      tkm_ShowNearestSurfaceVertex( Surf_tVertexSet_Original );
      break;
    case 'p':
      tkm_ShowNearestSurfaceVertex( Surf_tVertexSet_Pial );
      break;
    case 'x':
      eResult = DspA_SetOrientation( this, mri_tOrientation_Sagittal );
      break;
    case 'y':
      eResult = DspA_SetOrientation( this, mri_tOrientation_Horizontal );
      break;
    case 'z':
      eResult = DspA_SetOrientation( this, mri_tOrientation_Coronal );
      break;
    case '+':
      eResult = DspA_SetZoomLevel( this, DspA_knMaxZoomLevel );
      break;
    case '-':
      eResult = DspA_SetZoomLevel( this, 1 );
      break;
    }
  }

  if ( DspA_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleKeyDown_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetBrushInfoToDefault ( tkmDisplayAreaRef this,
                                       DspA_tBrush       iBrush ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make sure the brush is in bounds */
  if ( iBrush < 0 ||
       iBrush >= DspA_knNumBrushes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set the brush theshold info */
  switch ( iBrush ) {
  case DspA_tBrush_EditOne:
    sBrush.mInfo[iBrush].mLow = tkm_knEditToWhiteLow;
    sBrush.mInfo[iBrush].mHigh = tkm_knEditToWhiteHigh;
    sBrush.mInfo[iBrush].mNewValue = tkm_knEditToWhiteNewValue;
    break;
  case DspA_tBrush_EditTwo:
    sBrush.mInfo[iBrush].mLow = tkm_knEditToBlackLow;
    sBrush.mInfo[iBrush].mHigh = tkm_knEditToBlackHigh;
    sBrush.mInfo[iBrush].mNewValue = tkm_knEditToBlackNewValue;
    break;
  default:
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* if we're the currently focused display... */
  if ( sFocusedDisplay == this ) {

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %f %f %f %d %d",
              (int)iBrush,
              (float)sBrush.mInfo[iBrush].mLow,
              (float)sBrush.mInfo[iBrush].mHigh,
              (float)sBrush.mInfo[iBrush].mNewValue,
              (int)sBrush.mInfo[iBrush].mMode,
              (int)sBrush.mInfo[iBrush].mCloneSource );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushInfo, sTclArguments );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetBrushInfoToDefault: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_BrushVoxels_ ( tkmDisplayAreaRef this,
                              xVoxelRef         ipCenterVox,
                              void*             ipData,
                              void(*ipFunction)(xVoxelRef,int,void*) ) {

  DspA_tErr        eResult     = DspA_tErr_NoErr;
  int              nDimensionX = 0;
  int              nDimensionY = 0;
  int              nDimensionZ = 0;
  int              nXCenter    = 0;
  int              nYCenter    = 0;
  int              nZCenter    = 0;
  int              nXRadius    = 0;
  int              nYRadius    = 0;
  int              nZRadius    = 0;
  int              nXMin       = 0;
  int              nXMax       = 0;
  int              nYMin       = 0;
  int              nYMax       = 0;
  int              nZMin       = 0;
  int              nZMax       = 0;
  xVoxel           MRIIdx;
  xVoxelRef        paBrushedVoxs = NULL;
  int              nBrushedVoxs = 0;

  DebugEnterFunction( ("DspA_BrushVoxels_( this=%p, ipCenterVox=%p, "
                       "ipData=%p )", this, ipCenterVox, ipData) );

  /* Convert the center vox to an MRI index. Also get the dimensions
     of the right volume. */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
    Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Aux],
                             ipCenterVox, &MRIIdx );
    Volm_GetDimensions( this->mpVolume[tkm_tVolumeType_Aux],
                        &nDimensionX, &nDimensionY, &nDimensionZ );
  } else {
    Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                             ipCenterVox, &MRIIdx );
    Volm_GetDimensions( this->mpVolume[tkm_tVolumeType_Main],
                        &nDimensionX, &nDimensionY, &nDimensionZ );
  }

  /* get our center voxel. */
  nXCenter = xVoxl_GetX( &MRIIdx );
  nYCenter = xVoxl_GetY( &MRIIdx );
  nZCenter = xVoxl_GetZ( &MRIIdx );

  /* set all radii to the brush radius. we subtract one because of our
     looping bounds. */
  nXRadius = sBrush.mnRadius - 1;
  nYRadius = sBrush.mnRadius - 1;
  nZRadius = sBrush.mnRadius - 1;

  /* if we're not in 3d, set the same radius of the same plane as our
     current orientation to 0. */
  if ( !sBrush.mb3D ) {
    switch ( this->mOrientation ) {
    case mri_tOrientation_Coronal:
      nZRadius = 0;
      break;
    case mri_tOrientation_Horizontal:
      nYRadius = 0;
      break;
    case mri_tOrientation_Sagittal:
      nXRadius = 0;
      break;
    default:
      eResult = DspA_tErr_InvalidOrientation;
      goto error;
      break;
    }
  }


  /* Calc the limits and cap them. */
  nXMin = MAX( nXCenter - nXRadius, 0 );
  nXMax = MIN( nXCenter + nXRadius, nDimensionX );
  nYMin = MAX( nYCenter - nYRadius, 0 );
  nYMax = MIN( nYCenter + nYRadius, nDimensionY );
  nZMin = MAX( nZCenter - nZRadius, 0 );
  nZMax = MIN( nZCenter + nZRadius, nDimensionZ );

  /* set our voxel. */
  xVoxl_Set( &MRIIdx, nXMin, nYMin, nZMin );

  /* Allocate the array of voxels. Set it to the max it can be, with
     is the size of the cuboid in the bounds calc'd above. */
  paBrushedVoxs = malloc( sizeof(xVoxel) *
                          (nXMax - nXMin + 1) *
                          (nYMax - nYMin + 1) *
                          (nZMax - nZMin + 1));
  DebugAssertThrowX( (NULL != paBrushedVoxs),
                     eResult, DspA_tErr_AllocationFailed );

  /* loop around the starting point... */
  do {

    /* if we're circular, check the radius. if no good, continue. */
    if ( DspA_tBrushShape_Circle == sBrush.mShape &&
         ( (nXCenter - xVoxl_GetX(&MRIIdx)) *
           (nXCenter - xVoxl_GetX(&MRIIdx)) +
           (nYCenter - xVoxl_GetY(&MRIIdx)) *
           (nYCenter - xVoxl_GetY(&MRIIdx)) +
           (nZCenter - xVoxl_GetZ(&MRIIdx)) *
           (nZCenter - xVoxl_GetZ(&MRIIdx)) >
           (sBrush.mnRadius-1) * (sBrush.mnRadius-1) ) ) {
      continue;
    }

    /* Check the bounds */
    if( xVoxl_GetX(&MRIIdx) < 0 || xVoxl_GetX(&MRIIdx) >= nDimensionX ||
	xVoxl_GetY(&MRIIdx) < 0 || xVoxl_GetY(&MRIIdx) >= nDimensionY ||
	xVoxl_GetZ(&MRIIdx) < 0 || xVoxl_GetZ(&MRIIdx) >= nDimensionZ )
      continue;


    /* Add this voxel to the list and increment the count. */
    xVoxl_Copy( &(paBrushedVoxs[nBrushedVoxs]), &MRIIdx );
    nBrushedVoxs++;

  } while ( xVoxl_IncrementWithMinsUntilLimits( &MRIIdx, nXMin, nYMin,
            nXMax, nYMax, nZMax ));

  /* run the function on these voxels. */
  ipFunction( paBrushedVoxs, nBrushedVoxs, ipData );

  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_BrushVoxels_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  if ( NULL != paBrushedVoxs ) {
    free( paBrushedVoxs );
  }

  DebugExitFunction;

  return eResult;
}

void DspA_BrushVoxelsInThreshold_ ( xVoxelRef ipaVoxel, int inCount,
                                    void* ipData ) {

  DspA_tBrush brush  = DspA_tBrush_None;

  /* make sure the brush is in bounds */
  brush = *(DspA_tBrush*)ipData;
  if ( brush < 0 ||
       brush >= DspA_knNumBrushes ) {
    return;
  }

  /* edit the voxel according to what the brush target and mode is. */
  if ( DspA_tBrushTarget_Main == sBrush.mTarget ||
       DspA_tBrushTarget_MainAux == sBrush.mTarget ) {

    switch ( sBrush.mInfo[brush].mMode ) {
    case DspA_tBrushMode_Set:
      tkm_EditAnatomicalVolumeInRangeArray( tkm_tVolumeType_Main,
                                            ipaVoxel, inCount,
                                            sBrush.mInfo[brush].mLow,
                                            sBrush.mInfo[brush].mHigh,
                                            sBrush.mInfo[brush].mNewValue );
      break;
    case DspA_tBrushMode_Clone:
      tkm_CloneAnatomicalVolumeInRangeArray( tkm_tVolumeType_Main,
                                             sBrush.mInfo[brush].mCloneSource,
                                             ipaVoxel, inCount,
                                             sBrush.mInfo[brush].mLow,
                                             sBrush.mInfo[brush].mHigh );
      break;
    default:
      break;
    }
  }
  if ( DspA_tBrushTarget_Aux == sBrush.mTarget ||
       DspA_tBrushTarget_MainAux == sBrush.mTarget ) {

    switch ( sBrush.mInfo[brush].mMode ) {
    case DspA_tBrushMode_Set:
      tkm_EditAnatomicalVolumeInRangeArray( tkm_tVolumeType_Aux,
                                            ipaVoxel, inCount,
                                            sBrush.mInfo[brush].mLow,
                                            sBrush.mInfo[brush].mHigh,
                                            sBrush.mInfo[brush].mNewValue );
      break;
    case DspA_tBrushMode_Clone:
      tkm_CloneAnatomicalVolumeInRangeArray( tkm_tVolumeType_Aux,
                                             sBrush.mInfo[brush].mCloneSource,
                                             ipaVoxel, inCount,
                                             sBrush.mInfo[brush].mLow,
                                             sBrush.mInfo[brush].mHigh );
      break;
    default:
      break;
    }
  }
}

void DspA_SelectVoxels_ ( xVoxelRef ipaVoxel, int inCount, void* ipData ) {

  DspA_tSelectAction selectAction = DspA_tSelectAction_None;

  /* make sure the action is in bounds */
  selectAction = *(DspA_tSelectAction*)ipData;
  if ( selectAction < 0 ||
       selectAction >= DspA_knNumSelectActions ) {
    return;
  }

  /* select or deselect the voxel */
  switch ( selectAction ) {
  case DspA_tSelectAction_Select:
    tkm_SelectVoxelArray( ipaVoxel, inCount );
    break;
  case DspA_tSelectAction_Deselect:
    tkm_DeselectVoxelArray( ipaVoxel, inCount );
    break;
  default:
    break;
  }
}

void DspA_EditSegmentationVoxels_ ( xVoxelRef ipaVoxel, int inCount,
                                    void* ipData ) {

  DebugEnterFunction( ("DspA_EditSegmentationVoxels_( ipaVoxel=%p, "
                       "inCount=%d, ipData=%p", ipaVoxel, inCount, ipData) );

  tkm_EditSegmentationArray( sSegBrush.mDest,
                             ipaVoxel, inCount, sSegBrush.mNewValue );

  DebugExitFunction;
}

DspA_tErr DspA_SendMouseInfoToTcl ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xPoint2n  bufferPt;
  xVoxel    idx;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* translate current mouse location to volume idx. iof these fail, it's not
     fatal, since we're probably just switching and configuring the pane
     locations and the mouse is out of the new pane's location.*/
  eResult = DspA_ConvertScreenToBuffer_( this, &(this->mMouseLocation),
                                         &bufferPt );
  if ( DspA_tErr_NoErr != eResult ) {
    eResult = DspA_tErr_NoErr;
    goto cleanup;
  }
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, &idx );
  if ( DspA_tErr_NoErr != eResult ) {
    eResult = DspA_tErr_NoErr;
    goto cleanup;
  }

  /* if good, send to tcl */
  eResult = DspA_VerifyVolumeVoxel_( this, &idx );
  if ( DspA_tErr_NoErr == eResult ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Mouseover, &idx );
  }

  eResult = DspA_tErr_NoErr;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SendMouseInfoToTcl: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_Redraw_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult       = DspA_tErr_NoErr;
  MWin_tErr eWindowResult = MWin_tErr_NoErr;

  /* call our mom's redraw func. */
  eWindowResult = MWin_Redraw( this->mpWindow );
  if ( MWin_tErr_NoErr != eWindowResult ) {
    eResult = DspA_tErr_ErrorAccessingWindow;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_Redraw_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_HandleDraw_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* don't draw if our height and width are nothing */
  if ( this->mnHeight == 0
       || this->mnWidth == 0 )
    goto cleanup;

  /* if the slice changed, build the current frame buffer */
  if ( this->mbSliceChanged ) {

    /* Draw the anatomical or lack thereof. (If anatomical is not
       displayed, will just clear entire buffer to black.) */
    eResult = DspA_BuildCurrentFrame_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }

    /* Segmentation overlay */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_SegmentationVolumeOverlay] ) {
      eResult = DspA_DrawSegmentationOverlayToFrame_( this );
      if ( DspA_tErr_NoErr != eResult ) {
        DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
        eResult = DspA_tErr_NoErr;
      }
    }

    /* DTI overlay */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_DTIOverlay] ) {
      eResult = DspA_DrawDTIOverlayToFrame_( this );
      if ( DspA_tErr_NoErr != eResult ) {
        DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
        eResult = DspA_tErr_NoErr;
      }
    }

    /* Undoable voxels overlay */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_UndoVolume] ) {
      eResult = DspA_DrawUndoableVoxelsOverlayToFrame_( this );
      if ( DspA_tErr_NoErr != eResult ) {
        DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
        eResult = DspA_tErr_NoErr;
      }
    }

    /* Draw functional overlay */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ) {
      eResult = DspA_DrawFunctionalOverlayToFrame_( this );
      if ( DspA_tErr_NoErr != eResult ) {
        DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
        eResult = DspA_tErr_NoErr;
      }
    }

    /* Selection overlay */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_Selection] ) {
      eResult = DspA_DrawSelectionToFrame_( this );
      if ( DspA_tErr_NoErr != eResult ) {
        DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
        eResult = DspA_tErr_NoErr;
      }
    }

    /* Line tool lines. */
    eResult = DspA_DrawLinesToFrame_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  /* Draw the frame buffer to screen. */
  eResult = DspA_DrawFrameBuffer_ ( this );
  if ( DspA_tErr_NoErr != eResult ) {
    DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
    eResult = DspA_tErr_NoErr;
  }

  /* Now a couple overlays that draw with direct openGL commands. */
  /* Head points */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_HeadPoints] ) {
    eResult = DspA_DrawHeadPoints_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  /* Control points */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_ControlPoints] ) {
    eResult = DspA_DrawControlPoints_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

#if 0
  /* Vector field. Never really finished debugging this. */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_VectorField] ) {
    eResult = DspA_DrawVectorField_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
#endif

  /* Draw the surface */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_MainSurface] ||
       this->mabDisplayFlags[DspA_tDisplayFlag_OriginalSurface] ||
       this->mabDisplayFlags[DspA_tDisplayFlag_PialSurface] ) {
    eResult = DspA_DrawSurface_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  /* Draw the cursor */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_Cursor] ) {
    eResult = DspA_DrawCursor_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  /* if we're focused and our focus frame flag is on... */
  if ( this == sFocusedDisplay
       && TRUE == this->mabDisplayFlags[DspA_tDisplayFlag_FocusFrame] ) {
    /* draw a frame around us. */
    eResult = DspA_DrawFrameAroundDisplay_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  if ( this->mabDisplayFlags[DspA_tDisplayFlag_Axes] ) {
    eResult = DspA_DrawAxes_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  /* Line tool lines. */
  eResult = DspA_DrawLines_ ( this );
  if ( DspA_tErr_NoErr != eResult ) {
    DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
    eResult = DspA_tErr_NoErr;
  }

  /* Draw functional overlay color bar */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalColorScaleBar] ) {
    eResult = DspA_DrawFunctionalOverlay_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

  /* rebuilt all our slice changes, so we can clear this flag. */
  this->mbSliceChanged = FALSE;

  /* we have to send the cursor info again to update the stars around
     the currently active volume or segmentation (which may now be
     different) */
  if ( sFocusedDisplay == this ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
                                     this->mpCursor );
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Mouseover,
                                     this->mpMouseLocationAnaIdx );
  }

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleDraw_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_DrawFrameBuffer_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  DspA_SetUpOpenGLPort_( this );

  glDrawPixels ( this->mnVolumeSizeX, this->mnVolumeSizeY,
                 GL_RGBA, GL_UNSIGNED_BYTE, this->mpFrameBuffer );

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawFrameBuffer_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_DrawSurface_ ( tkmDisplayAreaRef this ) {

  DspA_tErr         eResult    = DspA_tErr_NoErr;
  int               nSurface   = 0;
  Surf_tVertexSet   vertexSet    = Surf_tVertexSet_Main;
  xGrowableArrayRef list       = NULL;
  float             faColor[3] = {0, 0, 0};

  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {

    /* only recalculate the draw list of the plane and orientation has
       changed. */
    if ( NULL != this->mpSurface[nSurface] ) {

      /* if the orienation or slice has changed... */
      if ( this->mbSliceChanged ) {

        /* rebuild the draw lists. */
        eResult = DspA_BuildSurfaceDrawLists_( this, nSurface );

        /* if we got an out of memory error.. */
        if ( DspA_tErr_OutOfMemory == eResult ) {

          /* delete draw lists and try again. */
          DspA_PurgeSurfaceLists_( this );
          eResult = DspA_BuildSurfaceDrawLists_ ( this, nSurface );

          /* still got an error, return it. */
          if ( DspA_tErr_NoErr != eResult )
            goto error;

        } else if ( DspA_tErr_NoErr != eResult ) {
          goto error;
        }
      }

      /* set up opengl */
      DspA_SetUpOpenGLPort_( this );

      /* for each surface type... */
      for ( vertexSet = Surf_tVertexSet_Main;
            vertexSet < Surf_knNumVertexSets; vertexSet++ ) {

        /* if this surface is visible... */
        if ( this->mabDisplayFlags[vertexSet + 
                                   DspA_tDisplayFlag_MainSurface] ) {

          /* get the list. */
          list = DspA_GetSurfaceList_( this, nSurface,
                                       this->mOrientation, vertexSet,
                                       DspA_GetCurrentSliceNumber_(this) );
          if ( NULL == list ) {
            eResult = DspA_tErr_ErrorAccessingSurfaceList;
            goto error;
          }

          /* set the color */
          xColr_PackFloatArray
            (&(this->maSurfaceLineColor[nSurface][vertexSet]), faColor );

          /* draw the points. */
          glLineWidth(  this->manSurfaceLineWidth[nSurface][vertexSet] );
          DspA_ParsePointList_( this, GL_LINES, list, faColor );

          /* if vertices are visible... */
          if (this->mabDisplayFlags[DspA_tDisplayFlag_DisplaySurfaceVertices]){

            /* invert color */
            faColor[0] = 1.0 - faColor[0];
            faColor[1] = 1.0 - faColor[1];
            faColor[2] = 1.0 - faColor[2];

            /* set point size. */
            glPointSize( DspA_knSurfaceVertexSize );

            /* draw the vertices. */
            DspA_ParsePointList_( this, GL_POINTS, list, faColor );
          }

          /* if we have a hilited vertex for this surface... */
          if ( vertexSet == this->mHilitedSurface
               && -1 != this->mnHilitedVertexIndex  ) {

            /* ?? */
          }
        }
      }
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawSurface_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_DrawCursor_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult      = DspA_tErr_NoErr;
  xPoint2n  bufferPt     = {0, 0};
  float     faColor[3];
  int       nWidth       = 0;
  int       nHeight      = 0;

  /* convert the voxel and see if it's on the screen. the screen is
     actually sized to the buffer size, i.e. volumeSize x volumeSize,
     so only convert to buffer space.*/
  eResult = DspA_ConvertVolumeToBuffer_ ( this, this->mpCursor, &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  //  bufferPt.mnY = GLDRAW_Y_FLIP(bufferPt.mnY);

  /* calculate width and height using scale */
  nWidth  = ((float) DspA_knCursorCrosshairSize / this->mfFrameBufferScaleX);
  nHeight = ((float) DspA_knCursorCrosshairSize / this->mfFrameBufferScaleY);

  nWidth = MAX( nWidth, 1 );
  nHeight = MAX( nHeight, 1 );

  /* draw the crosshair */
  faColor[0] = sCursorColor.mfRed;
  faColor[1] = sCursorColor.mfGreen;
  faColor[2] = sCursorColor.mfBlue;
  DspA_DrawMarker_( this, sCursorShape, faColor, &bufferPt,
                    DspA_knCursorCrosshairSize );

  DspA_SetUpOpenGLPort_( this );

  /* draw the edge markers */
  glLineWidth( 1 );
  glColor3fv( faColor );

  glBegin( GL_LINES );
  glVertex2d( bufferPt.mnX, 0 );
  glVertex2d( bufferPt.mnX, 3*nHeight );
  glEnd();

  glBegin( GL_LINES );
  glVertex2d( bufferPt.mnX, this->mnVolumeSizeY );
  glVertex2d( bufferPt.mnX, this->mnVolumeSizeY - 3*nHeight );
  glEnd();

  glBegin( GL_LINES );
  glVertex2d( 0, bufferPt.mnY );
  glVertex2d( 3*nWidth, bufferPt.mnY );
  glEnd();

  glBegin( GL_LINES );
  glVertex2d( this->mnVolumeSizeX, bufferPt.mnY );
  glVertex2d( this->mnVolumeSizeX - 3*nWidth, bufferPt.mnY );
  glEnd();

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawCursor_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawFrameAroundDisplay_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult      = DspA_tErr_NoErr;

  DspA_SetUpOpenGLPort_( this );

  /* draw a green box around us */
  glColor3f ( 0.0, 1.0, 0.0 );

  glBegin( GL_LINE_STRIP );
  glVertex2d( 1, 1 );
  glVertex2d( this->mnVolumeSizeX-1, 1 );
  glVertex2d( this->mnVolumeSizeX-1, this->mnVolumeSizeY-1 );
  glVertex2d( 1, this->mnVolumeSizeY-1 );
  glVertex2d( 1, 1 );
  glEnd ();

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawFrameAroundDisplay_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawAxes_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  xPoint2n  volumePt    = {0, 0};
  xVoxelRef  pVoxel      = NULL;

  DspA_SetUpOpenGLPort_( this );

  glColor3f ( 0.0, 1.0, 0.0 );

  /* draw arrows */
  switch ( this->mOrientation ) {

  case mri_tOrientation_Coronal:
    volumePt.mnX = 40;
    volumePt.mnY = this->mnVolumeSizeY - 20;
    DspA_DrawHorizontalArrow_( &volumePt, -20, "R" );
    volumePt.mnX = this->mnVolumeSizeX - 20;
    volumePt.mnY = 40;
    DspA_DrawVerticalArrow_( &volumePt, -20, "S" );
    break;
  case mri_tOrientation_Horizontal:
    volumePt.mnX = 40;
    volumePt.mnY = this->mnVolumeSizeY - 20;
    DspA_DrawHorizontalArrow_( &volumePt, -20, "R" );
    volumePt.mnX = this->mnVolumeSizeX - 20;
    volumePt.mnY = 40;
    DspA_DrawVerticalArrow_( &volumePt, -20, "A" );
    break;
  case mri_tOrientation_Sagittal:
    volumePt.mnX = this->mnVolumeSizeX - 40;
    volumePt.mnY = this->mnVolumeSizeY - 20;
    DspA_DrawHorizontalArrow_( &volumePt, 20, "A" );
    volumePt.mnX = 20;
    volumePt.mnY = 40;
    DspA_DrawVerticalArrow_( &volumePt, -20, "S" );
    break;
  default:
    break;
  }


  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawAxes_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pVoxel );

  return eResult;
}

DspA_tErr DspA_DrawLinesToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  xPoint2n  bufferPt    = {0,0};
  GLubyte*  pFrame      = NULL;
  xColor3f  color;
  xVoxel    anaIdx;
  int       nLineVox    = 0;
  int       yMin        = 0;
  int       yMax        = 0;
  int       yInc        = 0;

  /* If we don't have good vertices, bail. */
  if ( this->mLineVertex1.mnX < 0 || this->mLineVertex1.mnY < 0 ||
       this->mLineVertex2.mnX < 0 || this->mLineVertex2.mnY < 0 )
    goto cleanup;

  /* get a ptr to the frame buffer. */
  pFrame = this->mpFrameBuffer;

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY-1;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY-1;
    yMax = 0;
    yInc = -1;
  }

  /* Just loop through getting the anatomical color at this index. */
  for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
    for ( bufferPt.mnX = 0;
          bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

      /* get a volume voxel.*/
      eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
      if ( DspA_tErr_NoErr != eResult )
        goto error;

      for ( nLineVox = 0; nLineVox < this->mNumLineVoxels; nLineVox++ ) {

        if ( xVoxl_IsEqualInt( &(this->mLineVoxels[nLineVox]), &anaIdx ) ) {

          /* get the current color in the buffer */
          xColr_SetFloat( &color, (float)pFrame[DspA_knRedPixelCompIndex] /
                          (float)DspA_knMaxPixelValue,
                          (float)pFrame[DspA_knGreenPixelCompIndex] /
                          (float)DspA_knMaxPixelValue,
                          (float)pFrame[DspA_knBluePixelCompIndex] /
                          (float)DspA_knMaxPixelValue );

          /* make it redder */
          xColr_HilightComponent( &color, xColr_tComponent_Red );

          /* put it back */
          pFrame[DspA_knRedPixelCompIndex]   =
            (GLubyte)(color.mfRed * (float)DspA_knMaxPixelValue);
          pFrame[DspA_knGreenPixelCompIndex] =
            (GLubyte)(color.mfGreen * (float)DspA_knMaxPixelValue);
          pFrame[DspA_knBluePixelCompIndex]  =
            (GLubyte)(color.mfBlue * (float)DspA_knMaxPixelValue);
          pFrame[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;

          break;
        }
      }

      /* advance our pointer. */
      pFrame += DspA_knNumBytesPerPixel;
    }
  }

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawLinesToFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawLines_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  xPoint2f  drawPoint1;
  xPoint2f  drawPoint2;

  /* If we don't have good vertices, bail. */
  if ( this->mLineVertex1.mnX < 0 || this->mLineVertex1.mnY < 0 ||
       this->mLineVertex2.mnX < 0 || this->mLineVertex2.mnY < 0 )
    goto cleanup;

  /* convert to zoomed coords. */
  drawPoint1.mfX = 
    ((float)this->mnZoomLevel * 
     ((float)(this->mLineVertex1.mnX) - 
      xVoxl_GetFloatX(this->mpZoomCenter))) + 
    (float)(this->mnVolumeSizeX/2.0);
  drawPoint1.mfY = 
    ((float)this->mnZoomLevel * 
     ((float)(this->mLineVertex1.mnY) - 
      xVoxl_GetFloatY(this->mpZoomCenter))) + 
    (float)(this->mnVolumeSizeY/2.0);

  drawPoint2.mfX = 
    ((float)this->mnZoomLevel * 
     ((float)(this->mLineVertex2.mnX) - 
      xVoxl_GetFloatX(this->mpZoomCenter))) + 
    (float)(this->mnVolumeSizeX/2.0);
  drawPoint2.mfY = 
    ((float)this->mnZoomLevel * 
     ((float)(this->mLineVertex2.mnY) - 
      xVoxl_GetFloatY(this->mpZoomCenter))) + 
    (float)(this->mnVolumeSizeY/2.0);

  /* y flip */
  drawPoint1.mfY = GLDRAW_Y_FLIP(drawPoint1.mfY);
  drawPoint2.mfY = GLDRAW_Y_FLIP(drawPoint2.mfY);

  DspA_SetUpOpenGLPort_( this );

  /* draw a green line */
  glColor3f ( 0.0, 1.0, 0.0 );

  glBegin( GL_LINE_STRIP );
  glVertex2f( drawPoint1.mfX, drawPoint1.mfY );
  glVertex2f( drawPoint2.mfX, drawPoint2.mfY );
  glEnd ();

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawLines_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawFunctionalOverlay_ ( tkmDisplayAreaRef this ) {

  DspA_tErr             eResult      = DspA_tErr_NoErr;
  FunV_tErr             eFunctional  = FunV_tErr_NoError;
  FunV_tFunctionalValue min          = 0;
  FunV_tFunctionalValue mid          = 0;
  FunV_tFunctionalValue slope        = 0;
  FunV_tFunctionalValue max          = 0;
  int                   nY           = 0;
  FunV_tFunctionalValue funcValue    = 0.0;
  float                 absFuncValue = 0.0;
  xColor3f              color;
  xColor3f              newColor     = {0,0,0};
  int                   numDecimals  = 0;
  char                  sValue[1024] = "";
  char                  sFormat[1024]= "";
  int                   nChar        = 0;
  float                 funcPerLine  = 0;
  int                   posMinLine   = 0;
  int                   posMidLine   = 0;
  int                   negMinLine   = 0;
  int                   negMidLine   = 0;

  DspA_SetUpOpenGLPort_( this );


  /* draw the color scale bar. get threshold max and min. */
  eFunctional = FunV_GetThreshold( this->mpFunctionalVolume,
                                   &min, &mid, &slope );
  if ( FunV_tErr_NoError != eFunctional ) {
    eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
    goto error;
  }
  eFunctional = FunV_GetThresholdMax( this->mpFunctionalVolume, &max );
  if ( FunV_tErr_NoError != eFunctional ) {
    eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
    goto error;
  }

  funcPerLine = (float)(max*2.0) / (float)this->mnVolumeSizeY;
  posMinLine = this->mnVolumeSizeY/2 + (min / funcPerLine);
  posMidLine = this->mnVolumeSizeY/2 + (mid / funcPerLine);
  negMinLine = this->mnVolumeSizeY/2 - (min / funcPerLine);
  negMidLine = this->mnVolumeSizeY/2 - (mid / funcPerLine);

  color.mfRed = color.mfGreen = color.mfBlue;
  for ( nY = 0; nY < this->mnVolumeSizeY; nY++ ) {

    /* get an interpolated value within the range of -max to +max
       determined by the y value */
    funcValue = (FunV_tFunctionalValue)
                ( (float)(this->mnVolumeSizeY - nY) * 
                  (float)((max*2.0)/(float)this->mnVolumeSizeY) - max );

    /* get the functional color for this value */
    eFunctional = FunV_GetColorForValue( this->mpFunctionalVolume,
                                         funcValue, &color, &newColor );
    if ( FunV_tErr_NoError != eFunctional ) {
      eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
      goto error;
    }

    /* draw a colored line on the right side of the screen. */
    glColor3f ( newColor.mfRed, newColor.mfGreen, newColor.mfBlue );
    glBegin( GL_POLYGON );
    glVertex2i( this->mnVolumeSizeX - 10, GLDRAW_Y_FLIP(nY) );
    glVertex2i( this->mnVolumeSizeX - 1, GLDRAW_Y_FLIP(nY) );
    glVertex2i( this->mnVolumeSizeX - 1, GLDRAW_Y_FLIP(nY)+1 );
    glVertex2i( this->mnVolumeSizeX - 10, GLDRAW_Y_FLIP(nY)+1 );
    glEnd ();

    /* Draw a value marker at the top and every 50 pixels. */
    if ( nY == 0 ||
         nY == this->mnVolumeSizeY-1 ||
         nY == posMinLine ||
         nY == posMidLine ||
         nY == negMinLine ||
         nY == negMidLine ) {

      /* Draw an extra little line to our label. */
      glColor3f( 1.0, 1.0, 1.0 );
      glBegin (GL_LINES);
      glVertex2i( this->mnVolumeSizeX - 1, GLDRAW_Y_FLIP(nY) );
      glVertex2i( this->mnVolumeSizeX - 12, GLDRAW_Y_FLIP(nY) );
      glEnd ();

      /* Calc how many decimals our label should have. */
      absFuncValue = fabs( funcValue );
      if (absFuncValue > 1 || absFuncValue == 0) numDecimals = 2 ;
      else if (absFuncValue > .1) numDecimals = 3 ;
      else if (absFuncValue > .01) numDecimals = 4 ;
      else if (absFuncValue > 0.001) numDecimals = 5 ;
      else if (absFuncValue > 0.0001) numDecimals = 6 ;
      else if (absFuncValue > 0.00001) numDecimals = 7 ;
      else if (absFuncValue > 0.000001) numDecimals = 8 ;
      else if (absFuncValue > 0.0000001) numDecimals = 9 ;
      else numDecimals = 10 ;

      /* Create the label string. */
      sprintf( sFormat, "%%2.%df", numDecimals) ;
      sprintf( sValue, sFormat, funcValue );

      /* Draw it. */
      glColor3f( 1.0, 1.0, 1.0 );
      glRasterPos2i( this->mnVolumeSizeX - 10 - (strlen(sValue) * 4) - 2,
                     (nY==0 ? GLDRAW_Y_FLIP(nY+8) : GLDRAW_Y_FLIP(nY)) );
      for ( nChar = 0; nChar < strlen(sValue); nChar++ ) {
        glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sValue[nChar] );
      }
    }
  }

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawFunctionalOverlay_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

void DspA_DrawVerticalArrow_ ( xPoint2nRef iStart,
                               int         inLength,
                               char*       isLabel ) {

  int nChar = 0;
  int nStrLength = 0;
  int nHeadOffset = 0;

  if ( inLength > 0 ) {
    nHeadOffset = 3;
  } else {
    nHeadOffset = -3;
  }

  glBegin( GL_LINE_STRIP );
  glVertex2d( iStart->mnX, iStart->mnY );
  glVertex2d( iStart->mnX, iStart->mnY + inLength );
  glVertex2d( iStart->mnX - nHeadOffset,iStart->mnY + inLength - nHeadOffset );
  glVertex2d( iStart->mnX, iStart->mnY + inLength );
  glVertex2d( iStart->mnX + nHeadOffset,iStart->mnY + inLength - nHeadOffset );
  glEnd();

  glRasterPos2i( iStart->mnX + 2, iStart->mnY + 7 );
  nStrLength = strlen( isLabel );
  for ( nChar = 0; nChar < nStrLength; nChar++ )
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, isLabel[nChar] );

}

void DspA_DrawHorizontalArrow_ ( xPoint2nRef iStart,
                                 int         inLength,
                                 char*       isLabel ) {

  int nChar = 0;
  int nStrLength = 0;
  int nHeadOffset = 0;

  if ( inLength > 0 ) {
    nHeadOffset = 3;
  } else {
    nHeadOffset = -3;
  }

  glBegin( GL_LINE_STRIP );
  glVertex2d( iStart->mnX, iStart->mnY );
  glVertex2d( iStart->mnX + inLength, iStart->mnY );
  glVertex2d( iStart->mnX + inLength - nHeadOffset,iStart->mnY - nHeadOffset );
  glVertex2d( iStart->mnX + inLength, iStart->mnY );
  glVertex2d( iStart->mnX + inLength - nHeadOffset,iStart->mnY + nHeadOffset );
  glEnd();

  glRasterPos2i( iStart->mnX + 2, iStart->mnY + 7 );
  nStrLength = strlen( isLabel );
  for ( nChar = 0; nChar < nStrLength; nChar++ )
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, isLabel[nChar] );

}



DspA_tErr DspA_BuildCurrentFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr       eResult  = DspA_tErr_NoErr;
  xPoint2n        bufferPt = {0, 0};
  GLubyte*        pDest    = NULL;
  tkm_tVolumeType volume   = tkm_tVolumeType_Main;
  xVoxel          anaIdx;
  xColor3n        color    = {0,0,0};
  int             yMin     = 0;
  int             yMax     = 0;
  int             yInc     = 0;

  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;

  DisableDebuggingOutput;

  /* Decide which volume we're looking at. */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
    volume = tkm_tVolumeType_Aux;
  } else {
    volume = tkm_tVolumeType_Main;
  }

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY;
    yMax = 0;
    yInc = -1;
  }

  /* If the anatomical flag is not on, just paint the whole thing
     black. */
  if ( !this->mabDisplayFlags[DspA_tDisplayFlag_Anatomical] ) {

    for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
      for ( bufferPt.mnX = 0;
            bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

        /* set the pixel */
        pDest[DspA_knRedPixelCompIndex]   = 0;
        pDest[DspA_knGreenPixelCompIndex] = 0;
        pDest[DspA_knBluePixelCompIndex]  = 0;
        pDest[DspA_knAlphaPixelCompIndex] = 1;

        /* advance our pointer. */
        pDest += DspA_knNumBytesPerPixel;
      }
    }

  } else {

    /* If the maximum intensity projection flag is on, loop through
       getting the max int value. Otherwise loop through getting the
       normal color value. */
    if ( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] ) {

      for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
        for ( bufferPt.mnX = 0;
              bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

          /* get a volume voxel.*/

          eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
          if ( DspA_tErr_NoErr != eResult )
            goto error;

          Volm_GetMaxIntColorAtIdx( this->mpVolume[volume], &anaIdx,
                                    this->mOrientation, &color );

          /* set the pixel */
          pDest[DspA_knRedPixelCompIndex]   = (GLubyte)color.mnRed;
          pDest[DspA_knGreenPixelCompIndex] = (GLubyte)color.mnGreen;
          pDest[DspA_knBluePixelCompIndex]  = (GLubyte)color.mnBlue;
          pDest[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;

          /* advance our pointer. */
          pDest += DspA_knNumBytesPerPixel;
        }
      }

    } else {

      /* Just loop through getting the anatomical color at this index. */
      for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
        for ( bufferPt.mnX = 0;
              bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

          /* get a volume voxel.*/
          eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
          if ( DspA_tErr_NoErr != eResult )
            goto error;

          Volm_GetIntColorAtIdx( this->mpVolume[volume],
                                 &anaIdx, &color );


          /* set the pixel */
          pDest[DspA_knRedPixelCompIndex]   = (GLubyte)color.mnRed;
          pDest[DspA_knGreenPixelCompIndex] = (GLubyte)color.mnGreen;
          pDest[DspA_knBluePixelCompIndex]  = (GLubyte)color.mnBlue;
          pDest[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;

          /* advance our pointer. */
          pDest += DspA_knNumBytesPerPixel;
        }
      }
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_BuildCurrentFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  EnableDebuggingOutput;

  return eResult;
}

DspA_tErr DspA_DrawSegmentationOverlayToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr  eResult     = DspA_tErr_NoErr;
  Volm_tErr  eVolume     = Volm_tErr_NoErr;
  int        eCTAB       = NO_ERROR;
  tkm_tSegType seg       = tkm_tSegType_Main;
  xPoint2n   bufferPt    = {0, 0};
  GLubyte*   pDest       = NULL;
  xVoxel     anaIdx;
  xColor3n   color       = {0,0,0};
  float      fRed        = 0;
  float      fGreen      = 0;
  float      fBlue       = 0;
  xColor3n   newColor    = {0,0,0};
  int        yMin        = 0;
  int        yMax        = 0;
  int        yInc        = 0;
  float      value;
  float      fAlpha      = 0;
  float      fFinalAlpha = 0;

  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;

  DisableDebuggingOutput;

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY;
    yMax = 0;
    yInc = -1;
  }

  if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume] ) {
    seg = tkm_tSegType_Aux;
  } else {
    seg = tkm_tSegType_Main;
  }

  for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
    for ( bufferPt.mnX = 0;
          bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

      /* Get a volume voxel. */
      eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
      if ( DspA_tErr_NoErr != eResult )
        goto error;

      /* Get the value in the volume. */
      eVolume = Volm_GetValueAtIdxUnsafe( this->mSegmentationVolume[seg],
                                          &anaIdx, &value );
      if ( Volm_tErr_NoErr != eVolume )
        goto error;

      if ( 0 != value ) {

        /* Get a color from the look up table. If we got an error,
           draw this voxel as red with an alpha 1, overriding the
           user's overlay's alpha. Otherwise just use the user's
           alpha. */
        eCTAB = CTABrgbaAtIndexf( this->mSegmentationColorTable[seg],
                                  value, &fRed, &fGreen, &fBlue, &fAlpha );
        if ( NO_ERROR == eCTAB ) {
          fFinalAlpha = fAlpha * this->mfSegmentationAlpha;
        } else {
          fFinalAlpha = 1.0;
          fRed   = 1.0;
          fGreen = 0;
          fBlue  = 0;
        }

        /* Get the color at the dest. */
        color.mnRed   = pDest[DspA_knRedPixelCompIndex];
        color.mnGreen = pDest[DspA_knGreenPixelCompIndex];
        color.mnBlue  = pDest[DspA_knBluePixelCompIndex];

        newColor.mnRed =
          ((float)color.mnRed * (1.0 - fFinalAlpha)) +
          (fRed * fFinalAlpha * 255.0);
        newColor.mnGreen =
          ((float)color.mnGreen * (1.0 - fFinalAlpha)) +
          (fGreen * fFinalAlpha * 255.0);
        newColor.mnBlue =
          ((float)color.mnBlue * (1.0 - fFinalAlpha)) +
          (fBlue * fFinalAlpha * 255.0);

        /* set the pixel */
        pDest[DspA_knRedPixelCompIndex]   = (GLubyte)newColor.mnRed;
        pDest[DspA_knGreenPixelCompIndex] = (GLubyte)newColor.mnGreen;
        pDest[DspA_knBluePixelCompIndex]  = (GLubyte)newColor.mnBlue;
      }

      /* advance our pointer. */
      pDest += DspA_knNumBytesPerPixel;
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawSegmentationOverlayToFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  EnableDebuggingOutput;

  return eResult;
}

DspA_tErr DspA_DrawDTIOverlayToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr  eResult     = DspA_tErr_NoErr;
  xPoint2n   bufferPt    = {0, 0};
  GLubyte*   pDest       = NULL;
  xVoxel     anaIdx;
  xColr_tComponent comp;
  int   color;
  int   overlayColor;
  int   newColor;
  int        yMin        = 0;
  int        yMax        = 0;
  int        yInc        = 0;
  float      value       = 0;

  DisableDebuggingOutput;

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY;
    yMax = 0;
    yInc = -1;
  }

  for ( comp = xColr_tComponent_Red; comp <= xColr_tComponent_Blue; comp++ ) {

    /* get a ptr to the frame buffer. depending on the component, we
       want to stagger it so first we process all the red components,
       then the green, etc.*/
    pDest = this->mpFrameBuffer;
    switch ( comp ) {
    case xColr_tComponent_Red:
      pDest += DspA_knRedPixelCompIndex;
      break;
    case xColr_tComponent_Green:
      pDest += DspA_knGreenPixelCompIndex;
      break;
    case xColr_tComponent_Blue:
      pDest += DspA_knBluePixelCompIndex;
      break;
    default:
      break;
    }

    for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
      for ( bufferPt.mnX = 0;
            bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

        /* Get a volume voxel. */
        eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* Map to the colors. Get the value of the corresponding
           axis. x is in frame 0, y in 1, and z in two. */
        switch ( this->maDTIAxisForComponent[comp] ) {
        case tAxis_X:
          Volm_GetValueAtIdxFrameUnsafe( this->mpDTIVolume,
                                         &anaIdx, 0, &value );
          break;
        case tAxis_Y:
          Volm_GetValueAtIdxFrameUnsafe( this->mpDTIVolume,
                                         &anaIdx, 1, &value );
          break;
        case tAxis_Z:
          Volm_GetValueAtIdxFrameUnsafe( this->mpDTIVolume,
                                         &anaIdx, 2, &value );
          break;
        default:
          break;
        }

        if ( 0 != value ) {

          /* Convert float value to 0-255. Use fabs() here because the
             value in the volume could be negative. */
          overlayColor =
            (int)((float)fabs(value) * (float)DspA_knMaxPixelValue);

          /* Get the color at the dest. */
          color = *pDest;

          /* Do the blend. */
          newColor = ((float)color * (1.0 - this->mfDTIAlpha)) +
                     ((float)overlayColor * this->mfDTIAlpha);

          /* set the pixel */
          *pDest = (GLubyte)newColor;
        }

        /* advance our pointer. */
        pDest += DspA_knNumBytesPerPixel;
      }
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawDTIOverlayToFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  EnableDebuggingOutput;

  return eResult;
}

DspA_tErr DspA_DrawUndoableVoxelsOverlayToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr  eResult     = DspA_tErr_NoErr;
  xPoint2n   bufferPt    = {0, 0};
  GLubyte*   pDest       = NULL;
  xVoxel     anaIdx;
  xVoxel     MRIIdx;
  xColor3n   color;
  int        yMin        = 0;
  int        yMax        = 0;
  int        yInc        = 0;

  DisableDebuggingOutput;

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY;
    yMax = 0;
    yInc = -1;
  }

  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;

  for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
    for ( bufferPt.mnX = 0;
          bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

      /* Get a volume voxel. */
      eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
      if ( DspA_tErr_NoErr != eResult )
        goto error;

      Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                               &anaIdx, &MRIIdx );

      /* Is this voxel is in the undo volume? */
      if ( tkm_IsMRIIdxInUndoVolume ( &MRIIdx ) ) {

        /* Get the color at the dest. */
        color.mnRed   = pDest[DspA_knRedPixelCompIndex];
        color.mnGreen = pDest[DspA_knGreenPixelCompIndex];
        color.mnBlue  = pDest[DspA_knBluePixelCompIndex];

        /* Highlight the blue component. */
        color.mnBlue = MIN( DspA_knMaxPixelValue,
                            color.mnBlue + DspA_knMaxPixelValue / 2 );

        /* If the other components are two bright to let the
           highlight show, lower them a bit. */
        if ( color.mnBlue - color.mnGreen < DspA_knMaxPixelValue / 2 ) {
          color.mnGreen = color.mnBlue - DspA_knMaxPixelValue / 2;
        }
        if ( color.mnBlue - color.mnRed < DspA_knMaxPixelValue / 2 ) {
          color.mnRed = color.mnBlue - DspA_knMaxPixelValue / 2;
        }

        /* set the pixel */
        pDest[DspA_knRedPixelCompIndex]   = (GLubyte)color.mnRed;
        pDest[DspA_knGreenPixelCompIndex] = (GLubyte)color.mnGreen;
        pDest[DspA_knBluePixelCompIndex]  = (GLubyte)color.mnBlue;
      }

      /* advance our pointer. */
      pDest += DspA_knNumBytesPerPixel;
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawUndoableVoxelsOverlayToFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  EnableDebuggingOutput;

  return eResult;
}

DspA_tErr DspA_DrawFunctionalOverlayToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr             eResult     = DspA_tErr_NoErr;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  xPoint2n              bufferPt    = {0,0};
  FunV_tFunctionalValue funcValue   = 0.0;
  xColor3f              color;
  xColor3f              newColor    = {0,0,0};
  GLubyte*              pDest       = NULL;
  xVoxel                anaIdx;
  xVoxel                mriIdx;
  float                 anaValue    = 0;
  int                   yMin        = 0;
  int                   yMax        = 0;
  int                   yInc        = 0;

  DisableDebuggingOutput;

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY;
    yMax = 0;
    yInc = -1;
  }

  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;

  for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
    for ( bufferPt.mnX = 0;
          bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

      /* Get a volume voxel. */
      eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
      if ( DspA_tErr_NoErr != eResult )
        goto error;

      /* get a functional value. */
      Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                               &anaIdx, &mriIdx );
      eFunctional = FunV_GetValueAtMRIIdx( this->mpFunctionalVolume,
                                           &mriIdx, TRUE, &funcValue );

      /* If we're masking to the aux volume, check its value first. */
      if (this->mabDisplayFlags[DspA_tDisplayFlag_MaskFunctionalOverlayToAux]) 
      {

        Volm_GetValueAtIdx( this->mpVolume[tkm_tVolumeType_Aux],
                            &anaIdx, &anaValue );
      }

      /* if it was a valid voxel and if we're masking to aux the aux
      value is not 0 */
      if ( FunV_tErr_NoError == eFunctional &&
       ((this->mabDisplayFlags[DspA_tDisplayFlag_MaskFunctionalOverlayToAux] &&
         anaValue != 0)   ||
        !this->mabDisplayFlags[DspA_tDisplayFlag_MaskFunctionalOverlayToAux])) 
      {

        /* Get the current color and convert it to float. */
        color.mfRed   = (float)pDest[DspA_knRedPixelCompIndex] /
                        (float)DspA_knMaxPixelValue;
        color.mfGreen = (float)pDest[DspA_knGreenPixelCompIndex] /
                        (float)DspA_knMaxPixelValue;
        color.mfBlue  = (float)pDest[DspA_knBluePixelCompIndex] /
                        (float)DspA_knMaxPixelValue;


        /* get a color value. this function automatically does the
           appropriate blending. */
        eFunctional = FunV_GetColorForValue ( this->mpFunctionalVolume,
                                              funcValue, &color, &newColor );

        /* Do the blend and set the color back in the buffer,
           converting back to int. */
        pDest[DspA_knRedPixelCompIndex]   =
          (GLubyte)(((color.mfRed * (1.0 - this->mfFuncOverlayAlpha))+
                     (newColor.mfRed * this->mfFuncOverlayAlpha))*
                    DspA_knMaxPixelValue);
        pDest[DspA_knGreenPixelCompIndex] =
          (GLubyte)(((color.mfGreen * (1.0 - this->mfFuncOverlayAlpha))+
                     (newColor.mfGreen * this->mfFuncOverlayAlpha))*
                    DspA_knMaxPixelValue);
        pDest[DspA_knBluePixelCompIndex]  =
          (GLubyte)(((color.mfBlue * (1.0 - this->mfFuncOverlayAlpha))+
                     (newColor.mfBlue * this->mfFuncOverlayAlpha))*
                    DspA_knMaxPixelValue);

      } else if ( FunV_tErr_InvalidMRIIdx == eFunctional &&
                  this->mabDisplayFlags[DspA_tDisplayFlag_MaskToFunctionalOverlay]) {

        pDest[DspA_knRedPixelCompIndex]   = (GLubyte)0;
        pDest[DspA_knGreenPixelCompIndex] = (GLubyte)0;
        pDest[DspA_knBluePixelCompIndex]  = (GLubyte)0;
      }

      /* advance our pointer. */
      pDest += DspA_knNumBytesPerPixel;
    }
  }

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawFunctionalOverlayToFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_DrawControlPoints_ ( tkmDisplayAreaRef this ) {

  DspA_tErr    eResult    = DspA_tErr_NoErr;
  xListRef     list       = NULL;
  xList_tErr   eList      = xList_tErr_NoErr;
  x3Lst_tErr   e3DList    = x3Lst_tErr_NoErr;
  xVoxelRef    MRIIdx  = NULL;
  xVoxel       anaIdx;
  xPoint2n     bufferPt   = {0,0};
  float        faColor[3] = {0, 0, 0};

  /* decide which list we want out of the space. The space is in MRI
     Idx space, though, so */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    e3DList = x3Lst_GetItemsInZPlane
      ( this->mpControlPoints,
        DspA_GetCurrentSliceNumberInMRIIdx_(this),
        &list );
    break;
  case mri_tOrientation_Sagittal:
    e3DList = x3Lst_GetItemsInXPlane
      ( this->mpControlPoints,
        DspA_GetCurrentSliceNumberInMRIIdx_(this),
        &list );
    break;
  case mri_tOrientation_Horizontal:
    e3DList = x3Lst_GetItemsInYPlane
      ( this->mpControlPoints,
        DspA_GetCurrentSliceNumberInMRIIdx_(this),
        &list );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* check for error. */
  if ( x3Lst_tErr_NoErr != e3DList )
    goto error;

  /* if we got a list... */
  if ( NULL != list ) {

    /* set color */
    faColor[0] = 0;
    faColor[1] = 1;
    faColor[2] = 0;

    /* traverse the list */
    eList = xList_ResetPosition( list );
    void* pvoid = (void*) &MRIIdx;
    while ( (eList = xList_NextFromPos( list, (void**)pvoid ))
            != xList_tErr_EndOfList ) {

      if ( MRIIdx ) {

        /* Since our drawing function works in screen space, we need
           to convert it to the anaIdx space (which is screen
           space). */
        Volm_ConvertMRIIdxToIdx( this->mpVolume[tkm_tVolumeType_Main],
                                 MRIIdx, &anaIdx );

        /* convert the control point to be in the middle of voxel */
        xVoxl_SetFloatX( &anaIdx, floor(xVoxl_GetFloatX(&anaIdx)+0.5) );
        xVoxl_SetFloatY( &anaIdx, floor(xVoxl_GetFloatY(&anaIdx)+0.5) );
        xVoxl_SetFloatZ( &anaIdx, floor(xVoxl_GetFloatZ(&anaIdx)+0.5) );

        /* convert to buffer point. */
        eResult = DspA_ConvertVolumeToBuffer_( this, &anaIdx, &bufferPt );
        if ( DspA_tErr_NoErr != eResult )
          goto error;

        /* draw a point here. */
        DspA_DrawMarker_( this, DspA_tMarker_Crosshair, faColor, &bufferPt,
                          DspA_knControlPointCrosshairSize );
      }
    }

    if ( eList != xList_tErr_EndOfList )
      goto error;
  }

  goto cleanup;

error:

  if ( x3Lst_tErr_NoErr != e3DList )
    eResult = DspA_tErr_ErrorAccessingControlPoints;

  if ( xList_tErr_NoErr != eList )
    eResult = DspA_tErr_ErrorAccessingControlPoints;

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawControlPoints_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawVectorField_ ( tkmDisplayAreaRef this ) {

  DspA_tErr    eResult    = DspA_tErr_NoErr;
  float        faColor[3] = {0, 0, 0};
  xVoxel       start;
  xVoxel       direction;

  /* decide which list we want out of the space. */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    break;
  case mri_tOrientation_Sagittal:
    break;
  case mri_tOrientation_Horizontal:
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* set color */
  faColor[0] = 0;
  faColor[1] = 1;
  faColor[2] = 0;

  /* this was some test code. it still doesn't work quite right. */
  xVoxl_Set( &start, 120.5, 120.5, 120.5 );
  xVoxl_SetFloat( &direction, .5, .3, .1 );
  while ( xVoxl_IncrementWithMinUntilLimit( &start, 120.5, 121.5 ) ) {
    DspA_DrawVector_( this, faColor, &start, &direction );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawVectorField_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawSelectionToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr     eResult  = DspA_tErr_NoErr;
  xPoint2n      bufferPt = {0,0};
  GLubyte*      pFrame   = NULL;
  xColor3f      color;
  xVoxel        anaIdx;
  float         value    = 0;
  int           yMin     = 0;
  int           yMax     = 0;
  int           yInc     = 0;

  if ( !tkm_IsSelectionPresent() ) {
    goto cleanup;
  }

  /* get a ptr to the frame buffer. */
  pFrame = this->mpFrameBuffer;

  if ( mri_tOrientation_Horizontal == this->mOrientation ) {
    yMin = 0;
    yMax = this->mnVolumeSizeY;
    yInc = 1;
  } else {
    yMin = this->mnVolumeSizeY;
    yMax = 0;
    yInc = -1;
  }

  /* Just loop through getting the selection value at this index. */
  for ( bufferPt.mnY = yMin; bufferPt.mnY != yMax; bufferPt.mnY += yInc) {
    for ( bufferPt.mnX = 0;
          bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX ++ ) {

      /* get a volume voxel.*/
      eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
      if ( DspA_tErr_NoErr != eResult )
        goto error;

      if ( xVoxl_GetX(&anaIdx) != 0 &&
           xVoxl_GetY(&anaIdx) != 0 &&
           xVoxl_GetZ(&anaIdx) != 0 ) {

        Volm_GetValueAtIdxUnsafe( this->mpSelection, &anaIdx, &value );
        if ( value > 0 ) {

          /* get the current color in the buffer */
          xColr_SetFloat( &color, (float)pFrame[DspA_knRedPixelCompIndex] /
                          (float)DspA_knMaxPixelValue,
                          (float)pFrame[DspA_knGreenPixelCompIndex] /
                          (float)DspA_knMaxPixelValue,
                          (float)pFrame[DspA_knBluePixelCompIndex] /
                          (float)DspA_knMaxPixelValue );

          // LABEL drawing code - make it red instead of green (BRF)
          /* make it greener */
          xColr_HilightComponent( &color, xColr_tComponent_Green );
          xColr_HilightComponent( &color, xColr_tComponent_Red );

          /* put it back */
          pFrame[DspA_knRedPixelCompIndex]   =
            (GLubyte)(color.mfRed * (float)DspA_knMaxPixelValue);
          pFrame[DspA_knGreenPixelCompIndex] =
            (GLubyte)(color.mfGreen * (float)DspA_knMaxPixelValue);
          pFrame[DspA_knBluePixelCompIndex]  =
            (GLubyte)(color.mfBlue * (float)DspA_knMaxPixelValue);
          pFrame[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;
        }
      }

      /* advance our pointer. */
      pFrame += DspA_knNumBytesPerPixel;
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawSelectionToFrame_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawHeadPoints_ ( tkmDisplayAreaRef this ) {

  DspA_tErr            eResult     = DspA_tErr_NoErr;
  HPtL_tErr            eHeadPtList = HPtL_tErr_NoErr;
  HPtL_tHeadPointRef   pHeadPt     = NULL;
  xPoint2n             bufferPt    = {0,0};
  float                faColor[3]  = {0, 0, 0};
  xVoxel               anaIdx;
  int                  nSize       = 0;
  HPtL_tIterationPlane plane       = HPtL_tIterationPlane_X;

  /* find our iteration orientation */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] ) {
    plane = HPtL_tIterationPlane_All;
  } else {
    switch ( this->mOrientation ) {
    case mri_tOrientation_Horizontal:
      plane = HPtL_tIterationPlane_Y;
      break;
    case mri_tOrientation_Coronal:
      plane = HPtL_tIterationPlane_Z;
      break;
    case mri_tOrientation_Sagittal:
      plane = HPtL_tIterationPlane_X;
      break;
    default:
      eResult = DspA_tErr_InvalidOrientation;
      goto error;
      break;
    }
  }

  /* reset the iterator */
  eHeadPtList = HPtL_ResetIterator( this->mHeadPoints,
                                    plane,
                                    DspA_GetCurrentSliceNumber_( this ),
                                    1.0 );
  if ( HPtL_tErr_NoErr != eHeadPtList ) {
    goto error;
  }

  while ( (eHeadPtList = HPtL_NextPoint( this->mHeadPoints, &pHeadPt ))
          == HPtL_tErr_NoErr ) {

    /* color is red if hilited, purple if cardinal, green if not. */
    if ( pHeadPt == this->mpSelectedHeadPoint ) {
      faColor[0] = 1;
      faColor[1] = 0;
      faColor[2] = 0;
    } else if ( strcmp( pHeadPt->msLabel, "cardinal" ) == 0 ) {
      faColor[0] = 1;
      faColor[1] = 0;
      faColor[2] = 1;
    } else {
      faColor[0] = 0;
      faColor[1] = 1;
      faColor[2] = 0;
    }

    /* cardinal points are big */
    if ( strcmp( pHeadPt->msLabel, "cardinal" ) == 0 ) {
      nSize = DspA_knControlPointCrosshairSize * 2;
    } else {
      nSize = DspA_knControlPointCrosshairSize;
    }

    /* get our pt */
    xVoxl_Copy( &anaIdx, &(pHeadPt->mClientPoint) );

    /* convert to buffer point. */
    DisableDebuggingOutput;
    eResult = DspA_ConvertVolumeToBuffer_ ( this, &anaIdx, &bufferPt );
    EnableDebuggingOutput;
    if ( DspA_tErr_NoErr != eResult )
      continue;

    /* draw a point here. */
    DspA_DrawMarker_( this, DspA_tMarker_Diamond,
                      faColor, &bufferPt, nSize );
  }

  if ( eHeadPtList != HPtL_tErr_LastPoint )
    goto error;

  /* clean up error codes */
  eHeadPtList = HPtL_tErr_NoErr;
  eResult = DspA_tErr_NoErr;

  goto cleanup;

error:

  if ( HPtL_tErr_NoErr != eHeadPtList ) {
    eResult = DspA_tErr_ErrorAccessingHeadPointList;
    DebugPrint( ("HPtL error %d in DspA_DrawHeadPoints_: %s\n",
                 eHeadPtList, HPtL_GetErrorString(eHeadPtList) ) );
  }

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawHeadPoints_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_BuildSurfaceDrawLists_ ( tkmDisplayAreaRef this,
                                        tkm_tSurfaceType  iSurface) {

  DspA_tErr           eResult              = DspA_tErr_NoErr;
  Surf_tErr           eSurface             = Surf_tErr_NoErr;
  xGArr_tErr          eList                = xGArr_tErr_NoErr;
  xGrowableArrayRef   list                 = NULL;
  int                 nSlice               = 0;
  Surf_tVertexSet     vertexSet              = Surf_tVertexSet_Main;
  xPoint2f            zeroPoint            = {0,0};
  xVoxel              curPlane;
  DspA_tSurfaceListNode drawListNode;
  xVoxel              anaVertex;
  xVoxel              anaNeighborVertex;
  xVoxel              normAnaVertex;
  xVoxel              normAnaNeighborVertex;
  int                 nVertexIndex         = -1;
  int                 nNeighborVertexIndex = -1;
  xPoint2f            intersectionPt       = {0, 0};
  xPoint2f            interpIntersectionPt = {0, 0};
  tBoolean            bPointsOnThisFace    = FALSE;
  int                 annotation;
  int                 annotationRed;
  int                 annotationGreen;
  int                 annotationBlue;

  normAnaVertex.mfX=0;
  normAnaVertex.mfY=0;
  normAnaVertex.mfZ=0;
  normAnaNeighborVertex.mfX=0;
  normAnaNeighborVertex.mfY=0;
  normAnaNeighborVertex.mfZ=0;

  /* get the current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* set up for gl drawing */
  DspA_SetUpOpenGLPort_( this );

  /* for each vertex type... */
  for ( vertexSet = Surf_tVertexSet_Main;
        vertexSet < Surf_knNumVertexSets; vertexSet++ ) {

    /* only build if this surface is being displayed */
    if ( !this->mabDisplayFlags[vertexSet + DspA_tDisplayFlag_MainSurface] ) {
      continue;
    }

    /* check to see if we already have this list. if so, don't build
       it unless the slice changed flag is set.. */
    list = DspA_GetSurfaceList_( this, iSurface, this->mOrientation, vertexSet,
                                 DspA_GetCurrentSliceNumber_(this) );
    if ( NULL != list
         // rkt: why was this here? to update draw list when new surface
         // is loaded? but DspA_PurgeSurfaceLists_ should just be called
         // to do that.
         || !this->mbSliceChanged ) {
      continue;
    }

    /* make a new list. */
    DspA_NewSurfaceList_( this, iSurface, this->mOrientation, vertexSet,
                          DspA_GetCurrentSliceNumber_(this) );

    /* get the list. */
    list = DspA_GetSurfaceList_( this, iSurface, this->mOrientation, vertexSet,
                                 DspA_GetCurrentSliceNumber_(this) );
    if ( NULL == list ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }

    //    xUtil_StartTimer ();

    /* make a voxel with our current slice in it and
       set the iterator to our view */
    DspA_ConvertPlaneToVolume_ ( this, &zeroPoint, nSlice, this->mOrientation,
                                 &curPlane );
    eSurface = 
      Surf_SetIteratorPosition( this->mpSurface[iSurface], &curPlane );
    if ( Surf_tErr_NoErr != eSurface )
      goto error;

    bPointsOnThisFace = FALSE;

    /* while we have vertices to check.. */
    while ( (eSurface = 
             Surf_GetNextAndNeighborVertex( this->mpSurface[iSurface],
                                            vertexSet,
                                            &anaVertex,
                                            &nVertexIndex,
                                            &anaNeighborVertex,
                                            &nNeighborVertexIndex ))
            != Surf_tErr_LastFace ) {

      /* if the line between these two points intersects the
      current plane... */
      DspA_NormalizeVoxel_( &anaVertex,
                            this->mOrientation, &normAnaVertex );
      DspA_NormalizeVoxel_( &anaNeighborVertex,
                            this->mOrientation, &normAnaNeighborVertex );
      if ( xUtil_LineIntersectsPlane
           ( &normAnaVertex, &normAnaNeighborVertex,
             nSlice,
             &intersectionPt, &interpIntersectionPt)) {

        /* fill out a node. */
        drawListNode.mbVertex = TRUE;
        drawListNode.mnOriginalVertexIndex = nVertexIndex;
        drawListNode.mnNeighborVertexIndex = nNeighborVertexIndex;
        drawListNode.mIntersectionPoint = intersectionPt;
        drawListNode.mInterpIntersectionPoint = interpIntersectionPt;
        drawListNode.mOverrideColor = FALSE;

        /* Get the annotation value at this vertex. If it's non-zero,
           set the override flag to true and put the color of the
           annotation into the draw list node. */
        Surf_GetVertexAnnotationByIndex( this->mpSurface[iSurface],
                                         nVertexIndex, &annotation );
        if ( annotation ) {
          MRISAnnotToRGB( annotation,
                          annotationRed, annotationGreen, annotationBlue );
          xColr_SetFloat( &(drawListNode.mColor),
                          (float)annotationRed / 256.0,
                          (float)annotationGreen / 256.0,
                          (float)annotationBlue / 256.0 );

          drawListNode.mOverrideColor = TRUE;
        }

        /* the original vertex is just the (unnormalized) anatomical
           vertex. same with the neighbor vertex. */
        xVoxl_Copy( &drawListNode.mOriginalVertex, &anaVertex );
        xVoxl_Copy( &drawListNode.mNeighborVertex, &anaNeighborVertex );

        /* the interp vertex is a normalized point with x/y coords of
           the intersection point and the z coord of the cur slice,
           put through the DspA_UnnormalizeVoxel_ function. */
        xVoxl_SetFloat( &drawListNode.mInterpVertex,
                        interpIntersectionPt.mfX, interpIntersectionPt.mfY,
                        DspA_GetCurrentSliceNumber_( this ) );
        DspA_UnnormalizeVoxel_( &drawListNode.mInterpVertex,
                                this->mOrientation,
                                &drawListNode.mInterpVertex );

        /* add this node to the list. */
        eList = xGArr_Add( list, &drawListNode );
        if ( xGArr_tErr_NoErr != eList ) {
          DebugPrint( ("xGArr error %d in DspA_BuildSurfaceDrawLists_: %s\n",
                       eList, xGArr_GetErrorString( eList ) ) );
          eResult = DspA_tErr_ErrorAccessingSurfaceList;
          goto error;
        }

        bPointsOnThisFace = TRUE;
      }

      /* if we have a last vertex, and we drew points on this
      face add a face marker. */
      if ( eSurface == Surf_tErr_LastVertex
           && bPointsOnThisFace ) {

        /* fill out a node. */
        drawListNode.mbVertex = FALSE;

        eList = xGArr_Add( list, &drawListNode );
        if ( xGArr_tErr_NoErr != eList ) {
          DebugPrint( ("xGArr error %d in DspA_BuildSurfaceDrawLists_: %s\n",
                       eList, xGArr_GetErrorString( eList ) ) );
          eResult = DspA_tErr_ErrorAccessingSurfaceList;
          goto error;
        }

        bPointsOnThisFace = FALSE;
      }

    }

    /* if surface error is just on last face, clear it. */
    if ( Surf_tErr_LastFace == eSurface )
      eSurface = Surf_tErr_NoErr;

    /* check for other errors. */
    if ( Surf_tErr_NoErr != eSurface )
      goto error;

    //    xUtil_StopTimer( "build surface list" );

  } /* for each surface... */

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eSurface )
    eResult = DspA_tErr_ErrorAccessingSurface;

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_BuildSurfaceDrawLists_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_DrawMarker_ ( tkmDisplayAreaRef this,
                             DspA_tMarker      iMarker,
                             float*            ifaColor,
                             xPoint2nRef       ipWhere,
                             int               inSize ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  int       nWidth  = 0;
  int       nHeight = 0;

  DspA_SetUpOpenGLPort_( this );

  glLineWidth( 1 );

  switch ( iMarker ) {

  case DspA_tMarker_Crosshair:

    /* calculate width and height using scale */
    nWidth  = ((float)inSize / this->mfFrameBufferScaleX);
    nHeight = ((float)inSize / this->mfFrameBufferScaleY);

    nWidth = MAX( nWidth, 1 );
    nHeight = MAX( nHeight, 1 );

    DspA_SetUpOpenGLPort_( this );

    /* draw the crosshair */
    glColor3fv ( ifaColor );

    glBegin ( GL_LINES );
    glVertex2d ( ipWhere->mnX, ipWhere->mnY-nHeight );
    glVertex2d ( ipWhere->mnX, ipWhere->mnY+nHeight );
    glEnd ();

    glBegin ( GL_LINES );
    glVertex2d ( ipWhere->mnX-nWidth, ipWhere->mnY );
    glVertex2d ( ipWhere->mnX+nWidth, ipWhere->mnY );
    glEnd ();

    break;

  case DspA_tMarker_Diamond:

    /* calculate width and height using scale */
    nWidth  = ((float)inSize / this->mfFrameBufferScaleX);
    nHeight = ((float)inSize / this->mfFrameBufferScaleY);

    /* draw the diamond */
    glColor3fv ( ifaColor );

    glBegin ( GL_LINE_LOOP );
    glVertex2d ( ipWhere->mnX-nWidth, ipWhere->mnY );
    glVertex2d ( ipWhere->mnX, ipWhere->mnY-nHeight );
    glVertex2d ( ipWhere->mnX+nWidth, ipWhere->mnY );
    glVertex2d ( ipWhere->mnX, ipWhere->mnY+nHeight );
    glEnd ();

    break;

  default:
    break;
  }

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DspA_DrawMarker_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;

}

DspA_tErr DspA_DrawVector_ ( tkmDisplayAreaRef this,
                             float*            ifaColor,
                             xVoxelRef         ipVoxelStart,
                             xVoxelRef         ipVoxelDirection ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xPoint2f  bufferStart;
  xVoxel    voxelEnd;
  xPoint2f  bufferEnd;

  bufferStart.mfX=0;
  bufferStart.mfY=0;
  bufferEnd.mfX=0;
  bufferEnd.mfY=0;

  DspA_ConvertVolumeToBufferf_ ( this, ipVoxelStart, &bufferStart );
  xVoxl_SetFloatX( &voxelEnd,
                   xVoxl_GetFloatX( ipVoxelStart ) +
                   xVoxl_GetFloatX( ipVoxelDirection ));
  xVoxl_SetFloatY( &voxelEnd,
                   xVoxl_GetFloatY( ipVoxelStart ) +
                   xVoxl_GetFloatY( ipVoxelDirection ));
  xVoxl_SetFloatZ( &voxelEnd,
                   xVoxl_GetFloatZ( ipVoxelStart ) +
                   xVoxl_GetFloatZ( ipVoxelDirection ));
  DspA_ConvertVolumeToBufferf_ ( this, &voxelEnd, &bufferEnd );

  DspA_SetUpOpenGLPort_( this );

  glLineWidth( 1 );

  /* draw the vector */
  glBegin ( GL_LINES );
  glColor3f ( 1, 0, 0 );
  glVertex2f ( bufferStart.mfX, bufferStart.mfY );
  glColor3fv ( ifaColor );
  glVertex2f ( bufferEnd.mfX, bufferEnd.mfY );
  glEnd ();

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DrawVector_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;

}

DspA_tErr DspA_GetCursor ( tkmDisplayAreaRef this,
                           xVoxelRef          opCursor ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Return it. */
  xVoxl_Copy( opCursor, this->mpCursor );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetCursor: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_GetCursorInMRIIdx ( tkmDisplayAreaRef this,
                                   xVoxelRef         opMRIIdx ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xVoxel    MRIIdx;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Convert the cursor to MRI idx. */
  Volm_ConvertIdxToMRIIdx(this->mpVolume[tkm_tVolumeType_Main],
                          this->mpCursor, &MRIIdx);

  /* Return it. */
  xVoxl_Copy( opMRIIdx, &MRIIdx );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetCursorInMRIIdx: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_GetOrientation ( tkmDisplayAreaRef this,
                                mri_tOrientation* oOrientation ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* return the orientation */
  *oOrientation = this->mOrientation;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetOrientation: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_GetZoomLevel ( tkmDisplayAreaRef this,
                              int*              oZoomLevel ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* return the zoom level */
  *oZoomLevel = this->mnZoomLevel;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetZoomLevel: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

int DspA_GetCurrentSliceNumber_ ( tkmDisplayAreaRef this ) {

  int nSlice = 0;

  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    nSlice = xVoxl_GetZ( this->mpCursor );
    break;
  case mri_tOrientation_Horizontal:
    nSlice = xVoxl_GetY( this->mpCursor );
    break;
  case mri_tOrientation_Sagittal:
    nSlice = xVoxl_GetX( this->mpCursor );
    break;
  default:
    nSlice = 0;
    break;
  }

  return nSlice;
}

int DspA_GetCurrentSliceNumberInMRIIdx_ ( tkmDisplayAreaRef this ) {

  int nSlice = 0;
  xVoxel MRIIdx;

  Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                           this->mpCursor, &MRIIdx );

  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    nSlice = xVoxl_GetZ( &MRIIdx );
    break;
  case mri_tOrientation_Horizontal:
    nSlice = xVoxl_GetY( &MRIIdx );
    break;
  case mri_tOrientation_Sagittal:
    nSlice = xVoxl_GetX( &MRIIdx );
    break;
  default:
    nSlice = 0;
    break;
  }

  return nSlice;
}

DspA_tErr DspA_GetSelectedHeadPt ( tkmDisplayAreaRef   this,
                                   HPtL_tHeadPointRef* opHeadPoint ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* return the head point */
  *opHeadPoint = this->mpSelectedHeadPoint;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetSelectedHeadPt: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_GetClosestInterpSurfVoxel ( tkmDisplayAreaRef this,
    tkm_tSurfaceType  iSurface,
    Surf_tVertexSet   iSet,
    xVoxelRef         iAnaIdx,
    xVoxelRef         oOrigAnaIdx,
    xVoxelRef         oInterpAnaIdx,
    char*             osDescription ) {

  DspA_tErr             eResult               = DspA_tErr_NoErr;
  xGArr_tErr            eList                 = xGArr_tErr_NoErr;
  int                   nSlice                = 0;
  xGrowableArrayRef     list                  = NULL;
  DspA_tSurfaceListNode drawListNode;
  float                 dx                    = 0;
  float                 dy                    = 0;
  float                 dz                    = 0;
  float                 fDistance             = 0;
  float                 fLowestDistance       = 0;
  int                   nClosestIndex         = 0;
  int                   nClosestNeighborIndex = 0;
  xVoxel                closestAnaIdx;
  xVoxel                closestNeighborAnaIdx;
  xVoxel                closestInterpAnaIdx;
  xVoxel                closestRAS;
  xVoxel                closestNeighborRAS;
  xVoxel                closestInterpRAS;


  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  if ( iSurface < 0 || iSurface >= tkm_knNumSurfaceTypes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  if ( iSet < 0 || iSet >= Surf_knNumVertexSets ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  if ( NULL == iAnaIdx ||
       NULL == oInterpAnaIdx ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* walk through all surface lists. for each one... */
  fLowestDistance = 99999;
  nClosestIndex = -1;
  for ( nSlice = 0; nSlice < 256; nSlice++ ) {

    list = DspA_GetSurfaceList_( this, iSurface,
                                 this->mOrientation, iSet, nSlice );
    if ( NULL == list ) {
      continue;
    }

    /* walk through the list. */
    eList = xGArr_ResetIterator( list );
    if ( xGArr_tErr_NoErr != eList )
      goto error;
    while ( (eList = xGArr_NextItem( list, (void*)&drawListNode ))
            == xGArr_tErr_NoErr ) {

      /* skip if it's not a vertex. */
      if ( !drawListNode.mbVertex ) {
        continue;
      }

      /* calc the fake (unsquared) distance to the interp voxel. */
      dx = xVoxl_GetFloatX( iAnaIdx ) -
           xVoxl_GetFloatX( &drawListNode.mInterpVertex );
      dy = xVoxl_GetFloatY( iAnaIdx ) -
           xVoxl_GetFloatY( &drawListNode.mInterpVertex );
      dz = xVoxl_GetFloatZ( iAnaIdx ) -
           xVoxl_GetFloatZ( &drawListNode.mInterpVertex );
      fDistance = dx*dx + dy*dy + dz*dz;

      if ( fDistance < fLowestDistance ) {
        nClosestIndex = drawListNode.mnOriginalVertexIndex;
        nClosestNeighborIndex = drawListNode.mnNeighborVertexIndex;
        xVoxl_Copy( &closestInterpAnaIdx, &drawListNode.mInterpVertex );
        xVoxl_Copy( &closestAnaIdx, &drawListNode.mOriginalVertex );
        xVoxl_Copy( &closestNeighborAnaIdx, &drawListNode.mNeighborVertex );
        fLowestDistance = fDistance;
      }
    }
  }

  if ( -1 == nClosestIndex ) {
    eResult = DspA_tErr_CouldntFindClosestVoxel;
    goto error;
  }

  /* Copy the closest interp vertex out out. Copy the original if they
     want it. */
  xVoxl_Copy( oOrigAnaIdx, &closestAnaIdx );
  if ( NULL != oInterpAnaIdx ) {
    xVoxl_Copy( oInterpAnaIdx, &closestInterpAnaIdx );
  }

  /* calculate the actual distance. */
  fLowestDistance = sqrt( fLowestDistance );

  /* If the want one, make a string of info. */
  if ( NULL != osDescription ) {

    /* convert the ana idx coords to RAS. */
    Volm_ConvertIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                          &closestAnaIdx, &closestRAS );
    Volm_ConvertIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                          &closestNeighborAnaIdx, &closestNeighborRAS );
    Volm_ConvertIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                          &closestInterpAnaIdx, &closestInterpRAS );

    sprintf( osDescription, "Index: %d Distance: %.2f\n"
             "\t    Original vertex index: %d RAS: %.2f %.2f %.2f\n"
             "\t    Neighbor vertex index: %d RAS: %.2f %.2f %.2f\n"
             "\tIntersection vertex RAS: %.2f %.2f %.2f",
             nClosestIndex, fLowestDistance,
             nClosestIndex, xVoxl_ExpandFloat( &closestRAS ),
             nClosestNeighborIndex, xVoxl_ExpandFloat( &closestNeighborRAS ),
             xVoxl_ExpandFloat( &closestInterpRAS ) );
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GetClosestInterpSurfVoxel: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

  if ( xGArr_tErr_NoErr != eList ) {
    DebugPrint( ("Error %d in DspA_GetClosestInterpSurfVoxel: %s\n",
                 eResult, xGArr_GetErrorString(eList) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_BuildLineToolVoxelList_( tkmDisplayAreaRef this ) {

  DspA_tErr  eResult    = DspA_tErr_NoErr;
  xPoint2n   bufferPt;
  xVoxel     anaIdx;
  xVoxel     lineVox;
  int px, qx, tx, py, qy, ty;

  lineVox.mfX=0;
  lineVox.mfY=0;

  px = this->mLineVertex1.mnX;
  py = this->mLineVertex1.mnY;
  qx = this->mLineVertex2.mnX;
  qy = this->mLineVertex2.mnY;

  this->mNumLineVoxels = 0;
  for ( bufferPt.mnY = 0; 
        bufferPt.mnY < this->mnVolumeSizeY; 
        bufferPt.mnY += this->mnZoomLevel ) {
    for ( bufferPt.mnX = 0; 
          bufferPt.mnX < this->mnVolumeSizeX; 
          bufferPt.mnX += this->mnZoomLevel ) {

      eResult = DspA_ConvertBufferToVolume_ ( this, &bufferPt, &anaIdx );
      if ( DspA_tErr_NoErr == eResult ) {

        DspA_NormalizeVoxel_( &anaIdx, this->mOrientation, &lineVox );
        tx = xVoxl_GetX( &lineVox );
        ty = xVoxl_GetY( &lineVox );

        if ( ABS( (qy-py)*(tx-px) - (ty-py)*(qx-px) ) >=
             MAX( abs(qx-px), ABS(qy-py) ) )
          continue;
        if ( (qx<px && px<tx) || (qy<py && py<ty) )
          continue;
        if ( (tx<px && px<qx) || (ty<py && py<qy) )
          continue;
        if ( (px<qx && qx<tx) || (py<qy && qy<ty) )
          continue;
        if ( (tx<qx && qx<px) || (ty<qy && qy<py) )
          continue;

        if ( this->mNumLineVoxels < DspA_knMaxNumLineVoxels ) {
          xVoxl_Copy( &(this->mLineVoxels[this->mNumLineVoxels]), &anaIdx );
          this->mNumLineVoxels++;
        }
      }
    }
  }


  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_BuildLineToolVoxelList_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


DspA_tErr DspA_InitSurfaceLists_( tkmDisplayAreaRef this,
                                  int               inNumLists ) {

  DspA_tErr           eResult    = DspA_tErr_NoErr;
  int                 nSurface   = 0;
  int                 nDrawList  = 0;

  /* for each surface type */
  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {

    /* allocate surface point lists. */
    this->maSurfaceLists[nSurface] = 
      (xGrowableArrayRef*)
      malloc (sizeof(xGrowableArrayRef) * inNumLists );

    /* set all lists to null */
    for ( nDrawList = inNumLists - 1; nDrawList >= 0; nDrawList-- ) {
      this->maSurfaceLists[nSurface][nDrawList] = NULL;
    }
  }

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_InitSurfaceLists_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_PurgeSurfaceLists_ ( tkmDisplayAreaRef this ) {

  DspA_tErr           eResult    = DspA_tErr_NoErr;
  int                 nSurface   = 0;
  int                 nDrawList  = 0;
  xGArr_tErr          eList      = xGArr_tErr_NoErr;

  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {

    for ( nDrawList = DspA_GetNumSurfaceLists_( this ) - 1;
          nDrawList >= 0; nDrawList-- ) {

      /* if this is a list... */
      if ( NULL != this->maSurfaceLists[nSurface][nDrawList] ) {

        /* delete the list. */
        eList = xGArr_Delete( &this->maSurfaceLists[nSurface][nDrawList] );
        if ( xGArr_tErr_NoErr != eList ) {
          eResult = DspA_tErr_ErrorAccessingSurfaceList;
          goto error;
        }
      }
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_PurgeSurfaceLists_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_NewSurfaceList_ ( tkmDisplayAreaRef this,
                                 tkm_tSurfaceType  iSurface,
                                 mri_tOrientation  iOrientation,
                                 Surf_tVertexSet   iVertexSet,
                                 int               inSlice ) {

  DspA_tErr   eResult   = DspA_tErr_NoErr;
  int         nDrawList = 0;
  xGArr_tErr  eList     = xGArr_tErr_NoErr;

  /* get the list index. */
  nDrawList = DspA_GetSurfaceListIndex_( this, iOrientation,
                                         iVertexSet, inSlice );

  /* allocate a list. */
  xGArr_New( &this->maSurfaceLists[iSurface][nDrawList],
             sizeof( DspA_tSurfaceListNode ), 512 );
  if ( xGArr_tErr_NoErr != eList ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_NewSurfaceList_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

xGrowableArrayRef DspA_GetSurfaceList_ ( tkmDisplayAreaRef this,
    tkm_tSurfaceType  iSurface,
    mri_tOrientation  iOrientation,
    Surf_tVertexSet  iVertexSet,
    int               inSlice ) {

  xGrowableArrayRef pList = NULL;

  if ( NULL != this->maSurfaceLists ) {
    pList = this->maSurfaceLists [iSurface]
            [ DspA_GetSurfaceListIndex_( this, 
                                         iOrientation, 
                                         iVertexSet, 
                                         inSlice ) ];
  }

  return pList;
}

int DspA_GetNumSurfaceLists_ ( tkmDisplayAreaRef this ) {

  /* may be incorrect, assumes dimensions are all equal */
  return this->mnVolumeSizeZ * Surf_knNumVertexSets * mri_knNumOrientations;
}

int DspA_GetSurfaceListIndex_ ( tkmDisplayAreaRef this,
                                mri_tOrientation  iOrientation,
                                Surf_tVertexSet  iVertexSet,
                                int               inSlice ) {

  /* may be incorrect, assumes dimensions are all equal */
  return ((int)iOrientation * Surf_knNumVertexSets * this->mnVolumeSizeZ) +
         ((int)iVertexSet * this->mnVolumeSizeX) + inSlice;
}




DspA_tErr DspA_ConvertVolumeToBuffer_ ( tkmDisplayAreaRef this,
                                        xVoxelRef         ipVolumeVox,
                                        xPoint2nRef       opBufferPt ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fVolumeSize = (float) this->mnVolumeSizeX;

  /* verify the voxel */
  eResult = DspA_VerifyVolumeVoxel_( this, ipVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* first extract the points in the voxel that are on the same place as
     our orientation. */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    fX = xVoxl_GetFloatX( ipVolumeVox );
    fY = xVoxl_GetFloatY( ipVolumeVox );
    break;
  case mri_tOrientation_Horizontal:
    fX = xVoxl_GetFloatX( ipVolumeVox );
    fY = xVoxl_GetFloatZ( ipVolumeVox );
    break;
  case mri_tOrientation_Sagittal:
    fX = xVoxl_GetFloatZ( ipVolumeVox );
    fY = xVoxl_GetFloatY( ipVolumeVox );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* flip */
  fY = GLDRAW_Y_FLIP_FLOAT(fY);

  /* now zoom the coords to our zoomed buffer state */
  fX = (fZoomLevel * (fX - xVoxl_GetFloatX(this->mpZoomCenter))) +
       (fVolumeSize/2.0);
  fY = (fZoomLevel *
        (fY - GLDRAW_Y_FLIP_FLOAT(xVoxl_GetFloatY(this->mpZoomCenter)))) +
       (fVolumeSize/2.0);

  /* return the point */
  opBufferPt->mnX = (int) floor( fX );
  opBufferPt->mnY = (int) floor( fY );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DspA_ConvertVolumeToBuffer_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertVolumeToBufferf_ ( tkmDisplayAreaRef this,
    xVoxelRef         ipVolumeVox,
    xPoint2fRef       opBufferPt ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fVolumeSize = (float) this->mnVolumeSizeX;

  /* verify the voxel */
  eResult = DspA_VerifyVolumeVoxel_( this, ipVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* first extract the points in the voxel that are on the same place as
     our orientation. */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    fX = xVoxl_GetFloatX( ipVolumeVox );
    fY = xVoxl_GetFloatY( ipVolumeVox );
    break;
  case mri_tOrientation_Horizontal:
    fX = xVoxl_GetFloatX( ipVolumeVox );
    fY = xVoxl_GetFloatZ( ipVolumeVox );
    break;
  case mri_tOrientation_Sagittal:
    fX = xVoxl_GetFloatZ( ipVolumeVox );
    fY = xVoxl_GetFloatY( ipVolumeVox );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* flip */
  fY = GLDRAW_Y_FLIP_FLOAT(fY);

  /* now zoom the coords to our zoomed buffer state and return them */
  opBufferPt->mfX = (fZoomLevel * (fX - xVoxl_GetFloatX(this->mpZoomCenter))) +
                    (fVolumeSize/2.0);
  opBufferPt->mfY = 
    (fZoomLevel * 
     (fY - 
      GLDRAW_Y_FLIP_FLOAT(xVoxl_GetFloatY(this->mpZoomCenter)))) +
    (fVolumeSize/2.0);

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DspA_ConvertVolumeToBufferf_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertBufferToVolume_ ( tkmDisplayAreaRef this,
                                        xPoint2nRef       ipBufferPt,
                                        xVoxelRef         opVolumeVox ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
  xVoxel    anaIdx;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fBufferX    = (float) ipBufferPt->mnX;
  float     fBufferY    = (float) ipBufferPt->mnY;
  float     fVolumeSize = (float) this->mnVolumeSizeX;

  /* unzoom the coords */
  fX = fBufferX / fZoomLevel +
       ( xVoxl_GetFloatX(this->mpZoomCenter) - (fVolumeSize/2.0/fZoomLevel) );
  fY = fBufferY / fZoomLevel +
       ( xVoxl_GetFloatY(this->mpZoomCenter) - (fVolumeSize/2.0/fZoomLevel) );

  /* build a voxel out of these two coords and the current slice */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_SetFloat( &anaIdx, fX, fY, DspA_GetCurrentSliceNumber_( this ) );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_SetFloat( &anaIdx, fX, DspA_GetCurrentSliceNumber_( this ), fY );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetFloat( &anaIdx, DspA_GetCurrentSliceNumber_( this ), fY, fX );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* copy into the outgoing voxel to return it */
  xVoxl_Copy( opVolumeVox, &anaIdx );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DspA_ConvertBufferToVolume_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertPlaneToVolume_ ( tkmDisplayAreaRef this,
                                       xPoint2fRef       ipPlanePt,
                                       int               inSlice,
                                       mri_tOrientation  iOrientation,
                                       xVoxelRef         opVolumeVox ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  switch ( iOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_SetFloatX( opVolumeVox, ipPlanePt->mfX );
    xVoxl_SetFloatY( opVolumeVox, ipPlanePt->mfY );
    xVoxl_SetFloatZ( opVolumeVox, inSlice );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_SetFloatX( opVolumeVox, ipPlanePt->mfX );
    xVoxl_SetFloatY( opVolumeVox, inSlice );
    xVoxl_SetFloatZ( opVolumeVox, ipPlanePt->mfY );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetFloatX( opVolumeVox, inSlice );
    xVoxl_SetFloatY( opVolumeVox, ipPlanePt->mfY );
    xVoxl_SetFloatZ( opVolumeVox, ipPlanePt->mfX );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
  }

  return eResult;
}

DspA_tErr DspA_ConvertVolumeToPlane_ ( tkmDisplayAreaRef this,
                                       xVoxelRef         ipVolumeVox,
                                       mri_tOrientation  iOrientation,
                                       xPoint2fRef       opPlanePt,
                                       int*              onSlice ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  switch ( iOrientation ) {
  case mri_tOrientation_Coronal:
    opPlanePt->mfX = xVoxl_GetFloatX( ipVolumeVox );
    opPlanePt->mfY = xVoxl_GetFloatY( ipVolumeVox );
    *onSlice       = xVoxl_GetFloatZ( ipVolumeVox );
    break;
  case mri_tOrientation_Horizontal:
    opPlanePt->mfX = xVoxl_GetFloatX( ipVolumeVox );
    opPlanePt->mfY = xVoxl_GetFloatZ( ipVolumeVox );
    *onSlice       = xVoxl_GetFloatY( ipVolumeVox );
    break;
  case mri_tOrientation_Sagittal:
    opPlanePt->mfX = xVoxl_GetFloatZ( ipVolumeVox );
    opPlanePt->mfY = xVoxl_GetFloatY( ipVolumeVox );
    *onSlice       = xVoxl_GetFloatX( ipVolumeVox );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
  }

  return eResult;
}

DspA_tErr DspA_ConvertBufferToScreen_ ( tkmDisplayAreaRef this,
                                        xPoint2nRef       ipBufferPt,
                                        xPoint2nRef       opScreenPt ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  int       nX      = 0;
  int       nY      = 0;

  /* verify the buffer pt */
  eResult = DspA_VerifyBufferPoint_( this, ipBufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* use the frame scale factor to convert */
  nX = (float)ipBufferPt->mnX * this->mfFrameBufferScaleX;
  nY = (float)ipBufferPt->mnY * this->mfFrameBufferScaleY;

  /* return the point */
  opScreenPt->mnX = nX;
  opScreenPt->mnY = nY;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DspA_ConvertBufferToScreen_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertScreenToBuffer_ ( tkmDisplayAreaRef this,
                                        xPoint2nRef       ipScreenPt,
                                        xPoint2nRef       opBufferPt ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xPoint2n  localPt = {0, 0};
  float     fX      = 0;
  float     fY      = 0;
  int       nX      = 0;
  int       nY      = 0;

  /* first subtract our position from the point to get our local point. the
     y craziness has to do with the fact that not only is the y coordinate
     system flipped, but the displays are flipped in the medit window. this
     just kind of works, at least for 2x2. */
  localPt.mnX = ipScreenPt->mnX - this->mLocationInSuper.mnX;
  localPt.mnY = (((int)(ipScreenPt->mnY + this->mnHeight)) % this->mnHeight);

  /* verify the screen pt */
  eResult = DspA_VerifyScreenPoint_( this, &localPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* use the frame scale factor to convert */
  fX = (float)localPt.mnX / this->mfFrameBufferScaleX;
  fY = (float)localPt.mnY / this->mfFrameBufferScaleY;

  nX = (int) floor( fX );
  nY = (int) floor( fY );

  /* y flip */
  nY = GLDRAW_Y_FLIP(nY);

  /* return the point */
  opBufferPt->mnX = nX;
  opBufferPt->mnY = nY;

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_DspA_ConvertScreenToBuffer_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SendViewStateToTcl_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult                    = DspA_tErr_NoErr;
  char      sVolumeName[tkm_knNameLen] = "";
  char      sTclArguments[STRLEN]      = "";
  int       nFlag                      = 0;
  int       brush                      = 0;
  int       nSurface                   = 0;
  int       nVertexSet                 = 0;

  /* send the point info for our cursor */
  DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
                                   this->mpCursor );

  /* send volume name */
  Volm_CopyVolumeName( this->mpVolume[tkm_tVolumeType_Main],
                       sVolumeName, sizeof(sVolumeName) );
  sprintf( sTclArguments, "\"%s value\"", sVolumeName );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );

  /* send the aux volume name if it's loaded */
  if ( NULL != this->mpVolume[tkm_tVolumeType_Aux] ) {
    Volm_CopyVolumeName( this->mpVolume[tkm_tVolumeType_Aux],
                         sVolumeName, sizeof(sVolumeName) );
    sprintf( sTclArguments, "\"%s value\"", sVolumeName );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeName, sTclArguments );
  }

  /* send the orientation */
  sprintf ( sTclArguments, "%d", (int)this->mOrientation );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateOrientation, sTclArguments );

  /* send the zoom level */
  sprintf ( sTclArguments, "%d", (int)this->mnZoomLevel );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateZoomLevel, sTclArguments );

  /* send all flags */
  for ( nFlag = 0; nFlag < DspA_knNumDisplayFlags; nFlag++ ) {
    sprintf ( sTclArguments, "%d %d", (int)nFlag,
              (int)this->mabDisplayFlags[nFlag] );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateDisplayFlag, sTclArguments );
  }

  /* send tool update. */
  sprintf ( sTclArguments, "%d", (int)sTool );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateTool, sTclArguments );

  /* send brush info */
  sprintf ( sTclArguments, "%d %d %d",
            (int)sBrush.mnRadius, (int)sBrush.mShape,
            (int)sBrush.mb3D );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushShape, sTclArguments );

  /* send the threshold info. */
  for ( brush = 0; brush < DspA_knNumBrushes; brush++ ) {
    sprintf ( sTclArguments, "%d %f %f %f %d %d",
              (int)brush,
              (float)sBrush.mInfo[brush].mLow,
              (float)sBrush.mInfo[brush].mHigh,
              (float)sBrush.mInfo[brush].mNewValue,
              (int)sBrush.mInfo[brush].mMode,
              (int)sBrush.mInfo[brush].mCloneSource );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushInfo, sTclArguments );
  }

  /* send the cursor color */
  sprintf ( sTclArguments, "%f %f %f",
            sCursorColor.mfRed, sCursorColor.mfGreen, sCursorColor.mfBlue );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateCursorColor, sTclArguments );

  /* send the surface line info */
  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {
    for ( nVertexSet = 0; nVertexSet < Surf_knNumVertexSets; nVertexSet++ ) {
      sprintf ( sTclArguments, "%d %d %d", (int)nSurface, (int)nVertexSet,
                this->manSurfaceLineWidth[nSurface][nVertexSet] );
      tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceLineWidth,
                          sTclArguments );
      sprintf ( sTclArguments, "%d %d %f %f %f", 
                (int)nSurface,(int)nVertexSet,
                xColr_ExpandFloat
                ( &(this->maSurfaceLineColor[nSurface][nVertexSet])));
      tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceLineColor,
                          sTclArguments );
    }
  }

  return eResult;
}

DspA_tErr DspA_SendPointInformationToTcl_ ( tkmDisplayAreaRef this,
    DspA_tDisplaySet  iSet,
    xVoxelRef         iAnaIdx ) {

  xVoxel                voxel;
  int                   eCTAB              = NO_ERROR;
  char                  sTclArguments[STRLEN] = "";
  int                   nSlice             = 0;
  float                 fVolumeValue       = 0;
  HPtL_tHeadPointRef    pHeadPoint         = NULL;
  FunV_tErr             eFunctional        = FunV_tErr_NoError;
  xVoxel                funcIdx            = { 0, 0, 0 };
  xVoxel                funcRAS            = { 0, 0, 0 };
  FunV_tFunctionalValue funcValue          = 0;
  tBoolean              bFuncSelection     = FALSE;
  int                   nSegLabelIndex     = 0;
  float                 fSegLabelIndex     = 0;
  char                  sLabel[STRLEN]        = "";
  int                   nValue             = 0;
  DspA_tHistogramParams histoParams;
  float                 fDistance          = 0;

  xVoxel MRIIdx;
  xVoxel auxMRIIdx;

  /* send the anatomical index. */
  // translate the screen idx into the src Idx
  Volm_ConvertIdxToMRIIdx(this->mpVolume[tkm_tVolumeType_Main],
                          iAnaIdx, &MRIIdx);

  // **********************************************************************
  // To be implemented: need to tell whether these values are valid or not.
  sprintf(sTclArguments, "%s %d %d %d", 
          DspA_ksaDisplaySet[iSet], xVoxl_ExpandInt(&MRIIdx) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeCursor, sTclArguments );


  /* send the slice number */
  nSlice = DspA_GetCurrentSliceNumber_( this );
  sprintf( sTclArguments, "%d", nSlice );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeSlice, sTclArguments );


  /* also convert to RAS and send those coords along. for these we
     actually use the surface RAS coords unless it's a new-style
     surface. */
  if ( tkm_UseRealRAS() ) {
    Volm_ConvertMRIIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                             &MRIIdx, &voxel );
  } else {
    Volm_ConvertMRIIdxToSurfaceRAS( this->mpVolume[tkm_tVolumeType_Main],
                                    &MRIIdx, &voxel );
  }
  sprintf( sTclArguments, "%s %.1f %.1f %.1f",
           DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateRASCursor, sTclArguments );

  /* also convert to mni and send those coords along. */
  if (NULL != 
      this->mpVolume[tkm_tVolumeType_Main]->mpMriValues->linear_transform) {
    Volm_ConvertIdxToMNITal( this->mpVolume[tkm_tVolumeType_Main],
                             iAnaIdx, &voxel );
    sprintf( sTclArguments, "%s %.1f %.1f %.1f",
             DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateMNICursor, sTclArguments );


    /* and the tal coords */
    Volm_ConvertIdxToTal( this->mpVolume[tkm_tVolumeType_Main],
                          iAnaIdx, &voxel );
    sprintf( sTclArguments, "%s %.1f %.1f %.1f",
             DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateTalCursor, sTclArguments );
  }


  /* and the scanner coords. */
  Volm_ConvertIdxToScanner( this->mpVolume[tkm_tVolumeType_Main],
                            iAnaIdx, &voxel );
  sprintf( sTclArguments, "%s %.1f %.1f %.1f",
           DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateScannerCursor, sTclArguments );


  /* also get the volume value and send that along. */
  Volm_GetValueAtMRIIdx_( this->mpVolume[tkm_tVolumeType_Main],
                          &MRIIdx, &fVolumeValue );
  switch (this->mpVolume[tkm_tVolumeType_Main]->mpMriValues->type) {
  default:
    /* If there is an aux volume loaded and it's not being
    displayed, show the value in stars, otherwise no stars. */
    if ( (NULL == this->mpVolume[tkm_tVolumeType_Aux]) ||
         this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      sprintf( sTclArguments, "%s %d",
               DspA_ksaDisplaySet[iSet], (int)fVolumeValue );
    } else {
      sprintf( sTclArguments, "%s **%d**",
               DspA_ksaDisplaySet[iSet], (int)fVolumeValue );
    }
    break ;
  case MRI_FLOAT: {
    char fmt[STRLEN] ;
    float f ;
    int   decs ;

    f = fabs(fVolumeValue) ;
    if (f > 1) decs = 2 ;
    else if (f > .1) decs = 3 ;
    else if (f > .01) decs = 4 ;
    else if (f > 0.001) decs = 5 ;
    else if (f > 0.0001) decs = 6 ;
    else if (f > 0.00001) decs = 7 ;
    else if (f > 0.000001) decs = 8 ;
    else if (f > 0.0000001) decs = 9 ;
    else decs = 10 ;

    sprintf(fmt, "%%s %%2.%df", decs) ;
    sprintf( sTclArguments, fmt,
             DspA_ksaDisplaySet[iSet], fVolumeValue );
    break ;
  }
  }

  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValue, sTclArguments );

  /* send aux volume value if it's loaded. */
  if ( NULL != this->mpVolume[tkm_tVolumeType_Aux] ) {

    Volm_ConvertIdxToMRIIdx(this->mpVolume[tkm_tVolumeType_Aux],
                            iAnaIdx, &auxMRIIdx);

    Volm_GetValueAtMRIIdx_( this->mpVolume[tkm_tVolumeType_Aux],
                            &auxMRIIdx, &fVolumeValue );

    switch (this->mpVolume[tkm_tVolumeType_Aux]->mpMriValues->type) {
    default:
      if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
        sprintf( sTclArguments, "%s **%d**",
                 DspA_ksaDisplaySet[iSet], (int)fVolumeValue );
      } else {
        sprintf( sTclArguments, "%s %d",
                 DspA_ksaDisplaySet[iSet], (int)fVolumeValue );
      }
      break ;
    case MRI_FLOAT: {
      char fmt[STRLEN] ;
      float f ;
      int   decs ;

      f = fabs(fVolumeValue) ;
      if (f > 1)     decs = 2 ;
      else if (f > .1) decs = 3 ;
      else if (f > .01) decs = 4 ;
      else if (f > 0.001) decs = 5 ;
      else if (f > 0.0001) decs = 6 ;
      else if (f > 0.00001) decs = 7 ;
      else if (f > 0.000001) decs = 8 ;
      else if (f > 0.0000001) decs = 9 ;
      else decs = 10 ;

      sprintf(fmt, "%%s %%2.%df", decs) ;
      sprintf( sTclArguments, fmt,
               DspA_ksaDisplaySet[iSet], fVolumeValue );
      break ;
    }
    }
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeValue, sTclArguments );
  }

  /* also see if we have functional data and can send a value for that
     as well. */
  //  if ( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ) {
  if ( NULL != this->mpFunctionalVolume ) {

    DisableDebuggingOutput;
    if ( DspA_tDisplaySet_Cursor == iSet ) {

      /* if this is the cursor, use the voxel locations for the selection
      and the average value */
      eFunctional = FunV_GetAvgFunctionalValue( this->mpFunctionalVolume,
                    &funcValue, &funcIdx, &funcRAS,
                    &bFuncSelection );

    } else if ( DspA_tDisplaySet_Mouseover == iSet ) {

      /* if this is the mouseover, use the value at the current point,
      and convert the point that to the func idx and ras. */
      eFunctional = FunV_GetValueAtMRIIdx( this->mpFunctionalVolume,
                                           &MRIIdx, FALSE, &funcValue );
      if ( FunV_tErr_NoError == eFunctional ) {

        /* convert the points */
        FunV_ConvertMRIIdxToFuncIdx( this->mpFunctionalVolume,
                                     &MRIIdx, &funcIdx );
        FunV_ConvertMRIIdxToFuncRAS( this->mpFunctionalVolume,
                                     &MRIIdx, &funcRAS );
      }
    }
    EnableDebuggingOutput;

    /* send the value and the coords */
    sprintf( sTclArguments, "%s %f", DspA_ksaDisplaySet[iSet], funcValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateFunctionalValue,
                        sTclArguments );


    sprintf( sTclArguments, "%s %.1f %.1f %.1f",
             DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &funcIdx ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateFunctionalCoords,
                        sTclArguments );


    sprintf( sTclArguments, "%s %.1f %.1f %.1f",
             DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &funcRAS ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateFunctionalRASCoords,
                        sTclArguments );
  }


  /* and the seg label if we have one */
  if ( NULL != this->mSegmentationVolume[tkm_tSegType_Main] &&
       NULL != this->mSegmentationColorTable[tkm_tSegType_Main] ) {

    /* Get the value in ana idx space (this was messing up for some
       volumes with weird CRAS transforms) and get the corresponding
       label string. */
    Volm_GetValueAtIdx_( this->mSegmentationVolume[tkm_tSegType_Main],
                         iAnaIdx, &fSegLabelIndex );
    if ( 0 == fSegLabelIndex ) {
      strcpy( sLabel, "None" );
    } else {
      eCTAB =
        CTABcopyName( this->mSegmentationColorTable[tkm_tSegType_Main],
                      fSegLabelIndex, sLabel, sizeof(sLabel) );
      if ( NO_ERROR != eCTAB ) {
        strcpy( sLabel, "Out of bounds." );
      }
    }

    /* if this is a click, and this is the seg volume we're looking
       at, set the index */
    if ( DspA_tDisplaySet_Cursor == iSet &&
         !this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume]) {
      this->mnSegmentationVolumeIndex = nSegLabelIndex;
    }

    /* if this is a click with the edit tool and the volume count flag
       is on, and this is the volume we're looking at, calc and
       display the volume */
    if ( DspA_tDisplaySet_Cursor == iSet &&
         DspA_tTool_EditSegmentation == sTool &&
         this->mabDisplayFlags[DspA_tDisplayFlag_SegLabelVolumeCount] &&
         !this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume]) {

      tkm_CalcSegLabelVolume( tkm_tSegType_Main, &MRIIdx, &nValue );
      sprintf( sTclArguments, "%s \"%s (%d)\"",
               DspA_ksaDisplaySet[iSet], sLabel, nValue );

      /* else just the label */
    } else {
      /* If there is an aux volume loaded and it's not being
      displayed, show the value in stars, otherwise no stars. */
      if ( (NULL == this->mSegmentationVolume[tkm_tSegType_Aux]) ||
           this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume] ) {
        sprintf( sTclArguments, "%s \"%s\"",
                 DspA_ksaDisplaySet[iSet], sLabel );
      } else {
        sprintf( sTclArguments, "%s \"**%s**\"",
                 DspA_ksaDisplaySet[iSet], sLabel );
      }
    }
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSegLabel, sTclArguments );
  }


  /* and the aux seg label if we have one */
  if ( NULL != this->mSegmentationVolume[tkm_tSegType_Aux] ) {

    /* Get the value in ana idx space (this was messing up for some
       volumes with weird CRAS transforms) and get the corresponding
       label string. */
    Volm_GetValueAtIdx_( this->mSegmentationVolume[tkm_tSegType_Aux],
                         iAnaIdx, &fSegLabelIndex );
    if ( 0 == fSegLabelIndex ) {
      strcpy( sLabel, "None" );
    } else {
      eCTAB =
        CTABcopyName( this->mSegmentationColorTable[tkm_tSegType_Aux],
                      fSegLabelIndex, sLabel, sizeof(sLabel) );
      if ( NO_ERROR != eCTAB ) {
        strcpy( sLabel, "Out of bounds." );
      }
    }

    /* if this is a click, and this is the seg volume we're looking
       at, set the index */
    if ( DspA_tDisplaySet_Cursor == iSet &&
         this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume]) {
      this->mnSegmentationVolumeIndex = nSegLabelIndex;
    }

    /* if this is a click with the edit tool and the volume count flag
       is on, and this is the volume we're looking at, calc and
       display the volume */
    if ( DspA_tDisplaySet_Cursor == iSet &&
         DspA_tTool_EditSegmentation == sTool &&
         this->mabDisplayFlags[DspA_tDisplayFlag_SegLabelVolumeCount] &&
         this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume]) {

      tkm_CalcSegLabelVolume( tkm_tSegType_Aux, &MRIIdx, &nValue );
      sprintf( sTclArguments, "%s \"%s (%d)\"",
               DspA_ksaDisplaySet[iSet], sLabel, nValue );

      /* else just the label */
    } else {
      if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxSegmentationVolume] ) {
        sprintf( sTclArguments, "%s \"**%s**\"",
                 DspA_ksaDisplaySet[iSet], sLabel );
      } else {
        sprintf( sTclArguments, "%s \"%s\"",
                 DspA_ksaDisplaySet[iSet], sLabel );
      }
    }
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxSegLabel, sTclArguments );
  }


  /* and the head point label if it's on */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_HeadPoints]
       && NULL != this->mHeadPoints ) {

    /* get the closest head point */
    tkm_GetHeadPoint( iAnaIdx, this->mOrientation,
                      this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj],
                      &pHeadPoint );
    if ( NULL != pHeadPoint ) {

      /* get a nice label */
      sprintf( sTclArguments, "%s \"%s\"",
               DspA_ksaDisplaySet[iSet], pHeadPoint->msLabel );

    } else {

      sprintf( sTclArguments, "%s %s",
               DspA_ksaDisplaySet[iSet], "None" );
    }
    tkm_SendTclCommand( tkm_tTclCommand_UpdateHeadPointLabel, sTclArguments );
  }


  /* if we have gca data and this is the cursor, use gcadump to dump
     the info to the screen */
  if ( NULL != this->mGCAVolume &&
       DspA_tDisplaySet_Cursor == iSet ) {
    GCAdump( this->mGCAVolume,
             this->mpVolume[tkm_tVolumeType_Main]->mpMriValues,
             xVoxl_ExpandInt( iAnaIdx ), this->mGCATransform, stdout,
             this->mabDisplayFlags[DspA_tDisplayFlag_VerboseGCADump] );
  }

  if ( NULL != this->mVLI1 &&
       DspA_tDisplaySet_Cursor == iSet ) {
    int xn, yn, zn, label_counts_c1[256], label_counts_c2[256], index,
    inNumValues, n1, n2 ;
    VL *vl ;
    char name1[STRLEN], name2[STRLEN], *cp ;

    FileNameOnly(this->isVLI1_name, name1) ;
    cp = strrchr(name1, '.') ;
    if (cp) *cp = 0 ;
    FileNameOnly(this->isVLI2_name, name2) ;
    cp = strrchr(name2, '.') ;
    if (cp) *cp = 0 ;
    xn = nint(xVoxl_GetX(iAnaIdx) / this->mVLI1->resolution) ;
    yn = nint(xVoxl_GetY(iAnaIdx) / this->mVLI1->resolution) ;
    zn = nint(xVoxl_GetZ(iAnaIdx) / this->mVLI1->resolution) ;

    memset(label_counts_c1, 0, sizeof(label_counts_c1)) ;
    vl = &this->mVLI1->vl[xn][yn][zn] ;
    for (n1 = index = 0 ; index < vl->nlabels ; index++) {
      label_counts_c1[vl->labels[index]] += vl->counts[index] ;
      n1 += vl->counts[index] ;
    }

    memset(label_counts_c2, 0, sizeof(label_counts_c2)) ;
    vl = &this->mVLI2->vl[xn][yn][zn] ;
    for (n2 = index = 0 ; index < vl->nlabels ; index++) {
      label_counts_c2[vl->labels[index]] += vl->counts[index] ;
      n2 += vl->counts[index] ;
    }

    /* count # of different labels */
    for (inNumValues = index = 0 ; index < 256 ; index++)
      if ((label_counts_c1[index]) > 0 || label_counts_c2[index] > 0)
        inNumValues++ ;

    xUtil_snprintf( histoParams.msTitle, sizeof(histoParams.msTitle),
                    "labels of %s vs. %s at (%d,%d,%d)",
                    name1, name2,
                    xVoxl_ExpandInt( iAnaIdx ) );
    /* the titles of the axes */
    xUtil_strncpy( histoParams.msXAxisTitle, "Labels",
                   sizeof(histoParams.msXAxisTitle) );

    xUtil_strncpy( histoParams.msLabel1, name1,
                   sizeof(histoParams.msLabel1) );
    xUtil_strncpy( histoParams.msLabel2,  name2,
                   sizeof(histoParams.msLabel2) );

    /* declare arrays of values and labels */
    histoParams.mafValues1 = (float*) calloc( inNumValues, sizeof(float) );
    histoParams.mafValues2 = (float*) calloc( inNumValues, sizeof(float) );
    histoParams.masXAxisLabels = (char**) calloc( inNumValues, sizeof(char *));

    /* fill them up with some random nubmers and cma labels */
    histoParams.mnNumValues = inNumValues;
    for (nValue = index = 0 ; index < 256 ; index++) {
      if ((label_counts_c1[index]) == 0 && label_counts_c2[index] == 0)
        continue ;

      /* assign values for the elements */
      histoParams.mafValues1[nValue] = (float)label_counts_c1[index] ;
      histoParams.mafValues2[nValue] = (float)label_counts_c2[index] ;
      if (this->mabDisplayFlags[DspA_tDisplayFlag_HistogramPercentChange]) {
        xUtil_strncpy( histoParams.msYAxisTitle, "% of voxels with label",
                       sizeof(histoParams.msYAxisTitle) );
        histoParams.mafValues1[nValue] /= ((float)n1/100) ;
        histoParams.mafValues2[nValue] /= ((float)n2/100) ;
      } else
        xUtil_strncpy( histoParams.msYAxisTitle, "# of voxels with label",
                       sizeof(histoParams.msYAxisTitle) );


      /* allocate and set the label for this element */
      histoParams.masXAxisLabels[nValue] =
        (char*) malloc( sizeof(char) * DspA_knHistoTitleLength );
      xUtil_strncpy( histoParams.masXAxisLabels[nValue],
                     cma_label_to_name( index ), DspA_knHistoTitleLength );
      nValue++ ;
    }

    /* draw the thing */
    DspA_DrawHistogram( this, &histoParams );

    free( histoParams.mafValues1 );
    free( histoParams.mafValues2 );
    for ( nValue = 0; nValue < inNumValues; nValue++ )
      free( histoParams.masXAxisLabels[nValue] );
    free( histoParams.masXAxisLabels );
  }

  /* if we have a surfaec, update the surface distance. if this is the
     cursor, find the distance from the last cursor and print that. if
     this is the mouseover, print the distance from the current
     cursor.  */
  if ( NULL != this->mpSurface[tkm_tSurfaceType_Main] ) {

    if ( DspA_tDisplaySet_Cursor == iSet ) {
      Surf_GetDistance( this->mpSurface[tkm_tSurfaceType_Main],
                        iAnaIdx, this->mpLastCursor, &fDistance );
    } else {
      Surf_GetDistance( this->mpSurface[tkm_tSurfaceType_Main],
                        iAnaIdx, this->mpCursor, &fDistance );
    }


    /* send the update */
    sprintf( sTclArguments, "%s %f", DspA_ksaDisplaySet[iSet], fDistance );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateDistance, sTclArguments );
  }

  /* Send the line length update. */
  sprintf( sTclArguments, "%s %f", DspA_ksaDisplaySet[iSet],
           this->mLineDistance );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateLineLength, sTclArguments );

  return DspA_tErr_NoErr;
}

DspA_tErr DspA_DrawHistogram ( tkmDisplayAreaRef        this,
                               DspA_tHistogramParamsRef iParams ) {

  DspA_tErr eResult             = DspA_tErr_NoErr;
  char      sTclArgs[50000]     = "";
  char      sValues1List[10000] = "";
  char      sValues2List[10000] = "";
  char      sLabelsList[10000]  = "";
  int       nValue              = 0;

  DebugEnterFunction( ("DspA_ShowHistogram( this=%p, iParams=%p )",
                       this, iParams) );

  DebugNote( ("Verifying self") );
  eResult = DspA_Verify( this );
  DebugAssertThrow( (eResult == DspA_tErr_NoErr) );

  /* print all the titles */
  DebugNote( ("Printing titles") );
  xUtil_snprintf( sTclArgs, sizeof(sTclArgs),
                  "-title \"%s\" -xAxisTitle \"%s\" -yAxisTitle \"%s\" "
                  "-label1 \"%s\" -label2 \"%s\"", iParams->msTitle,
                  iParams->msXAxisTitle, iParams->msYAxisTitle,
                  iParams->msLabel1, iParams->msLabel2 );

  /* print the open brace for the list */
  DebugNote( ("Printing open braces") );
  xUtil_strncpy( sValues1List, "-values1 { ", sizeof(sValues1List) );
  xUtil_strncpy( sValues2List, "-values2 { ", sizeof(sValues2List) );
  xUtil_strncpy( sLabelsList, "-xAxisLabels { ", sizeof(sLabelsList) );

  /* for each value... */
  for ( nValue = 0; nValue < iParams->mnNumValues; nValue++ ) {

    /* build a list of values and labels */
    DebugNote( ("Printing elements %d/%d", nValue, iParams->mnNumValues) );
    sprintf( sValues1List, "%s %.2f",
             sValues1List, iParams->mafValues1[nValue] );
    sprintf( sValues2List, "%s %.2f",
             sValues2List, iParams->mafValues2[nValue] );
    sprintf( sLabelsList, "%s \"%s\"",
             sLabelsList, iParams->masXAxisLabels[nValue] );
  }

  /* print the closing brace for the list */
  DebugNote( ("Printing closing braces") );
  sprintf( sValues1List, "%s }", sValues1List );
  sprintf( sValues2List, "%s }", sValues2List );
  sprintf( sLabelsList, "%s }", sLabelsList );

  /* print the lists */
  DebugNote( ("Printing lists to main args") );
  sprintf( sTclArgs, "%s %s %s %s", sTclArgs,
           sValues1List, sValues2List, sLabelsList );

  tkm_SendTclCommand( tkm_tTclCommand_DrawHistogram, sTclArgs );

  DebugCatch;
  DebugCatchError( eResult, DspA_tErr_NoErr, DspA_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_SetSurfaceDistanceAtCursor ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fDistance   = 0;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make sure we have a surface. */
  if ( NULL != this->mpSurface[tkm_tSurfaceType_Main] ) {

    /* calc our distance. */
    Surf_GetDistance( this->mpSurface[tkm_tSurfaceType_Main],
                      this->mpCursor, this->mpLastCursor, &fDistance );

    /* set the distance for this ana idx */
    tkm_SetSurfaceDistance( this->mpCursor, fDistance );
  }


  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetSurfaceDistanceAtCursor: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_SetMRIValueAtCursorInSurface ( tkmDisplayAreaRef this,
    Surf_tVertexSet   iVertexSet) {

  DspA_tErr eResult  = DspA_tErr_NoErr;
  float     anaValue = 0;

  DebugEnterFunction( ("DspA_SetMRIValueAtCursorInSurface( this=%p, "
                       "iVertexSet=%d )", this, (int)iVertexSet) );

  /* verify us. */
  eResult = DspA_Verify ( this );
  DebugAssertThrow( (eResult == DspA_tErr_NoErr) );

  /* make sure we have a surface. */
  if ( NULL != this->mpSurface[tkm_tSurfaceType_Main] ) {

    /* Get the MRI value here. */
    Volm_GetValueAtIdx( this->mpVolume[tkm_tVolumeType_Main],
                        this->mpCursor, &anaValue );

    /* set the distance for this ana idx */
    tkm_SetMRIValueInSurface( this->mpCursor, iVertexSet, anaValue );
  }

  DebugCatch;
  DebugCatchError( eResult, DspA_tErr_NoErr, DspA_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_SmartCutAtCursor ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult      = DspA_tErr_NoErr;
  Volm_tErr eVolume      = Volm_tErr_NoErr;
  tkm_tVolumeType volume = tkm_tVolumeType_Main;
  int       curX         = 0;
  int       curY         = 0;
  int       curZ         = 0;
  int       minX         = 0;
  int       minY         = 0;
  int       minZ         = 0;
  int       maxX         = 0;
  int       maxY         = 0;
  int       maxZ         = 0;
  int       dMinX        = 0;
  int       dMinY        = 0;
  int       dMinZ        = 0;
  int       dMaxX        = 0;
  int       dMaxY        = 0;
  int       dMaxZ        = 0;
  int       dX           = 0;
  int       dY           = 0;
  int       dZ           = 0;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* Get the right volume based on which one we're viewing at the
     moment. Also gegt the dimensions from this volume. */
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
    volume = tkm_tVolumeType_Aux;
    eVolume = Volm_GetDimensions( this->mpVolume[tkm_tVolumeType_Aux],
                                  &maxX, &maxY, &maxZ );
  } else {
    volume = tkm_tVolumeType_Main;
    eVolume = Volm_GetDimensions( this->mpVolume[tkm_tVolumeType_Main],
                                  &maxX, &maxY, &maxZ );
  }

  /* Mins are zero. */
  minX = minY = minZ = 0;

  /* We actually want the max index, which is one less than the
     dimension. */
  maxX = maxX - 1;
  maxY = maxY - 1;
  maxZ = maxZ - 1;

  /* Get the cursor location. */
  curX = xVoxl_GetX( this->mpCursor );
  curY = xVoxl_GetY( this->mpCursor );
  curZ = xVoxl_GetZ( this->mpCursor );

  /* Find the distance between the cursor and the min and max
     indices. */
  dMinX = abs( minX - curX );
  dMaxX = abs( maxX - curX );
  dMinY = abs( minY - curY );
  dMaxY = abs( maxY - curY );
  dMinZ = abs( minZ - curZ );
  dMaxZ = abs( maxZ - curZ );

  /* Now find the lesser of the distance from min or max in a
     plane. */
  dX = MIN( dMinX, dMaxX );
  dY = MIN( dMinY, dMaxY );
  dZ = MIN( dMinZ, dMaxZ );

  /* Based on our orientation, set the distance that's parallel
     to the view to a big number so that it will not be the smallest
     distance. */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Sagittal:
    dX = 1000;
    break;
  case mri_tOrientation_Horizontal:
    dY = 1000;
    break;
  case mri_tOrientation_Coronal:
    dZ = 1000;
    break;
  default:
    break;
  }

  /* Now see in which dimension we have the smallest distance. Then,
     in that axis, find if we are closer to the min or max. Then call
     the set function with the proper region and set those values to
     0. */
  if ( dX < dY && dX < dZ ) {
    if ( dMinX < dMaxX ) {
      tkm_SetAnatomicalVolumeRegion( volume,
                                     minX, curX, minY, maxY, minZ, maxZ, 0 );
    } else {
      tkm_SetAnatomicalVolumeRegion( volume,
                                     curX, maxX, minY, maxY, minZ, maxZ, 0 );
    }
  } else if ( dY < dX && dY < dZ ) {
    if ( dMinY < dMaxY ) {
      tkm_SetAnatomicalVolumeRegion( volume,
                                     minX, maxX, minY, curY, minZ, maxZ, 0 );
    } else {
      tkm_SetAnatomicalVolumeRegion( volume,
                                     minX, maxX, curY, maxY, minZ, maxZ, 0 );
    }
  } else if ( dZ < dX && dZ < dY ) {
    if ( dMinZ < dMaxZ ) {
      tkm_SetAnatomicalVolumeRegion( volume,
                                     minX, maxX, minY, maxY, minZ, curZ, 0 );
    } else {
      tkm_SetAnatomicalVolumeRegion( volume,
                                     minX, maxX, minY, maxY, curZ, maxZ, 0 );
    }
  }

  /* Redraw ourselves. */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SmartCutAtCursor: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}



DspA_tErr DspA_AddLineToSelection ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult   = DspA_tErr_NoErr;
  xVoxel    MRIIdx;
  int       nLineVox  = 0;

  DebugEnterFunction( ("DspA_AddLineToSelection( this=%p )", this) );

  /* verify us. */
  eResult = DspA_Verify ( this );
  DebugAssertThrow( (eResult == DspA_tErr_NoErr) );

  /* Add each voxel to the selection. */
  for ( nLineVox = 0; nLineVox < this->mNumLineVoxels; nLineVox++ ) {

    /* Convert to mri idx. */
    Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                             &(this->mLineVoxels[nLineVox]),
                             &MRIIdx );

    tkm_SelectVoxel( &MRIIdx );
  }

  /* Redraw ourselves. */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  DebugCatch;
  DebugCatchError( eResult, DspA_tErr_NoErr, DspA_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

DspA_tErr DspA_WriteLineToLabel ( tkmDisplayAreaRef this,
                                  char*             isFileName) {

  DspA_tErr eResult   = DspA_tErr_NoErr;
  Volm_tErr eVolume   = Volm_tErr_NoErr;
  char          sFileName[tkm_knPathLen] = "";
  float         value                    = 0;
  xVoxel        MRIIdx;
  xVoxel        ras;
  LABEL*        pLabel                   = NULL;
  LABEL_VERTEX* pVertex                  = NULL;
  int           nVertex                   = 0;
  int           eLabel                   = 0;
  int       nLineVox  = 0;

  DebugEnterFunction( ("DspA_WriteLineToLabel( this=%p, isFileName=%s )",
                       this, isFileName) );

  /* verify us. */
  eResult = DspA_Verify ( this );
  DebugAssertThrow( (eResult == DspA_tErr_NoErr) );

  if ( this->mNumLineVoxels <= 0 )
    goto cleanup;

  /* make the file name */
  DebugNote( ("Making file name from %s", isFileName) );
  tkm_MakeFileName( isFileName, tkm_tFileName_Label,
                    sFileName, sizeof(sFileName) );

  /* allocate a label file with that number of voxels, and the passed
     in label file name. */
  DebugNote( ("Allocating label with %d voxels") );
  pLabel = LabelAlloc( this->mNumLineVoxels, NULL, sFileName );
  DebugAssertThrowX( (NULL != pLabel), eResult, tkm_tErr_CouldntAllocate );

  /* set the number of points in the label */
  pLabel->n_points = this->mNumLineVoxels;

  /* Copy subject name. */
  Volm_CopySubjectName( this->mpVolume[tkm_tVolumeType_Main],
                        pLabel->subject_name, 100 );

  /* Add each line voxel. */
  nVertex = 0;
  for ( nLineVox = 0; nLineVox < this->mNumLineVoxels; nLineVox++ ) {

    /* Get the value */
    Volm_GetValueAtIdx( this->mpVolume[tkm_tVolumeType_Main],
                        &(this->mLineVoxels[nLineVox]), &value );

    /* convert mri idx to surface ras. note we may use surface ras
       here because it ignores c_ras, which is what label files
       should to surfacebe comptaible with tksurfer.  */
    if ( tkm_UseRealRAS() ) {
      eVolume =
        Volm_ConvertIdxToRAS( this->mpVolume[tkm_tVolumeType_Main],
                              &(this->mLineVoxels[nLineVox]), &ras );
    } else {
      eVolume =
        Volm_ConvertIdxToMRIIdx( this->mpVolume[tkm_tVolumeType_Main],
                                 &(this->mLineVoxels[nLineVox]), &MRIIdx );
      eVolume =
        Volm_ConvertMRIIdxToSurfaceRAS( this->mpVolume[tkm_tVolumeType_Main],
                                        &MRIIdx, &ras );
    }
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                       eResult, tkm_tErr_ErrorAccessingVolume );

    /* get a ptr the vertex in the label file. */
    DebugNote( ("Geting ptr to vertex %d", nVertex) );
    pVertex = &(pLabel->lv[nVertex]);

    /* set the vertex */
    DebugNote( ("Setting values of vertex %d", nVertex) );
    pVertex->x = xVoxl_GetFloatX( &ras );
    pVertex->y = xVoxl_GetFloatY( &ras );
    pVertex->z = xVoxl_GetFloatZ( &ras );

    /* set the vno to -1, which is significant somewhere outside
       the realm of tkmedit. set stat value to the mri value. */
    pVertex->vno = -1;
    pVertex->stat = value;
    pVertex->deleted = FALSE;

    /* inc our global count. */
    nVertex++;
  }

  /* write the file */
  DebugNote( ("Writing label file to %s", sFileName) );
  eLabel = LabelWrite( pLabel, sFileName );
  DebugAssertThrowX( (NO_ERROR == eLabel),
                     eResult, tkm_tErr_CouldntWriteFile );

  DebugCatch;
  DebugCatchError( eResult, DspA_tErr_NoErr, DspA_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


DspA_tErr DspA_Verify ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this ) {
    eResult = DspA_tErr_InvalidPtr;
    goto cleanup;
  }

  /* check signature */
  if ( DspA_kSignature != this->mSignature ) {
    eResult = DspA_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_VerifyVolumeVoxel_  ( tkmDisplayAreaRef this,
                                     xVoxelRef          ipVoxel ) {
  DspA_tErr eResult = DspA_tErr_NoErr;

  /* make sure voxel is in bounds */
  if ( xVoxl_GetX( ipVoxel )    < this->mnMinVolumeIndexX
       || xVoxl_GetY( ipVoxel ) < this->mnMinVolumeIndexY
       || xVoxl_GetZ( ipVoxel ) < this->mnMinVolumeIndexZ
       || xVoxl_GetX( ipVoxel ) >= this->mnMaxVolumeIndexX
       || xVoxl_GetY( ipVoxel ) >= this->mnMaxVolumeIndexY
       || xVoxl_GetZ( ipVoxel ) >= this->mnMaxVolumeIndexZ ) {

    eResult = DspA_tErr_InvalidVolumeVoxel;
  }

  return eResult;
}

DspA_tErr DspA_VerifyScreenPoint_ ( tkmDisplayAreaRef this,
                                    xPoint2nRef       ipScreenPt ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* make sure screen point is in screen bounds */
  if ( ipScreenPt->mnX    < 0
       || ipScreenPt->mnY < 0
       || ipScreenPt->mnX >= this->mnWidth
       || ipScreenPt->mnY >= this->mnHeight ) {

    eResult = DspA_tErr_InvalidScreenPoint;
  }

  return eResult;
}

DspA_tErr DspA_VerifyBufferPoint_ ( tkmDisplayAreaRef this,
                                    xPoint2nRef       ipBufferPt ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* make sure the buffer point is in bounds */
  if ( ipBufferPt->mnX    < 0
       || ipBufferPt->mnY < 0
       || ipBufferPt->mnX >= this->mnVolumeSizeX
       || ipBufferPt->mnY >= this->mnVolumeSizeY ) {

    eResult = DspA_tErr_InvalidBufferPoint;
  }

  return eResult;
}

void DspA_SetUpOpenGLPort_ ( tkmDisplayAreaRef this ) {

  GLenum eGL;

  /* copy the frame buffer to the screen. */
  glMatrixMode   ( GL_PROJECTION );
  glLoadIdentity ();

  /* glOrtho ( left, right, bottom, top, near, far )
     ( left, bottom, -near ) is lower left corner of window
     ( right, top, -near ) is upper right */

  glOrtho        ( 0, this->mnVolumeSizeX, this->mnVolumeSizeY, 0, -1.0, 1.0 );
  //glOrtho ( 0, this->mnVolumeSizeX, 0, this->mnVolumeSizeX, -1.0, 1.0 );
  eGL = glGetError ();
  if ( GL_NO_ERROR != eGL )
    DebugPrint( ("glOrtho got error %d\n", eGL ) );

  glViewport( this->mLocationInSuper.mnX, this->mLocationInSuper.mnY,
              this->mnWidth, this->mnHeight );
  eGL = glGetError ();
  if ( GL_NO_ERROR != eGL )
    DebugPrint( ("glViewport got error %d\n", eGL ) );

  glRasterPos2i  ( 0, this->mnVolumeSizeX );
  //glRasterPos2i  ( 0, 0 );
  eGL = glGetError ();
  if ( GL_NO_ERROR != eGL )
    DebugPrint( ("glRasterPos2i got error %d\n", eGL ) );

  glPixelZoom ( this->mfFrameBufferScaleX, this->mfFrameBufferScaleY );
  eGL = glGetError ();
  if ( GL_NO_ERROR != eGL )
    DebugPrint( ("glPixelZoom got error %d\n", eGL ) );
}

void DspA_DebugPrint_ ( tkmDisplayAreaRef this ) {

  DebugPrint( ("tkmDisplayArea\n" ) );
  DebugPrint( ("\tx %d y %d w %d h %d\n",
               this->mLocationInSuper.mnX, this->mLocationInSuper.mnY,
               this->mnWidth, this->mnHeight ) );
  DebugPrint( ("\tvolume size %d, x scale %.2f y scale %.2f\n",
               this->mnVolumeSizeX, this->mfFrameBufferScaleX,
               this->mfFrameBufferScaleY ) );
  DebugPrint( ("\tzoom level %d center %d, %d, %d\n",
               this->mnZoomLevel, xVoxl_ExpandInt(this->mpZoomCenter) ) );
  DebugPrint( 
    ("\tcursor %d, %d, %d orientation %s slice %d tool %d\n",
     xVoxl_ExpandInt(this->mpCursor), DspA_ksaOrientation[this->mOrientation],
     DspA_GetCurrentSliceNumber_(this), sTool ) );
}

void DspA_Signal ( char* isFuncName, int inLineNum, DspA_tErr ieCode ) {

  DebugPrint( ("Signal in %s, line %d: %d, %s\n",
               isFuncName, inLineNum, ieCode, DspA_GetErrorString(ieCode) ) );
}

char* DspA_GetErrorString ( DspA_tErr ieCode ) {

  DspA_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= DspA_knNumErrorCodes ) {
    eCode = DspA_tErr_InvalidErrorCode;
  }

  return DspA_ksaErrorStrings [eCode];
}




DspA_tErr DspA_NormalizeVoxel_ ( xVoxelRef        ipAnaIdx,
                                 mri_tOrientation iOrientation,
                                 xVoxelRef        opNormIdx ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  switch ( iOrientation ) {
  case mri_tOrientation_Horizontal:
    xVoxl_SetFloat( opNormIdx,
                    xVoxl_GetFloatX(ipAnaIdx),
                    xVoxl_GetFloatZ(ipAnaIdx),
                    xVoxl_GetFloatY(ipAnaIdx) );
    break;
  case mri_tOrientation_Coronal:
    xVoxl_SetFloat( opNormIdx,
                    xVoxl_GetFloatX(ipAnaIdx),
                    xVoxl_GetFloatY(ipAnaIdx),
                    xVoxl_GetFloatZ(ipAnaIdx) );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetFloat( opNormIdx,
                    xVoxl_GetFloatZ(ipAnaIdx),
                    xVoxl_GetFloatY(ipAnaIdx),
                    xVoxl_GetFloatX(ipAnaIdx) );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_NormalizeVoxel: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_UnnormalizeVoxel_ ( xVoxelRef        ipNormIdx,
                                   mri_tOrientation iOrientation,
                                   xVoxelRef        opAnaIdx ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  switch ( iOrientation ) {
  case mri_tOrientation_Horizontal:
    xVoxl_SetFloat( opAnaIdx,
                    xVoxl_GetFloatX(ipNormIdx),
                    xVoxl_GetFloatZ(ipNormIdx),
                    xVoxl_GetFloatY(ipNormIdx) );
    break;
  case mri_tOrientation_Coronal:
    xVoxl_SetFloat( opAnaIdx,
                    xVoxl_GetFloatX(ipNormIdx),
                    xVoxl_GetFloatY(ipNormIdx),
                    xVoxl_GetFloatZ(ipNormIdx) );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetFloat( opAnaIdx,
                    xVoxl_GetFloatZ(ipNormIdx),
                    xVoxl_GetFloatY(ipNormIdx),
                    xVoxl_GetFloatX(ipNormIdx) );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_UnnormalizeVoxel: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

tBoolean xUtil_LineIntersectsPlane ( xVoxelRef         ipAnaIdxA,
                                     xVoxelRef         ipAnaIdxB,
                                     int               inPlane,
                                     xPoint2fRef       opIntersectionPt,
                                     xPoint2fRef     opInterpIntersectionPt ) {

  float    fPlane           = inPlane;
  float    fDistanceA       = 0;
  float    fDistanceB       = 0;
  float    fAlpha           = 1;
  xPoint2f intersectionPt   = {0, 0};
  xPoint2f interpIntersectionPt   = {0, 0};

  /* get distance from each to plane. */
  fDistanceA = xVoxl_GetFloatZ( ipAnaIdxA ) - fPlane;
  fDistanceB = xVoxl_GetFloatZ( ipAnaIdxB ) - fPlane;

  /* if product is negative or 0, they intersect the plane. */
  if ( fDistanceA * fDistanceB > 0.0 ) {
    return FALSE;
  }

  /* first find the uninterpolated point. the intersection is just the
     projection onto plane */
  intersectionPt.mfX = xVoxl_GetFloatX( ipAnaIdxA );
  intersectionPt.mfY = xVoxl_GetFloatY( ipAnaIdxA );


  /* now find the interpolated intersection. find an iterpolation
     factor, which is the intersection of the line the plane. */
  /* make sure they arn't on the same plane... */
  if ( xVoxl_GetFloatZ(ipAnaIdxB) - xVoxl_GetFloatZ(ipAnaIdxA) != 0.0 ) {

    fAlpha = (fPlane - xVoxl_GetFloatZ( ipAnaIdxA )) /
             (xVoxl_GetFloatZ( ipAnaIdxB ) - xVoxl_GetFloatZ( ipAnaIdxA ));

  } else {

    /* alpha is just 1. */
    fAlpha = 1.0;
  }

  /* interpolate to find the intersection. */
  interpIntersectionPt.mfX = 
    (xVoxl_GetFloatX( ipAnaIdxA ) +
     fAlpha * (xVoxl_GetFloatX(ipAnaIdxB) - xVoxl_GetFloatX(ipAnaIdxA)));
  interpIntersectionPt.mfY = 
    (xVoxl_GetFloatY( ipAnaIdxA ) +
     fAlpha * (xVoxl_GetFloatY(ipAnaIdxB) - xVoxl_GetFloatY(ipAnaIdxA)));

  /* return the points. */
  *opIntersectionPt = intersectionPt;
  *opInterpIntersectionPt = interpIntersectionPt;

  return TRUE;
}


DspA_tErr DspA_AdjustSurfaceAnaIdx ( tkmDisplayAreaRef this,
                                     xVoxelRef         iAnaIdx ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* these are .51 because tkmDisplayArea automatically moves cursors
     on .0 boundaries to .5 boundaries. when looking at an orig
     vertex, the bounds are on .0, and the cursor misses the
     vertex. .51 prevents that from happening. */
  xVoxl_SetFloatX( iAnaIdx, xVoxl_GetFloatX(iAnaIdx) + 0.51 );
  xVoxl_SetFloatY( iAnaIdx, xVoxl_GetFloatY(iAnaIdx) + 0.51 );
  xVoxl_SetFloatZ( iAnaIdx, xVoxl_GetFloatZ(iAnaIdx) + 0.51 );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_AdjustSurfaceAnaIdx: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_UnadjustSurfaceAnaIdx ( tkmDisplayAreaRef this,
                                       xVoxelRef         iAnaIdx ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  xVoxl_SetFloatX( iAnaIdx, xVoxl_GetFloatX(iAnaIdx) - 0.50 );
  xVoxl_SetFloatY( iAnaIdx, xVoxl_GetFloatY(iAnaIdx) - 0.50 );
  xVoxl_SetFloatZ( iAnaIdx, xVoxl_GetFloatZ(iAnaIdx) - 0.50 );

  goto cleanup;

error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_UnadjustSurfaceAnaIdx: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

DspA_tErr DspA_ParsePointList_( tkmDisplayAreaRef this,
                                GLenum            inMode,
                                xGrowableArrayRef iList,
                                float             ifaColor[3] ) {

  DspA_tErr  eResult        = DspA_tErr_NoErr;
  xGArr_tErr eList          = xGArr_tErr_NoErr;
  tBoolean   bOperationOpen = FALSE;
  DspA_tSurfaceListNode drawListNode;
  xPoint2f   drawPoint      = {0,0};
  float      faColor[3]     = {0, 0, 0};

  /* reset the list position. */
  eList = xGArr_ResetIterator( iList );
  if ( xGArr_tErr_NoErr != eList )
    goto error;

  /* start new operation. */
  glBegin( inMode );
  bOperationOpen = TRUE;

  /* each point has an uninterpolated point and an interpolated point
     in succession. there may be any number of these, and then a face
     marker, represented by the point -1, -1. */
  while ( (eList = xGArr_NextItem( iList, (void*)&drawListNode ))
          == xGArr_tErr_NoErr ) {

    /* if it's not a vertex... */
    if ( !drawListNode.mbVertex ) {

      /* if operation was still going, end it. */
      if ( bOperationOpen )
        glEnd();

      /* start new operation. */
      glBegin( inMode );
      bOperationOpen = TRUE;

    } else {

      /* switch on the draw flag to see which point to get. */
      if (this->mabDisplayFlags[DspA_tDisplayFlag_InterpolateSurfaceVertices]) 
      {
        drawPoint.mfX = drawListNode.mInterpIntersectionPoint.mfX;
        drawPoint.mfY = drawListNode.mInterpIntersectionPoint.mfY;
      } else {
        drawPoint.mfX = drawListNode.mIntersectionPoint.mfX;
        drawPoint.mfY = drawListNode.mIntersectionPoint.mfY;
      }

      /* convert to zoomed coords. */
      drawPoint.mfX = 
        ((float)this->mnZoomLevel * 
         (drawPoint.mfX - xVoxl_GetFloatX(this->mpZoomCenter))) + 
        (float)(this->mnVolumeSizeX/2.0);
      drawPoint.mfY = 
        ((float)this->mnZoomLevel * 
         (drawPoint.mfY - xVoxl_GetFloatY(this->mpZoomCenter))) + 
        (float)(this->mnVolumeSizeY/2.0);

      /* y flip */
      drawPoint.mfY = GLDRAW_Y_FLIP(drawPoint.mfY);

      /* If we have an override color, set it. */
      if ( drawListNode.mOverrideColor ) {
        xColr_PackFloatArray( &(drawListNode.mColor), faColor );
        glColor3fv( faColor );
      } else {
        glColor3fv( ifaColor );
      }

      /* and draw the pt. */
      glVertex2f( drawPoint.mfX, drawPoint.mfY );
    }
  }

  /* clear flag if last item */
  if ( eList == xGArr_tErr_LastItem )
    eList = xGArr_tErr_NoErr;

  /* check for other errors */
  if ( eList != xGArr_tErr_NoErr )
    goto error;

  /* end last operation */
  if ( bOperationOpen )
    glEnd();

  goto cleanup;

error:

  if ( xGArr_tErr_NoErr != eList )
    eResult = DspA_tErr_ErrorAccessingSurfaceList;

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_ParsePointList_: %s\n",
                 eResult, DspA_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}
