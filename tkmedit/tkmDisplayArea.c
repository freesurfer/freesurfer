#include "tkmDisplayArea.h"
#include "tkmMeditWindow.h"
#include "tkmFunctionalVolume.h"
#include "xUtilities.h"
#include "utils.h"

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
static DspA_tTool              sTool  = DspA_tTool_SelectVoxels;
static DspA_tBrushSettings     sBrush;
static DspA_tParcBrushSettings sParcBrush;

/* cursor info too */
static xColor3f     sCursorColor = { 1, 0, 0 };
static DspA_tMarker sCursorShape = DspA_tMarker_Crosshair;

/* static focused display. */
static tkmDisplayAreaRef sFocusedDisplay = NULL;

char DspA_ksaErrorStrings [DspA_knNumErrorCodes][STRLEN] = {
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
  "Invalid error code."
};

char DspA_ksaOrientation [mri_knNumOrientations][STRLEN] = {
  "Coronal",
  "Horizontal",
  "Sagittal"
};

char DspA_ksaSurface [Surf_knNumVertexSets][STRLEN] = {
  "Main",
  "Original",
  "Canonical"
};

char DspA_ksaDisplaySet [DspA_knNumDisplaySets][STRLEN] = {
  "cursor",
  "mouseover"
};

DspA_tErr DspA_New ( tkmDisplayAreaRef* oppWindow,
		     tkmMeditWindowRef  ipWindow ) {
  
  DspA_tErr         eResult      = DspA_tErr_NoErr;
  tkmDisplayAreaRef this         = NULL;
  int               nFlag        = 0;
  int               nSurface     = 0;
  int               nDTI         = 0;
  xColor3f          color;
  
  /* allocate us. */
  this = (tkmDisplayAreaRef) malloc( sizeof(tkmDisplayArea) );
  if( NULL == this ) {
    eResult = DspA_tErr_AllocationFailed;
    goto error;
  }
  
  /* set the signature */
  this->mSignature = DspA_kSignature;
  
  /* set the parent window. */
  this->mpWindow             = ipWindow;
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
  this->mnROIGroupIndex        = -1;
  for( nSurface = 0; nSurface < Surf_knNumVertexSets; nSurface++ )
    DspA_SetSurfaceLineWidth( this, nSurface, 1 );
  xColr_Set( &color, 1, 1, 0 );
  DspA_SetSurfaceLineColor( this, Surf_tVertexSet_Main, &color );
  xColr_Set( &color, 0, 1, 0 );
  DspA_SetSurfaceLineColor( this, Surf_tVertexSet_Original, &color );
  xColr_Set( &color, 1, 0, 0 );
  DspA_SetSurfaceLineColor( this, Surf_tVertexSet_Pial, &color );
  
  
  /* all our display flags start out false. */
  for( nFlag = 0; nFlag < DspA_knNumDisplayFlags; nFlag++ )
    this->mabDisplayFlags[nFlag] = FALSE;
  
  /* null ptrs for display data. */
  this->mpVolume                = NULL;
  this->mpAuxVolume             = NULL;
  this->mROIGroup               = NULL;
  for( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {
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
  for( nDTI = 0; nDTI < tkm_knNumDTIVolumeTypes; nDTI++ ) {
    this->mpDTIVolume[nDTI] = NULL;
  }
  
  /* set default brush info */
  sBrush.mnRadius = 1;
  sBrush.mShape   = DspA_tBrushShape_Circle;
  sBrush.mb3D     = FALSE;
  DspA_SetBrushInfoToDefault( this, DspA_tBrush_EditOne );
  DspA_SetBrushInfoToDefault( this, DspA_tBrush_EditTwo );
  
  /* default parc brush info */
  sParcBrush.mNewValue  = 0;
  sParcBrush.mb3D       = FALSE;
  sParcBrush.mSrc       = tkm_tVolumeType_Main;
  sParcBrush.mnFuzzy    = 0;
  sParcBrush.mnDistance = 0;
  
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
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* delete frame buffer */
  if( NULL != this->mpFrameBuffer )
    free( this->mpFrameBuffer );
  
  /* delete our voxels */
  xVoxl_Delete( &this->mpLastCursor );
  xVoxl_Delete( &this->mpCursor );
  xVoxl_Delete( &this->mpZoomCenter );
  
  /* trash the signature */
  this->mSignature = 0x1;
  
  /* delete us */
  free( this );
  
  /* return null */
  *ioppWindow = NULL;
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* set our location */
  this->mLocationInSuper = iLocation;
  
  /* set our size */
  this->mnWidth  = inWidth;
  this->mnHeight = inHeight;
  
  /* set our scale */
  if( this->mnVolumeSizeX > 0 ) {
    
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
  if( DspA_tErr_NoErr != eResult ) {
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
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* get the title information */
  Volm_CopySubjectName( this->mpVolume, sSubjectName, sizeof(sSubjectName) );
  Volm_CopyVolumeName( this->mpVolume, sVolumeName, sizeof(sVolumeName) );
  if( NULL != this->mpAuxVolume ) {
    Volm_CopyVolumeName( this->mpAuxVolume, 
			 sAuxVolumeName, sizeof(sAuxVolumeName) );
  }
  
  /* if we don't have an aux volume */
  if( NULL == this->mpAuxVolume ) {
    
    /* just use the subject and volume name */
    sprintf( sTitle, "%s: %s", sSubjectName, sVolumeName );
    
  } else {
    
    /* else see which one is displayed. use that name first and then
       the other name in parens. */
    if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      sprintf( sTitle, "%s: %s (** %s **)", 
	       sSubjectName, sVolumeName, sAuxVolumeName );
    } else {
      sprintf( sTitle, "%s: ** %s ** (%s)", 
	       sSubjectName, sVolumeName, sAuxVolumeName );
    }
  }
  
  /* set the window name. */
  MWin_SetWindowTitle( this->mpWindow, sTitle );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_UpdateWindowTitle: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
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
  
  DebugEnterFunction( ("DspA_SetVolume( this=%p, ipVolume=%p, inSize=%d,%d,%d )", 
		       this, ipVolume, inSizeX, inSizeY, inSizeZ) );
  
  xVoxl_New( &pCenter );
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the main volume */
  this->mpVolume = ipVolume;
  
  /* save the volume size */
  this->mnVolumeSizeX = inSizeX;
  this->mnVolumeSizeY = inSizeY;
  this->mnVolumeSizeZ = inSizeZ;
  
  /* if we alreayd have a frame buffer, delete it */
  if( NULL == this->mpFrameBuffer ) {
    free( this->mpFrameBuffer );
    this->mpFrameBuffer = NULL;
  }
  
  /* allocate a new one */
  nSize = MAX(MAX(this->mnVolumeSizeX, this->mnVolumeSizeY), this->mnVolumeSizeZ) ;
  this->mpFrameBuffer = (GLubyte*) malloc( nSize * nSize *
					   DspA_knNumBytesPerPixel );
  if( NULL == this->mpFrameBuffer ) {
    eResult = DspA_tErr_AllocationFailed;
    goto error;
  }
  
  /* set inital values for our buffer scale */
  this->mfFrameBufferScaleX = 
    (float)this->mnWidth  / (float)this->mnVolumeSizeX;
  this->mfFrameBufferScaleY = 
    (float)this->mnHeight  / (float)this->mnVolumeSizeY;
  
  /* initialize surface point lists. */
  eResult = DspA_InitSurfaceLists_( this, 
				    nSize * Surf_knNumVertexSets * mri_knNumOrientations );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* update window title */
  DspA_UpdateWindowTitle( this );
  
  /* send volume name */
  Volm_CopyVolumeName( this->mpVolume, sVolumeName, sizeof(sVolumeName) );
  sprintf( sTclArguments, "\"%s value\"", sVolumeName );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );
  
  /* get the center of the volume */
  xVoxl_Set( pCenter, floor(this->mnVolumeSizeX/2), 
	     floor(this->mnVolumeSizeY/2), floor(this->mnVolumeSizeZ/2) );
  
  /* set cursor and zoom center to middle of volume if not already set */
  if (xVoxl_GetX(this->mpCursor)==0 &&
      xVoxl_GetY(this->mpCursor)==0 &&
      xVoxl_GetZ(this->mpCursor)==0)
    {
      eResult = DspA_SetCursor( this, pCenter );
      if( DspA_tErr_NoErr != eResult )
	goto error;
      eResult = DspA_SetZoomCenter( this, pCenter );
      if( DspA_tErr_NoErr != eResult )
	goto error;
    }
  else
    {
      eResult = DspA_SetCursor( this, this->mpCursor );
      if( DspA_tErr_NoErr != eResult )
	goto error;
    }
  
  
  /* show cursor. */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Cursor, TRUE );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* show volume. */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Anatomical, TRUE );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* set dirty flag and redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* make sure this size is the same as the main volume size. */
  if( inSizeX != this->mnVolumeSizeX ||
      inSizeY != this->mnVolumeSizeY ||
      inSizeZ != this->mnVolumeSizeZ ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* save the aux volume */
  this->mpAuxVolume = ipVolume;
  
  /* if we got a volume... */
  if( NULL != this->mpAuxVolume ) {
    
    /* send volume name to tk window */
    Volm_CopyVolumeName( this->mpAuxVolume, sVolumeName, sizeof(sVolumeName) );
    xUtil_snprintf( sTclArguments, sizeof(sTclArguments),
		    "\"%s value\"", sVolumeName );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeName, sTclArguments );
    
    /* show the volume value */
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxValue, "1" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeOptions, "1" );
    
    /* if we're currently showing the aux volume, the slice has changed */
    if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      this->mbSliceChanged = TRUE;
    }
    
    /* if we're focused, send the new information for the cursor */
    if( sFocusedDisplay == this ) {
      DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
				       this->mpCursor );
    }
  } else {
    
    /* hide the volume value */
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxValue, "0" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeOptions, "0" );
    
    /* don't show the aux volume */
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_AuxVolume, FALSE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }    
  
  /* update window title */
  DspA_UpdateWindowTitle( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetAuxVolume: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

DspA_tErr DspA_SetROIGroup ( tkmDisplayAreaRef this,
			     mriVolumeRef      iGroup ) {
  
  DspA_tErr eResult            = DspA_tErr_NoErr;
  tBoolean  bHaveGroup         = FALSE;
  char      sTclArguments[STRLEN] = "";
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the group */
  this->mROIGroup = iGroup;
  
  /* turn stuff on or off based on if we have one. */
  if( this->mROIGroup != NULL ) {
    bHaveGroup = TRUE;
  }
  
  /* turn roi group on */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_ROIGroupOverlay,
				 bHaveGroup );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* show roi label */
  sprintf( sTclArguments, "%d", (int)bHaveGroup );
  tkm_SendTclCommand( tkm_tTclCommand_ShowROILabel, sTclArguments );
  tkm_SendTclCommand( tkm_tTclCommand_ShowROIGroupOptions, sTclArguments );
  
  /* if we're focused, send the new information for the cursor */
  if( sFocusedDisplay == this ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
				     this->mpCursor );
  }
  
  /* redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetROIGroup: %s\n",
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
  if( DspA_tErr_NoErr != eResult )
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
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* set the volume */
  this->mpFunctionalVolume = ipVolume;
  
  /* see if there's data... */
  eFunctional = FunV_IsOverlayPresent( this->mpFunctionalVolume, 
				       &bOverlayLoaded );
  if( FunV_tErr_NoError != eFunctional ) {
    DspA_Signal( "DspA_SetOverlayVolume", __LINE__, 
		 DspA_tErr_ErrorAccessingFunctionalVolume );
  }
  
  /* turn functional data and color scale bar on if there is */
  if( bOverlayLoaded ) {
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay,
				   TRUE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    eResult = DspA_SetDisplayFlag( this, 
				   DspA_tDisplayFlag_FunctionalColorScaleBar, TRUE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  /* if we're focused, send the new information for the cursor */
  if( sFocusedDisplay == this ) {
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
				     this->mpCursor );
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the control points space. */
  this->mpControlPoints = ipVoxels;
  
  /* turn control points on. */
  if( NULL != this->mpControlPoints ) {
    
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_ControlPoints, 
				   TRUE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }    
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetControlPointsSpace: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

DspA_tErr DspA_SetSelectionSpace( tkmDisplayAreaRef this, 
				  x3DListRef        ipVoxels ) {
  
  DspA_tErr eResult = DspA_tErr_NoErr;
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the selected voxels */
  this->mpSelection = ipVoxels;
  
  /* turn selection display on. */
  if( NULL != this->mpSelection ) {
    
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Selection, TRUE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  /* find the center of the selection and go there */
  eResult = DspA_SetCursorToCenterOfSpace( this, ipVoxels );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the selected voxels */
  this->mHeadPoints = iList;
  
  /* turn selection display on. */
  if( NULL != this->mHeadPoints ) {
    
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_HeadPoints, TRUE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    bHaveList = TRUE;
  }
  
  sprintf( sTclArguments, "%d", (int)bHaveList );
  tkm_SendTclCommand( tkm_tTclCommand_ShowHeadPointLabel, sTclArguments );
  
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the VLIs and TRANSFORM */
  this->mVLI1    = iVLI1;
  this->mVLI2    = iVLI2;
  strcpy(this->isVLI1_name, isVLI1_name) ;
  strcpy(this->isVLI2_name, isVLI2_name) ;
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetVLIs: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}


DspA_tErr DspA_SetDTIVolume  ( tkmDisplayAreaRef this,
			       tkm_tDTIVolumeType  iType,
			       mriVolumeRef        iVolume ) {
  
  
  DspA_tErr eResult             = DspA_tErr_NoErr;
  char       sTclArguments[tkm_knTclCmdLen] = "";
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  if( iType < 0 || iType >= tkm_knNumDTIVolumeTypes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* Save the volume. */
  this->mpDTIVolume[iType] = iVolume;
  
  /* Show the DTI options if we got a volume.  */
  sprintf( sTclArguments, "%d", (int)(NULL != iVolume) );
  tkm_SendTclCommand( tkm_tTclCommand_ShowDTIOptions, sTclArguments );
  
  /* turn DTI display on. */
  if( NULL != this->mpDTIVolume[iType] ) {
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_DTIOverlay, TRUE );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* save the GCA and LTA */
  this->mGCAVolume    = iVolume;
  this->mGCATransform = iTransform;
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* verify the cursor */
  eResult = DspA_VerifyVolumeVoxel_( this, ipCursor );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  
  /* allow the functional display to respond. */
  if( NULL != this->mpFunctionalVolume ) {
    eFunctional = FunV_AnatomicalVoxelClicked( this->mpFunctionalVolume,
					       ipCursor );
  }
  
  
  /* get our current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );
  
  /* Copy the current cursor into the last cursor. */
  xVoxl_Copy( this->mpLastCursor, this->mpCursor );
  
  /* set the cursor */
  xVoxl_Copy( this->mpCursor, ipCursor );
  
  /* if the new slice number is diffrent, set our dirty slice flag. */
  if( DspA_GetCurrentSliceNumber_( this ) != nSlice ) {
    this->mbSliceChanged = TRUE;
  }
  
  /* if cursor is .0 .0, change to .5 .5 so that it will draw in the
     center of a voxel on screen */
  DspA_ConvertVolumeToPlane_( this, this->mpCursor, this->mOrientation,
			      &planePt, &nSlice );
  if( planePt.mfX == (float)(int)planePt.mfX &&
      planePt.mfY == (float)(int)planePt.mfY ) {
    planePt.mfX += 0.5;
    planePt.mfY += 0.5;
    DspA_ConvertPlaneToVolume_( this, &planePt, nSlice,
				this->mOrientation, this->mpCursor );
  }
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* notify the window that the cursor has changed. */
    MWin_CursorChanged( this->mpWindow, this, this->mpCursor );
    
    /* send the information for this point */
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
				     this->mpCursor );
  }
  
  /* if we have head point data... */
  if(  NULL != this->mHeadPoints ) {
    
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
  if( DspA_tErr_NoErr != eResult ) {
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
  xVoxel     anaIdx;
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* convert the coord to the right space. */
  switch( iFromSpace ) {
  case mri_tCoordSpace_VolumeIdx:
    xVoxl_Copy( &anaIdx, ipCoord );
    break;
  case mri_tCoordSpace_RAS:
    Volm_ConvertRASToIdx( this->mpVolume, ipCoord, &anaIdx );
    break;
  case mri_tCoordSpace_Talairach:
    Volm_ConvertTalToIdx( this->mpVolume, ipCoord, &anaIdx );
    break;
  default:
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set the cursor */
  eResult = DspA_SetCursor( this, &anaIdx );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_ConvertAndSetCursor(%d, %d,%d,%d): %s\n",
		 eResult, (int)iFromSpace, xVoxl_ExpandInt(ipCoord),
		 DspA_GetErrorString(eResult) ) );
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* ignore if this is already our slice number */
  if( inSlice == DspA_GetCurrentSliceNumber_( this ) )
    goto cleanup;
  
  /* copy the cursor. */
  xVoxl_Copy( pCursor, this->mpCursor );
  
  /* change the slice */
  switch( this->mOrientation ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* update mouse over information if we're focused */
  if( sFocusedDisplay == this ) {
    eResult = DspA_SendMouseInfoToTcl( this );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* verify the orientation */
  if( iOrientation <= mri_tOrientation_None
      || iOrientation >= mri_knNumOrientations ) {
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
  }
  
  /* if the orientations are different, our slice will have changed. */
  if( this->mOrientation != iOrientation ) {
    this->mbSliceChanged = TRUE;
  }
  
  /* set the orientation */
  this->mOrientation = iOrientation;
  
  /* when zoomed, set the zoom center to the cursor so it appears to
     reorient around the cursor (zoom centers are only 2d and not valid
     when switching planes. */
  if( this->mnZoomLevel != DspA_knMinZoomLevel ) {
    eResult = DspA_SetZoomCenterToCursor( this );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
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
  if( DspA_tErr_NoErr != eResult ) {
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
  xVoxelRef pCenter            = NULL;
  int       nNewLevel          = 0;
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* verify the zoom level. */
  nNewLevel = this->mnZoomLevel;
  if( inLevel >= DspA_knMinZoomLevel 
      && inLevel <= DspA_knMaxZoomLevel ) {
    nNewLevel = inLevel;
  }
  
  DspA_SetZoomCenterToCursor( this );
  
  /* set our zoom level. */
  this->mnZoomLevel = nNewLevel;
  
  /* if this is zoom level one, set our center to the middle of
     this slice. */
  if( 1 == this->mnZoomLevel ) {
    
    /* set a voxel to 128,128,128, setting the zoom center will not change
       the current slice. */
    xVoxl_New( &pCenter );
    xVoxl_Set( pCenter, this->mnVolumeSizeX/2, 
	       this->mnVolumeSizeY/2, this->mnVolumeSizeZ/2 );
    eResult = DspA_SetZoomCenter( this, pCenter );
    xVoxl_Delete( &pCenter );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  /* this requires a rebuild of the frame buffer. */
  this->mbSliceChanged = TRUE;
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
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
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* verify the center */
  eResult = DspA_VerifyVolumeVoxel_( this, ipCenter );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  switch( this->mOrientation ) {
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
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* set the center to our cursor */
  DspA_SetZoomCenter( this, this->mpCursor );
  
  /* schedule a redraw */
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* set the values */
  this->mHilitedSurface        = inSurface;
  this->mnHilitedVertexIndex   = inVertex;
  
  /* schedule a redraw */
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HiliteSurfaceVertex: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}


DspA_tErr DspA_SetDisplayFlag ( tkmDisplayAreaRef this,
				DspA_tDisplayFlag iWhichFlag,
				tBoolean          ibNewValue ) {
  
  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  tBoolean  bNewValue          = FALSE;
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* check the flag */
  if( iWhichFlag <= DspA_tDisplayFlag_None
      || iWhichFlag >= DspA_knNumDisplayFlags ) {
    eResult = DspA_tErr_InvalidDisplayFlag;
    goto error;
  }
  
  /* get the new value. */
  bNewValue = ibNewValue;
  
  /* verify the flag. */
  switch( iWhichFlag ) {
    
  case DspA_tDisplayFlag_AuxVolume:
    
    /* if no aux volume, set to false. */
    if( NULL == this->mpAuxVolume )
      bNewValue = FALSE;
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    /* change the window title. */
    DspA_UpdateWindowTitle( this );
    
    break;
    
  case DspA_tDisplayFlag_Anatomical:
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_ROIGroupOverlay:
  case DspA_tDisplayFlag_ROIVolumeCount:
    
    /* if no roi group, set to false. */
    if( NULL == this->mROIGroup )
      bNewValue = FALSE;
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_MainSurface:
  case DspA_tDisplayFlag_OriginalSurface:
  case DspA_tDisplayFlag_CanonicalSurface:
  case DspA_tDisplayFlag_DisplaySurfaceVertices:
    
    /* if no surface, set to false. */
    if( NULL == this->mpSurface[tkm_tSurfaceType_Main] ) 
      bNewValue = FALSE;
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_InterpolateSurfaceVertices:
    
    /* if no surface, set to false. */
    if( NULL == this->mpSurface[tkm_tSurfaceType_Main] ) 
      bNewValue = FALSE;
    
    /* if the flag is different, purge lists and set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue ) {
      this->mbSliceChanged = TRUE;
      eResult = DspA_PurgeSurfaceLists_( this );
      if( DspA_tErr_NoErr != eResult )
	goto error;
    }
    
    break;
    
  case DspA_tDisplayFlag_FunctionalOverlay:
  case DspA_tDisplayFlag_FunctionalColorScaleBar:
  case DspA_tDisplayFlag_MaskToFunctionalOverlay:
    
    /* if no func data, set to false. */
    if( NULL == this->mpFunctionalVolume ) {
      bNewValue = FALSE;
    }
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_HistogramPercentChange:
    
    /* if no VLI data, set to false. */
    if( NULL == this->mVLI1 ) {
      bNewValue = FALSE;
    }
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_Selection:
  case DspA_tDisplayFlag_MaxIntProj:
  case DspA_tDisplayFlag_UndoVolume:
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_HeadPoints:
    
    /* check to see if we have a list */
    if( NULL == this->mHeadPoints )
      bNewValue = FALSE;
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  case DspA_tDisplayFlag_DTIOverlay:
    
    /* if no DTI data, set to false. */
    if( NULL == this->mpDTIVolume[tkm_tDTIVolumeType_X] ||
	NULL == this->mpDTIVolume[tkm_tDTIVolumeType_Y] ||
	NULL == this->mpDTIVolume[tkm_tDTIVolumeType_Z] ||
	NULL == this->mpDTIVolume[tkm_tDTIVolumeType_FA] ) {
      bNewValue = FALSE;
    }
    
    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;
    
    break;
    
  default:
    break;
  }
  
  /* set the value */
  this->mabDisplayFlags[iWhichFlag] = bNewValue;
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
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
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
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
  if( sFocusedDisplay == this ) {
    
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
  if( inRadius >= DspA_knMinBrushRadius 
      && inRadius <= DspA_knMaxBrushRadius ) {
    sBrush.mnRadius = inRadius;
  }
  
  if( iShape >= 0 
      && iShape < DspA_knNumBrushShapes ) {
    sBrush.mShape = iShape;
  }
  
  /* Set the brush info */
  sBrush.mb3D = ib3D;
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d",
	      (int)sBrush.mnRadius, (int)sBrush.mShape, (int)sBrush.mb3D );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushShape, sTclArguments );
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetBrush: %s\n",
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
  if( iBrush < 0 ||
      iBrush >= DspA_knNumBrushes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set the brush theshold info */
  sBrush.mInfo[ iBrush ] = *iInfo;
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d %d",
	      (int)iBrush,
	      (int)sBrush.mInfo[iBrush].mnLow,
	      (int)sBrush.mInfo[iBrush].mnHigh,
	      (int)sBrush.mInfo[iBrush].mnNewValue );
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
  
  if( NULL == iColor ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set cursor color and notify of change */
  sCursorColor = *iColor;
  DspA_Redraw_( this );
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
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
				     Surf_tVertexSet   iSurface,
				     int               inWidth ) {
  
  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  if( iSurface <= Surf_tVertexSet_None || 
      iSurface >= Surf_knNumVertexSets ||
      inWidth < 0 ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set surface width and notify of change */
  this->manSurfaceLineWidth[iSurface] = inWidth;
  DspA_Redraw_( this );
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d", iSurface,  
	      this->manSurfaceLineWidth[iSurface] );
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
				     Surf_tVertexSet   iSurface,
				     xColor3fRef       iColor ) {
  
  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  if( iSurface <= Surf_tVertexSet_None || 
      iSurface >= Surf_knNumVertexSets ||
      NULL == iColor ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set surface color and notify of change */
  this->maSurfaceLineColor[iSurface] = *iColor;
  DspA_Redraw_( this );
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %f %f %f", (int)iSurface,
	      xColr_ExpandFloat( &(this->maSurfaceLineColor[iSurface]) ) );
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

DspA_tErr DspA_SetCursorShape ( tkmDisplayAreaRef this,
				DspA_tMarker      iShape ) {
  
  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  if( iShape <= DspA_tMarker_None ||
      iShape >= DspA_knNumMarkers ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set cursor shape and notify of change */
  sCursorShape = iShape;
  DspA_Redraw_( this );
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
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

DspA_tErr DspA_SetParcBrushInfo ( tkmDisplayAreaRef        this,
				  DspA_tParcBrushSettings* iSettings ) {
  
  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[STRLEN] = "";
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  if( iSettings->mSrc < tkm_tVolumeType_Main ||
      iSettings->mSrc >= tkm_knNumVolumeTypes ||
      iSettings->mnFuzzy < 0 ||
      iSettings->mnFuzzy >= 256 ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set brush data */
  sParcBrush = *iSettings;
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf( sTclArguments, "%d %d %d %d %d", 
	     sParcBrush.mNewValue, sParcBrush.mb3D,
	     sParcBrush.mSrc, sParcBrush.mnFuzzy, sParcBrush.mnDistance );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateParcBrushInfo, sTclArguments );
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SetParcBrushInfo: %s\n",
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
  while( nSlice < 0 ) 
    nSlice += this->mnVolumeSizeZ;
  while( nSlice >= this->mnVolumeSizeZ )
    nSlice -= this->mnVolumeSizeZ;
  
  /* set the slice */
  eResult = DspA_SetSlice( this, nSlice );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* notify window of change if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
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

DspA_tErr DspA_Focus ( tkmDisplayAreaRef this ) {
  
  DspA_tErr eResult = DspA_tErr_NoErr;
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* if we're already focused, return. */
  if( this == sFocusedDisplay )
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

DspA_tErr DspA_SetCursorToCenterOfSpace ( tkmDisplayAreaRef this,
					  x3DListRef        ipList ) {
  
  DspA_tErr  eResult                 = DspA_tErr_NoErr;
  xList_tErr eList                   = xList_tErr_NoErr;
  x3Lst_tErr e3DList                 = x3Lst_tErr_NoErr;
  int        nSlice                  = 0;
  xListRef   list                    = NULL;
  int        nVoxels                 = 0;
  tBoolean   bInRun                  = FALSE;
  int        nNumSlicesInRun         = 0;
  int        nFirstSliceInRun        = 0;
  int        nNumSlicesInLongestRun  = 0;
  int        nFirstSliceInLongestRun = 0;
  int        nLastSliceInLongestRun  = 0;
  
  DebugEnterFunction( ("DspA_SetCursorToCenterOfSpace( this=%p, ipList=%p )",
		       this, ipList) );
  
  DebugNote( ("Verifying self") );
  eResult = DspA_Verify( this );
  DebugAssertThrow( (eResult == DspA_tErr_NoErr) );
  
  /* for each slice */
  for( nSlice = 0; nSlice < this->mnVolumeSizeZ; nSlice++ ) {
    
    /* depending on what orientation we are in, get a list of voxels sorted
       by our slice number */
    switch ( this->mOrientation ) {
    case mri_tOrientation_Coronal:
      DebugNote( ("Getting list in z plane for slice %d",nSlice) );
      e3DList = x3Lst_GetItemsInZPlane( this->mpSelection, nSlice, &list );
      break;
    case mri_tOrientation_Sagittal:
      DebugNote( ("Getting list in x plane for slice %d",nSlice) );
      e3DList = x3Lst_GetItemsInXPlane( this->mpSelection, nSlice, &list );
      break;
    case mri_tOrientation_Horizontal:
      DebugNote( ("Getting list in y plane for slice %d",nSlice) );
      e3DList = x3Lst_GetItemsInYPlane( this->mpSelection, nSlice, &list );
      break;
    default:
      eResult = DspA_tErr_InvalidOrientation;
      goto error;
      break;
    }
    
    /* see if we have any voxels here. */
    DebugNote( ("Getting number of voxels in list with xList_GetCount") );
    eList = xList_GetCount( list, &nVoxels );
    DebugAssertThrowX( (eList==xList_tErr_NoErr), 
		       eResult, DspA_tErr_ErrorAccessingList );
    
    /* if there's a voxel here... */
    if( nVoxels > 0 ) {
      
      /* if we're not in a run, start a run. set the ctr to one and save
	 the beginning slice of the run. */
      if( !bInRun ) {
	bInRun           = TRUE;
	nNumSlicesInRun  = 1;
	nFirstSliceInRun = nSlice;
	
	/* if we are in a run, inc the counter */
      } else {
	nNumSlicesInRun++;
      }
      
      /* if there's no voxel here, and we were in a run... */
    } else if( bInRun ) {
      
      /* end the run. if the run size is greater than our max, save 
	 the midpoint of the start and end of the run. */
      bInRun = FALSE;
      if( nNumSlicesInRun > nNumSlicesInLongestRun ) {
	nFirstSliceInLongestRun = nFirstSliceInRun;
	nLastSliceInLongestRun  = nSlice - 1;
	nNumSlicesInLongestRun  = nNumSlicesInRun;
      }
    }
  }
  
  /* if we found a run.. */
  if( nNumSlicesInLongestRun > 0 ) {
    
    /* set the slice to the middle of the run */
    DebugNote( ("Setting slice to middle of run") );
    DspA_SetSlice( this, 
		   ((nLastSliceInLongestRun - nFirstSliceInLongestRun) / 2) +
		   nFirstSliceInLongestRun);
  }
  
  DebugCatch;
  DebugCatchError( eResult, DspA_tErr_NoErr, DspA_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

DspA_tErr DspA_HandleEvent ( tkmDisplayAreaRef this, 
			     xGWin_tEventRef   ipEvent ){
  
  DspA_tErr   eResult     = DspA_tErr_NoErr;
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  if( xGWin_tEventType_MouseUp == ipEvent->mType
      || xGWin_tEventType_MouseDown == ipEvent->mType
      || xGWin_tEventType_MouseMoved == ipEvent->mType ) {
    
  }
  
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
    
    eResult = DspA_HandleMouseMoved_( this, ipEvent );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
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
  
  DspA_tErr   eResult     = DspA_tErr_NoErr;
  xPoint2n    bufferPt    = {0,0};
  xVoxelRef   pVolumeVox  = NULL;
  int         nParcIndex  = 0;
  DspA_tParcBrushSettings parcBrush;
  
  xVoxl_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
#if 0
  DebugPrint( ("Mouse up screen x %d y %d buffer x %d y %d volume %f %f %f\n",
	       ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY, 
	       xVoxl_ExpandFloat( pVolumeVox ) ) );
#endif
  
  /* if nav tool and not ctrl... */
  if( DspA_tTool_Navigate == sTool &&
      !ipEvent->mbCtrlKey ) {
    
    /* if there was little delta... */
    if( this->mTotalDelta.mfX > -1.0 && this->mTotalDelta.mfX < 1.0 &&
	this->mTotalDelta.mfY > -1.0 && this->mTotalDelta.mfY < 1.0 ) {
      
      /* do a single slice up or down or zoom in or out depending on what
	 vertical side of the screen we clicked on */
      switch( ipEvent->mButton ) {
      case 2:
	if( GLDRAW_Y_FLIP(bufferPt.mnY) > 128 ) {
	  DspA_SetSlice( this,this->mnOriginalSlice - 1 );
	} else {
	  DspA_SetSlice( this,this->mnOriginalSlice + 1 );
	}
	break;
      case 3:
	if( GLDRAW_Y_FLIP(bufferPt.mnY) > 128 ) {
	  DspA_SetZoomLevel( this,this->mnOriginalZoomLevel - 1 );
	} else {
	  DspA_SetZoomLevel( this,this->mnOriginalZoomLevel + 1 );
	}
	break;
      }
    }
  }
  
#if 0
  /* allow the functional display to respond. */
  if( NULL != this->mpFunctionalVolume ) {
    eFunctional = FunV_AnatomicalVoxelClicked( this->mpFunctionalVolume,
					       pVolumeVox );
    if( FunV_tErr_NoError != eFunctional ) {
      DebugPrint( ("DspA_HandleMouseUp_(): Error while passing clicked voxel to functional volume.\n" ) );
    }
  }
#endif   
  
  /* if this isn't the nav tool or is but ctrl is down... */
  if( !(DspA_tTool_Navigate == sTool &&
	!ipEvent->mbCtrlKey) ) {
    
    /* set the cursor. */
    eResult = DspA_SetCursor( this, pVolumeVox );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  /* if ctrl was down... */
  if ( ipEvent->mbCtrlKey ) {
    
    /* set zoom center to cursor */
    eResult = DspA_SetZoomCenterToCursor( this );
    if( DspA_tErr_NoErr != eResult )
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
  
  /* if edit parc tool... */
  if( DspA_tTool_EditSegmentation == sTool
      && NULL != this->mROIGroup
      && !ipEvent->mbCtrlKey) {
    
    switch( ipEvent->mButton ) {
      
      /* button two edits or color-sucks */
    case 2:
      
      /* shift-2 sucks the color */
      if( ipEvent->mbShiftKey ) {
	
	/* get the color and set our brush info with the same settings
	   except for the new color */
	tkm_GetROILabel( pVolumeVox, &nParcIndex, NULL );
	parcBrush = sParcBrush;
	parcBrush.mNewValue = nParcIndex;
	DspA_SetParcBrushInfo( this, &parcBrush );
	
      } else {
	
	eResult = DspA_BrushVoxels_( this, pVolumeVox,  
				     NULL, DspA_EditSegmentationVoxels_ );
	if( DspA_tErr_NoErr != eResult )
	  goto error;
	this->mbSliceChanged = TRUE;
      }
      
      break;
      
      /* button three does a flood fill */
    case 3:
      tkm_FloodFillSegmentation( pVolumeVox, sParcBrush.mNewValue, 
				 sParcBrush.mb3D, sParcBrush.mSrc, 
				 sParcBrush.mnFuzzy, sParcBrush.mnDistance );
      this->mbSliceChanged = TRUE;
      break;
    }
    
    /* send the cursor info again to update the new roi volume after 
       editing. */
    eResult = DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
					       pVolumeVox );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
  }
  
  /* if edit tool and shift */
  if( DspA_tTool_EditVoxels == sTool 
      && ipEvent->mbShiftKey ) {
    
    /* undo this voxel */
    tkm_RestoreUndoVolumeAroundAnaIdx( pVolumeVox );
  }
  
  /* if we're in control point selection mode... */
  if( DspA_tTool_EditCtrlPts == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* if button 2, make a new point here. */
    if ( 2 == ipEvent->mButton ) {
      
      tkm_MakeControlPoint( pVolumeVox );
      
      /* if button 3, delete the nearest point. */
    } else if ( 3 == ipEvent->mButton ) {
      
      tkm_RemoveControlPointWithinDist( pVolumeVox, this->mOrientation, 3 );
    }
  }
  
  /* most things want a mouse up after we're done, so schedule one and 
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
  
#if 0
  DebugPrint( ("Mouse down screen x %d y %d buffer x %d y %d volume %d %d %d\n",
	       ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
	       xVoxl_ExpandInt( pVolumeVox ) ) );
#endif
  
  /* if nav tool... */
  if( DspA_tTool_Navigate == sTool ) {
    
    this->mLastClick = ipEvent->mWhere;
    this->mTotalDelta.mfX = 0;
    this->mTotalDelta.mfY = 0;
    xVoxl_Copy( this->mpOriginalZoomCenter, this->mpZoomCenter );
    this->mnOriginalSlice = DspA_GetCurrentSliceNumber_( this );
    this->mnOriginalZoomLevel = this->mnZoomLevel;
  }
  
  /* button 2 or 3 edit tool with no modifiers: */
  if( ( 2 == ipEvent->mButton
	|| 3 == ipEvent->mButton )
      && DspA_tTool_EditVoxels == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* clear the undo list. */
    tkm_ClearUndoList();
    
    /* button determines the brush action */
    if( 2 == ipEvent->mButton ) {
      brushAction = DspA_tBrush_EditOne;
    } else if ( 3 == ipEvent->mButton ) {
      brushAction = DspA_tBrush_EditTwo;
    }
    
    /* brush the voxels */
    eResult = DspA_BrushVoxels_( this, pVolumeVox, (void*)&brushAction,
				 DspA_BrushVoxelsInThreshold_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* editing requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
    
  }
  
  /* select mode with no modfiers */
  if( ( 2 == ipEvent->mButton
	|| 3 == ipEvent->mButton )
      && DspA_tTool_SelectVoxels == sTool
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* button determines the select action */
    if( 2 == ipEvent->mButton ) {
      selectAction = DspA_tSelectAction_Select;
    } else if ( 3 == ipEvent->mButton ) {
      selectAction = DspA_tSelectAction_Deselect;
    }
    
    /* brush the voxels */
    eResult = DspA_BrushVoxels_( this, pVolumeVox, (void*)&selectAction,
				 DspA_SelectVoxels_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* selecting requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
    
  }
  
  /* if edit parc tool with button 1 or 3, clear the undo list */
  if( ( 1 == ipEvent->mButton
	|| 3 == ipEvent->mButton )
      && DspA_tTool_EditSegmentation == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    tkm_ClearUndoList();
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
  xVoxelRef          pVolumeVox   = NULL;
  DspA_tBrush        brushAction  = DspA_tBrush_None;
  DspA_tSelectAction selectAction = DspA_tSelectAction_None;
  xPoint2f           delta;
  xPoint2f           newCenterPt;
  xVoxel             newCenterIdx;
  int                nNewSlice    = 0;
  int                nNewZoomLevel= 0;
  
  xVoxl_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
#if 0
  DebugPrint( ("Mouse moved screen x %d y %d buffer x %d y %d volume %d %d %d\n",
	       ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
	       xVoxl_ExpandInt( pVolumeVox ) ) );
#endif
  
  eResult = DspA_VerifyVolumeVoxel_( this, pVolumeVox );
  if ( DspA_tErr_NoErr == eResult )
    DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Mouseover, 
				     pVolumeVox );
  
  /* save this mouse location */
  this->mMouseLocation = ipEvent->mWhere;
  
  /* if a button isn't down, skip this */
  if( ipEvent->mButton == 0 )
    goto cleanup;
  
  /* if nav tool and no ctrl key... */
  if( DspA_tTool_Navigate == sTool &&
      !ipEvent->mbCtrlKey ) {
    
    /* figure out how much we've moved */
    delta.mfX = (float)(this->mLastClick.mnX - ipEvent->mWhere.mnX) / (float)this->mnZoomLevel / 2.0;
    delta.mfY = (float)(this->mLastClick.mnY - ipEvent->mWhere.mnY) / (float)this->mnZoomLevel / 2.0;
    
    /* flip y if horizontal cuz our freaking screen is upside down */
    if( this->mOrientation == mri_tOrientation_Horizontal )
      delta.mfY = -delta.mfY;
    
    /* add to the total delta */
    this->mTotalDelta.mfX += delta.mfX;
    this->mTotalDelta.mfY += delta.mfY;
    
    /* save this mouse position */
    this->mLastClick = ipEvent->mWhere;
    
    /* what button? */
    switch( ipEvent->mButton ) {
      
      /* button one is 'drag' the volume view. add the rounded delta
	 to the original center and set the center */
    case 1:
      newCenterPt.mfX = this->mTotalDelta.mfX +
	xVoxl_GetX( this->mpOriginalZoomCenter );
      newCenterPt.mfY = this->mTotalDelta.mfY +
	xVoxl_GetY( this->mpOriginalZoomCenter );
      DspA_ConvertPlaneToVolume_( this, &newCenterPt, 0,
				  this->mOrientation, &newCenterIdx );
      DspA_SetZoomCenter( this, &newCenterIdx );
      break;
      
      /* button 2 is move the slice. add the rounded y delta to the 
	 slice */
    case 2:
      nNewSlice = this->mnOriginalSlice + rint(this->mTotalDelta.mfY);
      if( nNewSlice >= 0 && nNewSlice < this->mnVolumeSizeZ ) {
	DspA_SetSlice( this, nNewSlice );
      }
      break;
      
      /* button 3 is zoom. add the rounded delta to the zoom level */
    case 3:
      nNewZoomLevel = this->mnOriginalZoomLevel + rint(this->mTotalDelta.mfY);
      if ( nNewZoomLevel >= DspA_knMinZoomLevel 
	   && nNewZoomLevel <= DspA_knMaxZoomLevel ) {
	DspA_SetZoomLevel( this, nNewZoomLevel );
      }
      break;
    }
    
  }
  
  /* button 2 or 3 edit tool with no modifiers: */
  if( ( 2 == ipEvent->mButton
	|| 3 == ipEvent->mButton )
      && DspA_tTool_EditVoxels == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* button determines the brush action */
    if( 2 == ipEvent->mButton ) {
      brushAction = DspA_tBrush_EditOne;
    } else if ( 3 == ipEvent->mButton ) {
      brushAction = DspA_tBrush_EditTwo;
    }
    
    /* brush the voxels */
    eResult = DspA_BrushVoxels_( this, pVolumeVox, (void*)&brushAction,
				 DspA_BrushVoxelsInThreshold_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* editing requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
    
  }
  
  /* if edit parc tool button 1... */
  if( DspA_tTool_EditSegmentation == sTool
      && NULL != this->mROIGroup
      && ipEvent->mButton == 2
      && !(ipEvent->mbAltKey)) {
    
    /* edit the parc volume */
    eResult = DspA_BrushVoxels_( this, pVolumeVox, 
				 NULL, DspA_EditSegmentationVoxels_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* editing requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }
  
  /* select mode with no modfiers */
  if( ( 2 == ipEvent->mButton
	|| 3 == ipEvent->mButton )
      && DspA_tTool_SelectVoxels == sTool
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* button determines the select action */
    if( 2 == ipEvent->mButton ) {
      selectAction = DspA_tSelectAction_Select;
    } else if ( 3 == ipEvent->mButton ) {
      selectAction = DspA_tSelectAction_Deselect;
    }
    
    /* brush the voxels */
    eResult = DspA_BrushVoxels_( this, pVolumeVox,
				 (void*)&selectAction, DspA_SelectVoxels_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* selecting requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
    
  }
  
  eResult = DspA_tErr_NoErr;
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_HandleMouseMoved_: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  xVoxl_Delete( &pVolumeVox );
  
  return eResult;
}

DspA_tErr DspA_HandleKeyDown_ ( tkmDisplayAreaRef this, 
				xGWin_tEventRef   ipEvent ) {
  
  DspA_tErr eResult       = DspA_tErr_NoErr;
  MWin_tErr eWindowResult = MWin_tErr_NoErr;
  FunV_tErr eFunctional   = FunV_tErr_NoError;
  
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
    /* ctrl 1 goes to main volume */
    if ( ipEvent->mKey == '1' &&
	 ipEvent->mbCtrlKey ) {
      eResult = DspA_SetDisplayFlag( this, 
				     DspA_tDisplayFlag_AuxVolume, FALSE );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
      break;
    }
    
    /* ctrl 2 goes to aux volume */
    if ( ipEvent->mKey == '1' &&
	 ipEvent->mbCtrlKey ) {
      eResult = DspA_SetDisplayFlag( this, 
				     DspA_tDisplayFlag_AuxVolume, TRUE );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    
    /* any other number sets the radius of the brush. */
    eResult = DspA_SetBrushShape( this, (int)(ipEvent->mKey - '0'),
				  sBrush.mShape, sBrush.mb3D );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 'c':
    /* alt-c toggles main and aux view. */
    if( ipEvent->mbAltKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_AuxVolume );
      /* ctrl+c toggles cursor display */
    } else if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_Cursor );
    } else {
      /* c sets tool to control pts */
      eResult = DspA_SetTool( this, DspA_tTool_EditCtrlPts );
    }
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    break;
    
  case 'e':
    /* e sets tool to edit */
    eResult = DspA_SetTool( this, DspA_tTool_EditVoxels );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 'f':
    /* alt-f swaps the anatomical/overlay display */
    if( ipEvent->mbAltKey ) {
      if( NULL == this->mpFunctionalVolume )
	goto cleanup;
      
      if( 1 == this->mabDisplayFlags[DspA_tDisplayFlag_Anatomical] ) {
	DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Anatomical, 0 );
	DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay, 1 );
      } else {
	DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Anatomical, 1 );
	DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay, 0 );
      }
    }
    
    /* ctrl+f toggles functional overlay display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
					DspA_tDisplayFlag_FunctionalOverlay );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case 'g':
    
    /* g sets tool to edit segmentation */
    eResult = DspA_SetTool( this, DspA_tTool_EditSegmentation );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    break;
    
  case 'h':
    
    FunV_UseOverlayCache( this->mpFunctionalVolume, 
			  !this->mpFunctionalVolume->mbUseOverlayCache );
    if( this->mpFunctionalVolume->mbUseOverlayCache ) {
      DebugPrint( ("Cache enabled.\n" ) );
    } else {
      DebugPrint( ("Cache disabled.\n" ) );
    }
    break;
    
  case 'i':
    /* ctrl+i toggles interpolated vertices display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
					DspA_tDisplayFlag_InterpolateSurfaceVertices );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case 'k':
    /* ctrl+k toggles cursor linking */
    if( ipEvent->mbCtrlKey ) {
      eWindowResult = MWin_ToggleLinkedCursorFlag( this->mpWindow );
      if( MWin_tErr_NoErr != eWindowResult ) {
	eResult = DspA_tErr_ErrorAccessingWindow;
	goto error;
      }
    }
    break;
    
  case 'l':
    /* ctrl+l toggles selection display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_Selection );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case 'n':
    /* n sets tool to navigate */
    eResult = DspA_SetTool( this, DspA_tTool_Navigate );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 'm':
    /* m toggles the segmentation display */
    eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_ROIGroupOverlay);
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 'o':
    /* ctrl+o toggles original suface display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
					DspA_tDisplayFlag_OriginalSurface );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case 'p':
    
    /* ctrl+p toggles canonical/pial suface display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
					DspA_tDisplayFlag_CanonicalSurface );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    } else {
      /* p sets tool to edit */
      eResult = DspA_SetTool( this, DspA_tTool_EditSegmentation );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    
    break;
    
  case 's':
    /* ctrl+s toggles main suface display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_MainSurface );
    } else {
      /* s sets tool to select */
      eResult = DspA_SetTool( this, DspA_tTool_SelectVoxels );
    }
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 't':
    /* ctrl+s toggles ctrl pt display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
					DspA_tDisplayFlag_ControlPoints );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case 'v':
    /* ctrl-v toggles vertex display. */
    if ( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
					DspA_tDisplayFlag_DisplaySurfaceVertices );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case 'x': 
    /* x sets plane to sagittal */
    eResult = DspA_SetOrientation( this, mri_tOrientation_Sagittal );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 'y': 
    /* y sets plane to horizontal */
    eResult = DspA_SetOrientation( this, mri_tOrientation_Horizontal );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case 'z': 
    /* ctrl-z undos */
    if ( ipEvent->mbCtrlKey ) {
      tkm_RestoreUndoList();
    } else {
      /* z sets plane to coronal */
      eResult = DspA_SetOrientation( this, mri_tOrientation_Coronal );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
    }
    break;
    
  case xGWin_tKey_UpArrow:
  case xGWin_tKey_RightArrow:
    /* move up a slice */
    eResult = DspA_ChangeSliceBy( this, 1 );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    DspA_Redraw_( this );
    break;
    
  case xGWin_tKey_DownArrow:
  case xGWin_tKey_LeftArrow:
    /* move down a slice */
    eResult = DspA_ChangeSliceBy( this, -1 );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    DspA_Redraw_( this );
    break;
    
  case '+':
    if( NULL != this->mpFunctionalVolume ) {
      /* move the time point up */
      eFunctional = FunV_ChangeTimePointBy( this->mpFunctionalVolume, 1 );
      if( FunV_tErr_NoError != eFunctional ) {
	eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
	goto error;
      }
    }
    break;
    
  case '-':
    if( NULL != this->mpFunctionalVolume ) {
      /* move the time point down */
      eFunctional = FunV_ChangeTimePointBy( this->mpFunctionalVolume, -1 );
      if( FunV_tErr_NoError != eFunctional ) {
	eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
	goto error;
      }
    }
    break;
    
  }
  
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
  if( iBrush < 0 ||
      iBrush >= DspA_knNumBrushes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* set the brush theshold info */
  switch( iBrush ) {
  case DspA_tBrush_EditOne:
    sBrush.mInfo[iBrush].mnLow = tkm_knEditToWhiteLow;
    sBrush.mInfo[iBrush].mnHigh = tkm_knEditToWhiteHigh;
    sBrush.mInfo[iBrush].mnNewValue = tkm_knEditToWhiteNewValue;
    break;
  case DspA_tBrush_EditTwo:
    sBrush.mInfo[iBrush].mnLow = tkm_knEditToBlackLow;
    sBrush.mInfo[iBrush].mnHigh = tkm_knEditToBlackHigh;
    sBrush.mInfo[iBrush].mnNewValue = tkm_knEditToBlackNewValue;
    break;
  default:
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  
  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d %d",
	      (int)iBrush,
	      (int)sBrush.mInfo[iBrush].mnLow,
	      (int)sBrush.mInfo[iBrush].mnHigh,
	      (int)sBrush.mInfo[iBrush].mnNewValue );
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
			      void(*ipFunction)(xVoxelRef,void*) ) {
  
  DspA_tErr        eResult     = DspA_tErr_NoErr;
  int              nXCenter    = 0;
  int              nYCenter    = 0;
  int              nZCenter    = 0;
  int              nXRadius    = 0;
  int              nYRadius    = 0;
  int              nZRadius    = 0;
  int              nX          = 0;
  int              nY          = 0;
  int              nZ          = 0;
  xVoxelRef        pVolumeVox  = NULL;
  
  xVoxl_New( &pVolumeVox );
  
  /* get our center voxel. */
  nXCenter = xVoxl_GetX( ipCenterVox );
  nYCenter = xVoxl_GetY( ipCenterVox );
  nZCenter = xVoxl_GetZ( ipCenterVox );
  
  /* set all radii to the brush radius. we subtract one because of our 
     looping bounds. */
  nXRadius = sBrush.mnRadius - 1;
  nYRadius = sBrush.mnRadius - 1;
  nZRadius = sBrush.mnRadius - 1;
  
  /* if we're not in 3d, set the same radius of the same plane as our
     current orientation to 0. */
  if( !sBrush.mb3D ) {
    switch( this->mOrientation ) {
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
  
  /* loop around the starting point... */
  for ( nZ = nZCenter - nZRadius; nZ <= nZCenter + nZRadius; nZ++ ) {
    for ( nY = nYCenter - nYRadius; nY <= nYCenter + nYRadius; nY++ ) {
      for ( nX = nXCenter - nXRadius; nX <= nXCenter + nXRadius; nX++ ) {
	
	/* set our voxel. */
	xVoxl_Set( pVolumeVox, nX, nY, nZ );
	
	/* if it's valid... */
	if( DspA_tErr_InvalidVolumeVoxel != 
	    DspA_VerifyVolumeVoxel_( this, pVolumeVox ) ) {
	  
	  /* if we're circular, check the radius. if no good, continue. */
	  if( DspA_tBrushShape_Circle == sBrush.mShape 
	      && ( (nXCenter-nX)*(nXCenter-nX) +
		   (nYCenter-nY)*(nYCenter-nY) +
		   (nZCenter-nZ)*(nZCenter-nZ) >
		   (sBrush.mnRadius-1)*(sBrush.mnRadius-1) ) ) {
	    continue;
	  }
	  
	  /* run the function on this voxel. */
	  ipFunction( pVolumeVox, ipData );
	}
      }
    }
  }
  
  DspA_Redraw_( this );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_BrushVoxels_: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  xVoxl_Delete( &pVolumeVox );
  
  return eResult;
  
}

void DspA_BrushVoxelsInThreshold_ ( xVoxelRef ipVoxel, void* ipData ) {
  
  DspA_tBrush brush = DspA_tBrush_None;
  
  /* make sure the brush is in bounds */
  brush = *(DspA_tBrush*)ipData;
  if( brush < 0 ||
      brush >= DspA_knNumBrushes ) {
    return;
  }
  
  /* edit the voxel */
  tkm_EditVoxelInRange( ipVoxel, 
			sBrush.mInfo[brush].mnLow,
			sBrush.mInfo[brush].mnHigh,
			sBrush.mInfo[brush].mnNewValue );
}

void DspA_SelectVoxels_ ( xVoxelRef ipVoxel, void* ipData ) {
  
  DspA_tSelectAction selectAction = DspA_tSelectAction_None;
  
  /* make sure the action is in bounds */
  selectAction = *(DspA_tSelectAction*)ipData;
  if( selectAction < 0 ||
      selectAction >= DspA_knNumSelectActions ) {
    return;
  }
  
  /* select or deselect the voxel */
  switch( selectAction ) {
  case DspA_tSelectAction_Select:
    
    tkm_SelectVoxel( ipVoxel );
    break;
    
  case DspA_tSelectAction_Deselect:
    
    tkm_DeselectVoxel( ipVoxel );
    break;
    
  default:
    break;
  }
}

void DspA_EditSegmentationVoxels_ ( xVoxelRef ipVoxel, void* ipData ) {
  
  tkm_EditSegmentation( ipVoxel, sParcBrush.mNewValue );
}

DspA_tErr DspA_SelectCurrentROI ( tkmDisplayAreaRef this ) {
  
  DspA_tErr eResult = DspA_tErr_NoErr;
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* if we have the data and our index is good, tell tkmedit to select stuff */
  if( NULL != this->mROIGroup
      && -1 != this->mnROIGroupIndex ) {
    tkm_SelectCurrentROI( this->mnROIGroupIndex );
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_SelectCurrentROI: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

DspA_tErr DspA_GraphCurrentROIAvg ( tkmDisplayAreaRef this ) {
  
  DspA_tErr eResult = DspA_tErr_NoErr;
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* if we have the data and our index is good, tell tkmedit to handle it */
  if( NULL != this->mROIGroup
      && -1 != this->mnROIGroupIndex ) {
    tkm_GraphCurrentROIAvg( this->mnROIGroupIndex );
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_GraphCurrentROIAvg: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
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
  if( this->mbSliceChanged ) {
    eResult = DspA_BuildCurrentFrame_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
    
    /* draw the selection */
    if( this->mabDisplayFlags[DspA_tDisplayFlag_Selection] ) {
      eResult = DspA_DrawSelectionToFrame_( this );
      if ( DspA_tErr_NoErr != eResult ) {
	DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
	eResult = DspA_tErr_NoErr;
      }
    }
  }
  
  /* draw functional overlay. */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ) {
    eResult = DspA_DrawFunctionalOverlayToFrame_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
  
  /* draw the frame buffer */
  eResult = DspA_DrawFrameBuffer_ ( this );
  if ( DspA_tErr_NoErr != eResult ) {
    DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
    eResult = DspA_tErr_NoErr;
  }
  
  /* draw head points */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_HeadPoints] ) {
    eResult = DspA_DrawHeadPoints_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
  
  /* draw the control points */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_ControlPoints] ) {
    eResult = DspA_DrawControlPoints_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
  
  /* draw the control points */
  //  if( this->mabDisplayFlags[DspA_tDisplayFlag_VectorField] ) {
  eResult = DspA_DrawVectorField_( this );
  if ( DspA_tErr_NoErr != eResult ) {
    DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
    eResult = DspA_tErr_NoErr;
  }
  //  }
  
  /* draw the surface */
  eResult = DspA_DrawSurface_ ( this );
  if ( DspA_tErr_NoErr != eResult ) {
    DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
    eResult = DspA_tErr_NoErr;
  }
  
  /* draw the cursor */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_Cursor] ) {
    eResult = DspA_DrawCursor_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
  
  /* if we're focused and our focus frame flag is on... */
  if( this == sFocusedDisplay
      && TRUE == this->mabDisplayFlags[DspA_tDisplayFlag_FocusFrame] ) {
    /* draw a frame around us. */
    eResult = DspA_DrawFrameAroundDisplay_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
  
  if( this->mabDisplayFlags[DspA_tDisplayFlag_Axes] ) {
    eResult = DspA_DrawAxes_ ( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }
  
  /* rebuilt all our slice changes, so we can clear this flag. */
  this->mbSliceChanged = FALSE;
  
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
  
  glDrawPixels ( this->mnVolumeSizeX, this->mnVolumeSizeY,  /* almost certainly wrong (BRF) */
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
  Surf_tVertexSet   surface    = Surf_tVertexSet_Main;
  xGrowableArrayRef list       = NULL;
  float             faColor[3] = {0, 0, 0};
  
#if 0
  /* just draw surface and return. */
  eResult = DspA_DrawSurfaceDirect_( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  goto cleanup;
#endif
  
  for( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {
    
    /* only recalculate the draw list of the plane and orientation has
       changed. */
    if ( NULL != this->mpSurface[nSurface] ) {
      
      /* if the orienation or slice has changed... */
      if( this->mbSliceChanged ) {
	
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
      for ( surface = Surf_tVertexSet_Main;
	    surface < Surf_knNumVertexSets; surface++ ) {
	
	/* if this surface is visible... */
	if( this->mabDisplayFlags[surface + DspA_tDisplayFlag_MainSurface] ) {
	  
	  /* get the list. */
	  list = DspA_GetSurfaceList_( this, nSurface, 
				       this->mOrientation, surface,
				       DspA_GetCurrentSliceNumber_(this) );
	  if( NULL == list ) {
	    eResult = DspA_tErr_ErrorAccessingSurfaceList;
	    goto error;
	  }
	  
	  /* set the color */
	  xColr_PackFloatArray( &(this->maSurfaceLineColor[surface]), 
				faColor );
	  glColor3fv( faColor );
	  
	  /* draw the points. */
	  glLineWidth(  this->manSurfaceLineWidth[surface] );
	  DspA_ParsePointList_( this, GL_LINES, list );
	  
	  /* if vertices are visible... */
	  if( this->mabDisplayFlags[DspA_tDisplayFlag_DisplaySurfaceVertices] ) {
	    
	    /* invert color */
	    faColor[0] = 1.0 - faColor[0];
	    faColor[1] = 1.0 - faColor[1];
	    faColor[2] = 1.0 - faColor[2];
	    glColor3fv( faColor );
	    
	    /* set point size. */
	    glPointSize( DspA_knSurfaceVertexSize );
	    
	    /* draw the vertices. */
	    DspA_ParsePointList_( this, GL_POINTS, list );
	  }
	  
	  /* if we have a hilited vertex for this surface... */
	  if( surface == this->mHilitedSurface
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
  
#if 0
  
  int       nY          = 0;
  char      sLabel[STRLEN] = "";
  int       nLength     = 0;
  int       nChar       = 0;
  
  xVoxl_New( &pVoxel );
  
  DspA_SetUpOpenGLPort_( this );
  
  glColor3f ( 0.0, 1.0, 0.0 );
  
  volumePt.mnX = 0;
  
  /* all down the side, every size / 5 pixels */
  for ( nY = 0; nY < this->mnVolumeSizeY; nY += this->mnVolumeSizeY/5 ) {
    
    /* y flip the volume pt to flip the image over. */
    //volumePt.mnY = BUFFER_Y_FLIP(nY);
    volumePt.mnY = nY;
    
    /* get a volume voxel.*/
    eResult = DspA_ConvertBufferToVolume_ ( this, &volumePt, pVoxel );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* if this is good... */
    eResult = DspA_VerifyVolumeVoxel_( this, pVoxel );
    if( DspA_tErr_NoErr == eResult ) {
      
      /* convert the voxel coords to a string */
      sprintf( sLabel, "%d,%d,%d", xVoxl_ExpandInt(pVoxel) );
      
      /* draw a line here. */
      glBegin( GL_LINES );
      glVertex2d( 0, nY );
      glVertex2d( 10, nY );
      glEnd();
      
      /* get length of the string */
      nLength = strlen( sLabel );
      
      /* go to right place */
      glRasterPos2i( 0, nY + 5 );
      
      /* for every cahr in the string.. */
      for( nChar = 0; nChar < nLength; nChar++ ) {
	
	/* draw the char */
	glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
      }
    }
    
    eResult = DspA_tErr_NoErr;
  }
  
#endif
  
  DspA_SetUpOpenGLPort_( this );
  
  glColor3f ( 0.0, 1.0, 0.0 );
  
  /* draw arrows */
  switch( this->mOrientation ) {
    
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

void DspA_DrawVerticalArrow_ ( xPoint2nRef iStart,
			       int         inLength,
			       char*       isLabel ) {
  
  int nChar = 0;
  int nStrLength = 0;
  int nHeadOffset = 0;
  
  if( inLength > 0 ) {
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
  for( nChar = 0; nChar < nStrLength; nChar++ ) 
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, isLabel[nChar] );
  
}

void DspA_DrawHorizontalArrow_ ( xPoint2nRef iStart,
				 int         inLength,
				 char*       isLabel ) {
  
  int nChar = 0;
  int nStrLength = 0;
  int nHeadOffset = 0;
  
  if( inLength > 0 ) {
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
  for( nChar = 0; nChar < nStrLength; nChar++ ) 
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, isLabel[nChar] );
  
}



DspA_tErr DspA_BuildCurrentFrame_ ( tkmDisplayAreaRef this ) {
  
  DspA_tErr             eResult     = DspA_tErr_NoErr;
  xPoint2n              volumePt    = {0, 0};
  GLubyte*              pDest       = NULL;
  //  unsigned char         ucValue     = 0;
  xVoxelRef             pVoxel      = NULL;
  xVoxelRef             pFuncMin    = NULL;
  xVoxelRef             pFuncMax    = NULL;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  FunV_tFunctionalValue funcValue   = 0.0;
  xColor3f              color       = {0,0,0};
  xColor3f              funcColor   = {0,0,0};
  xColor3f              roiColor    = {0,0,0};
  xColor3f              dtiColor    = {0,0,0};
  int                   nY          = 0;
  
  //  xUtil_StartTimer();
  
  /* make our voxels */
  xVoxl_New ( &pVoxel );
  xVoxl_New ( &pFuncMin );
  xVoxl_New ( &pFuncMax );
  
  /* if we have and are showing functional data... */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] 
      && NULL != this->mpFunctionalVolume ) {
    
    /* get our overlay func bounds in antomical space */
    FunD_GetBoundsInAnatomical( this->mpFunctionalVolume->mpOverlayVolume, 
				pFuncMin, pFuncMax );
  }
  
  /* if we're in zoom level one, set zoom center to 128,128,128 */
  if( 1 == this->mnZoomLevel ) {
    
    /* set a voxel to 128,128,128, setting the zoom center will not change
       the current slice. */
    xVoxl_Set( pVoxel, this->mnVolumeSizeX/2, 
	       this->mnVolumeSizeY/2, this->mnVolumeSizeZ/2 );
    eResult = DspA_SetZoomCenter( this, pVoxel );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }
  
  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;
  
  DisableDebuggingOutput;
  
  /* go thru the buffer... */
  for ( nY = 0; nY < this->mnVolumeSizeY; nY ++ ) {
    for ( volumePt.mnX = 0; 
	  volumePt.mnX < this->mnVolumeSizeX; volumePt.mnX ++ ) {
      
      /* y flip the volume pt to flip the image over. */
      volumePt.mnY = BUFFER_Y_FLIP(nY);
      
      /* get a volume voxel.*/
      eResult = DspA_ConvertBufferToVolume_ ( this, &volumePt, pVoxel );
      if ( DspA_tErr_NoErr != eResult )
	goto error;
      
      /* check it. */
      eResult = DspA_VerifyVolumeVoxel_( this, pVoxel );
      if( DspA_tErr_NoErr == eResult ) {
	
	/* if we are drawing anatomical data... */
	if( this->mabDisplayFlags[DspA_tDisplayFlag_Anatomical] ) {
	  
	  /* get the normal or max color for the main or aux volume, 
	     depending on our display flags. */
	  if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
	    if( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] ) {
	      Volm_GetMaxColorAtIdx( this->mpAuxVolume, pVoxel,
				     this->mOrientation, &color );
	    } else {
	      Volm_GetColorAtIdx( this->mpAuxVolume, pVoxel, &color );
	    }
	  } else {
	    if( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] ) {
	      Volm_GetMaxColorAtIdx( this->mpVolume, pVoxel,
				     this->mOrientation, &color );
	    } else {
	      Volm_GetColorAtIdx( this->mpVolume, pVoxel, &color );
	    }
	  }
	  
	} else {
	  
	  /* color is just black */
	  color.mfRed   = 0;
	  color.mfGreen = 0;
	  color.mfBlue  = 0;
	}
	
	/* If we are showing the DTI volume, get a blended color at this
	   voxel and use this color. */
	if( this->mabDisplayFlags[DspA_tDisplayFlag_DTIOverlay] ) {
	  tkm_GetDTIColorAtVoxel( pVoxel, this->mOrientation, 
				  &color, &dtiColor );
	  color = dtiColor;
	}
	
	/* if we are showing roi... */
	if( this->mabDisplayFlags[DspA_tDisplayFlag_ROIGroupOverlay] ) {
	  
	  /* get roi color blended in. */
	  tkm_GetROIColorAtVoxel( pVoxel, &color, &roiColor );
	  color = roiColor;
	}
	
	/* if we are showing undoable voxels... */
	if( this->mabDisplayFlags[DspA_tDisplayFlag_UndoVolume] ) {
	  
	  /* is this is one of them, draw with a blue overlay */
	  if( tkm_IsAnaIdxInUndoVolume ( pVoxel ) ) {
	    xColr_HilightComponent( &color, xColr_tComponent_Blue );
	  }
	}
	
	/* if we have and are showing functional data, or our
	   mask arg is on... */
	if( ( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ||
	      this->mabDisplayFlags[DspA_tDisplayFlag_MaskToFunctionalOverlay] ) &&
	    NULL != this->mpFunctionalVolume ) {
	  
	  /* if this voxel is in our func bounds... */
	  if( xVoxl_GetX(pVoxel) >= xVoxl_GetX(pFuncMin)
	      && xVoxl_GetX(pVoxel) <= xVoxl_GetX(pFuncMax)
	      && xVoxl_GetY(pVoxel) >= xVoxl_GetY(pFuncMin)
	      && xVoxl_GetY(pVoxel) <= xVoxl_GetY(pFuncMax)
	      && xVoxl_GetZ(pVoxel) >= xVoxl_GetZ(pFuncMin)
	      && xVoxl_GetZ(pVoxel) <= xVoxl_GetZ(pFuncMax) ) {
	    
	    /* get a functional value. */
	    eFunctional = 
	      FunV_GetValueAtAnaIdx( this->mpFunctionalVolume,
				     pVoxel, &funcValue );
	    
	    /* if it was a valid voxel */
	    if( FunV_tErr_NoError == eFunctional ) {
	      
	      /* get a color value. use the red compoonent for
		 base value. */
	      eFunctional = FunV_GetColorForValue ( this->mpFunctionalVolume,
						    funcValue, 
						    &color, &funcColor );
	      
	      /* if the color is not all black.. */
	      if( funcColor.mfRed != 0 || 
		  funcColor.mfGreen != 0 || 
		  funcColor.mfBlue  != 0 ) {
		
		/* set the color to this func color */
		color = funcColor;
	      }
	    } 
	    
	  } else {
	    eFunctional = FunV_tErr_InvalidAnatomicalVoxel;
	  }
	  
	  /* if we are out of functional bounds and our mask arg is on, then 
	     this pixel is black. */
	  if( FunV_tErr_InvalidAnatomicalVoxel == eFunctional &&
	      this->mabDisplayFlags[DspA_tDisplayFlag_MaskToFunctionalOverlay]) {
	    color.mfRed   = 0;
	    color.mfGreen = 0;
	    color.mfBlue  = 0;
	  }
	}
	
      } else {
        
        /* voxel was out of bounds, set to out of bounds color. */
        color.mfRed   = 0;
        color.mfGreen = 0;
        color.mfBlue  = 0;
        
        /* clear error flag. */
        eResult = DspA_tErr_NoErr;
      }
      
      /* set the pixel */
      pDest[DspA_knRedPixelCompIndex]   = 
        (GLubyte)(color.mfRed * (float)DspA_knMaxPixelValue);
      pDest[DspA_knGreenPixelCompIndex] = 
        (GLubyte)(color.mfGreen * (float)DspA_knMaxPixelValue);
      pDest[DspA_knBluePixelCompIndex]  = 
        (GLubyte)(color.mfBlue * (float)DspA_knMaxPixelValue);
      pDest[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;
      
      /* advance our pointer. */
      pDest += DspA_knNumBytesPerPixel;
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
  
  /* delete the voxel. */
  xVoxl_Delete ( &pVoxel );
  xVoxl_Delete ( &pFuncMin );
  xVoxl_Delete ( &pFuncMax );
  
  //  xUtil_StopTimer( "build frame buffer" );
  
  return eResult;
}


DspA_tErr DspA_DrawFunctionalOverlayToFrame_ ( tkmDisplayAreaRef this ) {
  
  DspA_tErr             eResult     = DspA_tErr_NoErr;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  xPoint2n              bufferPt    = {0,0};
  FunV_tFunctionalValue max         = 0;
  FunV_tFunctionalValue funcValue   = 0.0;
  xColor3f              color       = {0,0,0};
  xColor3f              funcColor   = {0,0,0};
  GLubyte*              pFrame      = NULL;
  
  /* if we're drawing the color scale bar... */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalColorScaleBar] ) {
    
    /* draw the color scale bar. get threshold max. */
    eFunctional = FunV_GetThresholdMax( this->mpFunctionalVolume, &max );
    if( FunV_tErr_NoError != eFunctional ) {
      eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
      goto error;
    }
    
    /* down the buffer */
    for ( bufferPt.mnY = 0; 
	  bufferPt.mnY < this->mnVolumeSizeY; 
	  bufferPt.mnY++ ) {
      
      /* get an interpolated value within the range of -max to +max 
	 determined by the y value */
      funcValue = (FunV_tFunctionalValue) 
        ( (float)bufferPt.mnY * 
          (float)((max*2.0)/(float)this->mnVolumeSizeY) - max );
      
      /* get the functional color for this value */
      eFunctional = FunV_GetColorForValue( this->mpFunctionalVolume,
                                           funcValue, &color, &funcColor );
      if( FunV_tErr_NoError != eFunctional ) {
        eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
        goto error;
      }
      
      /* draw on the right side... */
      for ( bufferPt.mnX = this->mnVolumeSizeX - 10; 
            bufferPt.mnX < this->mnVolumeSizeX; bufferPt.mnX++ ) {
        
        /* write it back to the buffer. */
        pFrame = this->mpFrameBuffer + 
          ( (bufferPt.mnY * this->mnVolumeSizeY) + bufferPt.mnX ) * 
          DspA_knNumBytesPerPixel;
        
        pFrame[DspA_knRedPixelCompIndex]   = 
          (GLubyte)(funcColor.mfRed * DspA_knMaxPixelValue);
        pFrame[DspA_knGreenPixelCompIndex] = 
          (GLubyte)(funcColor.mfGreen * DspA_knMaxPixelValue);
        pFrame[DspA_knBluePixelCompIndex]  = 
          (GLubyte)(funcColor.mfBlue * DspA_knMaxPixelValue);
        pFrame[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;
      }
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
  xVoxelRef    controlPt  = NULL;
  xVoxel       convertedPt;
  xPoint2n     bufferPt   = {0,0};
  float        faColor[3] = {0, 0, 0};
  
  /* decide which list we want out of the space. */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    e3DList = x3Lst_GetItemsInZPlane( this->mpControlPoints, 
				      DspA_GetCurrentSliceNumber_(this),
				      &list );
    break;
  case mri_tOrientation_Sagittal:
    e3DList = x3Lst_GetItemsInXPlane( this->mpControlPoints, 
				      DspA_GetCurrentSliceNumber_(this),
				      &list );
    break;
  case mri_tOrientation_Horizontal:
    e3DList = x3Lst_GetItemsInYPlane( this->mpControlPoints, 
				      DspA_GetCurrentSliceNumber_(this),
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
    while( (eList = xList_NextFromPos( list, (void**)&controlPt )) 
	   != xList_tErr_EndOfList ) {
      
      if( controlPt ) {
	
	/* convert the control point to be in the middle of voxel */
	xVoxl_Copy( &convertedPt, controlPt );
	if( xVoxl_GetX( &convertedPt ) == xVoxl_GetFloatX( &convertedPt ) )
	  xVoxl_SetFloatX( &convertedPt, xVoxl_GetFloatX( &convertedPt ) + 0.5 );
	if( xVoxl_GetY( &convertedPt ) == xVoxl_GetFloatY( &convertedPt ) )
	  xVoxl_SetFloatY( &convertedPt, xVoxl_GetFloatY( &convertedPt ) + 0.5 );
	if( xVoxl_GetZ( &convertedPt ) == xVoxl_GetFloatZ( &convertedPt ) )
	  xVoxl_SetFloatZ( &convertedPt, xVoxl_GetFloatZ( &convertedPt ) + 0.5 );
	
	/* convert to buffer point. */
	eResult = DspA_ConvertVolumeToBuffer_ ( this, &convertedPt, &bufferPt );
	if ( DspA_tErr_NoErr != eResult )
	  goto error;
	
	/* draw a point here. */
	DspA_DrawMarker_( this, DspA_tMarker_Crosshair, faColor, &bufferPt, 
			  DspA_knControlPointCrosshairSize );
      }
    }
    
    if( eList != xList_tErr_EndOfList )
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
  //  xVoxel       start;
  //  xVoxel       direction;
  
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
#if 0  
  xVoxl_Set( &start, 120.5, 120.5, 120.5 );
  xVoxl_SetFloat( &direction, .5, .3, .1 );
  //  while( xVoxl_IncrementWithMinUntilLimit( &start, 120.5, 121.5 ) ) {
  
  DspA_DrawVector_( this, faColor, &start, &direction );
  //  }
#endif
  
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
  
  DspA_tErr     eResult           = DspA_tErr_NoErr;
  xListRef      list       = NULL;
  xList_tErr    eList      = xList_tErr_NoErr;
  x3Lst_tErr    e3DList    = x3Lst_tErr_NoErr;
  xVoxelRef     selection         = NULL;
  xPoint2n      bufferPt          = {0,0};
  GLubyte*      pFrame            = NULL;
  xColor3f      color;
  int           nX                = 0;
  int           nY                = 0;
  
  
  /* decide which list we want out of the space. */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    e3DList = x3Lst_GetItemsInZPlane( this->mpSelection,
				      DspA_GetCurrentSliceNumber_(this),
				      &list );
    break;
  case mri_tOrientation_Sagittal:
    e3DList = x3Lst_GetItemsInXPlane( this->mpSelection, 
				      DspA_GetCurrentSliceNumber_(this),
				      &list );
    break;
  case mri_tOrientation_Horizontal:
    e3DList = x3Lst_GetItemsInYPlane( this->mpSelection,
				      DspA_GetCurrentSliceNumber_(this),
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
    
    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&selection )) 
	   != xList_tErr_EndOfList ) {
      
      if( selection ) {
	
	/* covert to a buffer point. */
	eResult = DspA_ConvertVolumeToBuffer_ ( this, selection, &bufferPt );
	if ( DspA_tErr_NoErr != eResult )
	  goto error;
	
	/* y flip the volume pt to flip the image over. */
	bufferPt.mnY = Y_FLIP(bufferPt.mnY);
	
	/* write it back to the buffer */
	for( nY = bufferPt.mnY; 
	     BUFFER_Y_LT(nY,BUFFER_Y_INC(bufferPt.mnY,this->mnZoomLevel));
	     nY = BUFFER_Y_INC(nY,1) ) {
	  for( nX = bufferPt.mnX; nX < bufferPt.mnX + this->mnZoomLevel;nX++) {
	    
	    pFrame = this->mpFrameBuffer + 
	      ( (nY * this->mnVolumeSizeX) + nX ) * DspA_knNumBytesPerPixel;
	    
	    /* get the current color in the buffer */
	    xColr_Set( &color, (float)pFrame[DspA_knRedPixelCompIndex] /
		       (float)DspA_knMaxPixelValue,
		       (float)pFrame[DspA_knGreenPixelCompIndex] /
		       (float)DspA_knMaxPixelValue,
		       (float)pFrame[DspA_knBluePixelCompIndex] /
		       (float)DspA_knMaxPixelValue );
	    
	    /* make it greener */
	    xColr_HilightComponent( &color, xColr_tComponent_Green );
	    
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
      }
    }
    
    if( eList != xList_tErr_EndOfList )
      goto error;
  }
  
  goto cleanup;
  
 error:
  
  if ( x3Lst_tErr_NoErr != e3DList )
    eResult = DspA_tErr_ErrorAccessingSelection;
  
  if ( xList_tErr_NoErr != eList )
    eResult = DspA_tErr_ErrorAccessingSelection;
  
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
  if( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] ) {
    plane = HPtL_tIterationPlane_All;
  } else {
    switch( this->mOrientation ) {
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
  if( HPtL_tErr_NoErr != eHeadPtList ) {
    goto error;
  }
  
  while( (eHeadPtList = HPtL_NextPoint( this->mHeadPoints, &pHeadPt ))
	 == HPtL_tErr_NoErr ) {
    
    /* color is red if hilited, purple if cardinal, green if not. */
    if( pHeadPt == this->mpSelectedHeadPoint ) {
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
  
  if( eHeadPtList != HPtL_tErr_LastPoint )
    goto error;
  
  /* clean up error codes */
  eHeadPtList = HPtL_tErr_NoErr;
  eResult = DspA_tErr_NoErr;
  
  goto cleanup;
  
 error:
  
  if( HPtL_tErr_NoErr != eHeadPtList ) {
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
					tkm_tSurfaceType  iType) {
  
  DspA_tErr           eResult              = DspA_tErr_NoErr;
  Surf_tErr           eSurface             = Surf_tErr_NoErr;
  xGArr_tErr          eList                = xGArr_tErr_NoErr;
  xGrowableArrayRef   list                 = NULL;
  int                 nSlice               = 0;
  Surf_tVertexSet     surface              = Surf_tVertexSet_Main;
  xPoint2f            zeroPoint            = {0,0};
  xVoxelRef           curPlane             = NULL;
  DspA_tSurfaceListNode drawListNode;
  xVoxelRef           anaVertex            = NULL; 
  xVoxelRef           anaNeighborVertex    = NULL;
  int                 nVertexIndex         = -1;
  int                 nNeighborVertexIndex = -1;
  xPoint2f            intersectionPt       = {0, 0};
  xPoint2f            interpIntersectionPt = {0, 0};
  tBoolean            bPointsOnThisFace    = FALSE;
  
  xVoxl_New( &curPlane );
  xVoxl_New( &anaVertex );
  xVoxl_New( &anaNeighborVertex );
  
  /* get the current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );
  
  /* set up for gl drawing */
  DspA_SetUpOpenGLPort_( this );
  
  /* for each vertex type... */
  for ( surface = Surf_tVertexSet_Main;
	surface < Surf_knNumVertexSets; surface++ ) {
    
    /* only build if this surface is being displayed */
    if( !this->mabDisplayFlags[surface + DspA_tDisplayFlag_MainSurface] ) {
      continue;
    }
    
    /* check to see if we already have this list. if so, don't build it. */
    list = DspA_GetSurfaceList_( this, iType, this->mOrientation, surface,
				 DspA_GetCurrentSliceNumber_(this) );
    if( NULL != list ) {
      continue;
    }
    
    /* make a new list. */
    DspA_NewSurfaceList_( this, iType, this->mOrientation, surface, 
			  DspA_GetCurrentSliceNumber_(this) );
    
    /* get the list. */
    list = DspA_GetSurfaceList_( this, iType, this->mOrientation, surface,
				 DspA_GetCurrentSliceNumber_(this) );
    if( NULL == list ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }
    
    //    xUtil_StartTimer ();
    
    /* make a voxel with our current slice in it and
       set the iterator to our view */
    DspA_ConvertPlaneToVolume_ ( this, &zeroPoint, nSlice, this->mOrientation,
				 curPlane ); 
    eSurface = Surf_SetIteratorPosition( this->mpSurface[iType], curPlane );
    if( Surf_tErr_NoErr != eSurface )
      goto error;
    
    bPointsOnThisFace = FALSE;
    
    /* while we have vertices to check.. */
    while( (eSurface = Surf_GetNextAndNeighborVertex( this->mpSurface[iType],
						      surface,
						      anaVertex, 
						      &nVertexIndex,
						      anaNeighborVertex, 
						      &nNeighborVertexIndex ))
	   != Surf_tErr_LastFace ) {
      
      /* if the line between these two points intersects the
	 current plane... */
      DspA_NormalizeVoxel_( anaVertex, 
			    this->mOrientation, anaVertex );
      DspA_NormalizeVoxel_( anaNeighborVertex, 
			    this->mOrientation, anaNeighborVertex );
      if ( xUtil_LineIntersectsPlane( anaVertex, anaNeighborVertex, nSlice, 
				      &intersectionPt, &interpIntersectionPt)){
	/* fill out a node. */
	drawListNode.mbVertex = TRUE;
	drawListNode.mnOriginalVertexIndex = nVertexIndex;
	drawListNode.mIntersectionPoint = intersectionPt;
	drawListNode.mInterpIntersectionPoint = interpIntersectionPt;

	/* the original vertex is just the anatomical vertex. */
	xVoxl_Copy( &drawListNode.mOriginalVertex, anaVertex );

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
	if( xGArr_tErr_NoErr != eList ) {
	  DebugPrint( ("xGArr error %d in DspA_BuildSurfaceDrawLists_: %s\n",
		       eList, xGArr_GetErrorString( eList ) ) );
	  eResult = DspA_tErr_ErrorAccessingSurfaceList;
	  goto error;
	}

	bPointsOnThisFace = TRUE;
      }
      
      /* if we have a last vertex, and we drew points on this
	 face add a face marker. */
      if( eSurface == Surf_tErr_LastVertex
	  && bPointsOnThisFace ) {
	
	/* fill out a node. */
	drawListNode.mbVertex = FALSE;

	eList = xGArr_Add( list, &drawListNode );
	if( xGArr_tErr_NoErr != eList ) {
	  DebugPrint( ("xGArr error %d in DspA_BuildSurfaceDrawLists_: %s\n",
		       eList, xGArr_GetErrorString( eList ) ) );
	  eResult = DspA_tErr_ErrorAccessingSurfaceList;
	  goto error;
	}
	
	bPointsOnThisFace = FALSE;
      }
      
    }
    
    /* if surface error is just on last face, clear it. */
    if( Surf_tErr_LastFace == eSurface ) 
      eSurface = Surf_tErr_NoErr;
    
    /* check for other errors. */
    if( Surf_tErr_NoErr != eSurface )
      goto error;
    
    //    xUtil_StopTimer( "build surface list" );
    
  } /* for each surface... */
  
  goto cleanup;
  
 error:
  
  if( Surf_tErr_NoErr != eSurface )
    eResult = DspA_tErr_ErrorAccessingSurface;
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_BuildSurfaceDrawLists_: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  xVoxl_Delete( &curPlane );
  xVoxl_Delete( &anaVertex );
  xVoxl_Delete( &anaNeighborVertex );
  
  return eResult;
}

DspA_tErr DspA_DrawMarker_ ( tkmDisplayAreaRef this,
			     DspA_tMarker      iType,
			     float*            ifaColor,
			     xPoint2nRef       ipWhere,
			     int               inSize ) {
  
  DspA_tErr eResult = DspA_tErr_NoErr;
  int       nWidth  = 0;
  int       nHeight = 0;
  
  DspA_SetUpOpenGLPort_( this );
  
  glLineWidth( 1 );
  
  switch( iType ) {
    
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
  
#if 0
  fprintf( stderr, "vox start %.2f, %.2f, %.2f dir %.2f, %.2f, %.2f "
	   "bufst %.2f, %.2f bufend %.2f, %.2f\n",
	   xVoxl_ExpandFloat(ipVoxelStart),
	   xVoxl_ExpandFloat(ipVoxelDirection),  
	   bufferStart.mfX, bufferStart.mfY,
	   bufferEnd.mfX, bufferEnd.mfY );
#endif  
  
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
  
  /* return the cursor */
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
					   tkm_tSurfaceType  iType,
					   Surf_tVertexSet   iSet,
					   xVoxelRef         iAnaIdx,
					   xVoxelRef         oOrigAnaIdx,
					   xVoxelRef         oInterpAnaIdx,
					   char*             osDescription ) {

  DspA_tErr             eResult         = DspA_tErr_NoErr;
  xGArr_tErr            eList           = xGArr_tErr_NoErr;
  int                   nSlice          = 0;
  xGrowableArrayRef     list            = NULL;
  DspA_tSurfaceListNode drawListNode;
  float                 dx              = 0;
  float                 dy              = 0;
  float                 dz              = 0;
  float                 fDistance       = 0;
  float                 fLowestDistance = 0;
  int                   nClosestIndex   = 0;
  xVoxel                closestAnaIdx;
  xVoxel                closestInterpAnaIdx;

  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if( DspA_tErr_NoErr != eResult )
    goto error;
  if( iType < 0 || iType >= tkm_knNumSurfaceTypes ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  if( iSet < 0 || iSet >= Surf_knNumVertexSets ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }
  if( NULL == iAnaIdx ||
      NULL == oInterpAnaIdx ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* walk through all surface lists. for each one... */
  fLowestDistance = 99999;
  nClosestIndex = -1;
  for( nSlice = 0; nSlice < 256; nSlice++ ) {
    
    list = DspA_GetSurfaceList_( this, iType, 
				 this->mOrientation, iSet, nSlice );
    if( NULL == list ) {
      continue;
#if 0
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
#endif
    }
    
    /* walk through the list. */
    eList = xGArr_ResetIterator( list );
    if( xGArr_tErr_NoErr != eList )
      goto error;
    while( (eList = xGArr_NextItem( list, (void*)&drawListNode ))
	    == xGArr_tErr_NoErr ) {
      
      /* skip if it's not a vertex. */
      if( !drawListNode.mbVertex ) {
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
	xVoxl_Copy( &closestInterpAnaIdx, &drawListNode.mInterpVertex );
	xVoxl_Copy( &closestAnaIdx, &drawListNode.mOriginalVertex );
	fLowestDistance = fDistance;
      }
    }
  }  

  if( -1 == nClosestIndex ) {
    eResult = DspA_tErr_CouldntFindClosestVoxel;
    goto error;
  }
  
  /* Copy the closest interp vertex out out. Copy the original if they
     want it. */
  xVoxl_Copy( oOrigAnaIdx, &closestAnaIdx );
  if( NULL != oInterpAnaIdx ) {
    xVoxl_Copy( oInterpAnaIdx, &closestInterpAnaIdx );
  }

  /* calculate the actual distance. */
  fLowestDistance = sqrt( fLowestDistance );

  /* If the want one, make a string of info. */
  if( NULL != osDescription ) {
    sprintf( osDescription, "Index: %d Distance: %.2f\n"
	     "\t    Original vertex RAS Coords: %.2f %.2f %.2f\n"
	     "\tInterpolated vertex RAS Coords: %.2f %.2f %.2f", 
	     nClosestIndex, fLowestDistance, 
	     xVoxl_ExpandFloat( &closestAnaIdx ),
	     xVoxl_ExpandFloat( &closestInterpAnaIdx ) );
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


DspA_tErr DspA_InitSurfaceLists_( tkmDisplayAreaRef this,
				  int               inNumLists ) {
  
  DspA_tErr           eResult    = DspA_tErr_NoErr;
  int                 nSurface   = 0;
  int                 nDrawList  = 0;
  
  /* for each surface type */
  for ( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {
    
    /* allocate surface point lists. */
    this->maSurfaceLists[nSurface] = (xGrowableArrayRef*) 
      malloc (sizeof(xGrowableArrayRef) * inNumLists );
    
    /* set all lists to null */
    for( nDrawList = inNumLists - 1; nDrawList >= 0; nDrawList-- ) {
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
  
  for( nSurface = 0; nSurface < tkm_knNumSurfaceTypes; nSurface++ ) {
    
    for( nDrawList = DspA_GetNumSurfaceLists_( this ) - 1; 
	 nDrawList >= 0; nDrawList-- ) {
      
      /* if this is a list... */
      if( NULL != this->maSurfaceLists[nSurface][nDrawList] ) {
	
	/* delete the list. */
	eList = xGArr_Delete( &this->maSurfaceLists[nSurface][nDrawList] );
	if( xGArr_tErr_NoErr != eList ) {
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
				 tkm_tSurfaceType  iType,
				 mri_tOrientation  iOrientation,
				 Surf_tVertexSet   iSurface,
				 int               inSlice ) {
  
  DspA_tErr   eResult   = DspA_tErr_NoErr;
  int         nDrawList = 0;
  xGArr_tErr  eList     = xGArr_tErr_NoErr;
  
  /* get the list index. */
  nDrawList = DspA_GetSurfaceListIndex_( this, iOrientation, 
					 iSurface, inSlice );
  
  /* allocate a list. */
  xGArr_New( &this->maSurfaceLists[iType][nDrawList],
	     sizeof( DspA_tSurfaceListNode ), 512 );
  if( xGArr_tErr_NoErr != eList ) {
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
					 tkm_tSurfaceType  iType,
					 mri_tOrientation  iOrientation,
					 Surf_tVertexSet  iSurface,
					 int               inSlice ) {
  
  xGrowableArrayRef pList = NULL;
  
  if( NULL != this->maSurfaceLists ) {
    pList = this->maSurfaceLists [iType]
      [ DspA_GetSurfaceListIndex_( this, iOrientation, iSurface, inSlice ) ];
  }
  
  return pList;
}

int DspA_GetNumSurfaceLists_ ( tkmDisplayAreaRef this ) {
  
  /* may be incorrect, assumes dimensions are all equal */
  return this->mnVolumeSizeZ * Surf_knNumVertexSets * mri_knNumOrientations;
}

int DspA_GetSurfaceListIndex_ ( tkmDisplayAreaRef this,
				mri_tOrientation  iOrientation,
				Surf_tVertexSet  iSurface,
				int               inSlice ) {
  
  /* may be incorrect, assumes dimensions are all equal */
  return ((int)iOrientation * Surf_knNumVertexSets * this->mnVolumeSizeZ) +
    ((int)iSurface * this->mnVolumeSizeX) + inSlice;
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
  fY = (fZoomLevel * (fY - 
		      GLDRAW_Y_FLIP_FLOAT(xVoxl_GetFloatY(this->mpZoomCenter)))) +
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
  opBufferPt->mfY = (fZoomLevel * (fY - 
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
  
  switch( iOrientation ) {
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
  
  switch( iOrientation ) {
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
  char      sTclArguments[STRLEN]         = "";
  int       nFlag                      = 0;
  int       brush                      = 0;
  int       surface                    = 0;
  
  /* send the point info for our cursor */
  DspA_SendPointInformationToTcl_( this, DspA_tDisplaySet_Cursor,
				   this->mpCursor );
  
  /* send volume name */
  Volm_CopyVolumeName( this->mpVolume, sVolumeName, sizeof(sVolumeName) );
  sprintf( sTclArguments, "\"%s value\"", sVolumeName );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );
  
  /* send the aux volume name if it's loaded */
  if( NULL != this->mpAuxVolume ) {
    Volm_CopyVolumeName( this->mpAuxVolume, sVolumeName, sizeof(sVolumeName) );
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
  for( nFlag = 0; nFlag < DspA_knNumDisplayFlags; nFlag++ ) {
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
  for( brush = 0; brush < DspA_knNumBrushes; brush++ ) {
    sprintf ( sTclArguments, "%d %d %d %d",
	      (int)brush,
	      (int)sBrush.mInfo[brush].mnLow, 
	      (int)sBrush.mInfo[brush].mnHigh,
	      (int)sBrush.mInfo[brush].mnNewValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushInfo, sTclArguments );
  }
  
  /* send the cursor color */
  sprintf ( sTclArguments, "%f %f %f",
	    sCursorColor.mfRed, sCursorColor.mfGreen, sCursorColor.mfBlue );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateCursorColor, sTclArguments );
  
  /* send the surface line info */
  for( surface = 0; surface < Surf_knNumVertexSets; surface++ ) {
    sprintf ( sTclArguments, "%d %d", surface,  
	      this->manSurfaceLineWidth[surface] );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceLineWidth, 
			sTclArguments );
    sprintf ( sTclArguments, "%d %f %f %f", surface,
	      xColr_ExpandFloat( &(this->maSurfaceLineColor[surface]) ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceLineColor, 
			sTclArguments );
  }
  
  return eResult;
}

DspA_tErr DspA_SendPointInformationToTcl_ ( tkmDisplayAreaRef this,
					    DspA_tDisplaySet  iSet,
					    xVoxelRef         iAnaIdx ) {
  
  xVoxel                voxel;
  char                  sTclArguments[STRLEN] = "";
  int                   nSlice             = 0;
  float                 fVolumeValue       = 0;
  HPtL_tHeadPointRef    pHeadPoint         = NULL;
  FunV_tErr             eFunctional        = FunV_tErr_NoError;
  xVoxel                funcIdx;
  xVoxel                funcRAS;
  FunV_tFunctionalValue funcValue          = 0;
  tBoolean              bFuncSelection     = FALSE;
  int                   nROIIndex          = 0;
  char                  sLabel[STRLEN]        = "";
  int                   nValue             = 0;
  DspA_tHistogramParams histoParams;
  float                 fDistance          = 0;
  
  
  /* send the anatomical index. */
  sprintf( sTclArguments, "%s %d %d %d", 
	   DspA_ksaDisplaySet[iSet], xVoxl_ExpandInt( iAnaIdx ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeCursor, sTclArguments );
  
  /* send the slice number */
  nSlice = DspA_GetCurrentSliceNumber_( this );
  sprintf( sTclArguments, "%d", nSlice );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeSlice, sTclArguments );
  
  /* also convert to RAS and send those coords along. */
  Volm_ConvertIdxToRAS( this->mpVolume, iAnaIdx, &voxel );
  sprintf( sTclArguments, "%s %.1f %.1f %.1f", 
	   DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateRASCursor, sTclArguments );
  
  /* also convert to mni and send those coords along. */
  if (NULL != this->mpVolume->mpMriValues->linear_transform) {
    Volm_ConvertIdxToMNITal( this->mpVolume, iAnaIdx, &voxel );
    sprintf( sTclArguments, "%s %.1f %.1f %.1f", 
             DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateMNICursor, sTclArguments );
    
    /* and the tal coords */
    Volm_ConvertIdxToTal( this->mpVolume, iAnaIdx, &voxel );
    sprintf( sTclArguments, "%s %.1f %.1f %.1f", 
             DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateTalCursor, sTclArguments );
  }
  
  /* and the scanner coords */
  Volm_ConvertIdxToScanner( this->mpVolume, iAnaIdx, &voxel );
  sprintf( sTclArguments, "%s %.1f %.1f %.1f", 
	   DspA_ksaDisplaySet[iSet], xVoxl_ExpandFloat( &voxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateScannerCursor, sTclArguments );
  
  /* also get the volume value and send that along. */
  Volm_GetValueAtIdx( this->mpVolume, iAnaIdx, &fVolumeValue );
  switch (this->mpVolume->mpMriValues->type)
    {
    default:
      sprintf( sTclArguments, "%s %d", 
	       DspA_ksaDisplaySet[iSet], (int)fVolumeValue );
      break ;
    case MRI_FLOAT:
      {
	char fmt[STRLEN] ;
	float f ;
	int   decs ;
	
	f = fabs(fVolumeValue) ;
	if (f > 1)
	  decs = 2 ;
	else if (f > .1)
	  decs = 3 ;
	else if (f > .01)
	  decs = 4 ;
	else if (f > 0.001)
	  decs = 5 ;
	else
	  decs = 6 ;
	
	sprintf(fmt, "%%s %%2.%df", decs) ;
	sprintf( sTclArguments, fmt, 
		 DspA_ksaDisplaySet[iSet], fVolumeValue );
	break ;
      }
    }
  
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValue, sTclArguments );
  
  /* send aux volume value if it's loaded. */
  if( NULL != this->mpAuxVolume ) 
    {
      Volm_GetValueAtIdx( this->mpAuxVolume, iAnaIdx, &fVolumeValue );
      switch (this->mpAuxVolume->mpMriValues->type)
	{
	default:
	  sprintf( sTclArguments, "%s %d", 
		   DspA_ksaDisplaySet[iSet], (int)fVolumeValue );
	  break ;
	case MRI_FLOAT:
	  {
	    char fmt[STRLEN] ;
	    float f ;
	    int   decs ;
	    
	    f = fabs(fVolumeValue) ;
	    if (f > 1)
	      decs = 2 ;
	    else if (f > .1)
	      decs = 3 ;
	    else if (f > .01)
	      decs = 4 ;
	    else if (f > 0.001)
	      decs = 5 ;
	    else
	      decs = 6 ;
	    
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
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ) {
    
    DisableDebuggingOutput;
    if( DspA_tDisplaySet_Cursor == iSet ) {
      
      /* if this is the cursor, use the voxel locations for the selection
	 and the average value */
      eFunctional = FunV_GetAvgFunctionalValue( this->mpFunctionalVolume,
						&funcValue, &funcIdx, &funcRAS,
						&bFuncSelection );
      
    } else if( DspA_tDisplaySet_Mouseover == iSet ) {
      
      /* if this is the mouseover, use the value at the current point,
	 and convert the point that to the func idx and ras. */
      eFunctional = FunV_GetValueAtAnaIdx( this->mpFunctionalVolume,
					   iAnaIdx, &funcValue );
      if( FunV_tErr_NoError == eFunctional ) {
	
	/* convert the points */
	FunV_ConvertAnaIdxToFuncIdx( this->mpFunctionalVolume,
				     iAnaIdx, &funcIdx );
	FunV_ConvertAnaIdxToFuncRAS( this->mpFunctionalVolume,
				     iAnaIdx, &funcRAS );
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
  
  /* and the roi label if we have one */
  if( NULL != this->mROIGroup ) {
    tkm_GetROILabel( iAnaIdx, &nROIIndex, sLabel );
    
    /* if this is a click, set the index */
    if( DspA_tDisplaySet_Cursor == iSet ) {
      this->mnROIGroupIndex = nROIIndex;
    }
    
    /* if this is a click with the edit tool and the volume count flag is
       on, calc and display the volume */
    if( DspA_tDisplaySet_Cursor == iSet &&
        DspA_tTool_EditSegmentation == sTool &&
        this->mabDisplayFlags[DspA_tDisplayFlag_ROIVolumeCount] ) {
      tkm_CalcROIVolume( iAnaIdx, &nValue );
      sprintf( sTclArguments, "%s \"%s (%d)\"",
               DspA_ksaDisplaySet[iSet], sLabel, nValue );
      
      /* else just the label */
    } else {
      sprintf( sTclArguments, "%s \"%s\"",
               DspA_ksaDisplaySet[iSet], sLabel );
    }
    tkm_SendTclCommand( tkm_tTclCommand_UpdateROILabel, sTclArguments );
  }
  
  /* and the head point label if it's on */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_HeadPoints]
      && NULL != this->mHeadPoints ) {
    
    /* get the closest head point */
    tkm_GetHeadPoint( iAnaIdx, this->mOrientation, 
		      this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj],
		      &pHeadPoint );
    if( NULL != pHeadPoint ) {
      
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
  if( NULL != this->mGCAVolume &&
      DspA_tDisplaySet_Cursor == iSet ) {
    GCAdump( this->mGCAVolume, this->mpVolume->mpMriValues,
	     xVoxl_ExpandInt( iAnaIdx ), this->mGCATransform, stdout, 0 );
  }
  
  if( NULL != this->mVLI1 &&
      DspA_tDisplaySet_Cursor == iSet ) 
    {
      int xn, yn, zn, label_counts_c1[256], label_counts_c2[256], index,
        inNumValues, n1, n2 ;
      VL *vl ;
      char name1[STRLEN], name2[STRLEN], *cp ;
      
      FileNameOnly(this->isVLI1_name, name1) ;
      cp = strrchr(name1, '.') ; if (cp) *cp = 0 ;
      FileNameOnly(this->isVLI2_name, name2) ;
      cp = strrchr(name2, '.') ; if (cp) *cp = 0 ;
      xn = nint(xVoxl_GetX(iAnaIdx) / this->mVLI1->resolution) ;
      yn = nint(xVoxl_GetY(iAnaIdx) / this->mVLI1->resolution) ;
      zn = nint(xVoxl_GetZ(iAnaIdx) / this->mVLI1->resolution) ;
      
#if 0
      printf("voxel (%d, %d, %d) --> node (%d, %d, %d)\n",
	     xVoxl_ExpandInt( iAnaIdx ), xn, yn, zn) ;
      vl = &this->mVLI1->vl[xn][yn][zn] ;
      printf("%s:\n", name1) ;
      for (n = 0 ; n < vl->nlabels ; n++)
	{
	  if (vl->counts[n] > 0)
	    printf("%s (%d): %d\n", cma_label_to_name(vl->labels[n]),
		   vl->labels[n], vl->counts[n]) ;
	}
      vl = &this->mVLI2->vl[xn][yn][zn] ;
      printf("%s:\n", name2) ;
      for (n = 0 ; n < vl->nlabels ; n++)
	{
	  if (vl->counts[n] > 0)
	    printf("%s (%d): %d\n", cma_label_to_name(vl->labels[n]),
		   vl->labels[n], vl->counts[n]) ;
	}
#endif
      
      memset(label_counts_c1, 0, sizeof(label_counts_c1)) ;
      vl = &this->mVLI1->vl[xn][yn][zn] ;
      for (n1 = index = 0 ; index < vl->nlabels ; index++)
	{
	  label_counts_c1[vl->labels[index]] += vl->counts[index] ;
	  n1 += vl->counts[index] ;
	}
      
      memset(label_counts_c2, 0, sizeof(label_counts_c2)) ;
      vl = &this->mVLI2->vl[xn][yn][zn] ;
      for (n2 = index = 0 ; index < vl->nlabels ; index++)
	{
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
      for (nValue = index = 0 ; index < 256 ; index++)
	{
	  if ((label_counts_c1[index]) == 0 && label_counts_c2[index] == 0)
	    continue ;
	  
	  /* assign values for the elements */
	  histoParams.mafValues1[nValue] = (float)label_counts_c1[index] ;
	  histoParams.mafValues2[nValue] = (float)label_counts_c2[index] ;
	  if (this->mabDisplayFlags[DspA_tDisplayFlag_HistogramPercentChange])
	    {
	      xUtil_strncpy( histoParams.msYAxisTitle, "% of voxels with label",
			     sizeof(histoParams.msYAxisTitle) );
	      histoParams.mafValues1[nValue] /= ((float)n1/100) ;
	      histoParams.mafValues2[nValue] /= ((float)n2/100) ;
	    }
	  else
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
      for( nValue = 0; nValue < inNumValues; nValue++ )
	free( histoParams.masXAxisLabels[nValue] );
      free( histoParams.masXAxisLabels );
    }
  
#if 0
  /* histogram example */
  if( getenv("SEE_HISTO_DEMO") ) {
    if( DspA_tDisplaySet_Cursor == iSet ) {
      
      /* the title of the window and graph */
      DebugNote( ("Sprintfing coords") );
      xUtil_snprintf( histoParams.msTitle, sizeof(histoParams.msTitle),
		      "%d, %d, %d", xVoxl_ExpandInt( iAnaIdx ) );
      /* the titles of the axes */
      DebugNote( ("Copying axis titles") );
      xUtil_strncpy( histoParams.msXAxisTitle, "The Kinds of Stuff",
		     sizeof(histoParams.msXAxisTitle) );
      xUtil_strncpy( histoParams.msYAxisTitle, "Amount of Stuff",
		     sizeof(histoParams.msYAxisTitle) );
      /* the title of the elements in the legend */
      DebugNote( ("Copying element titles") );
      xUtil_strncpy( histoParams.msLabel1, "Some Stuff", 
		     sizeof(histoParams.msLabel1) );
      xUtil_strncpy( histoParams.msLabel2, "Other Stuff",
		     sizeof(histoParams.msLabel2) );
      
      /* declare arrays of values and labels */
#define knNumValues 20
      histoParams.mafValues1 = (float*) calloc( knNumValues, sizeof(float) );
      histoParams.mafValues2 = (float*) calloc( knNumValues, sizeof(float) );
      histoParams.masXAxisLabels = (char**) calloc( knNumValues, sizeof(float) );
      
      /* fill them up with some random nubmers and cma labels */
      histoParams.mnNumValues = knNumValues;
      for( nValue = 0; nValue < knNumValues; nValue++ ) {
	/* assign values for the elements */
	histoParams.mafValues1[nValue] = random()%100;
	histoParams.mafValues2[nValue] = random()%100;
	/* allocate and set the label for this element */
	histoParams.masXAxisLabels[nValue] = 
	  (char*) malloc( sizeof(char) * DspA_knHistoTitleLength );
	xUtil_strncpy( histoParams.masXAxisLabels[nValue], 
		       cma_label_to_name( nValue ), DspA_knHistoTitleLength );
      }
      
      /* draw the thing */
      DspA_DrawHistogram( this, &histoParams );
      
      free( histoParams.mafValues1 );
      free( histoParams.mafValues2 );
      for( nValue = 0; nValue < knNumValues; nValue++ )
	free( histoParams.masXAxisLabels[nValue] );
      free( histoParams.masXAxisLabels );
    }
  }
#endif
  
  
  /* if we have a surfaec, update the surface distance. if this is the
     cursor, find the distance from the last cursor and print that. if
     this is the mouseover, print the distance from the current
     cursor.  */
  if( NULL != this->mpSurface[tkm_tSurfaceType_Main] ) { 
    
    if( DspA_tDisplaySet_Cursor == iSet ) {
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
  for( nValue = 0; nValue < iParams->mnNumValues; nValue++ ) {
    
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
  if( NULL != this->mpSurface[tkm_tSurfaceType_Main] ) {
    
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
  if ( xVoxl_GetX( ipVoxel )    < 0
       || xVoxl_GetY( ipVoxel ) < 0
       || xVoxl_GetZ( ipVoxel ) < 0
       || xVoxl_GetX( ipVoxel ) >= this->mnVolumeSizeX
       || xVoxl_GetY( ipVoxel ) >= this->mnVolumeSizeY
       || xVoxl_GetZ( ipVoxel ) >= this->mnVolumeSizeZ ) {
    
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
  //glOrtho        ( 0, this->mnVolumeSizeX, 0, this->mnVolumeSizeX, -1.0, 1.0 );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint( ("glOrtho got error %d\n", eGL ) );
  
  glViewport( this->mLocationInSuper.mnX, this->mLocationInSuper.mnY,
	      this->mnWidth, this->mnHeight );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint( ("glViewport got error %d\n", eGL ) );
  
  glRasterPos2i  ( 0, this->mnVolumeSizeX );
  //glRasterPos2i  ( 0, 0 );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint( ("glRasterPos2i got error %d\n", eGL ) );
  
  glPixelZoom ( this->mfFrameBufferScaleX, this->mfFrameBufferScaleY );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
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
  DebugPrint( ("\tcursor %d, %d, %d orientation %s slice %d tool %d\n",
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

  switch( iOrientation ) {
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
  if( DspA_tErr_NoErr != eResult ) {
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

  switch( iOrientation ) {
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
  if( DspA_tErr_NoErr != eResult ) {
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
  if( xVoxl_GetFloatZ(ipAnaIdxB) - xVoxl_GetFloatZ(ipAnaIdxA) != 0.0 ) {
    
    fAlpha = (fPlane - xVoxl_GetFloatZ( ipAnaIdxA )) /
      (xVoxl_GetFloatZ( ipAnaIdxB ) - xVoxl_GetFloatZ( ipAnaIdxA ));
    
  } else {
    
    /* alpha is just 1. */
    fAlpha = 1.0;
  }
  
  /* interpolate to find the intersection. */
  interpIntersectionPt.mfX = (xVoxl_GetFloatX( ipAnaIdxA ) +
	 fAlpha * (xVoxl_GetFloatX(ipAnaIdxB) - xVoxl_GetFloatX(ipAnaIdxA)));
  interpIntersectionPt.mfY = (xVoxl_GetFloatY( ipAnaIdxA ) +
	 fAlpha * (xVoxl_GetFloatY(ipAnaIdxB) - xVoxl_GetFloatY(ipAnaIdxA)));
  
  /* return the points. */
  *opIntersectionPt = intersectionPt;
  *opInterpIntersectionPt = interpIntersectionPt;
  
  return TRUE;
}

DspA_tErr DspA_AdjustSurfaceDrawPoint_( tkmDisplayAreaRef this,
					xPoint2fRef       ipPoint ) {
  
  switch( this->mOrientation ) {
  case mri_tOrientation_Horizontal:
    ipPoint->mfX += 0.5;
    ipPoint->mfY += 0.5;
    break;
  case mri_tOrientation_Coronal:
    ipPoint->mfX += 0.5;
    ipPoint->mfY += 0.5;
    break;
  case mri_tOrientation_Sagittal:
    ipPoint->mfX += 0.5;
    ipPoint->mfY += 0.5;
    break;
  default:
    break;
  }
  
  return DspA_tErr_NoErr;
}


DspA_tErr DspA_AdjustSurfaceAnaIdx ( tkmDisplayAreaRef this,
				     xVoxelRef         iAnaIdx ) {
  
  DspA_tErr eResult = DspA_tErr_NoErr;
  
  /* verify us. */
  eResult = DspA_Verify( this );
  if( DspA_tErr_NoErr != eResult )
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
  if( DspA_tErr_NoErr != eResult ) {
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
  if( DspA_tErr_NoErr != eResult )
    goto error;
  
  xVoxl_SetFloatX( iAnaIdx, xVoxl_GetFloatX(iAnaIdx) - 0.50 );
  xVoxl_SetFloatY( iAnaIdx, xVoxl_GetFloatY(iAnaIdx) - 0.50 );
  xVoxl_SetFloatZ( iAnaIdx, xVoxl_GetFloatZ(iAnaIdx) - 0.50 );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_UnadjustSurfaceAnaIdx: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

DspA_tErr DspA_ParsePointList_( tkmDisplayAreaRef this,
				GLenum            inMode,
				xGrowableArrayRef iList ) {
  
  DspA_tErr  eResult        = DspA_tErr_NoErr;
  xGArr_tErr eList          = xGArr_tErr_NoErr;
  tBoolean   bOperationOpen = FALSE;
  DspA_tSurfaceListNode drawListNode;
  xPoint2f   drawPoint      = {0,0};
  
  /* reset the list position. */
  eList = xGArr_ResetIterator( iList );
  if( xGArr_tErr_NoErr != eList )
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
    if( !drawListNode.mbVertex ) {
      
      /* if operation was still going, end it. */ 
      if( bOperationOpen )
	glEnd();
      
      /* start new operation. */
      glBegin( inMode );
      bOperationOpen = TRUE;
      
    } else {
      
      /* switch on the draw flag to see which point to get. */
      if(this->mabDisplayFlags[DspA_tDisplayFlag_InterpolateSurfaceVertices]) {
	drawPoint.mfX = drawListNode.mInterpIntersectionPoint.mfX;
	drawPoint.mfY = drawListNode.mInterpIntersectionPoint.mfY;
      } else {
	drawPoint.mfX = drawListNode.mIntersectionPoint.mfX;
	drawPoint.mfY = drawListNode.mIntersectionPoint.mfY;
      }	
      
      /* adjust the point */
      DspA_AdjustSurfaceDrawPoint_( this, &drawPoint );
	
      /* convert to zoomed coords. */
      drawPoint.mfX = ((float)this->mnZoomLevel * (drawPoint.mfX - xVoxl_GetFloatX(this->mpZoomCenter))) + (float)(this->mnVolumeSizeX/2.0);
      drawPoint.mfY = ((float)this->mnZoomLevel * (drawPoint.mfY - xVoxl_GetFloatY(this->mpZoomCenter))) + (float)(this->mnVolumeSizeY/2.0);
      
      /* y flip */
      drawPoint.mfY = GLDRAW_Y_FLIP(drawPoint.mfY);
      
      /* and draw the pt. */
      glVertex2f( drawPoint.mfX, drawPoint.mfY );
    }
  }
  
  /* clear flag if last item */
  if( eList == xGArr_tErr_LastItem )
    eList = xGArr_tErr_NoErr;
  
  /* check for other errors */
  if( eList != xGArr_tErr_NoErr )
    goto error;
  
  /* end last operation */
  if( bOperationOpen )
    glEnd();
  
  goto cleanup;
  
 error:
  
  if( xGArr_tErr_NoErr != eList )
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in DspA_ParsePointList_: %s\n",
		 eResult, DspA_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}
