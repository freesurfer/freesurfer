#include "tkmDisplayArea.h"
#include "tkmMeditWindow.h"
#include "tkmFunctionalVolume.h"
#include "xUtilities.h"

/* i'm not sure what to do about these y flips. it seems that whenever we're
   using a point that's going to go into the buffer to be drawn to the screen,
   we should flip when not in horizontal view. when using regular gl drawing 
   commands to draw to the screen, only do it in horizontal orientation. */
#define BUFFER_Y_FLIP(y) ( this->mOrientation != mri_tOrientation_Horizontal? \
                          (this->mnVolumeSize - y) : y )
#define GLDRAW_Y_FLIP(y) ( this->mOrientation == mri_tOrientation_Horizontal? \
                          (this->mnVolumeSize - y) : y )
#define GLDRAW_Y_FLIP_FLOAT(y) \
                         ( this->mOrientation == mri_tOrientation_Horizontal? \
                          ((float)this->mnVolumeSize - y) : y )

//#define BUFFER_Y_FLIP(y) y
//#define GLDRAW_Y_FLIP(y) y

#define Y_FLIP(y)        (this->mnVolumeSize - y)


/* tool and brush info is static */
static DspA_tTool        sTool                = DspA_tTool_SelectVoxels;
static int               snBrushRadius        = 1;
static DspA_tBrushShape  sBrushShape          = DspA_tBrushShape_Circle;
static tBoolean          sbBrush3D            = FALSE;
static tVolumeValue      snBrushThresholdLow  = knMinVolumeValue;
static tVolumeValue      snBrushThresholdHigh = knMaxVolumeValue;
static tVolumeValue      snBrushNewValue      = knMaxVolumeValue/2;

/* static focused display. */
static tkmDisplayAreaRef sFocusedDisplay = NULL;

char DspA_ksaErrorStrings [DspA_knNumErrorCodes][256] = {
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
  "Couldn't access control points.",
  "Coulnd't access selection.",
  "Couldn't access parent window.",
  "Couldn't access functional volume.",
  "Couldn't access surface cache list.",
  "Invalid error code."
};

char DspA_ksaDisplayFlag [DspA_knNumDisplayFlags][256] = {

  "None",
  "AuxVolume",
  "Cursor",
  "MainSurface",
  "OriginalSurface",
  "CanonicalSurface",
  "InterpolateSurfaceVertices",
  "DisplaySurfaceVertices",
  "ControlPoints",
  "Selection",
  "FunctionalOverlay"
};

char DspA_ksaOrientation [mri_knNumOrientations][256] = {
  "Coronal",
  "Horizontal",
  "Sagittal"
};

char DspA_ksaSurface [Surf_knNumVertexSets][256] = {
  "Main",
  "Original",
  "Canonical"
};

DspA_tErr DspA_New ( tkmDisplayAreaRef* oppWindow,
         tkmMeditWindowRef  ipWindow ) {

  DspA_tErr         eResult      = DspA_tErr_NoErr;
  tkmDisplayAreaRef this         = NULL;
  int               nFlag        = 0;

  /* allocate us. */
  this = (tkmDisplayAreaRef) malloc ( sizeof(tkmDisplayArea) );
  if ( NULL == this ) {
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
  xVoxl_New ( &this->mpCursor );
  xVoxl_New ( &this->mpZoomCenter );

  /* stuff in default values for display states. */
  this->mOrientation           = mri_tOrientation_Coronal;
  this->mnZoomLevel            = 1;
  this->mnHilitedVertexIndex   = -1;
  this->mHilitedSurface        = Surf_tVertexSet_None;
  this->mbSliceChanged         = TRUE;
  this->mnVolumeSize           = 0;

  /* all our display flags start out false. */
  for ( nFlag = 0; nFlag < DspA_knNumDisplayFlags; nFlag++ )
    this->mabDisplayFlags[nFlag] = FALSE;

  /* null ptrs for display data. */
  this->mpVolume                = NULL;
  this->mpAuxVolume             = NULL;
  this->mpParcellationVolume    = NULL;
  this->mpSurface               = NULL;
  this->mpFunctionalVolume      = NULL;
  this->mpControlPoints         = NULL;
  this->mpSelectedControlPoints = NULL;
  this->mpSelection             = NULL;
  this->maSurfaceLists          = NULL;

  /* return window. */
  *oppWindow = this;

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_New: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* delete frame buffer */
  if( NULL != this->mpFrameBuffer )
    free( this->mpFrameBuffer );

  /* delete our voxels */
  xVoxl_Delete ( &this->mpCursor );
  xVoxl_Delete ( &this->mpZoomCenter );

  /* trash the signature */
  this->mSignature = 0x1;
 
  /* delete us */
  free ( this );
  
  /* return null */
  *ioppWindow = NULL;

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_Delete: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set our location */
  this->mLocationInSuper = iLocation;

  /* set our size */
  this->mnWidth  = inWidth;
  this->mnHeight = inHeight;

  /* set our scale */
  if ( this->mnVolumeSize > 0 ) {

    this->mfFrameBufferScaleX = 
      (float)this->mnWidth  / (float)this->mnVolumeSize;
    this->mfFrameBufferScaleY = 
      (float)this->mnHeight / (float)this->mnVolumeSize;
  }

  /* rebuild our current frame and redraw. */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetLocationInSuper: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_UpdateWindowTitle ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  char      sTitle[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make a window title. */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
    sprintf( sTitle, "Tkmedit: %s %s",
       tkm_GetSubjectName(), tkm_GetAuxVolumeName() );
  } else {
    sprintf( sTitle, "Tkmedit: %s %s",
       tkm_GetSubjectName(), tkm_GetVolumeName() );
  }
  
  /* set the window name. */
  MWin_SetWindowTitle( this->mpWindow, sTitle );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_UpdateWindowTitle: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetVolume ( tkmDisplayAreaRef this,
         tVolumeRef        ipVolume,
         int               inSize ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
 xVoxelRef  pCenter            = NULL;
  char      sTclArguments[256] = "";

  xVoxl_New( &pCenter );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the main volume */
  this->mpVolume = ipVolume;

  /* save the volume size */
  this->mnVolumeSize = inSize;

  /* if we alreayd have a frame buffer, delete it */
  if ( NULL == this->mpFrameBuffer ) {
    free ( this->mpFrameBuffer );
    this->mpFrameBuffer = NULL;
  }

  /* allocate a new one */
  this->mpFrameBuffer = (GLubyte*) malloc ( this->mnVolumeSize *
              this->mnVolumeSize * 
              DspA_knNumBytesPerPixel );
  if( NULL == this->mpFrameBuffer ) {
    eResult = DspA_tErr_AllocationFailed;
    goto error;
  }

  /* set inital values for our buffer scale */
  this->mfFrameBufferScaleX = 
    (float)this->mnWidth  / (float)this->mnVolumeSize;
  this->mfFrameBufferScaleY = 
    (float)this->mnHeight  / (float)this->mnVolumeSize;

  /* initialize surface point lists. */
  eResult = DspA_InitSurfaceLists_( this, 
      this->mnVolumeSize * Surf_knNumVertexSets * mri_knNumOrientations );
  if( DspA_tErr_NoErr != eResult )
    goto error;

  /* send volume name to tk window */
  sprintf( sTclArguments, "\"%s value\"", tkm_GetVolumeName() );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );
   
  /* get the center of the volume */
  xVoxl_Set( pCenter, this->mnVolumeSize/2, 
       this->mnVolumeSize/2, this->mnVolumeSize/2 );

  /* set cursor and zoom center to middle of volume */
  eResult = DspA_SetCursor( this, pCenter );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  eResult = DspA_SetZoomCenter( this, pCenter );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* show cursor. */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_Cursor, TRUE );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set dirty flag and redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetVolume: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete( &pCenter );

  return eResult;
}

DspA_tErr DspA_SetAuxVolume ( tkmDisplayAreaRef this,
            tVolumeRef        ipVolume,
            int               inSize ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make sure this size is the same as the main volume size. */
  if( inSize != this->mnVolumeSize ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* save the aux volume */
  this->mpAuxVolume = ipVolume;

  /* if we got a volume... */
  if( NULL != this->mpAuxVolume ) {
    
    /* send volume name to tk window */
    sprintf( sTclArguments, "\"%s value\"", tkm_GetAuxVolumeName() );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeName, sTclArguments );

    /* show the volume value */
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxValue, "1" );

  } else {

    /* hide the volume value */
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxValue, "0" );
  }    

  /* need to redraw if we're currently showing the aux volume */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetAuxVolume: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetParcellationVolume ( tkmDisplayAreaRef this,
               tVolumeRef        ipVolume,
               int               inSize ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  tBoolean  bHaveVolume        = FALSE;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the parcellation volume */
  this->mpParcellationVolume = ipVolume;

  /* turn stuff on or off based on if we have one. */
  if( this->mpParcellationVolume != NULL ) {
    bHaveVolume = TRUE;
  }

  /* turn parcellation on */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_ParcellationOverlay,
         bHaveVolume );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* show parc label */
  sprintf( sTclArguments, "%d", (int)bHaveVolume );
  tkm_SendTclCommand( tkm_tTclCommand_ShowParcellationLabel, sTclArguments );

  /* redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetAuxVolume: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetSurface ( tkmDisplayAreaRef this, 
          mriSurfaceRef     ipSurface ) {

  DspA_tErr        eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the surface */
  this->mpSurface = ipSurface;

  /* purge all the surface lists. */
  DspA_PurgeSurfaceLists_( this );

  /* set slice dirty flag and redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* set surface to null */
  this->mpSurface = NULL;

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetSurface: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
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

  /* turn functional data on if there is */
  if( bOverlayLoaded ) {
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_FunctionalOverlay,
           TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetOverlayVolume: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetControlPointsSpace ( tkmDisplayAreaRef this,
               x3DListRef        ipVoxels ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the control points space. */
  this->mpControlPoints = ipVoxels;

  /* turn control points on. */
  if( NULL != this->mpControlPoints ) {
    
    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_ControlPoints, 
           TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }    

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetControlPointsSpace: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetControlPointsSelectionList ( tkmDisplayAreaRef this,
                 xListRef          ipVoxels ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the list of selected control points */
  this->mpSelectedControlPoints = ipVoxels;

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetControlPointsSelectionList: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetSelectionSpace ( tkmDisplayAreaRef this, 
           x3DListRef        ipVoxels ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the selected voxels */
  this->mpSelection = ipVoxels;

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
    DebugPrint "Error %d in DspA_SetSelectionSpace: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetCursor ( tkmDisplayAreaRef this, 
        xVoxelRef          ipCursor ) {

  DspA_tErr             eResult            = DspA_tErr_NoErr;
 xVoxelRef              pVoxel             = NULL;
  int                   nSlice             = 0;
  unsigned char         ucVolumeValue      = 0;
  FunV_tErr             eFunctional        = FunV_tErr_NoError;
  FunV_tFunctionalValue functionalValue   = 0;
  char                  sTclArguments[256] = "";
 
  xVoxl_New( &pVoxel );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the cursor */
  eResult = DspA_VerifyVolumexVoxl_ ( this, ipCursor );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* get our current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* set the cursor */
  xVoxl_Copy( this->mpCursor, ipCursor );
  
  /* if the new slice number is diffrent, set our dirty slice flag. */
  if( DspA_GetCurrentSliceNumber_( this ) != nSlice ) {
    this->mbSliceChanged = TRUE;
  }

  /* if cursor is .0 .0 .0, change to .5 .5 .5 so that it will draw in the
     center of a voxel on screen. */
  if( xVoxl_GetFloatX(this->mpCursor) == 
      (float)(int)xVoxl_GetFloatX(this->mpCursor)
      && xVoxl_GetFloatY(this->mpCursor) == 
      (float)(int)xVoxl_GetFloatY(this->mpCursor)
      && xVoxl_GetFloatZ(this->mpCursor) == 
      (float)(int)xVoxl_GetFloatZ(this->mpCursor) ) {
    
    xVoxl_SetFloat( this->mpCursor,
        xVoxl_GetFloatX(this->mpCursor) + 0.5,
        xVoxl_GetFloatY(this->mpCursor) + 0.5,
        xVoxl_GetFloatZ(this->mpCursor) + 0.5 );
  }

  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* notify the window that the cursor has changed. */
    MWin_CursorChanged( this->mpWindow, this, this->mpCursor );

    /* send the cursor. */
    sprintf( sTclArguments, "%d %d %d",
       xVoxl_ExpandInt( this->mpCursor ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeCursor, sTclArguments );
    
    /* also convert to RAS and send those coords along. */
    tkm_ConvertVolumeToRAS( this->mpCursor, pVoxel );
    sprintf( sTclArguments, "%.1f %.1f %.1f", xVoxl_ExpandFloat( pVoxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateRASCursor, sTclArguments );
    
    /* also convert to Tal and send those coords along. */
    tkm_ConvertVolumeToTal( this->mpCursor, pVoxel );
    sprintf( sTclArguments, "%.1f %.1f %.1f", xVoxl_ExpandFloat( pVoxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateTalCursor, sTclArguments );
    
    /* also get the volume value and send that along. */
    ucVolumeValue = tkm_GetVolumeValue( this->mpVolume, this->mpCursor );
    sprintf( sTclArguments, "%d", ucVolumeValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValue, sTclArguments );
    
    /* and the aux volume value if there is one. */
    if( NULL != this->mpAuxVolume ) {
      ucVolumeValue = tkm_GetVolumeValue( this->mpAuxVolume, this->mpCursor );
      sprintf( sTclArguments, "%d", ucVolumeValue );
      tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeValue, 
        sTclArguments );
    }

    /* and the parcellation label if it's on */
    if( this->mabDisplayFlags[DspA_tDisplayFlag_ParcellationOverlay] 
  && NULL != this->mpParcellationVolume ) {
      tkm_GetParcellationLabel( this->mpCursor, sTclArguments );
      tkm_SendTclCommand( tkm_tTclCommand_UpdateParcellationLabel, 
        sTclArguments );
    }

    /* also see if we have functional data and can send a value for that
       as well. */     
    if ( NULL != this->mpFunctionalVolume ) {
      DisableDebuggingOutput;
      eFunctional = FunV_GetValueAtAnaIdx( this->mpFunctionalVolume,
             this->mpCursor, 
             &functionalValue );
      EnableDebuggingOutput;
      if( FunV_tErr_NoError == eFunctional ) {
  sprintf( sTclArguments, "%f", functionalValue );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateFunctionalValue, sTclArguments );
      }   
    }
  }

  /* schedule a redraw */
  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetCursor(%d,%d,%d): %s\n",
      eResult, xVoxl_ExpandInt(ipCursor),
      DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete( &pVoxel );

  return eResult;
}

DspA_tErr DspA_SetOrientation ( tkmDisplayAreaRef this, 
        mri_tOrientation  iOrientation ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the orientation */
  if ( iOrientation <= mri_tOrientation_None
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
    
    /* send the orientation */
    sprintf ( sTclArguments, "%d", (int)this->mOrientation );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateOrientation, sTclArguments );
  }

  /* schedule a redraw */
  DspA_Redraw_ ( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetOrientation: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetZoomLevel ( tkmDisplayAreaRef this,
            int               inLevel ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";
 xVoxelRef  pCenter            = NULL;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the zoom level. */
  if ( inLevel < DspA_knMinZoomLevel 
       || inLevel > DspA_knMaxZoomLevel ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* set our zoom level. */
  this->mnZoomLevel = inLevel;

  /* if this is zoom level one, set our center to the middle of
     this slice. */
  if( 1 == this->mnZoomLevel ) {

    /* set a voxel to 128,128,128, setting the zoom center will not change
       the current slice. */
    xVoxl_New( &pCenter );
    xVoxl_Set( pCenter, this->mnVolumeSize/2, 
         this->mnVolumeSize/2, this->mnVolumeSize/2 );
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
    sprintf ( sTclArguments, "%d", (int)this->mnZoomLevel );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateZoomLevel, sTclArguments );
  }

  /* schedule a redraw */
  DspA_Redraw_ ( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetZoomLevel: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the center */
  eResult = DspA_VerifyVolumexVoxl_ ( this, ipCenter );
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
  DspA_Redraw_ ( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetZoomCenter: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetZoomCenterToCursor ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the center to our cursor */
  DspA_SetZoomCenter( this, this->mpCursor );

  /* schedule a redraw */
  DspA_Redraw_ ( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetZoomCenterToCursor: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_HiliteSurfaceVertex ( tkmDisplayAreaRef this,
             Surf_tVertexSet  inSurface,
             int               inVertex ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the values */
  this->mHilitedSurface        = inSurface;
  this->mnHilitedVertexIndex   = inVertex;

  /* schedule a redraw */
  DspA_Redraw_ ( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_HiliteSurfaceVertex: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_SetDisplayFlag ( tkmDisplayAreaRef this,
        DspA_tDisplayFlag iWhichFlag,
        tBoolean          ibNewValue ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";
  tBoolean  bNewValue          = FALSE;

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

    break;

  case DspA_tDisplayFlag_ParcellationOverlay:

    /* if no parcellation volume, set to false. */
    if( NULL == this->mpParcellationVolume )
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
    if( NULL == this->mpSurface ) 
      bNewValue = FALSE;

    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;
    
  case DspA_tDisplayFlag_InterpolateSurfaceVertices:

    /* if no surface, set to false. */
    if( NULL == this->mpSurface ) 
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

    /* if no func data, set to false. */
    if( NULL == this->mpFunctionalVolume ) {
      DebugPrint "DspA_SetDisplayFlag ( DspA_tDisplayFlag_FunctionalOverlay, true ): func data wasn't loaded.\n" EndDebugPrint;
      bNewValue = FALSE;
    }

    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

  case DspA_tDisplayFlag_MaxIntProj:

    /* if the flag is different, set dirty flag */
    if( this->mabDisplayFlags[iWhichFlag] != bNewValue )
      this->mbSliceChanged = TRUE;

    break;

  default:
    break;
  }

  /* set the value */
  this->mabDisplayFlags[iWhichFlag] = bNewValue;

  /* if it was the aux volume flag, change the window title. */
  if( DspA_tDisplayFlag_AuxVolume == iWhichFlag ) {

    DspA_UpdateWindowTitle( this );
  }

  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* notify the window that the flag has changed. */
    MWin_DisplayFlagChanged( this->mpWindow, this, iWhichFlag,
           this->mabDisplayFlags[iWhichFlag] );

    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d", (int)iWhichFlag, (int)bNewValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateDisplayFlag, sTclArguments );
  }

  /* schedule a redraw */
  DspA_Redraw_ ( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_ToggleDisplayFlag: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetTool ( tkmDisplayAreaRef this,
       DspA_tTool        iTool ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";

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
    DebugPrint "Error %d in DspA_SetTool: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetBrush ( tkmDisplayAreaRef this,
        int               inRadius,
        DspA_tBrushShape  iShape,
        tBoolean          ib3D ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* check the radius */
  if ( inRadius < DspA_knMinBrushRadius 
       || inRadius > DspA_knMaxBrushRadius ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* Set the brush info */
  snBrushRadius = inRadius;
  sBrushShape   = iShape;
  sbBrush3D     = ib3D;

  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d",
        (int)snBrushRadius, (int)sBrushShape, 
        (int)sbBrush3D );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrush, sTclArguments );
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetBrush: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SetBrushThreshold ( tkmDisplayAreaRef this,
           tVolumeValue      inLow,
           tVolumeValue      inHigh,
           tVolumeValue      inNewValue ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* set the brush theshold info */
  snBrushThresholdLow  = inLow;
  snBrushThresholdHigh = inHigh;
  snBrushNewValue      = inNewValue;

  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* send the tcl update. */
    sprintf ( sTclArguments, "%d %d %d",
        (int)snBrushThresholdLow, (int)snBrushThresholdHigh, 
        (int)snBrushNewValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushThreshold, sTclArguments );
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_SetBrushThreshold: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_ChangeSliceBy_ ( tkmDisplayAreaRef this,
         int               inDelta ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
 xVoxelRef  pCursor = NULL;

  xVoxl_New( &pCursor );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* copy the cursor. */
  xVoxl_Copy( pCursor, this->mpCursor );

  /* apply the increment to the proper part of the cursor */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_SetZ( pCursor, xVoxl_GetZ( pCursor ) + inDelta );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_SetY( pCursor, xVoxl_GetY( pCursor ) + inDelta );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetX( pCursor, xVoxl_GetX( pCursor ) + inDelta );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
  }

  /* set the cursor. */
  eResult = DspA_SetCursor( this, pCursor );
   if ( DspA_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_ChangeSliceBy_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_Focus: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_Blur: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_GetPosition: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

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

    /* schedule a redraw - right now just call the draw function */
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
    DebugPrint "Error %d in DspA_HandleEvent: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_HandleMouseUp_ ( tkmDisplayAreaRef this, 
        xGWin_tEventRef   ipEvent ) {

  DspA_tErr   eResult     = DspA_tErr_NoErr;
  xPoint2n    bufferPt    = {0,0};
 xVoxelRef    pVolumeVox  = NULL;
  FunV_tErr   eFunctional = FunV_tErr_NoError;

  xVoxl_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  DebugPrint "Mouse up screen x %d y %d buffer x %d y %d volume %d %d %d\n",
    ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
    xVoxl_ExpandInt( pVolumeVox ) EndDebugPrint;

  /* set the cursor. */
  eResult = DspA_SetCursor( this, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* allow the functional display to respond. */
  if( NULL != this->mpFunctionalVolume ) {
    eFunctional = FunV_AnatomicalVoxelClicked( this->mpFunctionalVolume,
                 pVolumeVox );
    if( FunV_tErr_NoError != eFunctional ) {
      DebugPrint "DspA_HandleMouseUp_(): Error while passing clicked voxel to functional volume.\n" EndDebugPrint;
    }
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
  
  /* if we're in control point selection mode... */
  if( DspA_tTool_SelectCtrlPts == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
      
    /* if button 2, do a ctrl pt select. */
    if ( 2 == ipEvent->mButton ) {

      /* if shift is down... */
      if ( ipEvent->mbShiftKey ) {
  
  /* add nearest ctrl pt to selection */
  tkm_AddNearestCtrlPtToSelection( pVolumeVox, this->mOrientation );

  /* no shift key...*/
      } else {

  /* clear selection and select nearest. */
  tkm_DeselectAllCtrlPts();
  tkm_AddNearestCtrlPtToSelection( pVolumeVox, this->mOrientation );
      }

      /* if button 3, do a ctrl pt unselect. */
    } else if ( 3 == ipEvent->mButton ) {

      /* remove nearest point from selection */
      tkm_RemoveNearestCtrlPtFromSelection( pVolumeVox, this->mOrientation );
    }
  }

  /* most things want a mouse up after we're done, so schedule one and 
     let the draw function figure out if we don't need it. */
  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_HandleMouseUp_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleMouseDown_ ( tkmDisplayAreaRef this, 
          xGWin_tEventRef   ipEvent ) {

  DspA_tErr   eResult = DspA_tErr_NoErr;
  xPoint2n    bufferPt    = {0,0};
 xVoxelRef    pVolumeVox  = NULL;

  xVoxl_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
#if 0
  DebugPrint "Mouse down screen x %d y %d buffer x %d y %d volume %d %d %d\n",
    ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
    xVoxl_ExpandInt( pVolumeVox ) EndDebugPrint;
#endif

  /* edit mode with no modifiers: */
  if( DspA_tTool_EditVoxels == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* clear the undo list. */
    tkm_ClearUndoList();
    
    /* if button 2, do a brush to white. */
    if( 2 == ipEvent->mButton ) {
      
      /* set our brush up for this to white stuff */
      DspA_SetBrushThreshold( this, tkm_knEditToWhiteLow, 
            tkm_knEditToWhiteHigh, 
            tkm_knEditToWhiteNewValue );

      eResult = DspA_BrushVoxels_( this, pVolumeVox, 
           DspA_BrushVoxelsInThreshold_ );
      if( DspA_tErr_NoErr != eResult )
  goto error;
    
      /* editing requires us to rebuild buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );

    /* if button 3, do a brush to black. */
    } else if ( 3 == ipEvent->mButton ) {

      /* set our brush up for this to black stuff */
      DspA_SetBrushThreshold( this, tkm_knEditToBlackLow, 
            tkm_knEditToBlackHigh, 
            tkm_knEditToBlackNewValue );

      eResult = DspA_BrushVoxels_( this, pVolumeVox, 
           DspA_BrushVoxelsInThreshold_ );
      if( DspA_tErr_NoErr != eResult )
  goto error;

      /* editing requires us to rebuild buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );
    }
  }

  /* custom edit mode with no modifiers */
  if( DspA_tTool_CustomEditVoxels == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey
      && (2 == ipEvent->mButton || 3 == ipEvent->mButton) ) {

    /* clear the undo list. */
    tkm_ClearUndoList();
    
    /* brush the voxels with the current threshold setting */
    eResult = DspA_BrushVoxels_( this, pVolumeVox, 
         DspA_BrushVoxelsInThreshold_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* editing requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  /* select mode with no modfiers */
  if( DspA_tTool_SelectVoxels == sTool
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* if button 2, do a brush select. */
    if ( 2 == ipEvent->mButton ) {
      eResult = DspA_BrushVoxels_( this, pVolumeVox, tkm_SelectVoxel );
      if( DspA_tErr_NoErr != eResult )
  goto error;
      
      /* selecting requires us to rebuild buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );

      /* if button 3, do a brush unselect. */
    } else if ( 3 == ipEvent->mButton ) {
      eResult = DspA_BrushVoxels_( this, pVolumeVox, tkm_DeselectVoxel );
      if( DspA_tErr_NoErr != eResult )
  goto error;

      /* selecting requires us to rebuild buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_HandleMouseDown_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleMouseMoved_ ( tkmDisplayAreaRef this, 
           xGWin_tEventRef   ipEvent ) {
  DspA_tErr   eResult = DspA_tErr_NoErr;
  xPoint2n    bufferPt    = {0,0};
 xVoxelRef    pVolumeVox  = NULL;

  xVoxl_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

#if 0
  DebugPrint "Mouse moved screen x %d y %d buffer x %d y %d volume %d %d %d\n",
    ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
    xVoxl_ExpandInt( pVolumeVox ) EndDebugPrint;
#endif

  /* edit mode with no modifiers: */
  if( DspA_tTool_EditVoxels == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* if button 2, do a brush to white. */
    if ( 2 == ipEvent->mButton ) {
      eResult = DspA_BrushVoxels_( this, pVolumeVox, 
           DspA_BrushVoxelsInThreshold_ );
      if( DspA_tErr_NoErr != eResult )
  goto error;
      
      /* editing requires us to rebuild the buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );
    
      /* if button 3, do a brush to black. */
    } else if ( 3 == ipEvent->mButton ) {
      eResult = DspA_BrushVoxels_( this, pVolumeVox, 
           DspA_BrushVoxelsInThreshold_ );
      if( DspA_tErr_NoErr != eResult )
  goto error;

      /* editing requires us to rebuild the buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );
    }
  }

  /* custom edit mode with no modifiers and button 2 or 3 */
  if( DspA_tTool_CustomEditVoxels == sTool 
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey
      && (2 == ipEvent->mButton || 3 == ipEvent->mButton) ) {

    /* brush the voxels with the current threshold setting */
    eResult = DspA_BrushVoxels_( this, pVolumeVox, 
         DspA_BrushVoxelsInThreshold_ );
    if( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* editing requires us to rebuild buffer. */
    this->mbSliceChanged = TRUE;
    DspA_Redraw_( this );
  }

  /* select mode with no modfiers */
  if( DspA_tTool_SelectVoxels == sTool
      && !ipEvent->mbCtrlKey
      && !ipEvent->mbAltKey ) {
    
    /* if button 2, do a brush select. */
    if ( 2 == ipEvent->mButton ) {
      eResult = DspA_BrushVoxels_( this, pVolumeVox, tkm_SelectVoxel );
      if( DspA_tErr_NoErr != eResult )
  goto error;
      
      /* selecting requires us to rebuild the buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );
      
      /* if button 3, do a brush unselect. */
    } else if ( 3 == ipEvent->mButton ) {
      eResult = DspA_BrushVoxels_( this, pVolumeVox, tkm_DeselectVoxel );
      if( DspA_tErr_NoErr != eResult )
  goto error;
      
      /* selecting requires us to rebuild the buffer. */
      this->mbSliceChanged = TRUE;
      DspA_Redraw_( this );
    }
  }
  
  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_HandleMouseMoved_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    /* ctrl 1 goes to main volume */
    if ( ipEvent->mbCtrlKey ) {
      eResult = DspA_SetDisplayFlag( this, 
             DspA_tDisplayFlag_AuxVolume, FALSE );
      if ( DspA_tErr_NoErr != eResult )
  goto error;
    }
    break;

  case '2':
    /* ctrl 2 goes to aux volume */
    if ( ipEvent->mbCtrlKey ) {
      eResult = DspA_SetDisplayFlag( this, 
             DspA_tDisplayFlag_AuxVolume, TRUE );
      if ( DspA_tErr_NoErr != eResult )
  goto error;
    }
    break;

  case 'c':
    /* alt-c toggles main and aux view. */
    if( ipEvent->mbAltKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_AuxVolume );
      if ( DspA_tErr_NoErr != eResult )
  goto error;
    }
    /* ctrl+c toggles cursor display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_Cursor );
      if ( DspA_tErr_NoErr != eResult )
  goto error;
    }
    break;

  case 'f':
    /* ctrl+f toggles functional overlay display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, 
          DspA_tDisplayFlag_FunctionalOverlay );
      if ( DspA_tErr_NoErr != eResult )
  goto error;
    }
    break;

  case 'h':

    FunV_UseOverlayCache( this->mpFunctionalVolume, 
        !this->mpFunctionalVolume->mbUseOverlayCache );
    if( this->mpFunctionalVolume->mbUseOverlayCache ) {
      DebugPrint "Cache enabled.\n" EndDebugPrint;
    } else {
      DebugPrint "Cache disabled.\n" EndDebugPrint;
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
    }
    break;

  case 's':
    /* ctrl+s toggles main suface display */
    if( ipEvent->mbCtrlKey ) {
      eResult = DspA_ToggleDisplayFlag( this, DspA_tDisplayFlag_MainSurface );
      if ( DspA_tErr_NoErr != eResult )
  goto error;
    }
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
    /* z sets plane to coronal */
    eResult = DspA_SetOrientation( this, mri_tOrientation_Coronal );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
    
  case xGWin_tKey_UpArrow:
  case xGWin_tKey_RightArrow:
    /* move up a slice */
    eResult = DspA_ChangeSliceBy_( this, 1 );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    DspA_Redraw_( this );
    break;
    
  case xGWin_tKey_DownArrow:
  case xGWin_tKey_LeftArrow:
    /* move down a slice */
    eResult = DspA_ChangeSliceBy_( this, -1 );
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
    DebugPrint "Error %d in DspA_HandleKeyDown_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_BrushVoxels_ ( tkmDisplayAreaRef this,
           xVoxelRef          ipCenterVox,
            void(*ipFunction)(xVoxelRef) ) {

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
  nXRadius = snBrushRadius - 1;
  nYRadius = snBrushRadius - 1;
  nZRadius = snBrushRadius - 1;

  /* if we're not in 3d, set the same radius of the same plane as our
     current orientation to 0. */
  if( !sbBrush3D ) {
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
      DspA_VerifyVolumexVoxl_( this, pVolumeVox ) ) {

    /* if we're circular, check the radius. if no good, continue. */
    if( DspA_tBrushShape_Circle == sBrushShape 
        && ( (nXCenter-nX)*(nXCenter-nX) +
       (nYCenter-nY)*(nYCenter-nY) +
       (nZCenter-nZ)*(nZCenter-nZ) >
       (snBrushRadius-1)*(snBrushRadius-1) ) ) {
      continue;
    }

    /* run the function on this voxel. */
    ipFunction( pVolumeVox );
  }
      }
    }
  }

  DspA_Redraw_( this );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_BrushVoxels_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:
  xVoxl_Delete( &pVolumeVox );

  return eResult;

}

void DspA_BrushVoxelsInThreshold_ (xVoxelRef ipVoxel ) {

  tkm_EditVoxelInRange( ipVoxel, 
      snBrushThresholdLow,
      snBrushThresholdHigh,
      snBrushNewValue );
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
    DebugPrint "Error %d in DspA_Redraw_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  }

  /* draw the selection */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_Selection] ) {
    eResult = DspA_DrawSelectionToFrame_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
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

  /* draw the control points */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_ControlPoints] ) {
    eResult = DspA_DrawControlPointsToFrame_( this );
    if ( DspA_tErr_NoErr != eResult ) {
      DspA_Signal( "DspA_HandleDraw_", __LINE__, eResult );
      eResult = DspA_tErr_NoErr;
    }
  }

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
    DebugPrint "Error %d in DspA_HandleDraw_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_DrawFrameBuffer_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  DspA_SetUpOpenGLPort_( this );

  glDrawPixels ( this->mnVolumeSize, this->mnVolumeSize,
     GL_RGBA, GL_UNSIGNED_BYTE, this->mpFrameBuffer );
  
  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawFrameBuffer_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_DrawSurface_ ( tkmDisplayAreaRef this ) {

  DspA_tErr         eResult    = DspA_tErr_NoErr;
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

  /* only recalculate the draw list of the plane and orientation has
     changed. */
  if ( NULL != this->mpSurface ) {
    
    /* if the orienation or slice has changed... */
    if( this->mbSliceChanged ) {
      
      /* rebuild the draw lists. */
      eResult = DspA_BuildSurfaceDrawLists_( this );

      /* if we got an out of memory error.. */
      if ( DspA_tErr_OutOfMemory == eResult ) {
  
  /* delete draw lists and try again. */
  DspA_PurgeSurfaceLists_( this );
  eResult = DspA_BuildSurfaceDrawLists_ ( this );
  
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
  list = DspA_GetSurfaceList_( this, this->mOrientation, surface,
             DspA_GetCurrentSliceNumber_(this) );
  if( NULL == list ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }
  
  /* choose and set the color */
  switch( surface ) {
  case Surf_tVertexSet_Main:
    faColor[0] = 1.0; faColor[1] = 1.0; faColor[2] = 0.0;
    break;
  case Surf_tVertexSet_Original:
    faColor[0] = 0.0; faColor[1] = 1.0; faColor[2] = 0.0;
    break;
  case Surf_tVertexSet_Pial:
    faColor[0] = 1.0; faColor[1] = 0.0; faColor[2] = 0.0;
    break;
  default:
    faColor[0] = 0.5; faColor[1] = 0.5; faColor[2] = 0.5;
    break;
  }

  glColor3fv( faColor );
  
  /* draw the points. */
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
  
  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawSurface_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_DrawCursor_ ( tkmDisplayAreaRef this ) {

  DspA_tErr eResult      = DspA_tErr_NoErr;
  xPoint2n  bufferPt     = {0, 0};
  int       nWidth       = 0;
  int       nHeight      = 0;

  /* convert the voxel and see if it's on the screen. the screen is 
     actually sized to the buffer size, i.e. volumeSize x volumeSize, 
     so only convert to buffer space.*/
  eResult = DspA_ConvertVolumeToBuffer_ ( this, this->mpCursor, &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  bufferPt.mnY = GLDRAW_Y_FLIP(bufferPt.mnY);

  /* calculate width and height using scale */
  nWidth  = ((float) DspA_knCursorCrosshairSize / this->mfFrameBufferScaleX);
  nHeight = ((float) DspA_knCursorCrosshairSize / this->mfFrameBufferScaleY);

  DspA_SetUpOpenGLPort_( this );

  /* draw the cursor */
  glColor3f ( 1.0, 0.0, 0.0 );
  
  glBegin ( GL_LINES );
  glVertex2d ( bufferPt.mnX, bufferPt.mnY-nHeight );
  glVertex2d ( bufferPt.mnX, bufferPt.mnY+nHeight );
  glEnd ();
  
  glBegin ( GL_LINES );
  glVertex2d ( bufferPt.mnX-nWidth, bufferPt.mnY );
  glVertex2d ( bufferPt.mnX+nWidth, bufferPt.mnY );
  glEnd ();

  glBegin( GL_LINES );
  glVertex2d( bufferPt.mnX, 0 );
  glVertex2d( bufferPt.mnX, 3*nHeight );
  glEnd();

  glBegin( GL_LINES );
  glVertex2d( bufferPt.mnX, this->mnVolumeSize );
  glVertex2d( bufferPt.mnX, this->mnVolumeSize - 3*nHeight );
  glEnd();

  glBegin( GL_LINES );
  glVertex2d( 0, bufferPt.mnY );
  glVertex2d( 3*nWidth, bufferPt.mnY );
  glEnd();

  glBegin( GL_LINES );
  glVertex2d( this->mnVolumeSize, bufferPt.mnY );
  glVertex2d( this->mnVolumeSize - 3*nWidth, bufferPt.mnY );
  glEnd();

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawCursor_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  glVertex2d( this->mnVolumeSize-1, 1 );
  glVertex2d( this->mnVolumeSize-1, this->mnVolumeSize-1 );
  glVertex2d( 1, this->mnVolumeSize-1 );
  glVertex2d( 1, 1 );
  glEnd ();
  
  goto cleanup;

  goto error;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawFrameAroundDisplay_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  char      sLabel[256] = "";
  int       nLength     = 0;
  int       nChar       = 0;

  xVoxl_New( &pVoxel );

  DspA_SetUpOpenGLPort_( this );

  glColor3f ( 0.0, 1.0, 0.0 );

  volumePt.mnX = 0;

  /* all down the side, every size / 5 pixels */
  for ( nY = 0; nY < this->mnVolumeSize; nY += this->mnVolumeSize/5 ) {
      
    /* y flip the volume pt to flip the image over. */
    //volumePt.mnY = BUFFER_Y_FLIP(nY);
    volumePt.mnY = nY;

    /* get a volume voxel.*/
    eResult = DspA_ConvertBufferToVolume_ ( this, &volumePt, pVoxel );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* if this is good... */
    eResult = DspA_VerifyVolumexVoxl_( this, pVoxel );
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

  /* draw arrows */
  switch( this->mOrientation ) {

  case mri_tOrientation_Coronal:
    volumePt.mnX = 40; 
    volumePt.mnY = this->mnVolumeSize - 20;
    DspA_DrawHorizontalArrow_( &volumePt, -20, "R" );
    volumePt.mnX = this->mnVolumeSize - 20;
    volumePt.mnY = 40; 
    DspA_DrawVerticalArrow_( &volumePt, -20, "S" );
    break;
  case mri_tOrientation_Horizontal:
    volumePt.mnX = 40; 
    volumePt.mnY = this->mnVolumeSize - 20;
    DspA_DrawHorizontalArrow_( &volumePt, -20, "R" );
    volumePt.mnX = this->mnVolumeSize - 20;
    volumePt.mnY = 40; 
    DspA_DrawVerticalArrow_( &volumePt, -20, "A" );
    break;
  case mri_tOrientation_Sagittal:
   volumePt.mnX = this->mnVolumeSize - 40; 
    volumePt.mnY = this->mnVolumeSize - 20;
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
    DebugPrint "Error %d in DspA_DrawAxes_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  unsigned char         ucValue     = 0;
 xVoxelRef              pVoxel      = NULL;
 xVoxelRef              pFuncMin    = NULL;
 xVoxelRef              pFuncMax    = NULL;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  FunV_tFunctionalValue funcValue   = 0.0;
  xColor3f              color       = {0,0,0};
  xColor3f              funcColor   = {0,0,0};
  tBoolean              bPixelSet   = FALSE;
  int                   nY          = 0;

  xUtil_StartTimer();

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
    xVoxl_Set( pVoxel, this->mnVolumeSize/2, 
         this->mnVolumeSize/2, this->mnVolumeSize/2 );
    eResult = DspA_SetZoomCenter( this, pVoxel );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;

  DisableDebuggingOutput;

  /* go thru the buffer... */
  for ( nY = 0; nY < this->mnVolumeSize; nY ++ ) {
    for ( volumePt.mnX = 0; 
    volumePt.mnX < this->mnVolumeSize; volumePt.mnX ++ ) {

      /* haven't set this pixel yet. */
      bPixelSet = FALSE;

      /* y flip the volume pt to flip the image over. */
      volumePt.mnY = BUFFER_Y_FLIP(nY);

      /* get a volume voxel.*/
      eResult = DspA_ConvertBufferToVolume_ ( this, &volumePt, pVoxel );
      if ( DspA_tErr_NoErr != eResult )
  goto error;

      /* check it. */
      eResult = DspA_VerifyVolumexVoxl_( this, pVoxel );
      if( DspA_tErr_NoErr == eResult ) {
  
  /* if we are showing parcellation... */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_ParcellationOverlay] ) {
    
    /* get parcellation color. */
    tkm_GetParcellationColor( pVoxel, &color );

    /* if it's not zero, pixel is set. */
    if ( color.mfRed != 0 || color.mfBlue != 0 || color.mfGreen != 0 ) {
      bPixelSet = TRUE;
    }

  }

  /* if we didn't set the pixel color from the parcellation... */
  if ( !bPixelSet ) {
    
    /* get the plain anatomical value from the main or aux
       volume. */
    if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      if( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] )
        ucValue = tkm_GetMaxIntProjValue( this->mpAuxVolume,
            this->mOrientation,
            pVoxel );
        else
    ucValue = tkm_GetVolumeValue ( this->mpAuxVolume, pVoxel );
    } else {
      if( this->mabDisplayFlags[DspA_tDisplayFlag_MaxIntProj] )
        ucValue = tkm_GetMaxIntProjValue( this->mpVolume,
            this->mOrientation,
            pVoxel );
      else
        ucValue = tkm_GetVolumeValue ( this->mpVolume, pVoxel );
    }

    /* get the color */
    tkm_GetAnatomicalVolumeColor( ucValue, &color );
  }


  /* if we have and are showing functional data... */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] 
      && NULL != this->mpFunctionalVolume ) {
           
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
    DebugPrint "Error %d in DspA_BuildCurrentFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  EnableDebuggingOutput;

  /* delete the voxel. */
  xVoxl_Delete ( &pVoxel );
  xVoxl_Delete ( &pFuncMin );
  xVoxl_Delete ( &pFuncMax );

  xUtil_StopTimer( "build frame buffer" );

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
  
  /* draw the color scale bar. get threshold max. */
  eFunctional = FunV_GetThresholdMax( this->mpFunctionalVolume, &max );
  if( FunV_tErr_NoError != eFunctional ) {
    eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
    goto error;
  }

  /* down the buffer */
  for ( bufferPt.mnY = 0; bufferPt.mnY < this->mnVolumeSize; bufferPt.mnY++ ) {
      
    /* get an interpolated value within the range of -max to +max 
       determined by the y value */
    funcValue = (FunV_tFunctionalValue) 
      ( (float)bufferPt.mnY * 
  (float)((max*2.0)/(float)this->mnVolumeSize) - max );
      
      /* get the functional color for this value */
    eFunctional = FunV_GetColorForValue( this->mpFunctionalVolume,
           funcValue, &color, &funcColor );
    if( FunV_tErr_NoError != eFunctional ) {
      eResult = DspA_tErr_ErrorAccessingFunctionalVolume;
      goto error;
    }
  
    /* draw on the right side... */
    for ( bufferPt.mnX = this->mnVolumeSize - 10; 
    bufferPt.mnX < this->mnVolumeSize; bufferPt.mnX++ ) {
      
      /* write it back to the buffer. */
      pFrame = this->mpFrameBuffer + 
  ( (bufferPt.mnY * this->mnVolumeSize) + bufferPt.mnX ) * 
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
  
  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawFunctionalOverlayToFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


DspA_tErr DspA_DrawControlPointsToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr    eResult    = DspA_tErr_NoErr;
  xListRef     list       = NULL;
  xList_tErr   eList      = xList_tErr_NoErr;
  x3Lst_tErr   e3DList    = x3Lst_tErr_NoErr;
  xVoxelRef    controlPt  = NULL;
  tBoolean     bSelected  = FALSE;
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

    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&controlPt )) 
     != xList_tErr_EndOfList ) {

      if( controlPt ) {

  /* see if it's selected. */
  eList = xList_IsInList( this->mpSelectedControlPoints, 
        controlPt, &bSelected );
  if ( xList_tErr_NoErr != eList )
    goto error;

  /* color is green if not selected, yellow if so. */
  if ( bSelected ) {
    faColor[0] = 1;
    faColor[1] = 1;
    faColor[2] = 0;
  } else {
    faColor[0] = 0;
    faColor[1] = 1;
    faColor[2] = 0;
  }      
  
  /* convert to buffer point. */
  eResult = DspA_ConvertVolumeToBuffer_ ( this, controlPt, &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* draw a point here. */
  DspA_DrawCrosshairIntoFrame_( this, faColor, &bufferPt, 
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
    DebugPrint "Error %d in DspA_DrawControlPointsToFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
  int           nValue            = 0;
  int           nIntensifiedValue = 0;
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
  
  /* get the value at this point. */
  nValue = tkm_GetVolumeValue( this->mpVolume, selection );
    
  /* calculate the intensified green component. */
  nIntensifiedValue = nValue + DspA_knSelectionIntensityIncrement;
  if ( nIntensifiedValue > DspA_knMaxPixelValue ) {
    nIntensifiedValue = DspA_knMaxPixelValue;
  }
  
  /* lower the other components if they are too high. */
  if ( nIntensifiedValue - nValue <= 
       DspA_knMinSelectionIntensityDiff ) {  
    nValue -= DspA_knSelectionIntensityIncrement;
    if ( nValue < 0 ) 
      nValue = 0;
  }
  
  /* covert to a buffer point. */
  eResult = DspA_ConvertVolumeToBuffer_ ( this, selection, &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
    
  /* y flip the volume pt to flip the image over. */
  bufferPt.mnY = BUFFER_Y_FLIP(bufferPt.mnY);
  
  /* write it back to the buffer. */
  for( nY = bufferPt.mnY; nY < bufferPt.mnY + this->mnZoomLevel; nY++ ) {
    for( nX = bufferPt.mnX; nX < bufferPt.mnX + this->mnZoomLevel;nX++) {
      
      pFrame = this->mpFrameBuffer + 
        ( (nY * this->mnVolumeSize) + nX ) * DspA_knNumBytesPerPixel;
      pFrame[DspA_knRedPixelCompIndex]   = (GLubyte)nValue;
      pFrame[DspA_knGreenPixelCompIndex] = (GLubyte)nIntensifiedValue;
      pFrame[DspA_knBluePixelCompIndex]  = (GLubyte)nValue;
      pFrame[DspA_knAlphaPixelCompIndex] = (GLubyte)DspA_knMaxPixelValue;
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
    DebugPrint "Error %d in DspA_DrawSelectionToFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_BuildSurfaceDrawLists_ ( tkmDisplayAreaRef this ) {

  DspA_tErr           eResult           = DspA_tErr_NoErr;
  Surf_tErr           eSurface          = Surf_tErr_NoErr;
  xGArr_tErr          eList             = xGArr_tErr_NoErr;
  xGrowableArrayRef   list              = NULL;
  int                 nSlice            = 0;
  Surf_tVertexSet     surface           = Surf_tVertexSet_Main;
  xPoint2n            zeroPoint         = {0,0};
  xVoxelRef           curPlane          = NULL;
  xVoxelRef           anaVertex         = NULL; 
  xVoxelRef           anaNeighborVertex = NULL;
  xPoint2f            intersectionPt    = {0, 0};
  tBoolean            bPointsOnThisFace = FALSE;

  xVoxl_New( &curPlane );
  xVoxl_New( &anaVertex );
  xVoxl_New( &anaNeighborVertex );

  /* get the current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* set up for gl drawing */
  DspA_SetUpOpenGLPort_( this );

  /* for each surface type... */
  for ( surface = Surf_tVertexSet_Main;
  surface < Surf_knNumVertexSets; surface++ ) {
    
    /* only build if this surface is being displayed */
    if( !this->mabDisplayFlags[surface + DspA_tDisplayFlag_MainSurface] ) {
      continue;
    }

    /* check to see if we already have this list. if so, don't build it. */
    list = DspA_GetSurfaceList_( this, this->mOrientation, surface,
         DspA_GetCurrentSliceNumber_(this) );
    if( NULL != list ) {
      continue;
    }
    
    /* make a new list. */
    DspA_NewSurfaceList_( this, this->mOrientation, surface, 
        DspA_GetCurrentSliceNumber_(this) );
    
    /* get the list. */
    list = DspA_GetSurfaceList_( this, this->mOrientation, surface,
         DspA_GetCurrentSliceNumber_(this) );
    if( NULL == list ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }

    xUtil_StartTimer ();

    /* make a voxel with our current slice in it and
       set the iterator to our view */
    DspA_ConvertPlaneToVolume_ ( this, &zeroPoint, nSlice, this->mOrientation,
         curPlane );
    eSurface = Surf_SetIteratorPosition( this->mpSurface, curPlane );
    if( Surf_tErr_NoErr != eSurface )
      goto error;

    bPointsOnThisFace = FALSE;

    /* while we have vertices to check.. */
    while( (eSurface = Surf_GetNextAndNeighborVertex( this->mpSurface,
                  surface,
                  anaVertex, 
                  anaNeighborVertex ))
     != Surf_tErr_LastFace ) {

      /* convert coords to that z is the current plane */
      switch( this->mOrientation ) {
      case mri_tOrientation_Horizontal:
  xVoxl_SetFloat( anaVertex, 
      xVoxl_GetFloatX(anaVertex), 
      xVoxl_GetFloatZ(anaVertex),
      xVoxl_GetFloatY(anaVertex) );
  xVoxl_SetFloat( anaNeighborVertex,
      xVoxl_GetFloatX(anaNeighborVertex),
      xVoxl_GetFloatZ(anaNeighborVertex), 
      xVoxl_GetFloatY(anaNeighborVertex) );
  break;
      case mri_tOrientation_Coronal:
  xVoxl_SetFloat( anaVertex, 
      xVoxl_GetFloatX(anaVertex),
      xVoxl_GetFloatY(anaVertex),
      xVoxl_GetFloatZ(anaVertex) );
  xVoxl_SetFloat( anaNeighborVertex, 
      xVoxl_GetFloatX(anaNeighborVertex), 
      xVoxl_GetFloatY(anaNeighborVertex),
      xVoxl_GetFloatZ(anaNeighborVertex) );
  break;
      case mri_tOrientation_Sagittal:
  xVoxl_SetFloat( anaVertex,
      xVoxl_GetFloatZ(anaVertex),
      xVoxl_GetFloatY(anaVertex),
      xVoxl_GetFloatX(anaVertex) );
  xVoxl_SetFloat( anaNeighborVertex, 
      xVoxl_GetFloatZ(anaNeighborVertex),
      xVoxl_GetFloatY(anaNeighborVertex), 
      xVoxl_GetFloatX(anaNeighborVertex) );
  break;
      default:
  break;
      }

      /* if the line between these two points intersects the
   current plane... */
      if ( xUtil_LineIntersectsPlane( anaVertex, anaNeighborVertex, nSlice, 
   this->mabDisplayFlags[DspA_tDisplayFlag_InterpolateSurfaceVertices],
              &intersectionPt ) ) {
  /* add this point */
  eList = xGArr_Add( list, &intersectionPt );
  if( xGArr_tErr_NoErr != eList ) {
    DebugPrint "xGArr error %d in DspA_BuildSurfaceDrawLists_: %s\n",
      eList, xGArr_GetErrorString( eList ) EndDebugPrint;
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }

  bPointsOnThisFace = TRUE;
      }

      /* if we have a last vertex, and we drew points on this
   face add a face marker. */
      if( eSurface == Surf_tErr_LastVertex
    && bPointsOnThisFace ) {

  intersectionPt.mfX = -1;
  intersectionPt.mfY = -1;
  eList = xGArr_Add( list, &intersectionPt );
  if( xGArr_tErr_NoErr != eList ) {
    DebugPrint "xGArr error %d in DspA_BuildSurfaceDrawLists_: %s\n",
      eList, xGArr_GetErrorString( eList ) EndDebugPrint;
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
    
    xUtil_StopTimer( "build surface list" );

  } /* for each surface... */

  goto cleanup;

 error:

  if( Surf_tErr_NoErr != eSurface )
    eResult = DspA_tErr_ErrorAccessingSurface;
 
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_BuildSurfaceDrawLists_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete( &curPlane );
  xVoxl_Delete( &anaVertex );
  xVoxl_Delete( &anaNeighborVertex );

  return eResult;
  }

DspA_tErr DspA_DrawCrosshairIntoFrame_ ( tkmDisplayAreaRef this,
           float*            ifaColor,
           xPoint2nRef       ipWhere,
           int               inSize ) {
  DspA_tErr eResult = DspA_tErr_NoErr;
  int       nWidth  = 0;
  int       nHeight = 0;
  
  /* calculate width and height using scale */
  nWidth  = ((float)inSize / this->mfFrameBufferScaleX);
  nHeight = ((float)inSize / this->mfFrameBufferScaleY);

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

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DspA_DrawCrosshairIntoFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_GetCursor: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_GetOrientation: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in DspA_GetZoomLevel: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
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

DspA_tErr DspA_InitSurfaceLists_( tkmDisplayAreaRef this,
          int               inNumLists ) {

  DspA_tErr           eResult    = DspA_tErr_NoErr;
  int                 nDrawList  = 0;

  /* allocate surface point lists. */
  this->maSurfaceLists = (xGrowableArrayRef*) 
    malloc (sizeof(xGrowableArrayRef) * inNumLists );

  /* set all lists to null */
  for( nDrawList = inNumLists - 1; nDrawList >= 0; nDrawList-- ) {
    this->maSurfaceLists[nDrawList] = NULL;
  }

  goto cleanup;

  goto error;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_InitSurfaceLists_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_PurgeSurfaceLists_ ( tkmDisplayAreaRef this ) {

  DspA_tErr           eResult    = DspA_tErr_NoErr;
  int                 nDrawList  = 0;
  xGArr_tErr          eList      = xGArr_tErr_NoErr;

  for( nDrawList = DspA_GetNumSurfaceLists_( this ) - 1; 
       nDrawList >= 0;
       nDrawList-- ) {

    /* if this is a list... */
    if( NULL != this->maSurfaceLists[nDrawList] ) {
      
      /* delete the list. */
      eList = xGArr_Delete( &this->maSurfaceLists[nDrawList] );
      if( xGArr_tErr_NoErr != eList ) {
  eResult = DspA_tErr_ErrorAccessingSurfaceList;
  goto error;
      }
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_PurgeSurfaceLists_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_NewSurfaceList_ ( tkmDisplayAreaRef this,
         mri_tOrientation  iOrientation,
         Surf_tVertexSet  iSurface,
         int               inSlice ) {
  
  DspA_tErr   eResult   = DspA_tErr_NoErr;
  int         nDrawList = 0;
  xGArr_tErr  eList     = xGArr_tErr_NoErr;

  /* get the list index. */
  nDrawList = DspA_GetSurfaceListIndex_( this, iOrientation, 
           iSurface, inSlice );

  /* allocate a list. */
  xGArr_New( &this->maSurfaceLists[nDrawList],
       sizeof( xPoint2f ), 512 );
  if( xGArr_tErr_NoErr != eList ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_NewSurfaceList_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

xGrowableArrayRef DspA_GetSurfaceList_ ( tkmDisplayAreaRef this,
           mri_tOrientation  iOrientation,
           Surf_tVertexSet  iSurface,
           int               inSlice ) {

  xGrowableArrayRef pList = NULL;

  if( NULL != this->maSurfaceLists ) {
    pList = this->maSurfaceLists 
     [ DspA_GetSurfaceListIndex_( this, iOrientation, iSurface, inSlice ) ];
  }

  return pList;
}

int DspA_GetNumSurfaceLists_ ( tkmDisplayAreaRef this ) {

  return this->mnVolumeSize * Surf_knNumVertexSets * mri_knNumOrientations;
}

int DspA_GetSurfaceListIndex_ ( tkmDisplayAreaRef this,
        mri_tOrientation  iOrientation,
        Surf_tVertexSet  iSurface,
        int               inSlice ) {

  return ((int)iOrientation * Surf_knNumVertexSets * this->mnVolumeSize) +
   ((int)iSurface * this->mnVolumeSize) + inSlice;
}




DspA_tErr DspA_ConvertVolumeToBuffer_ ( tkmDisplayAreaRef this,
        xVoxelRef          ipVolumeVox,
          xPoint2nRef       opBufferPt ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fVolumeSize = (float) this->mnVolumeSize;

  /* verify the voxel */
  eResult = DspA_VerifyVolumexVoxl_( this, ipVolumeVox );
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

  /* now zoom the coords to our zoomed buffer state */
  fX = (fZoomLevel * (fX - xVoxl_GetFloatX(this->mpZoomCenter))) +
   (fVolumeSize/2.0);
  fY = (fZoomLevel * (fY - xVoxl_GetFloatY(this->mpZoomCenter))) +
    (fVolumeSize/2.0);

  /* return the point */
  opBufferPt->mnX = (int) rint( fX );
  opBufferPt->mnY = (int) rint( fY );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DspA_ConvertVolumeToBuffer_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertBufferToVolume_ ( tkmDisplayAreaRef this,
          xPoint2nRef       ipBufferPt,
        xVoxelRef          opVolumeVox ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
 xVoxelRef  pVolumeVox  = NULL;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fBufferX    = (float) ipBufferPt->mnX;
  float     fBufferY    = (float) ipBufferPt->mnY;
  float     fVolumeSize = (float) this->mnVolumeSize;

  xVoxl_New ( &pVolumeVox );

  /* unzoom the coords */
  fX = fBufferX / fZoomLevel + 
    ( xVoxl_GetFloatX(this->mpZoomCenter) - (fVolumeSize/2.0/fZoomLevel) );
  fY = fBufferY / fZoomLevel + 
    ( xVoxl_GetFloatY(this->mpZoomCenter) - (fVolumeSize/2.0/fZoomLevel) );

  /* build a voxel out of these two coords and the current slice */
  switch ( this->mOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_SetFloat( pVolumeVox, fX, fY, DspA_GetCurrentSliceNumber_( this ) );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_SetFloat( pVolumeVox, fX, DspA_GetCurrentSliceNumber_( this ), fY );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetFloat( pVolumeVox, DspA_GetCurrentSliceNumber_( this ), fY, fX );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* copy into the outgoing voxel to return it */
  xVoxl_Copy( opVolumeVox, pVolumeVox );
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DspA_ConvertBufferToVolume_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete ( &pVolumeVox );
  
  return eResult;
}

DspA_tErr DspA_ConvertPlaneToVolume_ ( tkmDisplayAreaRef this,
               xPoint2nRef       ipPlanePt,
               int               inSlice,
               mri_tOrientation  iOrientation,
              xVoxelRef          opVolumeVox ) {
  
  DspA_tErr eResult     = DspA_tErr_NoErr;
 xVoxelRef  pVolumeVox  = NULL;

  xVoxl_New ( &pVolumeVox );

  /* build a voxel out of these two coords and the slice */
  switch ( iOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_Set( pVolumeVox, 
         ipPlanePt->mnX, ipPlanePt->mnY, inSlice );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_Set( pVolumeVox, 
         ipPlanePt->mnX, inSlice, ipPlanePt->mnY );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_Set( pVolumeVox, 
         inSlice, ipPlanePt->mnY, ipPlanePt->mnX );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* copy into the outgoing voxel to return it */
  xVoxl_Copy( opVolumeVox, pVolumeVox );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_ConvertPlaneToVolume_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete ( &pVolumeVox );

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
    DebugPrint "Error %d in DspA_DspA_ConvertBufferToScreen_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_ConvertScreenToBuffer_ ( tkmDisplayAreaRef this,
          xPoint2nRef       ipScreenPt,
          xPoint2nRef       opBufferPt ) {

  DspA_tErr eResult = DspA_tErr_NoErr;
  xPoint2n  localPt = {0, 0};
  int       nX      = 0;
  int       nY      = 0;

  /* first subtract our position from the point to get our local point. the
     y craziness has to do with the fact that not only is the y coordinate
     system flipped, but the displays are flipped in the medit window. this
     just kind of works, at least for 2x2. */
  localPt.mnX = ipScreenPt->mnX - this->mLocationInSuper.mnX;
  localPt.mnY = ((ipScreenPt->mnY + this->mnHeight) % this->mnHeight);

  /* verify the screen pt */
  eResult = DspA_VerifyScreenPoint_( this, &localPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  /* use the frame scale factor to convert */
  nX = (float)localPt.mnX / this->mfFrameBufferScaleX;
  nY = (float)localPt.mnY / this->mfFrameBufferScaleY;

  /* y flip */
  nY = GLDRAW_Y_FLIP(nY);

  /* return the point */
  opBufferPt->mnX = nX;
  opBufferPt->mnY = nY;

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DspA_ConvertScreenToBuffer_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

DspA_tErr DspA_SendViewStateToTcl_ ( tkmDisplayAreaRef this ) {

  DspA_tErr             eResult            = DspA_tErr_NoErr;
  char                  sTclArguments[256] = "";
 xVoxelRef              pVoxel             = NULL;
  unsigned char         ucVolumeValue      = 0;
  FunV_tErr             eFunctional        = FunV_tErr_NoError;
  FunV_tFunctionalValue functionalValue    = 0;
  int                   nFlag              = 0;

  xVoxl_New( &pVoxel );

  /* send the cursor. */
  sprintf( sTclArguments, "%d %d %d",
     xVoxl_ExpandInt( this->mpCursor ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeCursor, sTclArguments );

  /* also convert to RAS and send those coords along. */
  tkm_ConvertVolumeToRAS( this->mpCursor, pVoxel );
  sprintf( sTclArguments, "%.1f %.1f %.1f", xVoxl_ExpandFloat( pVoxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateRASCursor, sTclArguments );
 
  /* also convert to Tal and send those coords along. */
  tkm_ConvertVolumeToTal( this->mpCursor, pVoxel );
  sprintf( sTclArguments, "%.1f %.1f %.1f", xVoxl_ExpandFloat( pVoxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateTalCursor, sTclArguments );
 
  /* send volume name */
  sprintf( sTclArguments, "\"%s value\"", tkm_GetVolumeName() );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );

  /* also get the volume value and send that along. */
  ucVolumeValue = tkm_GetVolumeValue( this->mpVolume, this->mpCursor );
  sprintf( sTclArguments, "%d", ucVolumeValue );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValue, sTclArguments );

  /* if aux volume is loaded... */
  if( NULL != this->mpAuxVolume ) {
    
    /* send name */
    sprintf( sTclArguments, "\"%s value\"", tkm_GetAuxVolumeName() );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeName, sTclArguments );
    
    /* and value */
    ucVolumeValue = tkm_GetVolumeValue( this->mpAuxVolume, this->mpCursor );
    sprintf( sTclArguments, "%d", ucVolumeValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeValue, sTclArguments );
  }
  
  /* also get the volume value and send that along. */
  ucVolumeValue = tkm_GetVolumeValue( this->mpVolume, this->mpCursor );
  sprintf( sTclArguments, "%d", ucVolumeValue );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValue, sTclArguments );

  /* and the parcellation label if it's on */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_ParcellationOverlay] 
      && NULL != this->mpParcellationVolume ) {
    tkm_GetParcellationLabel( this->mpCursor, sTclArguments );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateParcellationLabel, 
      sTclArguments );
  }

  /* also see if we have functional data and can send a value for that
     as well. */     
  if ( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ) {
    eFunctional = FunV_GetValueAtAnaIdx( this->mpFunctionalVolume,
           this->mpCursor, 
           &functionalValue );
    if( FunV_tErr_NoError == eFunctional ) {
      sprintf( sTclArguments, "%f", functionalValue );
      tkm_SendTclCommand( tkm_tTclCommand_UpdateFunctionalValue, sTclArguments );
    }
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
      (int)snBrushRadius, (int)sBrushShape, 
      (int)sbBrush3D );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateBrush, sTclArguments );

  /* send the threshold info. */
  sprintf ( sTclArguments, "%d %d %d",
      (int)snBrushThresholdLow, (int)snBrushThresholdHigh, 
      (int)snBrushNewValue );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateBrushThreshold, sTclArguments );

  xVoxl_Delete( &pVoxel );

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

DspA_tErr DspA_VerifyVolumexVoxl_  ( tkmDisplayAreaRef this,
            xVoxelRef          ipVoxel ) {
  DspA_tErr eResult = DspA_tErr_NoErr;

  /* make sure voxel is in bounds */
  if ( xVoxl_GetX( ipVoxel )    < 0
       || xVoxl_GetY( ipVoxel ) < 0
       || xVoxl_GetZ( ipVoxel ) < 0
       || xVoxl_GetX( ipVoxel ) >= this->mnVolumeSize
       || xVoxl_GetY( ipVoxel ) >= this->mnVolumeSize
       || xVoxl_GetZ( ipVoxel ) >= this->mnVolumeSize ) {

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
       || ipBufferPt->mnX >= this->mnVolumeSize
       || ipBufferPt->mnY >= this->mnVolumeSize ) {

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

  glOrtho        ( 0, this->mnVolumeSize, this->mnVolumeSize, 0, -1.0, 1.0 );
  //glOrtho        ( 0, this->mnVolumeSize, 0, this->mnVolumeSize, -1.0, 1.0 );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint "glOrtho got error %d\n", eGL EndDebugPrint;

  glViewport( this->mLocationInSuper.mnX, this->mLocationInSuper.mnY,
        this->mnWidth, this->mnHeight );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint "glViewport got error %d\n", eGL EndDebugPrint;

  glRasterPos2i  ( 0, this->mnVolumeSize );
  //glRasterPos2i  ( 0, 0 );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint "glRasterPos2i got error %d\n", eGL EndDebugPrint;

  glPixelZoom ( this->mfFrameBufferScaleX, this->mfFrameBufferScaleY );
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL )
    DebugPrint "glPixelZoom got error %d\n", eGL EndDebugPrint;
}

void DspA_DebugPrint_ ( tkmDisplayAreaRef this ) {

  DebugPrint "tkmDisplayArea\n" EndDebugPrint;
  DebugPrint "\tx %d y %d w %d h %d\n", 
    this->mLocationInSuper.mnX, this->mLocationInSuper.mnY,
    this->mnWidth, this->mnHeight EndDebugPrint;
  DebugPrint "\tvolume size %d, x scale %.2f y scale %.2f\n",
    this->mnVolumeSize, this->mfFrameBufferScaleX, 
    this->mfFrameBufferScaleY EndDebugPrint;
  DebugPrint "\tzoom level %d center %d, %d, %d\n",
    this->mnZoomLevel, xVoxl_ExpandInt(this->mpZoomCenter) EndDebugPrint;
  DebugPrint "\tcursor %d, %d, %d orientation %s slice %d tool %d\n",
    xVoxl_ExpandInt(this->mpCursor), DspA_ksaOrientation[this->mOrientation],
    DspA_GetCurrentSliceNumber_(this), sTool EndDebugPrint;
}

void DspA_Signal ( char* isFuncName, int inLineNum, DspA_tErr ieCode ) {

  DebugPrint "Signal in %s, line %d: %d, %s\n", 
    isFuncName, inLineNum, ieCode, DspA_GetErrorString(ieCode) EndDebugPrint;
}

char* DspA_GetErrorString ( DspA_tErr ieCode ) {

  DspA_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= DspA_knNumErrorCodes ) {
    eCode = DspA_tErr_InvalidErrorCode;
  }

  return DspA_ksaErrorStrings [eCode];
}

tBoolean xUtil_FaceIntersectsPlane( MRI_SURFACE*     ipSurface,
            face_type*       ipFace,
            int              inPlane,
            Surf_tVertexSet iSurface,
            mri_tOrientation iOrientation ) {

  int          nVertexA    = 0;
  int          nVertexB    = 0;
  vertex_type* pVertexA    = NULL;
  vertex_type* pVertexB    = NULL;
 xVoxelRef     pVoxelA     = NULL;
 xVoxelRef     pVoxelB     = NULL;
  float        fDistanceA  = 0;
  float        fDistanceB  = 0;
  tBoolean     bIntersects = FALSE;

  xVoxl_New( &pVoxelA );
  xVoxl_New( &pVoxelB );

  /* for every vertex... */
  for ( nVertexA = 0; nVertexA < VERTICES_PER_FACE; nVertexA++ ) {
    
    /* get the neigboring vertex. */
    nVertexB = nVertexA - 1;
    if ( nVertexB < 0 )
      nVertexB = VERTICES_PER_FACE - 1;
    
    /* get the vertices. */
    pVertexA = &(ipSurface->vertices[ipFace->v[nVertexA]]);
    pVertexB = &(ipSurface->vertices[ipFace->v[nVertexB]]);

    /* convert them to normalized voxels. */
    xUtil_NormalizeVertexToVoxel( pVertexA, iSurface, 
          iOrientation, pVoxelA );
    xUtil_NormalizeVertexToVoxel( pVertexB, iSurface, 
          iOrientation, pVoxelB );

    /* find the distance between the two points and the plane. */
    fDistanceA = xVoxl_GetZ( pVoxelA ) - inPlane;
    fDistanceB = xVoxl_GetZ( pVoxelB ) - inPlane;

    /* if product is negative, they cross the plane. if it is 0,
       they both lie on the plane. */
    if ( fDistanceA * fDistanceB <= 0 ) {
      bIntersects = TRUE;
      goto cleanup;
    }
    
  }

  goto cleanup;

  goto error;
 error:

  /* print error message */

 cleanup:

  xVoxl_Delete( &pVoxelA );
  xVoxl_Delete( &pVoxelB );

  return bIntersects;
}

void xUtil_NormalizeVertexToVoxel( vertex_type*     ipVertex,
           Surf_tVertexSet iSurface,
           mri_tOrientation iOrientation,
          xVoxelRef         opVoxel ) {

 xVoxelRef pRASVox = NULL;
 xVoxelRef pAnatomicalVox = NULL;
  Real rXVox = 0;
  Real rYVox = 0;
  Real rZVox = 0;

  xVoxl_New( &pRASVox );
  xVoxl_New( &pAnatomicalVox );
  
  switch( iSurface ) {
  case Surf_tVertexSet_Main:
    xVoxl_SetFloat( pRASVox, ipVertex->x, ipVertex->y, ipVertex->z );
    break;
  case Surf_tVertexSet_Pial:
    xVoxl_SetFloat( pRASVox, ipVertex->cx, ipVertex->cy, ipVertex->cz );
    break;
  case Surf_tVertexSet_Original:
    xVoxl_SetFloat( pRASVox, ipVertex->origx, ipVertex->origy, ipVertex->origz );
    break;
  default:
    break;
  }

  tkm_ConvertRASToVolume( pRASVox, pAnatomicalVox );
  rXVox = xVoxl_GetFloatX( pAnatomicalVox );
  rYVox = xVoxl_GetFloatY( pAnatomicalVox );
  rZVox = xVoxl_GetFloatZ( pAnatomicalVox );


  switch( iOrientation ) {
  case mri_tOrientation_Horizontal:
    xVoxl_SetFloat( opVoxel, (float)rXVox, (float)rZVox, (float)rYVox );
    break;
  case mri_tOrientation_Coronal:
    xVoxl_SetFloat( opVoxel, (float)rXVox, (float)rYVox, (float)rZVox );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_SetFloat( opVoxel, (float)rZVox, (float)rYVox, (float)rXVox );
    break;
  default:
    break;
  }

  xVoxl_Delete( &pRASVox );
  xVoxl_Delete( &pAnatomicalVox );
}

tBoolean xUtil_LineIntersectsPlane(xVoxelRef         ipLineVoxA,
           xVoxelRef         ipLineVoxB,
            int              inPlane,
            tBoolean         ipInterpolate,
            xPoint2fRef      opIntersectionPt ) {
  
  float    fPlane           = inPlane;
  float    fDistanceA       = 0;
  float    fDistanceB       = 0;
  float    fAlpha           = 1;
  xPoint2f intersectionPt   = {0, 0};

  /* get distance from each to plane. */
  fDistanceA = xVoxl_GetFloatZ( ipLineVoxA ) - fPlane;
  fDistanceB = xVoxl_GetFloatZ( ipLineVoxB ) - fPlane;

  /* if product is negative or 0, they intersect the plane. */
  if ( fDistanceA * fDistanceB > 0.0 ) {
    return FALSE;
  }

  /* if averaging the vertices, find an iterpolation factor, which is the
     intersection of the line the plane. */
  if( ipInterpolate ) {
    
    /* make sure they arn't on the same plane... */
    if( xVoxl_GetFloatZ(ipLineVoxB) - xVoxl_GetFloatZ(ipLineVoxA) != 0.0 ) {

      fAlpha = (fPlane - xVoxl_GetFloatZ( ipLineVoxA )) /
  (xVoxl_GetFloatZ( ipLineVoxB ) - xVoxl_GetFloatZ( ipLineVoxA ));
      
    } else {

      /* alpha is just 1. */
      fAlpha = 1.0;
    }

    /* interpolate to find the intersection. */
    intersectionPt.mfX = (xVoxl_GetFloatX( ipLineVoxA ) +
      fAlpha * (xVoxl_GetFloatX(ipLineVoxB) - xVoxl_GetFloatX(ipLineVoxA)));
    intersectionPt.mfY = (xVoxl_GetFloatY( ipLineVoxA ) +
      fAlpha * (xVoxl_GetFloatY(ipLineVoxB) - xVoxl_GetFloatY(ipLineVoxA)));

  } else {
    
    /* if no interpolationg, intersection is projection onto plane */
    intersectionPt.mfX = xVoxl_GetFloatX( ipLineVoxB );
    intersectionPt.mfY = xVoxl_GetFloatY( ipLineVoxB );
  }

  /*
  DebugPrint "intersecting %.2f,%.2f,%.2f to %.2f,%.2f,%.2f on %d",
    xVoxl_ExpandFloat( ipLineVoxA ), xVoxl_ExpandFloat( ipLineVoxB ),
    inPlane EndDebugPrint;

  DebugPrint "hit %d,%d\n",
    intersectionPt.mnX, intersectionPt.mnY EndDebugPrint;
  */

  /* return the point. */
  *opIntersectionPt = intersectionPt;

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

DspA_tErr DspA_ParsePointList_( tkmDisplayAreaRef this,
        GLenum            inMode,
        xGrowableArrayRef iList ) {

  DspA_tErr  eResult        = DspA_tErr_NoErr;
  xGArr_tErr eList          = xGArr_tErr_NoErr;
  tBoolean   bOperationOpen = FALSE;
  xPoint2f   intersectionPt = {0,0};
  xPoint2f   drawPoint      = {0,0};

  /* reset the list position. */
  eList = xGArr_ResetIterator( iList );
  if( xGArr_tErr_NoErr != eList )
    goto error;
  
  /* start new operation. */
  glBegin( inMode );
  bOperationOpen = TRUE;
  
  while ( (eList = xGArr_NextItem( iList, (void*)&intersectionPt ))
    == xGArr_tErr_NoErr ) {
    
    /* if it is -1,-1, it's a face marker... */
    if( -1.0 == intersectionPt.mfX 
  && -1.0 == intersectionPt.mfY ) {
      
      /* if operation was still going, end it. */ 
      if( bOperationOpen )
  glEnd();
      
      /* start new operation. */
      glBegin( inMode );
      bOperationOpen = TRUE;
      
    } else {
      
      drawPoint.mfX = intersectionPt.mfX;
      drawPoint.mfY = intersectionPt.mfY;
      
      /* adjust the point */
      DspA_AdjustSurfaceDrawPoint_( this, &drawPoint );
      
      /* convert to zoomed coords. */
      drawPoint.mfX = ((float)this->mnZoomLevel * (drawPoint.mfX - xVoxl_GetFloatX(this->mpZoomCenter))) + (float)(this->mnVolumeSize/2.0);
      drawPoint.mfY = ((float)this->mnZoomLevel * (drawPoint.mfY - xVoxl_GetFloatY(this->mpZoomCenter))) + (float)(this->mnVolumeSize/2.0);
    
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
    DebugPrint "Error %d in DspA_ParsePointList_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }
    
 cleanup:
  
  return eResult;
}
