#include "tkmDisplayArea.h"
#include "tkmMeditWindow.h"
#include "tkmFunctionalVolume.h"
#include "mritransform.h"
#include "xUtilities.h"

/* i'm not sure what to do about these y flips. it seems that whenever we're
   using a point that's going to go into the buffer to be drawn to the screen,
   we should flip when not in horizontal view. when using regular gl drawing 
   commands to draw to the screen, only do it in horizontal orientation. */
#define BUFFER_Y_FLIP(y) ( this->mOrientation != tkm_tOrientation_Horizontal? \
                          (this->mnVolumeSize - y) : y )
#define GLDRAW_Y_FLIP(y) ( this->mOrientation == tkm_tOrientation_Horizontal? \
                          (this->mnVolumeSize - y) : y )
#define GLDRAW_Y_FLIP_FLOAT(y) \
                         ( this->mOrientation == tkm_tOrientation_Horizontal? \
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
  "Error creating surface hash table.",
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

char DspA_ksaOrientation [tkm_knNumOrientations][256] = {
  "Coronal",
  "Horizontal",
  "Sagittal"
};

char DspA_ksaSurface [tkm_knNumSurfaceTypes][256] = {
  "Main",
  "Original",
  "Canonical"
};

DspA_tErr DspA_New ( tkmDisplayAreaRef* oppWindow,
         tkmMeditWindowRef  ipWindow ) {

  DspA_tErr         eResult   = DspA_tErr_NoErr;
  tkmDisplayAreaRef this      = NULL;
  int               nFlag     = 0;

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
  Voxel_New ( &this->mpCursor );
  Voxel_New ( &this->mpZoomCenter );

  /* stuff in default values for display states. */
  this->mOrientation           = tkm_tOrientation_Coronal;
  this->mnZoomLevel            = 1;
  this->mnHilitedVertexIndex   = -1;
  this->mHilitedSurface        = tkm_tSurfaceType_None;
  this->mbSliceChanged         = TRUE;

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
  
  DspA_tErr         eResult = DspA_tErr_NoErr;
  tkmDisplayAreaRef this    = NULL;

  /* get us */
  this = *ioppWindow;
    
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* delete our voxels */
  Voxel_Delete ( &this->mpCursor );
  Voxel_Delete ( &this->mpZoomCenter );

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
  VoxelRef  pCenter            = NULL;
  char      sTclArguments[256] = "";

  Voxel_New( &pCenter );

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

  /* set inital values for our buffer scale */
  this->mfFrameBufferScaleX = 
    (float)this->mnWidth  / (float)this->mnVolumeSize;
  this->mfFrameBufferScaleY = 
    (float)this->mnHeight  / (float)this->mnVolumeSize;

  /* initialize surface point lists. */
  eResult = DspA_InitSurfaceLists_( this, 
      this->mnVolumeSize * tkm_knNumSurfaceTypes * tkm_knNumOrientations );
  if( DspA_tErr_NoErr != eResult )
    goto error;

  /* send volume name to tk window */
  sprintf( sTclArguments, "\"%s value\"", tkm_GetVolumeName() );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeName, sTclArguments );
   
  /* get the center of the volume */
  Voxel_Set( pCenter, this->mnVolumeSize/2, 
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

  Voxel_Delete( &pCenter );

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

  /* send volume name to tk window */
  sprintf( sTclArguments, "\"%s value\"", tkm_GetAuxVolumeName() );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeName, sTclArguments );

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

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* make sure this size is the same as the main volume size. */
  if( inSize != this->mnVolumeSize ) {
    eResult = DspA_tErr_InvalidParameter;
    goto error;
  }

  /* save the parcellation volume */
  this->mpParcellationVolume = ipVolume;

  /* turn parcellation on */
  eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_ParcellationOverlay,
         TRUE );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

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
          MRI_SURFACE*      ipSurface ) {

  DspA_tErr        eResult = DspA_tErr_NoErr;
  tkm_tSurfaceType surface = tkm_tSurfaceType_None;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* save the surface */
  this->mpSurface = ipSurface;

#if 0
  /* if we have a surface... */
  if( NULL != this->mpSurface ) {

    /* for each surface... */
    for( surface = tkm_tSurfaceType_Current; 
   surface < tkm_knNumSurfaceTypes; surface++ ) {

      /* if we have a hash table for this surface, delete it. */
      if( NULL != (this->mpSurfaceHashTable[surface]) )
  MHTfree( &(this->mpSurfaceHashTable[surface]) );

      /* create the hash table for this surface */
      DebugPrint "DspA_SetSurface(): building table for %d...",
  (int) surface EndDebugPrint;
      this->mpSurfaceHashTable[surface] = 
  MHTfillTableAtResolution( this->mpSurface, NULL, surface+1, 
          DspA_kfHashTableResolution );
      DebugPrint " done!\n" EndDebugPrint;
      if( NULL == (this->mpSurfaceHashTable[surface]) ) {
  eResult = DspA_tErr_ErrorCreatingSurfaceHashTable;
  goto error;
      }
    }
  }
#endif

  /* if not null, turn surfaces and interpolation on. */
  if( NULL != this->mpSurface ) {

    eResult = DspA_SetDisplayFlag( this, DspA_tDisplayFlag_MainSurface, TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    eResult = DspA_SetDisplayFlag( this, 
      DspA_tDisplayFlag_InterpolateSurfaceVertices, TRUE );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* purge all the surface lists. */
  DspA_PurgeSurfaceLists_( this );

  /* set slice dirty flag and redraw */
  this->mbSliceChanged = TRUE;
  DspA_Redraw_( this );

  goto cleanup;

 error:

#if 0
  /* dispose anything that was allocated */
  for( surface = tkm_tSurfaceType_None; 
       surface < tkm_knNumSurfaceTypes; surface++ ) {
    
    /* if we have a hash table for this surface, delete it. */
    if( NULL != (this->mpSurfaceHashTable[surface]) )
      MHTfree( &(this->mpSurfaceHashTable[surface]) );
  }
#endif

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
               VoxelSpaceRef     ipVoxels ) {

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
                 VoxelListRef      ipVoxels ) {

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
           VoxelSpaceRef     ipVoxels ) {

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
         VoxelRef          ipCursor ) {

  DspA_tErr             eResult            = DspA_tErr_NoErr;
  VoxelRef              pVoxel             = NULL;
  int                   nSlice             = 0;
  unsigned char         ucVolumeValue      = 0;
  FunV_tErr             eFunctional        = FunV_tErr_NoError;
  FunV_tFunctionalValue functionalValue   = 0;
  char                  sTclArguments[256] = "";
 
  Voxel_New( &pVoxel );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the cursor */
  eResult = DspA_VerifyVolumeVoxel_ ( this, ipCursor );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* get our current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* set the cursor */
  Voxel_Copy( this->mpCursor, ipCursor );
  
  /* if the new slice number is diffrent, set our dirty slice flag. */
  if( DspA_GetCurrentSliceNumber_( this ) != nSlice ) {
    this->mbSliceChanged = TRUE;
  }

  /* if cursor is .0 .0 .0, change to .5 .5 .5 so that it will draw in the
     center of a voxel on screen. */
  if( Voxel_GetFloatX(this->mpCursor) == 
      (float)(int)Voxel_GetFloatX(this->mpCursor)
      && Voxel_GetFloatY(this->mpCursor) == 
      (float)(int)Voxel_GetFloatY(this->mpCursor)
      && Voxel_GetFloatZ(this->mpCursor) == 
      (float)(int)Voxel_GetFloatZ(this->mpCursor) ) {
    
    Voxel_SetFloat( this->mpCursor,
        Voxel_GetFloatX(this->mpCursor) + 0.5,
        Voxel_GetFloatY(this->mpCursor) + 0.5,
        Voxel_GetFloatZ(this->mpCursor) + 0.5 );
  }

  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
    /* notify the window that the cursor has changed. */
    MWin_CursorChanged( this->mpWindow, this, this->mpCursor );

    /* send the cursor. */
    sprintf( sTclArguments, "%d %d %d",
       EXPAND_VOXEL_INT( this->mpCursor ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeCursor, sTclArguments );
    
    /* also convert to RAS and send those coords along. */
    tkm_ConvertVolumeToRAS( this->mpCursor, pVoxel );
    sprintf( sTclArguments, "%.1f %.1f %.1f", EXPAND_VOXEL_FLOAT( pVoxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateRASCursor, sTclArguments );
    
    /* also convert to Tal and send those coords along. */
    tkm_ConvertVolumeToTal( this->mpCursor, pVoxel );
    sprintf( sTclArguments, "%.1f %.1f %.1f", EXPAND_VOXEL_FLOAT( pVoxel ) );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateTalCursor, sTclArguments );
    
    /* also get the volume value and send that along. */
    ucVolumeValue = tkm_GetVolumeValue( this->mpVolume, this->mpCursor );
    sprintf( sTclArguments, "%d", ucVolumeValue );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValue, sTclArguments );
    
    /* and the aux volume value if there is one. */
    if( NULL != this->mpAuxVolume ) {
      ucVolumeValue = tkm_GetVolumeValue( this->mpAuxVolume, this->mpCursor );
      sprintf( sTclArguments, "%d", ucVolumeValue );
      tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeValue, sTclArguments );
    }

    /* also see if we have functional data and can send a value for that
       as well. */     
    if ( NULL != this->mpFunctionalVolume ) {
      DisableDebuggingOutput;
      eFunctional = FunV_GetValueAtAnatomicalVoxel( this->mpFunctionalVolume,
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
      eResult, EXPAND_VOXEL_INT(ipCursor),
      DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  Voxel_Delete( &pVoxel );

  return eResult;
}

DspA_tErr DspA_SetOrientation ( tkmDisplayAreaRef this, 
        tkm_tOrientation  iOrientation ) {

  DspA_tErr eResult            = DspA_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the orientation */
  if ( iOrientation <= tkm_tOrientation_None
       || iOrientation >= tkm_knNumOrientations ) {
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
  VoxelRef  pCenter            = NULL;

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
    Voxel_New( &pCenter );
    Voxel_Set( pCenter, this->mnVolumeSize/2, 
         this->mnVolumeSize/2, this->mnVolumeSize/2 );
    eResult = DspA_SetZoomCenter( this, pCenter );
    Voxel_Delete( &pCenter );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* this requires a rebuild of the frame buffer. */
  this->mbSliceChanged = TRUE;

  /* if we're the currently focused display... */
  if( sFocusedDisplay == this ) {
    
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
             VoxelRef          ipCenter ) {

  DspA_tErr eResult  = DspA_tErr_NoErr;
  int       nX       = 0;
  int       nY       = 0;
  
  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* verify the center */
  eResult = DspA_VerifyVolumeVoxel_ ( this, ipCenter );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    nX = Voxel_GetX( ipCenter );
    nY = Voxel_GetY( ipCenter );
    break;
  case tkm_tOrientation_Horizontal:
    nX = Voxel_GetX( ipCenter );
    nY = Voxel_GetZ( ipCenter );
    break;
  case tkm_tOrientation_Sagittal:
    nX = Voxel_GetZ( ipCenter );
    nY = Voxel_GetY( ipCenter );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }
  
  /* set the center */
  Voxel_SetX( this->mpZoomCenter, nX );
  Voxel_SetY( this->mpZoomCenter, nY );
  
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
             tkm_tSurfaceType  inSurface,
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
  VoxelRef  pCursor = NULL;

  Voxel_New( &pCursor );

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* copy the cursor. */
  Voxel_Copy( pCursor, this->mpCursor );

  /* apply the increment to the proper part of the cursor */
  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    Voxel_SetZ( pCursor, Voxel_GetZ( pCursor ) + inDelta );
    break;
  case tkm_tOrientation_Horizontal:
    Voxel_SetY( pCursor, Voxel_GetY( pCursor ) + inDelta );
    break;
  case tkm_tOrientation_Sagittal:
    Voxel_SetX( pCursor, Voxel_GetX( pCursor ) + inDelta );
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

  Voxel_Delete( &pCursor );

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
  VoxelRef    pVolumeVox  = NULL;
  FunV_tErr   eFunctional = FunV_tErr_NoError;

  Voxel_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  DebugPrint "Mouse up screen x %d y %d buffer x %d y %d volume %d %d %d\n",
    ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
    EXPAND_VOXEL_INT( pVolumeVox ) EndDebugPrint;

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

  Voxel_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleMouseDown_ ( tkmDisplayAreaRef this, 
          xGWin_tEventRef   ipEvent ) {

  DspA_tErr   eResult = DspA_tErr_NoErr;
  xPoint2n    bufferPt    = {0,0};
  VoxelRef    pVolumeVox  = NULL;

  Voxel_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
#if 0
  DebugPrint "Mouse down screen x %d y %d buffer x %d y %d volume %d %d %d\n",
    ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
    EXPAND_VOXEL_INT( pVolumeVox ) EndDebugPrint;
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

  Voxel_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleMouseMoved_ ( tkmDisplayAreaRef this, 
           xGWin_tEventRef   ipEvent ) {
  DspA_tErr   eResult = DspA_tErr_NoErr;
  xPoint2n    bufferPt    = {0,0};
  VoxelRef    pVolumeVox  = NULL;

  Voxel_New( &pVolumeVox );
  
  eResult = DspA_ConvertScreenToBuffer_( this, &(ipEvent->mWhere), &bufferPt );
  if ( DspA_tErr_NoErr != eResult )
    goto error;
  
  eResult = DspA_ConvertBufferToVolume_( this, &bufferPt, pVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

#if 0
  DebugPrint "Mouse moved screen x %d y %d buffer x %d y %d volume %d %d %d\n",
    ipEvent->mWhere.mnX, ipEvent->mWhere.mnY, bufferPt.mnX, bufferPt.mnY,
    EXPAND_VOXEL_INT( pVolumeVox ) EndDebugPrint;
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

  Voxel_Delete( &pVolumeVox );

  return eResult;
}

DspA_tErr DspA_HandleKeyDown_ ( tkmDisplayAreaRef this, 
        xGWin_tEventRef   ipEvent ) {

  DspA_tErr eResult       = DspA_tErr_NoErr;
  MWin_tErr eWindowResult = MWin_tErr_NoErr;

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
    eResult = DspA_SetOrientation( this, tkm_tOrientation_Sagittal );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;
        
  case 'y': 
    /* y sets plane to horizontal */
    eResult = DspA_SetOrientation( this, tkm_tOrientation_Horizontal );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    break;

  case 'z': 
    /* z sets plane to coronal */
    eResult = DspA_SetOrientation( this, tkm_tOrientation_Coronal );
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
            VoxelRef          ipCenterVox,
            void(*ipFunction)(VoxelRef) ) {

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
  VoxelRef         pVolumeVox  = NULL;

  Voxel_New( &pVolumeVox );
  
  /* get our center voxel. */
  nXCenter = Voxel_GetX( ipCenterVox );
  nYCenter = Voxel_GetY( ipCenterVox );
  nZCenter = Voxel_GetZ( ipCenterVox );

  /* set all radii to the brush radius. we subtract one because of our 
     looping bounds. */
  nXRadius = snBrushRadius - 1;
  nYRadius = snBrushRadius - 1;
  nZRadius = snBrushRadius - 1;

  /* if we're not in 3d, set the same radius of the same plane as our
     current orientation to 0. */
  if( !sbBrush3D ) {
    switch( this->mOrientation ) {
    case tkm_tOrientation_Coronal:
      nZRadius = 0;
      break;
    case tkm_tOrientation_Horizontal:
      nYRadius = 0;
      break;
    case tkm_tOrientation_Sagittal:
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
  Voxel_Set( pVolumeVox, nX, nY, nZ );

  /* if it's valid... */
  if( DspA_tErr_InvalidVolumeVoxel != 
      DspA_VerifyVolumeVoxel_( this, pVolumeVox ) ) {

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
  Voxel_Delete( &pVolumeVox );

  return eResult;

}

void DspA_BrushVoxelsInThreshold_ ( VoxelRef ipVoxel ) {

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
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* draw the selection */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_Selection] ) {
    eResult = DspA_DrawSelectionToFrame_( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* draw functional overlay. */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] ) {
    eResult = DspA_DrawFunctionalOverlayToFrame_( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* draw the frame buffer */
  eResult = DspA_DrawFrameBuffer_ ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* draw the control points */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_ControlPoints] ) {
    eResult = DspA_DrawControlPointsToFrame_( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* draw the surface */
  eResult = DspA_DrawSurface_ ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* draw the cursor */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_Cursor] ) {
    eResult = DspA_DrawCursor_ ( this );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* rebuilt all our slice changes, so we can clear this flag. */
  this->mbSliceChanged = FALSE;

  goto cleanup;

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

  DspA_tErr        eResult   = DspA_tErr_NoErr;
  tkm_tSurfaceType surface   = tkm_tSurfaceType_Current;
  xListRef         pList     = NULL;
  float            faColor[3] = {0, 0, 0};

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
    for ( surface = tkm_tSurfaceType_Current;
    surface < tkm_knNumSurfaceTypes; surface++ ) {
      
      /* if this surface is visible... */
      if( this->mabDisplayFlags[surface + DspA_tDisplayFlag_MainSurface] ) {
  
  /* get the list. */
  pList = DspA_GetSurfaceList_( this, this->mOrientation, surface,
              DspA_GetCurrentSliceNumber_(this) );
  if( NULL == pList ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }
  
  /* choose and set the color */
  switch( surface ) {
  case tkm_tSurfaceType_Current:
    faColor[0] = 1.0; faColor[1] = 1.0; faColor[2] = 0.0;
    break;
  case tkm_tSurfaceType_Original:
    faColor[0] = 0.0; faColor[1] = 1.0; faColor[2] = 0.0;
    break;
  case tkm_tSurfaceType_Canonical:
    faColor[0] = 1.0; faColor[1] = 0.0; faColor[2] = 0.0;
    break;
  default:
    faColor[0] = 0.5; faColor[1] = 0.5; faColor[2] = 0.5;
    break;
  }

  glColor3fv( faColor );
  
  /* draw the points. */
  DspA_ParsePointList_( this, GL_LINES, pList );

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
    DspA_ParsePointList_( this, GL_POINTS, pList );
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


DspA_tErr DspA_BuildCurrentFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr             eResult     = DspA_tErr_NoErr;
  xPoint2n              volumePt    = {0, 0};
  GLubyte*              pDest       = NULL;
  unsigned char         ucValue     = 0;
  VoxelRef              pVoxel      = NULL;
  VoxelRef              pFuncRASVox = NULL;
  Real                  rRASX       = 0;
  Real                  rRASY       = 0;
  Real                  rRASZ       = 0;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  FunV_tFunctionalValue funcValue   = 0.0;
  float                 fRed        = 0;
  float                 fGreen      = 0;
  float                 fBlue      = 0;
  tBoolean              bPixelSet   = FALSE;
  int                   nY          = 0;

  xUtil_StartTimer();

  /* make our voxel */
  Voxel_New ( &pVoxel );
  Voxel_New ( &pFuncRASVox );

  /* if we're in zoom level one, set zoom center to 128,128,128 */
  if( 1 == this->mnZoomLevel ) {

    /* set a voxel to 128,128,128, setting the zoom center will not change
       the current slice. */
    Voxel_Set( pVoxel, this->mnVolumeSize/2, 
         this->mnVolumeSize/2, this->mnVolumeSize/2 );
    eResult = DspA_SetZoomCenter( this, pVoxel );
    if( DspA_tErr_NoErr != eResult )
      goto error;
  }

  /* get a ptr to the frame buffer. */
  pDest = this->mpFrameBuffer;

  DisableDebuggingOutput;

  /* go thru the buffer... */
  //  for ( volumePt.mnY = 0; volumePt.mnY < this->mnVolumeSize; volumePt.mnY++ ) {
  for ( nY = 0; nY < this->mnVolumeSize; nY ++ ) {
    for ( volumePt.mnX = 0; volumePt.mnX < this->mnVolumeSize; volumePt.mnX ++ ) {

      /* haven't set this pixel yet. */
      bPixelSet = FALSE;

      /* y flip the volume pt to flip the image over. */
      volumePt.mnY = BUFFER_Y_FLIP(nY);

      /* get a volume voxel.*/
      eResult = DspA_ConvertBufferToVolume_ ( this, &volumePt, pVoxel );
      if ( DspA_tErr_NoErr != eResult )
  goto error;

      /* check it. */
      eResult = DspA_VerifyVolumeVoxel_( this, pVoxel );
      if( DspA_tErr_NoErr == eResult ) {
  
  /* if we have and are showing functional data... */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_FunctionalOverlay] 
      && NULL != this->mpFunctionalVolume ) {
           
    /* convert to ras */
    trans_VoxelIndexToRAS( EXPAND_VOXEL_INT(pVoxel), 
         &rRASX, &rRASY, &rRASZ );
    Voxel_Set( pFuncRASVox, rRASX, rRASY, rRASZ );

    /* get a functional value. */
    eFunctional = FunV_GetValueAtRASVoxel( this->mpFunctionalVolume,
             pFuncRASVox, &funcValue );

    /* if it was a valid voxel */
    if( FunV_tErr_NoError == eFunctional ) {

      /* get the anatomical value for the offset color */
      ucValue = tkm_GetVolumeValue ( this->mpVolume, pVoxel );
      tkm_GetAnatomicalVolumeColor( ucValue, &fRed, &fGreen, &fBlue );
      
      /* get a color value. */
      eFunctional = FunV_GetColorForValue ( this->mpFunctionalVolume,
              funcValue, fRed,
              &fRed, &fGreen, &fBlue );

      /* if the color is not all black.. */
      if( fRed != 0 || fBlue != 0 || fGreen != 0 ) {
        
        /* pixel is set. */
        bPixelSet = TRUE;
      }
    } 
  } 
  
  /* if we are showing parcellation... */
  if( this->mabDisplayFlags[DspA_tDisplayFlag_ParcellationOverlay] ) {
    
    /* get parcellation color. */
    tkm_GetParcellationColor( pVoxel, &fRed, &fGreen, &fBlue );

    /* if it's not zero... */
    if ( fRed != 0 || fBlue != 0 || fGreen != 0 ) {
      
      /* pixel is set. */
      bPixelSet = TRUE;
    }
  }

  /* if we didn't set the pixel color from functional data... */
  if ( !bPixelSet ) {
    
    /* get the plain anatomical value from the main or aux
       volume. */
    if( this->mabDisplayFlags[DspA_tDisplayFlag_AuxVolume] ) {
      ucValue = tkm_GetVolumeValue ( this->mpAuxVolume, pVoxel );
    } else {
      ucValue = tkm_GetVolumeValue ( this->mpVolume, pVoxel );
    }

    /* get the color */
    tkm_GetAnatomicalVolumeColor( ucValue, &fRed, &fGreen, &fBlue );

    /* pixel is set. */
    bPixelSet = TRUE;
  }

      } else {

  /* voxel was out of bounds, set to out of bounds color. */
  fRed   = (GLubyte)0;
  fGreen = (GLubyte)0;
  fBlue  = (GLubyte)0;

  /* clear error flag. */
  eResult = DspA_tErr_NoErr;
      }  

      /* set the pixel */
      pDest[DspA_knRedPixelCompIndex]   = 
  (GLubyte)(fRed * (float)DspA_knMaxPixelValue);
      pDest[DspA_knGreenPixelCompIndex] = 
  (GLubyte)(fGreen * (float)DspA_knMaxPixelValue);
      pDest[DspA_knBluePixelCompIndex]  = 
  (GLubyte)(fBlue * (float)DspA_knMaxPixelValue);
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
  Voxel_Delete ( &pVoxel );
  Voxel_Delete ( &pFuncRASVox );

  xUtil_StopTimer( "build frame buffer" );

  return eResult;
}


DspA_tErr DspA_DrawFunctionalOverlayToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr             eResult     = DspA_tErr_NoErr;
  FunV_tErr             eFunctional = FunV_tErr_NoError;
  xPoint2n              bufferPt    = {0,0};
  FunV_tFunctionalValue max         = 0;
  FunV_tFunctionalValue funcValue   = 0.0;
  float                 fRed        = 0;
  float                 fGreen      = 0;
  float                 fBlue       = 0;
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
           funcValue, 0,
           &fRed, &fGreen, &fBlue );
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
   (GLubyte)(fRed * DspA_knMaxPixelValue);
      pFrame[DspA_knGreenPixelCompIndex] = 
  (GLubyte)(fGreen * DspA_knMaxPixelValue);
      pFrame[DspA_knBluePixelCompIndex]  = 
  (GLubyte)(fBlue * DspA_knMaxPixelValue);
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
  VoxelListRef pList      = NULL;
  char         eList      = 0;
  char         eSpace     = 0;
  int          nNumPoints = 0;
  int          nListIndex = 0;
  VoxelRef     pControlPt = NULL;
  tBoolean     bSelected  = FALSE;
  xPoint2n     bufferPt   = {0,0};
  float        faColor[3] = {0, 0, 0};

  /* new voxel */
  Voxel_New( &pControlPt );

  /* decide which list we want out of the space. */
  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    eSpace = VSpace_GetVoxelsInZPlane( this->mpControlPoints, 
               DspA_GetCurrentSliceNumber_(this), 
               &pList );
    break;
  case tkm_tOrientation_Sagittal:
    eSpace = VSpace_GetVoxelsInXPlane( this->mpControlPoints, 
               DspA_GetCurrentSliceNumber_(this), 
               &pList );
    break;
  case tkm_tOrientation_Horizontal:
    eSpace = VSpace_GetVoxelsInYPlane( this->mpControlPoints, 
               DspA_GetCurrentSliceNumber_(this), 
               &pList );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* check for error. */
  if ( kVSpaceErr_NoErr != eSpace ) {
    eResult = DspA_tErr_ErrorAccessingControlPoints;
    goto error;
  }

  /* get the number of points. */
  eList = VList_GetCount( pList, &nNumPoints );
  if ( kVListErr_NoErr != eList ) {
    eResult = DspA_tErr_ErrorAccessingControlPoints;
    goto error;
  }

  /* for each one... */
  for ( nListIndex = 0; nListIndex < nNumPoints; nListIndex++ ) {
    
    /* get the point. */
    eList = VList_GetNthVoxel( pList, nListIndex, pControlPt );
    if ( kVListErr_NoErr != eList ) {
      eResult = DspA_tErr_ErrorAccessingControlPoints;
      goto error;
    }

    /* see if it's selected. */
    eList = VList_IsInList( this->mpSelectedControlPoints, 
          pControlPt, &bSelected );
    if ( kVListErr_NoErr != eList ) {
      eResult = DspA_tErr_ErrorAccessingControlPoints;
      goto error;
    }

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
    eResult = DspA_ConvertVolumeToBuffer_ ( this, pControlPt, &bufferPt );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* draw a point here. */
    DspA_DrawCrosshairIntoFrame_( this, faColor, &bufferPt, 
          DspA_knControlPointCrosshairSize );
  }

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawControlPointsToFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  /* delete the voxel */
  Voxel_Delete( &pControlPt );

  return eResult;
}

DspA_tErr DspA_DrawSelectionToFrame_ ( tkmDisplayAreaRef this ) {

  DspA_tErr     eResult           = DspA_tErr_NoErr;
  VoxelListRef  pList             = NULL;
  char          eList             = 0;
  char          eSpace            = 0;
  int           nNumPoints        = 0;
  int           nListIndex        = 0;
  VoxelRef      pSelectedPt       = NULL;
  xPoint2n      bufferPt          = {0,0};
  GLubyte*      pFrame            = NULL;
  int           nValue            = 0;
  int           nIntensifiedValue = 0;
  int           nX                = 0;
  int           nY                = 0;
  
  /* new voxel */
  Voxel_New( &pSelectedPt );

  /* decide which list we want out of the space. */
  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    eSpace = VSpace_GetVoxelsInZPlane( this->mpSelection, 
               DspA_GetCurrentSliceNumber_(this),
               &pList );
    break;
  case tkm_tOrientation_Sagittal:
    eSpace = VSpace_GetVoxelsInXPlane( this->mpSelection, 
               DspA_GetCurrentSliceNumber_(this), 
               &pList );
    break;
  case tkm_tOrientation_Horizontal:
    eSpace = VSpace_GetVoxelsInYPlane( this->mpSelection, 
               DspA_GetCurrentSliceNumber_(this),
               &pList );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* check for error. */
  if ( kVSpaceErr_NoErr != eSpace ) {
    DebugPrint "Error in VSpace_GetVoxelsInXYZPlane: %s\n",
      VSpace_GetErrorString( eSpace ) EndDebugPrint;
    eResult = DspA_tErr_ErrorAccessingSelection;
    goto error;
  }

  /* get the number of points. */
  eList = VList_GetCount( pList, &nNumPoints );
  if ( kVListErr_NoErr != eList ) {
    DebugPrint "Error in VSpace_GetCount: %s\n",
      VSpace_GetErrorString( eSpace ) EndDebugPrint;
    eResult = DspA_tErr_ErrorAccessingSelection;
    goto error;
  }

  /* for each one... */
  for ( nListIndex = 0; nListIndex < nNumPoints; nListIndex++ ) {

    /* get the point. */
    eList = VList_GetNthVoxel( pList, nListIndex, pSelectedPt );
    if ( kVListErr_NoErr != eList ) {
      DebugPrint "Error in VSpace_GetNthVoxel(%d): %s\n",
  nListIndex, VSpace_GetErrorString( eSpace ) EndDebugPrint;
      eResult = DspA_tErr_ErrorAccessingSelection;
      goto error;
    }

    /* get the value at this point. */
    nValue = tkm_GetVolumeValue( this->mpVolume, pSelectedPt );
    
    /* calculate the intensified green component. */
    nIntensifiedValue = nValue + DspA_knSelectionIntensityIncrement;
    if ( nIntensifiedValue > DspA_knMaxPixelValue ) {
      nIntensifiedValue = DspA_knMaxPixelValue;
    }

    /* lower the other components if they are too high. */
    if ( nIntensifiedValue - nValue <= DspA_knMinSelectionIntensityDiff ) {  
  nValue -= DspA_knSelectionIntensityIncrement;
  if ( nValue < 0 ) 
    nValue = 0;
    }

    /* covert to a buffer point. */
    eResult = DspA_ConvertVolumeToBuffer_ ( this, pSelectedPt, &bufferPt );
    if ( DspA_tErr_NoErr != eResult )
      goto error;
    
    /* y flip the volume pt to flip the image over. */
    bufferPt.mnY = BUFFER_Y_FLIP(bufferPt.mnY);

    /* write it back to the buffer. */
    for( nY = bufferPt.mnY; nY < bufferPt.mnY + this->mnZoomLevel; nY++ ) {
      for( nX = bufferPt.mnX; nX < bufferPt.mnX + this->mnZoomLevel; nX++ ) {

  pFrame = this->mpFrameBuffer + ( (nY * this->mnVolumeSize) + nX ) * 
    DspA_knNumBytesPerPixel;
  pFrame[DspA_knRedPixelCompIndex]   = (unsigned char)nValue;
  pFrame[DspA_knGreenPixelCompIndex] = (unsigned char)nIntensifiedValue;
  pFrame[DspA_knBluePixelCompIndex]  = (unsigned char)nValue;
  pFrame[DspA_knAlphaPixelCompIndex] = DspA_knMaxPixelValue;
      }
    }
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawSelectionToFrame_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  /* delete the voxel */
  Voxel_Delete( &pSelectedPt );

  return eResult;
}

DspA_tErr DspA_BuildSurfaceDrawLists_ ( tkmDisplayAreaRef this ) {

  DspA_tErr        eResult         = DspA_tErr_NoErr;
  xListRef         pList           = NULL;
  xList_tErr       eList           = xList_tErr_NoErr;
  xPoint2n         volumePt        = {0,0};
  VoxelRef         pAnatomicalVox  = NULL;
  VoxelRef         pRASVox         = NULL;
  Real             rRASX           = 0;
  Real             rRASY           = 0;
  Real             rRASZ           = 0;
  MHBT*            pBucket         = NULL;
  int              nBin            = 0;
  int              nSlice          = 0;
  tkm_tSurfaceType surface         = tkm_tSurfaceType_Current;
  int              nFace           = 0;
  int              nVertex         = 0;
  int              nNeighbor       = 0;
  face_type*       pFace           = NULL;
  vertex_type*     pVertex         = NULL;
  vertex_type*     pNeighbor       = NULL;
  VoxelRef         pVertexVox      = NULL;
  VoxelRef         pNeighborVox    = NULL;
  xPoint2f         intersectionPt  = {0, 0};
  tkmPointListNodeRef pPointNode = NULL;
  tBoolean         bPointsInThisFace = FALSE;

  Voxel_New( &pAnatomicalVox );
  Voxel_New( &pRASVox );
  Voxel_New( &pVertexVox );
  Voxel_New( &pNeighborVox );

  /* get the current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* set up for gl drawing */
  DspA_SetUpOpenGLPort_( this );

  /* for each surface type... */
  for ( surface = tkm_tSurfaceType_Current;
  surface < tkm_knNumSurfaceTypes; surface++ ) {
    
    /* only build if this surface is being displayed */
    if( !this->mabDisplayFlags[surface + DspA_tDisplayFlag_MainSurface] ) {
      continue;
    }


    /* check to see if we already have this list. if so, don't build it. */
    pList = DspA_GetSurfaceList_( this, this->mOrientation, surface,
          DspA_GetCurrentSliceNumber_(this) );
    if( NULL != pList ) {
      continue;
    }
    
    /* make a new list. */
    DspA_NewSurfaceList_( this, this->mOrientation, surface, 
        DspA_GetCurrentSliceNumber_(this) );
    
    /* get the list. */
    pList = DspA_GetSurfaceList_( this, this->mOrientation, surface,
          DspA_GetCurrentSliceNumber_(this) );
    if( NULL == pList ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }

    /* create a point list node */
    pPointNode = (tkmPointListNodeRef) malloc (sizeof(tkmPointListNode));
    if( NULL == pPointNode ) {
      eResult = DspA_tErr_OutOfMemory;
      goto cleanup;
    }
    pPointNode->mnNumPoints = 0;

    /* insert it into the list. */
    eList = xList_InsertItem ( pList, pPointNode );
    if ( xList_tErr_NoErr != eList ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto cleanup;
    }

#if 0
    DebugPrint "Building list %d for surface %d... ",
      nSlice, (int)surface EndDebugPrint;

    /* for each point in a slice of the surface.. */
    for( volumePt.mnY = 0; volumePt.mnY < this->mnVolumeSize; volumePt.mnY++ ){
      for( volumePt.mnX = 0; volumePt.mnX < this->mnVolumeSize; volumePt.mnX++ ){
    
  /* make a volume voxel. */
  eResult = DspA_ConvertBufferToVolume_ ( this, &volumePt, 
            pAnatomicalVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* trans to ras */
  trans_VoxelIndexToRAS( EXPAND_VOXEL_INT(pAnatomicalVox), 
             &rRASX, &rRASY, &rRASZ );
  Voxel_Set( pRASVox, rRASX, rRASY, rRASZ );


  /* get the bucket at this point */
  pBucket = MHTgetBucket( this->mpSurfaceHashTable[surface],
        EXPAND_VOXEL_FLOAT(pRASVox) );

  /* if we got one... */
  if( NULL != pBucket ) {

    /* for each bin... */
    for( nBin = 0; nBin < pBucket->nused; nBin++ ) {

      /* get the face number. */
      nFace = (pBucket->bins[nBin]).fno;
#endif
      
    /* for each face... */
    for ( nFace = 0; nFace < this->mpSurface->nfaces; nFace++ ) {
        
      /* get the face. */
      pFace = &(this->mpSurface->faces[nFace]);
      if( NULL == pFace ) {
  eResult = DspA_tErr_ErrorAccessingSurfaceList;
  DebugPrint "Face %d was null!\n", nFace EndDebugPrint;
  goto error;
      }
      
      /* haven't drawn any points in this face */
      bPointsInThisFace = FALSE;
      
      /* get the last neighrboring vertex. */
      nNeighbor = VERTICES_PER_FACE - 1;
      
      /* for each point in the face.. */
      for ( nVertex = 0; nVertex < VERTICES_PER_FACE; nVertex++ ) {
  
  /* get the vertices. */
  pVertex   = &(this->mpSurface->vertices[pFace->v[nVertex]]);
  pNeighbor = &(this->mpSurface->vertices[pFace->v[nNeighbor]]);
  
  if( NULL == pVertex ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    DebugPrint "Vertex %d was null!\n", pFace->v[nVertex] EndDebugPrint;
    goto error;
  }
  
  if( NULL == pNeighbor ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    DebugPrint "Vertex %d was null!\n", pFace->v[nNeighbor] EndDebugPrint;
    goto error;
  }
  
  /* convert them to normalized voxels. */
  xUtil_NormalizeVertexToVoxel( pVertex, surface, 
              this->mOrientation, pVertexVox );
  xUtil_NormalizeVertexToVoxel( pNeighbor, surface, 
              this->mOrientation, pNeighborVox );
  
  /* if the line between these two points intersects the
     current plane... */
  if ( xUtil_LineIntersectsPlane( pNeighborVox, pVertexVox, nSlice, 
       this->mabDisplayFlags[DspA_tDisplayFlag_InterpolateSurfaceVertices],
          &intersectionPt ) ) {
    
    /* add this point to the list node. */
    pPointNode->mafPoints[pPointNode->mnNumPoints][0] = 
      intersectionPt.mfX;
    pPointNode->mafPoints[pPointNode->mnNumPoints][1] = 
      intersectionPt.mfY;
    pPointNode->mnNumPoints++;
    
    /* if we are at the limit of points... */
    if( pPointNode->mnNumPoints >= 
        DspA_knMaxPointsPerPointListNode ) {
      
      /* create a new node. */
      pPointNode = 
        (tkmPointListNodeRef) malloc (sizeof(tkmPointListNode));
      if( NULL == pPointNode ) {
        eResult = DspA_tErr_OutOfMemory;
        goto cleanup;
      }
      
      pPointNode->mnNumPoints = 0;
      
      /* insert it into the list. */
      eList = xList_InsertItem ( pList, pPointNode );
      if ( xList_tErr_NoErr != eList ) {
        eResult = DspA_tErr_ErrorAccessingSurfaceList;
        goto cleanup;
      }
    }
    
    /* we're drawing a point this face, so the our flag to true. */
    bPointsInThisFace = TRUE;
  }
  
  /* use this vertex as the next neighbor. */
  nNeighbor = nVertex;      
      }
    
      /* if we drew any points for this face... */
      if( bPointsInThisFace ) {
  
  /* add a point -1 -1 to signify end of face. */
  pPointNode->mafPoints[pPointNode->mnNumPoints][0] = -1.0;
  pPointNode->mafPoints[pPointNode->mnNumPoints][1] = -1.0;
  pPointNode->mnNumPoints++;
  
  /* if we are at the limit of points... */
  if( pPointNode->mnNumPoints >= DspA_knMaxPointsPerPointListNode ) {
    
    /* create a new node. */
    pPointNode = 
      (tkmPointListNodeRef) malloc (sizeof(tkmPointListNode));
    if( NULL == pPointNode ) {
      eResult = DspA_tErr_OutOfMemory;
      goto cleanup;
    }
    
    pPointNode->mnNumPoints = 0;
    
    /* insert it into the list. */
    eList = xList_InsertItem ( pList, pPointNode );
    if ( xList_tErr_NoErr != eList ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto cleanup;
    }
  }
      } /* if points in this face.. */
    } /* for each face... */
#if 0
  } /* if bucket is not null... */
      } /* x */
    } /* y */
    DebugPrint " done.\n" EndDebugPrint;
#endif
   
  } /* for each surface... */

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_BuildSurfaceDrawLists__: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  Voxel_Delete( &pAnatomicalVox );
  Voxel_Delete( &pRASVox );
  Voxel_Delete( &pVertexVox );
  Voxel_Delete( &pNeighborVox );

  return eResult;
}

DspA_tErr DspA_DrawSurfaceDirect_ ( tkmDisplayAreaRef this ) {

  DspA_tErr        eResult         = DspA_tErr_NoErr;
  int              nSlice          = 0;
  tkm_tSurfaceType surface         = tkm_tSurfaceType_Current;
  int              nFace           = 0;
  int              nVertex         = 0;
  int              nNeighbor       = 0;
  face_type*       pFace           = NULL;
  vertex_type*     pVertex         = NULL;
  vertex_type*     pNeighbor       = NULL;
  VoxelRef         pVertexVox      = NULL;
  VoxelRef         pNeighborVox    = NULL;
  xPoint2f         intersectionPt  = {0, 0};
  xPoint2n         drawPt          = {0,0};
  float            faColor[3]      = {0,0,0};
  

  Voxel_New( &pVertexVox );
  Voxel_New( &pNeighborVox );

  /* get the current slice. */
  nSlice = DspA_GetCurrentSliceNumber_( this );

  /* set up for gl drawing */
  DspA_SetUpOpenGLPort_( this );

  /* for each surface type... */
  for ( surface = tkm_tSurfaceType_Current;
  surface < tkm_knNumSurfaceTypes; surface++ ) {
    
    /* if this surface is visible... */
    if( this->mabDisplayFlags[surface+DspA_tDisplayFlag_MainSurface] ) {

      /* choose and set the color */
      switch( surface ) {
      case tkm_tSurfaceType_Current:
  faColor[0] = 1.0; faColor[1] = 1.0; faColor[2] = 0.0;
  break;
      case tkm_tSurfaceType_Original:
  faColor[0] = 0.0; faColor[1] = 1.0; faColor[2] = 0.0;
  break;
      case tkm_tSurfaceType_Canonical:
  faColor[0] = 1.0; faColor[1] = 0.0; faColor[2] = 0.0;
  break;
      default:
  faColor[0] = 0.5; faColor[1] = 0.5; faColor[2] = 0.5;
  break;
      }
      
      glColor3fv( faColor );
      
      /* for each face... */
      for ( nFace = 0; nFace < this->mpSurface->nfaces; nFace++ ) {
  
  /* get the face. */
  pFace = &(this->mpSurface->faces[nFace]);
  
  if( NULL == pFace ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }
  
  /* start a line */
  glBegin( GL_LINES );
  
  /* get the last neighrboring vertex. */
  nNeighbor = VERTICES_PER_FACE - 1;
  
  /* for each point in the face.. */
  for ( nVertex = 0; nVertex < VERTICES_PER_FACE; nVertex++ ) {
  
    /* get the vertices. */
    pVertex   = &(this->mpSurface->vertices[pFace->v[nVertex]]);
    pNeighbor = &(this->mpSurface->vertices[pFace->v[nNeighbor]]);
    
    if( NULL == pVertex ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }
    
    if( NULL == pNeighbor ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }
    
    /* convert them to normalized voxels. */
    xUtil_NormalizeVertexToVoxel( pVertex, surface, 
          this->mOrientation, pVertexVox );
    xUtil_NormalizeVertexToVoxel( pNeighbor, surface, 
          this->mOrientation, pNeighborVox );
    
    /* if the line between these two points intersects the
       current plane... */
    if ( xUtil_LineIntersectsPlane( pNeighborVox, pVertexVox, nSlice, 
            this->mabDisplayFlags[DspA_tDisplayFlag_InterpolateSurfaceVertices],
            &intersectionPt ) ) {
      
      /* adjust the point */
      //DspA_AdjustSurfaceDrawPoint_( this, &intersectionPt );

      /* convert to zoomed coords. */
      drawPt.mnX = (int) ((float)this->mnZoomLevel * (intersectionPt.mfX - Voxel_GetFloatX(this->mpZoomCenter))) + (float)(this->mnVolumeSize/2);
      drawPt.mnY = (int) ((float)this->mnZoomLevel * (intersectionPt.mfY - Voxel_GetFloatY(this->mpZoomCenter))) + (float)(this->mnVolumeSize/2);
      
      /* y flip */
      drawPt.mnY = GLDRAW_Y_FLIP(drawPt.mnY);

      /* and draw the pt. */
      glVertex2i( drawPt.mnX, drawPt.mnY );
    }
    
    /* use this vertex as the next neighbor. */
    nNeighbor = nVertex;      
  }
  
  /* end the line. */
  glEnd();
      
      }
    }    
  }

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DrawSurfaceDirect_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  Voxel_Delete( &pVertexVox );
  Voxel_Delete( &pNeighborVox );

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
         VoxelRef          opCursor ) {

  DspA_tErr eResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = DspA_Verify ( this );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* return the cursor */
  Voxel_Copy( opCursor, this->mpCursor );

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
        tkm_tOrientation* oOrientation ) {

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

int DspA_GetCurrentSliceNumber_ ( tkmDisplayAreaRef this ) {

  int nSlice = 0;

  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    nSlice = Voxel_GetZ( this->mpCursor );
    break;
  case tkm_tOrientation_Horizontal:
    nSlice = Voxel_GetY( this->mpCursor );
    break;
  case tkm_tOrientation_Sagittal:
    nSlice = Voxel_GetX( this->mpCursor );
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
  this->maSurfaceLists = (xListRef*) malloc (sizeof(xListRef) * inNumLists );

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
  xList_tErr          eList      = xList_tErr_NoErr;
  tkmPointListNodeRef pPointNode = NULL;

  for( nDrawList = DspA_GetNumSurfaceLists_( this ) - 1; 
       nDrawList >= 0;
       nDrawList-- ) {

    /* if this is a list... */
    if( NULL != this->maSurfaceLists[nDrawList] ) {
      
      /* pop each item and delete the point node. */
      eList = xList_tErr_NoErr;
      while( eList != xList_tErr_EndOfList ) {

  eList = xList_PopItem( this->maSurfaceLists[nDrawList],
             (void**)&pPointNode );
  if( xList_tErr_NoErr != eList
      && xList_tErr_EndOfList != eList ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }
  
  if( NULL != pPointNode )
    free( pPointNode );
      }
      
      /* delete the list. */
      eList = xList_Delete( &this->maSurfaceLists[nDrawList] );
      if( xList_tErr_NoErr != eList ) {
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
         tkm_tOrientation  iOrientation,
         tkm_tSurfaceType  iSurface,
         int               inSlice ) {
  
  DspA_tErr   eResult   = DspA_tErr_NoErr;
  int         nDrawList = 0;
  xList_tErr  eList     = xList_tErr_NoErr;

  /* get the list index. */
  nDrawList = DspA_GetSurfaceListIndex_( this, iOrientation, 
           iSurface, inSlice );

  /* allocate a list. */
  xList_New( &this->maSurfaceLists[nDrawList] );
  if( xList_tErr_NoErr != eList ) {
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

xListRef DspA_GetSurfaceList_ ( tkmDisplayAreaRef this,
        tkm_tOrientation  iOrientation,
        tkm_tSurfaceType  iSurface,
        int               inSlice ) {

  xListRef pList = NULL;

  if( NULL != this->maSurfaceLists ) {
    pList = this->maSurfaceLists 
     [ DspA_GetSurfaceListIndex_( this, iOrientation, iSurface, inSlice ) ];
  }

  return pList;
}

int DspA_GetNumSurfaceLists_ ( tkmDisplayAreaRef this ) {

  return this->mnVolumeSize * tkm_knNumSurfaceTypes * tkm_knNumOrientations;
}

int DspA_GetSurfaceListIndex_ ( tkmDisplayAreaRef this,
        tkm_tOrientation  iOrientation,
        tkm_tSurfaceType  iSurface,
        int               inSlice ) {

  return ((int)iOrientation * tkm_knNumSurfaceTypes * this->mnVolumeSize) +
   ((int)iSurface * this->mnVolumeSize) + inSlice;
}




DspA_tErr DspA_ConvertVolumeToBuffer_ ( tkmDisplayAreaRef this,
          VoxelRef          ipVolumeVox,
          xPoint2nRef       opBufferPt ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fVolumeSize = (float) this->mnVolumeSize;

  /* verify the voxel */
  eResult = DspA_VerifyVolumeVoxel_( this, ipVolumeVox );
  if ( DspA_tErr_NoErr != eResult )
    goto error;

  /* first extract the points in the voxel that are on the same place as
     our orientation. */
  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    fX = Voxel_GetFloatX( ipVolumeVox );
    fY = Voxel_GetFloatY( ipVolumeVox );
    break;
  case tkm_tOrientation_Horizontal:
    fX = Voxel_GetFloatX( ipVolumeVox );
    fY = Voxel_GetFloatZ( ipVolumeVox );
    break;
  case tkm_tOrientation_Sagittal:
    fX = Voxel_GetFloatZ( ipVolumeVox );
    fY = Voxel_GetFloatY( ipVolumeVox );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* now zoom the coords to our zoomed buffer state */
  fX = (fZoomLevel * (fX - Voxel_GetFloatX(this->mpZoomCenter))) +
   (fVolumeSize/2.0);
  fY = (fZoomLevel * (fY - Voxel_GetFloatY(this->mpZoomCenter))) +
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
          VoxelRef          opVolumeVox ) {

  DspA_tErr eResult     = DspA_tErr_NoErr;
  float     fX          = 0;
  float     fY          = 0;
  VoxelRef  pVolumeVox  = NULL;
  float     fZoomLevel  = (float) this->mnZoomLevel;
  float     fBufferX    = (float) ipBufferPt->mnX;
  float     fBufferY    = (float) ipBufferPt->mnY;
  float     fVolumeSize = (float) this->mnVolumeSize;

  Voxel_New ( &pVolumeVox );

  /* unzoom the coords */
  fX = fBufferX / fZoomLevel + 
    ( Voxel_GetFloatX(this->mpZoomCenter) - (fVolumeSize/2.0/fZoomLevel) );
  fY = fBufferY / fZoomLevel + 
    ( Voxel_GetFloatY(this->mpZoomCenter) - (fVolumeSize/2.0/fZoomLevel) );

  /* build a voxel out of these two coords and the current slice */
  switch ( this->mOrientation ) {
  case tkm_tOrientation_Coronal:
    Voxel_SetFloat( pVolumeVox, fX, fY, DspA_GetCurrentSliceNumber_( this ) );
    break;
  case tkm_tOrientation_Horizontal:
    Voxel_SetFloat( pVolumeVox, fX, DspA_GetCurrentSliceNumber_( this ), fY );
    break;
  case tkm_tOrientation_Sagittal:
    Voxel_SetFloat( pVolumeVox, DspA_GetCurrentSliceNumber_( this ), fY, fX );
    break;
  default:
    eResult = DspA_tErr_InvalidOrientation;
    goto error;
    break;
  }

  /* copy into the outgoing voxel to return it */
  Voxel_Copy( opVolumeVox, pVolumeVox );

  goto cleanup;

 error:

  /* print error message */
  if ( DspA_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in DspA_DspA_ConvertBufferToVolume_: %s\n",
      eResult, DspA_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  Voxel_Delete ( &pVolumeVox );

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
  VoxelRef              pVoxel             = NULL;
  unsigned char         ucVolumeValue      = 0;
  FunV_tErr             eFunctional        = FunV_tErr_NoError;
  FunV_tFunctionalValue functionalValue    = 0;
  int                   nFlag              = 0;

  Voxel_New( &pVoxel );

  /* send the cursor. */
  sprintf( sTclArguments, "%d %d %d",
     EXPAND_VOXEL_INT( this->mpCursor ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeCursor, sTclArguments );

  /* also convert to RAS and send those coords along. */
  tkm_ConvertVolumeToRAS( this->mpCursor, pVoxel );
  sprintf( sTclArguments, "%.1f %.1f %.1f", EXPAND_VOXEL_FLOAT( pVoxel ) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateRASCursor, sTclArguments );
 
  /* also convert to Tal and send those coords along. */
  tkm_ConvertVolumeToTal( this->mpCursor, pVoxel );
  sprintf( sTclArguments, "%.1f %.1f %.1f", EXPAND_VOXEL_FLOAT( pVoxel ) );
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

  /* also see if we have functional data and can send a value for that
     as well. */     
  if ( NULL != this->mpFunctionalVolume ) {
    eFunctional = FunV_GetValueAtAnatomicalVoxel( this->mpFunctionalVolume,
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

  Voxel_Delete( &pVoxel );

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
             VoxelRef          ipVoxel ) {
  DspA_tErr eResult = DspA_tErr_NoErr;

  /* make sure voxel is in bounds */
  if ( Voxel_GetX( ipVoxel )    < 0
       || Voxel_GetY( ipVoxel ) < 0
       || Voxel_GetZ( ipVoxel ) < 0
       || Voxel_GetX( ipVoxel ) >= this->mnVolumeSize
       || Voxel_GetY( ipVoxel ) >= this->mnVolumeSize
       || Voxel_GetZ( ipVoxel ) >= this->mnVolumeSize ) {

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
    this->mnZoomLevel, EXPAND_VOXEL_INT(this->mpZoomCenter) EndDebugPrint;
  DebugPrint "\tcursor %d, %d, %d orientation %s slice %d tool %d\n",
    EXPAND_VOXEL_INT(this->mpCursor), DspA_ksaOrientation[this->mOrientation],
    DspA_GetCurrentSliceNumber_(this), sTool EndDebugPrint;
}

void DspA_Signal ( char* isFuncName, int inLineNum, DspA_tErr ieCode ) {

  DebugPrint "Signal in %s, line %d: %d, %s", 
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
            tkm_tSurfaceType iSurface,
            tkm_tOrientation iOrientation ) {

  int          nVertexA    = 0;
  int          nVertexB    = 0;
  vertex_type* pVertexA    = NULL;
  vertex_type* pVertexB    = NULL;
  VoxelRef     pVoxelA     = NULL;
  VoxelRef     pVoxelB     = NULL;
  float        fDistanceA  = 0;
  float        fDistanceB  = 0;
  tBoolean     bIntersects = FALSE;

  Voxel_New( &pVoxelA );
  Voxel_New( &pVoxelB );

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
    fDistanceA = Voxel_GetZ( pVoxelA ) - inPlane;
    fDistanceB = Voxel_GetZ( pVoxelB ) - inPlane;

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

  Voxel_Delete( &pVoxelA );
  Voxel_Delete( &pVoxelB );

  return bIntersects;
}

void xUtil_NormalizeVertexToVoxel( vertex_type*     ipVertex,
           tkm_tSurfaceType iSurface,
           tkm_tOrientation iOrientation,
           VoxelRef         opVoxel ) {

  Real rXVox = 0;
  Real rYVox = 0;
  Real rZVox = 0;
  
  switch( iSurface ) {
  case tkm_tSurfaceType_Current:
    trans_RASToVoxel( ipVertex->x, ipVertex->y, ipVertex->z,
          &rXVox,      &rYVox,      &rZVox );
    break;
  case tkm_tSurfaceType_Canonical:
    trans_RASToVoxel( ipVertex->cx, ipVertex->cy, ipVertex->cz,
          &rXVox,       &rYVox,       &rZVox );
    break;
  case tkm_tSurfaceType_Original:
    trans_RASToVoxel( ipVertex->origx, ipVertex->origy, ipVertex->origz,
          &rXVox,          &rYVox,          &rZVox );
    break;
  default:
    break;
  }

  switch( iOrientation ) {
  case tkm_tOrientation_Horizontal:
    Voxel_SetFloat( opVoxel, (float)rXVox, (float)rZVox, (float)rYVox );
    break;
  case tkm_tOrientation_Coronal:
    Voxel_SetFloat( opVoxel, (float)rXVox, (float)rYVox, (float)rZVox );
    break;
  case tkm_tOrientation_Sagittal:
    Voxel_SetFloat( opVoxel, (float)rZVox, (float)rYVox, (float)rXVox );
    break;
  default:
    break;
  }
}

tBoolean xUtil_LineIntersectsPlane( VoxelRef         ipLineVoxA,
            VoxelRef         ipLineVoxB,
            int              inPlane,
            tBoolean         ipInterpolate,
            xPoint2fRef      opIntersectionPt ) {
  
  float    fPlane           = inPlane;
  float    fDistanceA       = 0;
  float    fDistanceB       = 0;
  float    fAlpha           = 1;
  xPoint2f intersectionPt   = {0, 0};

  /* get distance from each to plane. */
  fDistanceA = Voxel_GetFloatZ( ipLineVoxA ) - fPlane;
  fDistanceB = Voxel_GetFloatZ( ipLineVoxB ) - fPlane;

  /* if product is negative or 0, they intersect the plane. */
  if ( fDistanceA * fDistanceB > 0.0 ) {
    return FALSE;
  }

  /* if averaging the vertices, find an iterpolation factor, which is the
     intersection of the line the plane. */
  if( ipInterpolate ) {
    
    /* make sure they arn't on the same plane... */
    if( Voxel_GetFloatZ(ipLineVoxB) - Voxel_GetFloatZ(ipLineVoxA) != 0.0 ) {

      fAlpha = (fPlane - Voxel_GetFloatZ( ipLineVoxA )) /
  (Voxel_GetFloatZ( ipLineVoxB ) - Voxel_GetFloatZ( ipLineVoxA ));
      
    } else {

      /* alpha is just 1. */
      fAlpha = 1.0;
    }

    /* interpolate to find the intersection. */
    intersectionPt.mfX = (Voxel_GetFloatX( ipLineVoxA ) +
      fAlpha * (Voxel_GetFloatX(ipLineVoxB) - Voxel_GetFloatX(ipLineVoxA)));
    intersectionPt.mfY = (Voxel_GetFloatY( ipLineVoxA ) +
      fAlpha * (Voxel_GetFloatY(ipLineVoxB) - Voxel_GetFloatY(ipLineVoxA)));

  } else {
    
    /* if no interpolationg, intersection is projection onto plane */
    intersectionPt.mfX = Voxel_GetFloatX( ipLineVoxB );
    intersectionPt.mfY = Voxel_GetFloatY( ipLineVoxB );
  }

  /*
  DebugPrint "intersecting %.2f,%.2f,%.2f to %.2f,%.2f,%.2f on %d",
    EXPAND_VOXEL_FLOAT( ipLineVoxA ), EXPAND_VOXEL_FLOAT( ipLineVoxB ),
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
  case tkm_tOrientation_Horizontal:
    ipPoint->mfX += 0.5;
    ipPoint->mfY += 0.5;
    break;
  case tkm_tOrientation_Coronal:
    ipPoint->mfX += 0.5;
    ipPoint->mfY += 0.5;
    break;
  case tkm_tOrientation_Sagittal:
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
        xListRef          ipList ) {

  DspA_tErr        eResult   = DspA_tErr_NoErr;
  tBoolean bOperationOpen = FALSE;
  xPoint2f drawPoint      = {0,0};
  int              nNodePoint = 0;
  tkmPointListNodeRef pPointNode = NULL;
  xList_tErr       eList     = xList_tErr_NoErr;

  /* reset the list position. */
  eList = xList_ResetPosition ( ipList );
  if ( xList_tErr_NoErr != eList ) {
    eResult = DspA_tErr_ErrorAccessingSurfaceList;
    goto error;
  }
  
  /* start new operation. */
  glBegin( inMode );
  bOperationOpen = TRUE;

  pPointNode = NULL;
  eList = xList_tErr_NoErr;
  while ( eList != xList_tErr_EndOfList ) {
    
    /* get the next point node in the list. */
    eList = xList_GetNextItemFromPosition ( ipList, (void**)&pPointNode );
    if ( xList_tErr_NoErr != eList
   && xList_tErr_EndOfList != eList ) {
      eResult = DspA_tErr_ErrorAccessingSurfaceList;
      goto error;
    }
    
    /* if we got a node... */
    if ( NULL != pPointNode ) {
      
      /* for each point in the node... */
      for( nNodePoint = 0; 
     nNodePoint < pPointNode->mnNumPoints; 
     nNodePoint++ ) {
  
  /* if it is -1,-1, it's a face marker... */
  if( -1.0 == pPointNode->mafPoints[nNodePoint][0]
      && -1.0 == pPointNode->mafPoints[nNodePoint][1] ) {
    
    /* if operation was still going, end it. */ 
    if( bOperationOpen )
      glEnd();
    
    /* start new operation. */
    glBegin( inMode );
    bOperationOpen = TRUE;
    
  } else {
    
    drawPoint.mfX = pPointNode->mafPoints[nNodePoint][0];
    drawPoint.mfY = pPointNode->mafPoints[nNodePoint][1];

    /* adjust the point */
    DspA_AdjustSurfaceDrawPoint_( this, &drawPoint );

    /* convert to zoomed coords. */
    drawPoint.mfX = ((float)this->mnZoomLevel * (drawPoint.mfX - Voxel_GetFloatX(this->mpZoomCenter))) + (float)(this->mnVolumeSize/2.0);
    drawPoint.mfY = ((float)this->mnZoomLevel * (drawPoint.mfY - Voxel_GetFloatY(this->mpZoomCenter))) + (float)(this->mnVolumeSize/2.0);
    
    /* y flip */
    drawPoint.mfY = GLDRAW_Y_FLIP(drawPoint.mfY);
    
    /* and draw the pt. */
    glVertex2f( drawPoint.mfX, drawPoint.mfY );
  }
      }
    }
  }

  /* end last operation */
  if( bOperationOpen ) {
    glEnd();
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

  
