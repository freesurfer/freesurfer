#include "tkmMeditWindow.h"
#include "tkmDisplayArea.h"

char *MWin_ksaErrorStrings [MWin_knNumErrorCodes] = {

  "No error",
  "Allocation failed.",
  "Allocation of display area failed.",
  "Error accessing GLUT window.",
  "Invalid ptr to object (was probably NULL).",
  "Invalid signature found while verifying object.",
  "Invalid display index.",
  "Invalid display area.",
  "Invalid display configuration",
  "Invalid coordinate space.",
  "Error accessing display.",
  "Wrong number of arguments.",
  "Invalid event type.",
  "Invalid error code."
};

MWin_tErr MWin_New ( tkmMeditWindowRef* oppWindow,
		     char*              isTitle,
		     int                inWidth, 
		     int                inHeight ) {

  MWin_tErr         eResult     = MWin_tErr_NoErr;
  tkmMeditWindowRef this        = NULL;
  char              sTitle[256] = "";
  DspA_tErr         eDispResult = DspA_tErr_NoErr;
  int               nDisplay    = 0;

  /* allocate us. */
  this = (tkmMeditWindowRef) malloc ( sizeof(tkmMeditWindow) );
  if ( NULL == this ) {
    eResult = MWin_tErr_AllocationFailed;
    goto error;
  }

  /* set the signature */
  this->mSignature = MWin_kSignature;

  /* create the title. */
  sprintf ( sTitle, "Medit: %s", isTitle );

  /* create the window */
  xGWin_New ( &(this->mpGLutWindow), inWidth, inHeight, sTitle );
  xGWin_SetEventHandlerFunc ( this->mpGLutWindow, 
            MWin_EventCallback, this );
  xGWin_ActivateIdleEvents ( this->mpGLutWindow );
  xGWin_ActivatePassiveMotionEvents ( this->mpGLutWindow );

  /* set the size */
  this->mnWidth  = inWidth;
  this->mnHeight = inHeight;

  /* set cursor to crosshair */
  //  glutSetCursor( GLUT_CURSOR_CROSSHAIR );

  /* set display pointers to null */
  for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ )
    this->mapDisplays[nDisplay] = NULL;

  /* for each display */
  for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {

    /* create a new display */
    eDispResult = DspA_New ( &(this->mapDisplays[nDisplay]), this );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_DisplayAreaAllocationFailed;
      goto error;
    }

    DspA_SetID( this->mapDisplays[nDisplay], nDisplay );
  }

  /* default last clicked area. */
  this->mnLastClickedArea = 0;

  /* set the default configuration and position the displays */
  eResult =
    MWin_SetDisplayConfiguration ( this, 1, 1, MWin_tLinkPolicy_None );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* start out linking the cursor. */
  this->mbLinkedCursor = TRUE;

  /* not accepting tcl commands yet. */
  this->mbAcceptingTclCommands = FALSE;

  /* return window. */
  *oppWindow = this;

  goto cleanup;

 error:
  
  /* if we allocated some stuff... */
  if ( NULL != this ) {

    /* kill the displays if we allocated them */
    for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      if ( NULL != this->mapDisplays[nDisplay] )
  free ( this->mapDisplays[nDisplay] );
    }
    
    /* kill the main storage */
    free ( this );
  } 
   
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_New: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_Delete ( tkmMeditWindowRef* ioppWindow ) {
  
  MWin_tErr         eResult     = MWin_tErr_NoErr;
  tkmMeditWindowRef this        = NULL;
  DspA_tErr         eDispResult = DspA_tErr_NoErr;
  int               nDisplay    = 0;

  /* get us */
  this = *ioppWindow;
    
  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* for each display */
  for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {

    /* delete display */
    eDispResult = DspA_Delete ( &(this->mapDisplays[nDisplay]) );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_DisplayAreaAllocationFailed;
      goto error;
    }
  }

  /* delete the window */
  xGWin_Delete ( &(this->mpGLutWindow) );

  /* trash the signature */
  this->mSignature = 0x1;
 
  /* delete us */
  free ( this );
  
  /* return null */
  *ioppWindow = NULL;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_Delete: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetWindowTitle ( tkmMeditWindowRef this,
        char*             isTitle ) {

  MWin_tErr  eResult  = MWin_tErr_NoErr;
  xGWin_tErr eGWin    = xGWin_tErr_NoErr;

  DebugEnterFunction( ("MWin_SetWindowTitle( %p, isTitle=%s )",
           this, isTitle) );

  /* verify us. */
  eResult = MWin_Verify ( this );
  DebugAssertThrow( MWin_tErr_NoErr == eResult );

  /* set the gl window title */
  eGWin = xGWin_SetWindowTitle( this->mpGLutWindow, isTitle );
  DebugAssertThrowX( eGWin == xGWin_tErr_NoErr, eResult, 
         MWin_tErr_ErrorAccessingWindow );

  DebugCatch;
  DebugCatchError( eResult, MWin_tErr_NoErr, MWin_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

MWin_tErr MWin_GetWindowSize ( tkmMeditWindowRef this,
             int*              onX,
             int*              onY,
             int*              onWidth,
             int*              onHeight ) {

  MWin_tErr  eResult  = MWin_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* set stuff */
  *onX      = glutGet( GLUT_WINDOW_X );
  *onY      = glutGet( GLUT_WINDOW_Y );
  *onWidth  = this->mnWidth;
  *onHeight = this->mnHeight;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetWindowTitle: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetDisplayConfiguration ( tkmMeditWindowRef this,
           int               inCols, 
           int               inRows,
           MWin_tLinkPolicy  iPolicy ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  if( inCols < 0 || inRows < 0 ||
      inCols * inRows > MWin_knMaxNumAreas ) {
    eResult = MWin_tErr_InvalidDisplayConfiguration;
    goto error;
  }

  /* save the configuration */
  this->mnCols      = inCols;
  this->mnRows      = inRows;
  this->mLinkPolicy = iPolicy;

  /* position the subpanes */
  eResult = MWin_PositionDisplays_ ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if config is 1x1, turn focus frame off, otherwise turn it on. */
  eResult = MWin_SetDisplayFlag( this, -1, DspA_tDisplayFlag_FocusFrame,
         (this->mnCols == 1 && this->mnRows == 1)?
         FALSE:TRUE );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetDisplayConfiguration: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_PositionDisplays_ ( tkmMeditWindowRef this ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDisplay    = 0;
  int       nMosaicHackedIndex    = 0;
  xPoint2n  location    = {0, 0};
  int       nZoomLevel  = 0;
  xVoxel    cursor;
  int       nX          = 0;
  int       nY          = 0;
  int       nWidth      = 0;
  int       nHeight     = 0;
  mri_tOrientation orientation = mri_tOrientation_Coronal;

  /* if we're on 1x1, just make the focused one big and the rest
     little */
  if( this->mnCols == 1 && this->mnRows == 1 ) {
    
    for( nDisplay = nDisplay; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      location.mnX = location.mnY = 0;
      if( nDisplay == this->mnLastClickedArea ) {
	nWidth = this->mnWidth;
	nHeight = this->mnHeight;
      } else {
	nWidth = nHeight = 0;
      }
      eDispResult = DspA_SetPosition ( this->mapDisplays[nDisplay],
				       location, nWidth, nHeight );
      if ( DspA_tErr_NoErr != eDispResult ) {
	eResult = MWin_tErr_ErrorAccessingDisplay;
	goto error;
      }
    }
    
    goto cleanup;
  }
  
  /* get the zoom level, cursor, and orientation of the last clicked pane */
  eDispResult = DspA_GetZoomLevel( this->mapDisplays[this->mnLastClickedArea], 
				   &nZoomLevel );
  if ( DspA_tErr_NoErr != eDispResult )
    goto error;
  
  eDispResult = DspA_GetCursor( this->mapDisplays[this->mnLastClickedArea],
				&cursor );
  if ( DspA_tErr_NoErr != eDispResult )
    goto error;
  
  eDispResult =DspA_GetOrientation( this->mapDisplays[this->mnLastClickedArea],
				    &orientation );
  if ( DspA_tErr_NoErr != eDispResult )
    goto error;
  
  nDisplay = 0;
  for( nY = 0; nY < this->mnRows; nY++ ) {
    for( nX = 0; nX < this->mnCols; nX++ ) {
      
      /* if mult orienations, set different orientation on each */
      if( MWin_tLinkPolicy_MultipleOrientations == this->mLinkPolicy ) {
	
	/* set its position and size */
	location.mnX = (nX * this->mnWidth / this->mnCols);
	location.mnY = (nY * this->mnHeight / this->mnRows);
	eDispResult = DspA_SetPosition ( this->mapDisplays[nDisplay],
					 location, 
					 this->mnWidth / this->mnCols, 
					 this->mnHeight / this->mnRows);
	if ( DspA_tErr_NoErr != eDispResult )
	  goto error;
	
	/* set an orientation */
	eDispResult = 
	  DspA_SetOrientation ( this->mapDisplays[nDisplay], 
				nDisplay % mri_knNumOrientations );
	if ( DspA_tErr_NoErr != eDispResult )
	  goto error;
	
	/* set them all to the same cursor */
	eDispResult = DspA_SetCursor( this->mapDisplays[nDisplay], &cursor );
	if ( DspA_tErr_NoErr != eDispResult )
	  goto error;
	
	
      } else if( MWin_tLinkPolicy_Mosaic == this->mLinkPolicy ) {
	
  /* this is a yflip */
  nMosaicHackedIndex = ((this->mnRows - 1) - nY) * this->mnCols + nX;

  /* set its position and size */
  location.mnX = (nX * this->mnWidth / this->mnCols);
  location.mnY = (nY * this->mnHeight / this->mnRows);
  eDispResult = DspA_SetPosition ( this->mapDisplays[nMosaicHackedIndex],
           location, 
           this->mnWidth / this->mnCols, 
           this->mnHeight / this->mnRows);
  if ( DspA_tErr_NoErr != eDispResult )
    goto error;
      
  /* if mosaic, set same orientation and then set the slice */
  eDispResult = DspA_SetOrientation( this->mapDisplays[nMosaicHackedIndex], 
              orientation );
  if ( DspA_tErr_NoErr != eDispResult )
    goto error;

  eDispResult = DspA_SetSlice ( this->mapDisplays[nMosaicHackedIndex], 
              nDisplay *
              (255 / (this->mnRows*this->mnCols)) );
  if ( DspA_tErr_NoErr != eDispResult )
    goto error;
      }

      /* set zoom level to that of the last clicked one */
      eDispResult = DspA_SetZoomLevel( this->mapDisplays[nDisplay], 
               nZoomLevel );
      if ( DspA_tErr_NoErr != eDispResult )
    goto error;
      
      nDisplay++;
    }
  }
  
  /* set all other displays to 0 */
  for( nDisplay = nDisplay; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
    
    location.mnX = location.mnY = 0;
    eDispResult = DspA_SetPosition ( this->mapDisplays[nDisplay],
             location, 0, 0 );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }
  
  goto cleanup;
  
 error:
  
  if ( DspA_tErr_NoErr != eDispResult )
    eResult = MWin_tErr_ErrorAccessingDisplay;
 
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_PositionDisplays_: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetVolume ( tkmMeditWindowRef this,
                           int               inDispIndex,
                           mriVolumeRef      ipVolume,
                           int               inSizeX, 
                           int               inSizeY,
                           int               inSizeZ ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;
  
  DebugEnterFunction( ("MWin_SetVolume( this=%p, inDispIndex=%d, "
           "ipVolume=%p, inSize=%d,%d,%d )", this, inDispIndex,
           ipVolume, inSizeX, inSizeY, inSizeZ) );

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the volume */
  for ( nDispIndex = nDispIndexMin; 
	nDispIndex < nDispIndexMax; 
	nDispIndex++ ) {
    
    eDispResult = DspA_SetVolume ( this->mapDisplays[nDispIndex],
				   ipVolume, inSizeX, inSizeY, inSizeZ );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  DebugExitFunction;

  return eResult;
}

MWin_tErr MWin_SetAuxVolume ( tkmMeditWindowRef this,
                              int               inDispIndex,
                              mriVolumeRef      ipVolume,
                              int               inSizeX,
                              int               inSizeY,
                              int               inSizeZ ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the aux volume */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetAuxVolume ( this->mapDisplays[nDispIndex],
              ipVolume, inSizeX, inSizeY, inSizeZ );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetAuxVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetSegmentationVolume ( tkmMeditWindowRef this,
				       tkm_tSegType      iType,
				       int               inDispIndex,
				       mriVolumeRef      iVolume ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;
  
  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the volume */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = 
      DspA_SetSegmentationVolume ( this->mapDisplays[nDispIndex], 
				   iType, iVolume );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetSegmentationVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetSegmentationColorTable ( tkmMeditWindowRef this,
					   tkm_tSegType      iType,
					   int               inDispIndex,
					   mriColorLookupTableRef iCLUT ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;
  
  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the volume */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = 
      DspA_SetSegmentationColorTable ( this->mapDisplays[nDispIndex], 
				       iType, iCLUT );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetSegmentationVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetSurface ( tkmMeditWindowRef this, 
          int               inDispIndex,
          tkm_tSurfaceType  iType,
          mriSurfaceRef     ipSurface ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the surface */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {
    
    eDispResult = DspA_SetSurface ( this->mapDisplays[nDispIndex],
            iType, ipSurface );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetSurface: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetOverlayVolume ( tkmMeditWindowRef      this,
          int                    inDispIndex,
          tkmFunctionalVolumeRef ipVolume ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the overlay volume */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetOverlayVolume ( this->mapDisplays[nDispIndex],
            ipVolume );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetOverlayVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetControlPointsSpace ( tkmMeditWindowRef this,
               int               inDispIndex,
               x3DListRef        ipVoxels ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the control points */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetControlPointsSpace ( this->mapDisplays[nDispIndex],
                 ipVoxels );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetControlPointsSpace: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}


MWin_tErr MWin_SetSelectionSpace ( tkmMeditWindowRef this, 
				   int               inDispIndex,
				   mriVolumeRef      ipVolume ) {
  
  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the selection */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetSelectionSpace ( this->mapDisplays[nDispIndex],
					   ipVolume );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetSelectionSpace: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetHeadPointList ( tkmMeditWindowRef   this, 
          int                 inDispIndex,
          mriHeadPointListRef iList ) {
  
  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the list */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetHeadPointList ( this->mapDisplays[nDispIndex],
            iList );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetHeadPointList: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetVLIs  ( tkmMeditWindowRef this,
        int               inDispIndex,
        VLI*              iVLI1,
        VLI*              iVLI2,
        char*             isVLI1_name,
        char*             isVLI2_name ) {
  
  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the VLIs */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetVLIs ( this->mapDisplays[nDispIndex],
        iVLI1, iVLI2, isVLI1_name, isVLI2_name );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetVLIs: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetDTIVolume  ( tkmMeditWindowRef  this,
			       int                inDispIndex,
			       mriVolumeRef       iVolume ) {
  
  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the DTI volume */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetDTIVolume ( this->mapDisplays[nDispIndex], iVolume );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetDTIVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}


MWin_tErr MWin_SetGCA ( tkmMeditWindowRef   this, 
      int                 inDispIndex,
      GCA*                iVolume,
      TRANSFORM*          iTransform ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the gca */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetGCA ( this->mapDisplays[nDispIndex],
        iVolume, iTransform );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetGCA: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetLinkedCursorFlag ( tkmMeditWindowRef this, 
             tBoolean          ibLinkedCursor ) {

  MWin_tErr eResult            = MWin_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* set cursor link status. */
  this->mbLinkedCursor = ibLinkedCursor;

  /* send update to tk window. */
  sprintf( sTclArguments, "%d", (int)this->mbLinkedCursor );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateLinkedCursorFlag, sTclArguments );

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetLinkedCursorFlag: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_ToggleLinkedCursorFlag ( tkmMeditWindowRef this ) {

  MWin_tErr eResult            = MWin_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* set cursor link status. */
  eResult = MWin_SetLinkedCursorFlag( this, !(this->mbLinkedCursor) );
  if( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_ToggleLinkedCursorFlag: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetCursor ( tkmMeditWindowRef this, 
         int               inDispIndex,
        xVoxelRef          ipCursor ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the cursor */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetCursor ( this->mapDisplays[nDispIndex],
           ipCursor );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_ConvertAndSetCursor ( tkmMeditWindowRef this, 
             int               inDispIndex,
             mri_tCoordSpace   iFromSpace,
             xVoxelRef         ipCursor ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the cursor */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_ConvertAndSetCursor ( this->mapDisplays[nDispIndex],
               iFromSpace, ipCursor );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_ConvertAndSetCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}
MWin_tErr MWin_SetOrientation ( tkmMeditWindowRef this, 
        int               inDispIndex,
        mri_tOrientation  iOrientation ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;
  int       nSlice        = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the orientation */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    /* if in mosiac mode, save the slice */
    if( MWin_tLinkPolicy_Mosaic == this->mLinkPolicy ) {
      DspA_GetSlice( this->mapDisplays[nDispIndex], &nSlice );
    }

    eDispResult = DspA_SetOrientation ( this->mapDisplays[nDispIndex],
          iOrientation );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }

    /* if in mosiac mode, restore the slice */
    if( MWin_tLinkPolicy_Mosaic == this->mLinkPolicy ) {
      DspA_SetSlice( this->mapDisplays[nDispIndex], nSlice );
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetOrientation: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetSlice ( tkmMeditWindowRef this, 
        int               inDispIndex,
        int               inSlice ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the slice number */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetSlice ( this->mapDisplays[nDispIndex],
          inSlice );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetSlice: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetZoomCenter ( tkmMeditWindowRef this, 
             int               inDispIndex,
            xVoxelRef          ipCenter ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the center of the view */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetZoomCenter ( this->mapDisplays[nDispIndex],
               ipCenter );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetZoomCenter: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetZoomCenterToCursor ( tkmMeditWindowRef this,
               int               inDispIndex ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /*  set the zoom cetner to the cursor */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetZoomCenterToCursor ( this->mapDisplays[nDispIndex] );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetZoomCenterToCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetZoomLevel ( tkmMeditWindowRef this, 
            int               inDispIndex,
            int               inLevel ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /*  set the zoom cetner to the cursor */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetZoomLevel ( this->mapDisplays[nDispIndex], inLevel );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetZoomLevel: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_HiliteSurfaceVertex ( tkmMeditWindowRef this,
             int               inDispIndex,
             Surf_tVertexSet  inSurface,
             int               inVertex ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* hilite a vertex */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_HiliteSurfaceVertex ( this->mapDisplays[nDispIndex],
               inSurface, inVertex );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_HiliteSurfaceVertex: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetDisplayFlag ( tkmMeditWindowRef this,
        int               inDispIndex,
        DspA_tDisplayFlag iWhichFlag,
        tBoolean          ibNewValue ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set a display flag */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetDisplayFlag ( this->mapDisplays[nDispIndex],
          iWhichFlag, ibNewValue );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetDisplayFlag: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetSegmentationAlpha ( tkmMeditWindowRef this,
				      int               inDispIndex,
				      float             ifAlpha ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the alpha */
  for ( nDispIndex = nDispIndexMin; 
	nDispIndex < nDispIndexMax; 
	nDispIndex++ ) {
    
    eDispResult = DspA_SetSegmentationAlpha ( this->mapDisplays[nDispIndex],
					      ifAlpha );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetSegmentationAlpha: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetDTIAlpha ( tkmMeditWindowRef this,
				      int               inDispIndex,
				      float             ifAlpha ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the alpha */
  for ( nDispIndex = nDispIndexMin; 
	nDispIndex < nDispIndexMax; 
	nDispIndex++ ) {
    
    eDispResult = DspA_SetDTIAlpha ( this->mapDisplays[nDispIndex], ifAlpha );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetDTIAlpha: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetDTIAxisForComponent ( tkmMeditWindowRef this,
					int               inDispIndex,
					tkm_tAxis         iAxis,
					xColr_tComponent  iComponent ) {

  MWin_tErr eResult       = MWin_tErr_NoErr;
  DspA_tErr eDispResult   = DspA_tErr_NoErr;
  int       nDispIndex    = 0;
  int       nDispIndexMin = inDispIndex;
  int       nDispIndexMax = inDispIndex+1;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the display index. */
  eResult = MWin_VerifyDisplayIndex ( this, inDispIndex );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if working on all displays, set the iteration bounds. */
  if ( MWin_kAllDisplayAreas == inDispIndex ) {
    nDispIndexMin = 0;
    nDispIndexMax = MWin_knMaxNumAreas;
  }

  /* set the axis */
  for ( nDispIndex = nDispIndexMin; 
	nDispIndex < nDispIndexMax; 
	nDispIndex++ ) {
    
    eDispResult = 
      DspA_SetDTIAxisForComponent ( this->mapDisplays[nDispIndex], 
				    iAxis, iComponent );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SetDTIAxisForComponent: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_GetCursor ( tkmMeditWindowRef this,
        xVoxelRef          opCursor ) {
  
  MWin_tErr eResult      = MWin_tErr_NoErr;
  DspA_tErr eDispResult  = DspA_tErr_NoErr;
 xVoxelRef  pCursor      = NULL;

  /* new cursor. */
  xVoxl_New( &pCursor );

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* got the cursor from the last clicked display. */
  eDispResult = DspA_GetCursor ( this->mapDisplays[this->mnLastClickedArea],
         pCursor );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
  /* return the cursor */
  xVoxl_Copy( opCursor, pCursor );

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_GetCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  /* delete cursor */
  xVoxl_Delete( &pCursor );
  
  return eResult;
}

MWin_tErr MWin_GetOrientation ( tkmMeditWindowRef this,
        mri_tOrientation*   oOrientation ) {

  MWin_tErr        eResult      = MWin_tErr_NoErr;
  DspA_tErr        eDispResult  = DspA_tErr_NoErr;
  mri_tOrientation orientation  = mri_tOrientation_None;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* got the orientation from the last clicked display. */
  eDispResult = 
    DspA_GetOrientation( this->mapDisplays[this->mnLastClickedArea],
       &orientation );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
  /* return the orientation */
  *oOrientation = orientation;

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_GetOrientation: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_GetSelectedHeadPt ( tkmMeditWindowRef   this,
           HPtL_tHeadPointRef* opHeadPoint ) {

  MWin_tErr          eResult      = MWin_tErr_NoErr;
  DspA_tErr          eDispResult  = DspA_tErr_NoErr;
  HPtL_tHeadPointRef pHeadPoint   = NULL;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* got the head point from the last clicked display. */
  eDispResult = 
    DspA_GetSelectedHeadPt( this->mapDisplays[this->mnLastClickedArea],
          &pHeadPoint );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
  /* return the head point */
  *opHeadPoint = pHeadPoint;

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_GetSelectedHeadPt: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_GetClosestInterpSurfVoxel ( tkmMeditWindowRef this,
					   tkm_tSurfaceType  iType,
					   Surf_tVertexSet   iSet,
					   xVoxelRef         iAnaIdx,
					   xVoxelRef         oOrigAnaIdx,
					   xVoxelRef         oInterpAnaIdx,
					   char*             osDescription ) {

  MWin_tErr          eResult      = MWin_tErr_NoErr;
  DspA_tErr          eDispResult  = DspA_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* Pass the function to the last clicked display. */
  eDispResult = 
    DspA_GetClosestInterpSurfVoxel( this->mapDisplays[this->mnLastClickedArea],
				    iType, iSet, iAnaIdx, 
				    oOrigAnaIdx, oInterpAnaIdx,
				    osDescription);
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_GetClosestInterpSurfVoxel: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_AdjustSurfaceAnaIdx ( tkmMeditWindowRef   this,
             xVoxelRef           iAnaIdx ) {

  MWin_tErr          eResult      = MWin_tErr_NoErr;
  DspA_tErr          eDispResult  = DspA_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* tell the display area to adjust it */
  eDispResult = 
    DspA_AdjustSurfaceAnaIdx( this->mapDisplays[this->mnLastClickedArea],
            iAnaIdx );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_AdjustSurfaceAnaIdx: %s\n",
     eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_UnadjustSurfaceAnaIdx ( tkmMeditWindowRef   this,
               xVoxelRef           iAnaIdx ) {

  MWin_tErr          eResult      = MWin_tErr_NoErr;
  DspA_tErr          eDispResult  = DspA_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* tell the display area to adjust it */
  eDispResult = 
    DspA_UnadjustSurfaceAnaIdx( this->mapDisplays[this->mnLastClickedArea],
        iAnaIdx );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_UnadjustSurfaceAnaIdx: %s\n",
     eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_CursorChanged  ( tkmMeditWindowRef this,
        tkmDisplayAreaRef ipDisplay,
        xVoxelRef         ipCursor ) {

  MWin_tErr        eResult      = MWin_tErr_NoErr;
  DspA_tErr        eDispResult  = DspA_tErr_NoErr;
  int              nDisplay     = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* if we're not linking cursors, return. */
  if( !this->mbLinkedCursor )
    goto cleanup;

  /* if we are in mult mode.. */
  if( MWin_tLinkPolicy_MultipleOrientations == this->mLinkPolicy ) {
    
    /* for every display... */
    for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      /* if this is not the display that called us... */
      if( ipDisplay != (this->mapDisplays[nDisplay]) ) {
	
	/* set the cursor. */
	eDispResult = DspA_SetCursor ( this->mapDisplays[nDisplay], ipCursor );
	if ( DspA_tErr_NoErr != eDispResult ) {
	  eResult = MWin_tErr_ErrorAccessingDisplay;
	  goto error;
	}
  
      }
    }
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_CursorChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_ZoomLevelChanged  ( tkmMeditWindowRef this,
           tkmDisplayAreaRef ipDisplay,
           int               inZoomLevel ) {

  MWin_tErr        eResult      = MWin_tErr_NoErr;
  DspA_tErr        eDispResult  = DspA_tErr_NoErr;
  int              nDisplay     = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* if we're not linking cursors, return. */
  if( !this->mbLinkedCursor )
    goto cleanup;

  /* if we're not in single or mult orientation mode.. */
  if( MWin_tLinkPolicy_None == this->mLinkPolicy ||
      MWin_tLinkPolicy_MultipleOrientations == this->mLinkPolicy ) {
    
    /* for every display... */
    for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      /* if this is not the display that called us... */
      if( ipDisplay != (this->mapDisplays[nDisplay]) ) {
  
  /* set the zoom level */
  eDispResult = DspA_SetZoomLevel ( this->mapDisplays[nDisplay],
            inZoomLevel );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
      }
    }
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_ZoomLevelChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_DisplayFlagChanged  ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             DspA_tDisplayFlag iWhichFlag,
             tBoolean          ibNewValue ) {

  MWin_tErr        eResult      = MWin_tErr_NoErr;
  DspA_tErr        eDispResult  = DspA_tErr_NoErr;
  int              nDisplay     = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* if we're not linking cursors, return. */
  if( !this->mbLinkedCursor )
    goto cleanup;

  /* for every display... */
  for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {

    /* if this is not the display that called us... */
    if( ipDisplay != (this->mapDisplays[nDisplay]) ) {

      /* set the flag. */
      eDispResult = DspA_SetDisplayFlag ( this->mapDisplays[nDisplay], 
            iWhichFlag, ibNewValue );
      if ( DspA_tErr_NoErr != eDispResult ) {
  eResult = MWin_tErr_ErrorAccessingDisplay;
  goto error;
      }
    }
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_DisplayFlagChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_OrientationChanged  ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             mri_tOrientation  iOrientation ) {

  MWin_tErr  eResult      = MWin_tErr_NoErr;
  DspA_tErr  eDispResult  = DspA_tErr_NoErr;
  int        nDisplay     = 0;
  int        nSlice       = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* if we're not linking cursors, return. */
  if( !this->mbLinkedCursor )
    goto cleanup;

  /* if we are in mosaic mode.. */
  if( MWin_tLinkPolicy_Mosaic == this->mLinkPolicy ) {
    
    /* for every display... */
    for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      /* if this is not the display that called us... */
      if( ipDisplay != (this->mapDisplays[nDisplay]) ) {
  
  /* get current slice */
  DspA_GetSlice( this->mapDisplays[nDisplay], &nSlice );

  /* set the orientation. */
  eDispResult = DspA_SetOrientation ( this->mapDisplays[nDisplay], 
              iOrientation );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
      }
  
  /* set slice */
  DspA_SetSlice( this->mapDisplays[nDisplay], nSlice );
  
      }
    }
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_OrientationChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_SliceChanged  ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             int               inDelta ) {

  MWin_tErr  eResult      = MWin_tErr_NoErr;
  DspA_tErr  eDispResult  = DspA_tErr_NoErr;
  int        nDisplay     = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* if we're not linking cursors, return. */
  if( !this->mbLinkedCursor )
    goto cleanup;

  /* if we are in mosaic mode.. */
  if( MWin_tLinkPolicy_Mosaic == this->mLinkPolicy ) {
    
    /* for every display... */
    for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      /* if this is not the display that called us... */
      if( ipDisplay != (this->mapDisplays[nDisplay]) ) {
  
  /* change the slice */
  eDispResult = DspA_ChangeSliceBy ( this->mapDisplays[nDisplay], 
             inDelta );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
      }
    }
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_SliceChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;
}

void MWin_EventCallback ( void* ipWindow,
        xGWin_tEventRef ipEvent ) {

  /* typecast it and send it along */
  MWin_HandleEvent( (tkmMeditWindowRef)ipWindow, ipEvent );
}

void MWin_HandleEvent ( tkmMeditWindowRef   this, 
      xGWin_tEventRef     ipEvent ) {

  MWin_tErr eResult            = MWin_tErr_NoErr;
  DspA_tErr eDispResult        = DspA_tErr_NoErr;
  tBoolean  bWasHit            = FALSE;
  int       nDispIndex         = 0;
  int       nFlippedY          = 0;
  xPoint2n  location           = {0,0};
  int       nWidth             = 0;
  int       nHeight            = 0;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  switch ( ipEvent->mType ) {

  case xGWin_tEventType_KeyDown:
  case xGWin_tEventType_MouseDown:
  case xGWin_tEventType_MouseUp:
  case xGWin_tEventType_MouseMoved:

    /*
    if ( xGWin_tEventType_Idle != ipEvent->mType ) {
      DebugPrint( ("MWin_HandleEvent: Got\n" ) );
      xGWin_DebugPrintEvent( ipEvent );
    } 
    */

    /* hack! if this is a key-down and it's the home key, move the 
       tool window to right under us. */
    if( xGWin_tEventType_KeyDown == ipEvent->mType 
  && xGWin_tKey_Home == ipEvent->mKey ) {
      MWin_PlaceToolWindow_( this );
      goto cleanup;
    }

    /* look for tab. if so, change focus. */
    if( xGWin_tEventType_KeyDown == ipEvent->mType
  && xGWin_tKey_Tab == ipEvent->mKey ) {
      if( ipEvent->mbShiftKey ) {
  MWin_ChangeFocusedDisplayAreaBy_( this, -1 );
      } else {
  MWin_ChangeFocusedDisplayAreaBy_( this, 1 );
      }
      goto cleanup;
    }

    /* flip the y. i hate this. */
    nFlippedY = this->mnHeight - ipEvent->mWhere.mnY;
    
    /* these events have valid mouse points. only send the events to areas */
    /* contain the point. */
    /* go thru the display areas... */
    for ( nDispIndex = 0; 
    nDispIndex < MWin_knMaxNumAreas; 
    nDispIndex++ ) {
      
      /* assume not hit */
      bWasHit = FALSE;

      /* if this was a mouse event... */
      if( xGWin_tEventType_MouseMoved == ipEvent->mType
    || xGWin_tEventType_MouseDown == ipEvent->mType
    || xGWin_tEventType_MouseUp == ipEvent->mType ) {

  /* get the size and location */
  eDispResult = DspA_GetPosition ( this->mapDisplays[nDispIndex],
           &location, &nWidth, &nHeight );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
      
  /* find if it was hit. note we're using the flipped y here. */
  if ( ( ipEvent->mWhere.mnX    >= location.mnX 
         && ipEvent->mWhere.mnX <= location.mnX + nWidth ) &&
       ( nFlippedY              >= location.mnY
         && nFlippedY           <= location.mnY + nHeight ) ) {

    /* pass event to this pane. */
    bWasHit = TRUE;
  }

  /* else if it was a key event.. */
      } else if( xGWin_tEventType_KeyDown == ipEvent->mType ) {

  /* was hit if this is the focused pane */
  if( nDispIndex == this->mnLastClickedArea ) {
    bWasHit = TRUE;
  }
      }

      /* if this was hit... */
      if( bWasHit ) {
  
  /* if this was not the last clicked area and this was not
      a moved event... */
  if( nDispIndex != this->mnLastClickedArea &&
      ipEvent->mType != xGWin_tEventType_MouseMoved ) {

    /* save new focused display. HACK: we need to do this first, because
       calling DspA_Focus makes the display panel send all its view
       state information to tcl (DspA_SendViewStateToTcl_) which
       calls stuff like UpdateOrientation. tcl updates its orientation
       variable which triggers some of its update methods
       (UpdateOrientationWrapper), which calls SetOrientation, which
       sets the orientation of the pane[this->mnLastClickedArea]. arg */
    this->mnLastClickedArea = nDispIndex;
    
    /* focus on this display. */
    eDispResult = DspA_Focus( this->mapDisplays[nDispIndex] );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }

  }

  /* pass the event along. if it's a mouse moved event, only pass it
   if this is the focused pane. */
  if( (ipEvent->mType != xGWin_tEventType_MouseMoved) ||
      (nDispIndex == this->mnLastClickedArea &&
       ipEvent->mType == xGWin_tEventType_MouseMoved) ) {
    
    eDispResult = DspA_HandleEvent ( this->mapDisplays[nDispIndex],
             ipEvent );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }
      }
    }
    
    break;

  case xGWin_tEventType_Resize:

    /* set our width and height */
    this->mnWidth  = ipEvent->mWhere.mnX;
    this->mnHeight = ipEvent->mWhere.mnY;
    
    /* resize the opengl port */
    glutReshapeWindow ( this->mnWidth, this->mnHeight );
    
    /* recalc the display sizes. */
    MWin_PositionDisplays_( this );

    /* move the tool window to our position */
    MWin_PlaceToolWindow_( this );

    /* redraw */
    MWin_Redraw( this );
    break;

  case xGWin_tEventType_Draw:

    /* all display areas get and draw events. */
    for ( nDispIndex = 0; 
    nDispIndex < MWin_knMaxNumAreas; 
    nDispIndex++ ) {

      eDispResult = DspA_HandleEvent ( this->mapDisplays[nDispIndex],
               ipEvent );
      if ( DspA_tErr_NoErr != eDispResult ) {
  eResult = MWin_tErr_ErrorAccessingDisplay;
  goto error;
      }
    }

    /* call our draw handler. */
    MWin_HandleDraw_( this );

    break;

  case xGWin_tEventType_Idle:
    
    /* if idle, call tkmedit idle event handler. */
    if ( xGWin_tEventType_Idle == ipEvent->mType ) {
      tkm_HandleIdle();
    }
    
    break;

  default:

    eResult = MWin_tErr_InvalidEvent;
    goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_Event: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  /* we're responsbile for deleting the event */
  xGWin_DeleteEvent( &ipEvent );

  //  return eResult;
}

MWin_tErr MWin_Redraw ( tkmMeditWindowRef this ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;

  /* post a redisplay. */
  glutPostRedisplay();

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_Redraw: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_RedrawAll ( tkmMeditWindowRef this ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  int       nDispIndex  = 0;

  /* go thru each display area and set its rebuild slice flag. */
  for ( nDispIndex = 0; 
	nDispIndex < MWin_knMaxNumAreas; 
	nDispIndex++ ) {
    
    this->mapDisplays[nDispIndex]->mbSliceChanged = TRUE;
  }

  /* now draw ourselves. */
  MWin_Redraw( this );

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_RedrawAll: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_ForceRedraw ( tkmMeditWindowRef this ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDisplay    = DspA_tErr_NoErr;
  int       nDispIndex  = 0;

  /* go thru each display area and set its rebuild slice flag
     and force it to draw. */
  for ( nDispIndex = 0; 
  nDispIndex < MWin_knMaxNumAreas; 
  nDispIndex++ ) {
    
    this->mapDisplays[nDispIndex]->mbSliceChanged = TRUE;
    eDisplay = DspA_HandleDraw_( this->mapDisplays[nDispIndex] );

  }

  /* now force a redraw */
  MWin_HandleDraw_( this );

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_ForceRedraw: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_HandleDraw_ ( tkmMeditWindowRef this ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;

  /* swap the buffers */
  glutSwapBuffers();

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_HandleDraw_: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_PlaceToolWindow_ ( tkmMeditWindowRef this ) {

  MWin_tErr eResult            = MWin_tErr_NoErr;
  char      sTclArguments[256] = "";

  /* move tool window to just under our position */
  sprintf( sTclArguments, "+%d+%d",
     glutGet( GLUT_WINDOW_X ),
     glutGet( GLUT_WINDOW_Y ) + this->mnHeight + 
     MWin_knSpaceBetweenWindowAndPanel );
  tkm_SendTclCommand( tkm_tTclCommand_MoveToolWindow, sTclArguments );

  /* raise it */
  tkm_SendTclCommand( tkm_tTclCommand_RaiseToolWindow, "" );
    
  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_PlaceToolWindow_: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;

}

MWin_tErr MWin_ChangeFocusedDisplayAreaBy_ ( tkmMeditWindowRef this,
               int               inDelta ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDispIndex  = 0;
  int       nNumAreas   = 0;

  /* get the current display index */
  nDispIndex = this->mnLastClickedArea;
  
  /* how many areas do we have? */
  nNumAreas = this->mnCols * this->mnRows;

  /* inc. if more than max, set to 0. */
  nDispIndex += inDelta;
  while( nDispIndex >= nNumAreas )
    nDispIndex -= nNumAreas;
  while( nDispIndex < 0 ) 
    nDispIndex += MWin_knMaxNumAreas;

  /* if we changed... */
  if( this->mnLastClickedArea != nDispIndex ) {
    
    /* focus on this display. */
    this->mnLastClickedArea = nDispIndex;
    eDispResult = DspA_Focus( this->mapDisplays[this->mnLastClickedArea] );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
    
    /* redraw */
    MWin_Redraw( this );
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_ChangeFocusedDisplayAreaBy_: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;

}



MWin_tErr MWin_RegisterTclCommands ( tkmMeditWindowRef this,
             Tcl_Interp*       ipInterp ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* register our tcl commands */
  Tcl_CreateCommand ( ipInterp, "SetLinkedCursorFlag",
		      MWin_TclSetLinkedCursorFlag,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetCursor",
		      MWin_TclSetCursor,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetSlice",
		      MWin_TclSetSlice,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetOrientation",
		      MWin_TclSetOrientation,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetZoomLevel",
		      MWin_TclSetZoomLevel,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetZoomCenter",
		      MWin_TclSetZoomCenter,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetDisplayConfig",
		      MWin_TclSetDisplayConfig,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetDisplayFlag",
		      MWin_TclSetDisplayFlag,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetTool",
		      MWin_TclSetTool,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetBrushTarget",
		      MWin_TclSetBrushTarget,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetBrushShape",
		      MWin_TclSetBrushShape,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetBrushInfo",
		      MWin_TclSetBrushInfo,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetBrushInfoToDefaults",
		      MWin_TclSetBrushInfoToDefaults,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetCursorColor",
		      MWin_TclSetCursorColor,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetCursorShape",
		      MWin_TclSetCursorShape,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetSurfaceLineWidth",
		      MWin_TclSetSurfaceLineWidth,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetSurfaceLineColor",
		      MWin_TclSetSurfaceLineColor,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetFloodSelectParams",
		      MWin_TclSetFloodSelectParams,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetSegBrushInfo",
		      MWin_TclSetSegBrushInfo,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SelectCurrentSegLabel",
		      MWin_TclSelectCurrentSegLabel,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "GraphCurrentSegLabelAvg",
		      MWin_TclGraphCurrentSegLabelAvg,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetSurfaceDistanceAtCursor",
		      MWin_TclSetSurfaceDistanceAtCursor,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SmartCutAtCursor",
		      MWin_TclSmartCutAtCursor,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "RedrawAll",
		      MWin_TclRedrawAll,
		      (ClientData) this, (Tcl_CmdDeleteProc*) NULL );

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_RegisterTclCommands: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_AcceptTclCommands ( tkmMeditWindowRef this,
           tBoolean          ibAccept ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  /* set accept status. */
  this->mbAcceptingTclCommands = ibAccept;

  /* if going to true... */
  if( this->mbAcceptingTclCommands ) {

    /* verify the last clicked display area index. */
    eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
    if ( MWin_tErr_NoErr != eResult )
      goto error;
    
    /* focus on this pane. this will tell it to send its current display
       state to the tk window. */
    eDispResult = DspA_Focus( this->mapDisplays[this->mnLastClickedArea] );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in MWin_GetOrientation: %s\n",
      eResult, MWin_GetErrorString(eResult) ) );
  }
  
 cleanup:
  
  return eResult;


}

MWin_tErr MWin_Verify ( tkmMeditWindowRef this ) {

  MWin_tErr eResult = MWin_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this ) {
    eResult = MWin_tErr_InvalidPtr;
    goto cleanup;
  }
  
  /* check signature */
  if ( MWin_kSignature != this->mSignature ) {
    eResult = MWin_tErr_InvalidSignature;
    goto cleanup;
  }

 cleanup:

  return eResult;
}


MWin_tErr MWin_VerifyDisplayIndex ( tkmMeditWindowRef this,
            int               inDisplayIndex ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;

  /* check the number bounds */
  if ( ( inDisplayIndex < 0 
   || inDisplayIndex >= MWin_knMaxNumAreas )
       && MWin_kAllDisplayAreas != inDisplayIndex ) {
    eResult = MWin_tErr_InvalidDisplayIndex;
    goto cleanup;
  }

  /* check the actual display area */
  if ( MWin_kAllDisplayAreas != inDisplayIndex ) {
    eDispResult = DspA_Verify ( this->mapDisplays[inDisplayIndex] );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_InvalidDisplayArea;
      goto cleanup;
    }
  }

 cleanup:

  return eResult;
}



char* MWin_GetErrorString ( MWin_tErr ieCode ) {

  MWin_tErr eCode = ieCode;

  if ( ieCode < 0
       || ieCode >= MWin_knNumErrorCodes ) {
    eCode = MWin_tErr_InvalidErrorCode;
  }

  return MWin_ksaErrorStrings [eCode];
}

int MWin_TclSetLinkedCursorFlag ( ClientData  iClientData, 
          Tcl_Interp* ipInterp,
          int         argc,
          char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  tBoolean          bFlag        = FALSE;
  char              sError[256]  = "";       
 xVoxelRef          pCursor      = NULL;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
   
  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }
  
  /* parse the args and set the cursor. */
  bFlag = (tBoolean) atoi( argv[1] );

  /* set our flag. */
  eResult = MWin_SetLinkedCursorFlag( this, bFlag );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
    
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    
    sprintf ( sError, "Error %d in MWin_TclSetLinkedCursorFlag: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );
    
    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }
  
  eTclResult = TCL_ERROR;
  
 cleanup:

  /* delete the voxel */
  xVoxl_Delete( &pCursor );

  return eTclResult;
}

int MWin_TclSetCursor ( ClientData  iClientData, 
      Tcl_Interp* ipInterp,
      int         argc,
      char*       argv[] ) {
  
  tkmMeditWindowRef this               = NULL;
  int               eTclResult         = TCL_OK;
  MWin_tErr         eResult            = MWin_tErr_NoErr;
  char              sError[256]        = "";       
  xVoxelRef         pCursor            = NULL;
  mri_tCoordSpace   coordSpace         = mri_tCoordSpace_None;

  /* new the voxel */
  xVoxl_New( &pCursor );

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
   
  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }
  
  /* parse the args and set the cursor. */
  coordSpace = (mri_tCoordSpace) atoi( argv[1] );
  xVoxl_SetFloat( pCursor, atof( argv[2] ), atof( argv[3] ), atof( argv[4] ) );

  /* convert and set the cursor of the last clicked display. */
  eResult = MWin_ConvertAndSetCursor( this, this->mnLastClickedArea, 
              coordSpace, pCursor );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* zoom around new cursor. */
  eResult = MWin_SetZoomCenterToCursor ( this, -1 );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    
    sprintf ( sError, "Error %d in MWin_TclSetCursor: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );
    
    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }
  
  eTclResult = TCL_ERROR;
  
 cleanup:

  /* delete the voxel */
  xVoxl_Delete( &pCursor );

  return eTclResult;
}

int MWin_TclSetSlice ( ClientData  ipClientData, 
           Tcl_Interp* ipInterp,
           int         argc,
           char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  char              sError[256]  = "";       
  int               nSlice       = 0;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a slice. */
  nSlice = atoi( argv[1] );

  /* set the slice of the last clicked display. */
  eResult = MWin_SetSlice ( this, this->mnLastClickedArea, nSlice );
  if ( MWin_tErr_NoErr != eResult ) {
    goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetSlice: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}
int MWin_TclSetOrientation ( ClientData  ipClientData, 
           Tcl_Interp* ipInterp,
           int         argc,
           char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  char              sError[256]  = "";       
  mri_tOrientation  orientation  = mri_tOrientation_None;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get an orientation. */
  orientation = (mri_tOrientation) atoi( argv[1] );

  /* set the orientation of the last clicked display. */
  eResult = MWin_SetOrientation ( this, this->mnLastClickedArea, orientation );
  if ( MWin_tErr_NoErr != eResult ) {
    goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetOrientation: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}
int MWin_TclSetZoomLevel ( ClientData  ipClientData, 
         Tcl_Interp* ipInterp,
         int         argc,
         char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  char              sError[256]  = "";       
  int               nZoomLevel   = 0;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a zoom level. */
  nZoomLevel = atoi( argv[1] );

  /* set the zoom level of the last clicked display. */
  eResult = MWin_SetZoomLevel( this, this->mnLastClickedArea, nZoomLevel );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetZoomLevel: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetZoomCenter ( ClientData  iClientData, 
          Tcl_Interp* ipInterp,
          int         argc,
          char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  char              sError[256]  = "";       
 xVoxelRef          pCenter      = NULL;

  /* new the voxel */
  xVoxl_New( &pCenter );

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
   
  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }
  
  /* parse the args and set the Center. */
  xVoxl_Set( pCenter, atoi( argv[1] ), atoi( argv[2] ), atoi( argv[3] ) );

  /* set the Center of the last clicked display. */
  eResult = MWin_SetZoomCenter( this, this->mnLastClickedArea, pCenter );
  if ( MWin_tErr_NoErr != eResult )
    goto error;
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    
    sprintf ( sError, "Error %d in MWin_TclSetZoomCenter: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );
    
    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }
  
  eTclResult = TCL_ERROR;
  
 cleanup:

  /* delete the voxel */
  xVoxl_Delete( &pCenter );

  return eTclResult;
}


int MWin_TclSetDisplayConfig ( ClientData  ipClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] ) {

  tkmMeditWindowRef  this          = NULL;
  int                eTclResult    = TCL_OK;
  MWin_tErr          eResult       = MWin_tErr_NoErr;
  char               sError[256]   = "";       
  int                nCols         = 0;
  int                nRows         = 0;
  MWin_tLinkPolicy   policy        = MWin_tLinkPolicy_None;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get an configuration. */
  nCols = atoi( argv[1] );
  nRows = atoi( argv[2] );
  policy = (MWin_tLinkPolicy) atoi( argv[3] );

  /* set our configuration */
  eResult = MWin_SetDisplayConfiguration( this, nCols, nRows, policy );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetDisplayConfig: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetDisplayFlag ( ClientData  ipClientData, 
           Tcl_Interp* ipInterp,
           int         argc,
           char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  char              sError[256]  = "";       
  DspA_tDisplayFlag flag         = DspA_tDisplayFlag_None;
  tBoolean          bValue       = FALSE;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the number of arguments. */
  if ( argc < 3 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a flag and a value. */
  flag   = (DspA_tDisplayFlag) atoi( argv[1] );
  bValue = (tBoolean) atoi( argv[2] );

  /* set the flag of the last clicked display. */
  eResult = MWin_SetDisplayFlag( this, this->mnLastClickedArea, flag, bValue );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetDisplayFlag: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetTool ( ClientData  ipClientData, 
          Tcl_Interp* ipInterp,
          int         argc,
          char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  DspA_tTool        tool         = DspA_tTool_None;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a tool. */
  tool = (DspA_tTool) atoi( argv[1] );

  /* set the tool of the last clicked display. */
  eDispResult = 
    DspA_SetTool ( this->mapDisplays[this->mnLastClickedArea],
       tool );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetTool: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetBrushTarget ( ClientData  ipClientData, 
			     Tcl_Interp* ipInterp,
			     int         argc,
			     char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  DspA_tBrushTarget target       = DspA_tBrushTarget_None;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a target. */
  target  = (DspA_tBrushTarget) atoi( argv[1] );

  /* set the brush of the last clicked display. */
  eDispResult = 
    DspA_SetBrushTarget ( this->mapDisplays[this->mnLastClickedArea], target );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetBrushTarget: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetBrushShape ( ClientData  ipClientData, 
          Tcl_Interp* ipInterp,
          int         argc,
          char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  int               nRadius      = 0;
  DspA_tBrushShape  shape        = DspA_tBrushShape_None;
  tBoolean          b3D          = FALSE;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a radius, shape, and 3d flag. */
  nRadius = (int) atoi( argv[1] );
  shape   = (DspA_tBrushShape) atoi( argv[2] );
  b3D     = (tBoolean) atoi( argv[3] );

  /* set the brush of the last clicked display. */
  eDispResult = 
    DspA_SetBrushShape ( this->mapDisplays[this->mnLastClickedArea],
       nRadius, shape, b3D );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetBrushShape: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetBrushInfo ( ClientData  ipClientData, 
         Tcl_Interp* ipInterp,
         int         argc,
         char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  int               nBrush       = 0;
  int               nLow         = 0;
  int               nHigh        = 0;
  int               nNewValue    = 0;
  DspA_tBrushInfo   brushInfo;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 5 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a low, high, and new value */
  nBrush    = (int) atoi( argv[1] );
  nLow      = (int) atoi( argv[2] );
  nHigh     = (int) atoi( argv[3] );
  nNewValue = (int) atoi( argv[4] );

  /* make a struct */
  brushInfo.mnLow      = nLow;
  brushInfo.mnHigh     = nHigh;
  brushInfo.mnNewValue = nNewValue;

  /* set the brush of the last clicked display. */
  eDispResult = 
    DspA_SetBrushInfo ( this->mapDisplays[this->mnLastClickedArea],
      nBrush, &brushInfo );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetBrushInfo: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetBrushInfoToDefaults ( ClientData  ipClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  int               nBrush       = 0;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a brush */
  nBrush    = (int) atoi( argv[1] );

  /* call on the last clicked display. */
  eDispResult = 
    DspA_SetBrushInfoToDefault ( this->mapDisplays[this->mnLastClickedArea],
         nBrush  );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetBrushInfoToDefaults: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetCursorColor   ( ClientData  ipClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  xColor3f          cursorColor;           

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a color */
  cursorColor.mfRed = atof(argv[1]);
  cursorColor.mfGreen = atof(argv[2]);
  cursorColor.mfBlue = atof(argv[3]);

  /* call on the last clicked display. */
  eDispResult = 
    DspA_SetCursorColor ( this->mapDisplays[this->mnLastClickedArea],
        &cursorColor  );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetCursorColor: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetCursorShape ( ClientData  ipClientData, 
           Tcl_Interp* ipInterp,
           int         argc,
           char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  DspA_tMarker      shape        = DspA_tMarker_None;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a shape */
  shape = (DspA_tMarker) atoi( argv[1] );

  /* call on the last clicked display. */
  eDispResult = 
    DspA_SetCursorShape ( this->mapDisplays[this->mnLastClickedArea], shape );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetCursorShape: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetSurfaceLineWidth ( ClientData  ipClientData, 
				  Tcl_Interp* ipInterp,
				  int         argc,
				  char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  Surf_tVertexSet   set          = Surf_tVertexSet_None;
  tkm_tSurfaceType  surface      = tkm_tSurfaceType_Main;
  int               nWidth       = 0;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args */
  surface = (tkm_tSurfaceType) atoi( argv[1] );
  set     = (Surf_tVertexSet) atoi( argv[2] );
  nWidth  = atoi( argv[3] );

  /* call the on the last clicked display */
  eDispResult = 
    DspA_SetSurfaceLineWidth ( this->mapDisplays[this->mnLastClickedArea], 
			       surface, set, nWidth );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetSurfaceLineWidth: %s\n",
        eResult, MWin_GetErrorString(eResult) );
    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
} 

int MWin_TclSetSurfaceLineColor ( ClientData  ipClientData, 
          Tcl_Interp* ipInterp,
          int         argc,
          char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  Surf_tVertexSet   set          = Surf_tVertexSet_None;
  tkm_tSurfaceType  surface      = tkm_tSurfaceType_Main;
  xColor3f          color;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 6 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args */
  surface = (tkm_tSurfaceType) atoi( argv[1] );
  set     = (Surf_tVertexSet) atoi( argv[2] );
  color.mfRed   = atof( argv[3] );
  color.mfGreen = atof( argv[4] );
  color.mfBlue  = atof( argv[5] );

  /* call the on the last clicked display */
   eDispResult = 
     DspA_SetSurfaceLineColor ( this->mapDisplays[this->mnLastClickedArea], 
				surface, set, &color );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetSurfaceLineColor: %s\n",
        eResult, MWin_GetErrorString(eResult) );
    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
} 


int MWin_TclSetFloodSelectParams ( ClientData  ipClientData, 
				   Tcl_Interp* ipInterp,
				   int         argc,
				   char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  DspA_tFloodSelectSettings settings;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 5 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args */
  settings.mb3D        = (int) atoi( argv[1] );
  settings.mSrc        = (tkm_tVolumeTarget) atoi( argv[2] );
  settings.mnFuzzy     = (int) atoi( argv[3] );
  settings.mnDistance  = (int) atoi( argv[4] );

  /* call on the last clicked display. */
  eDispResult = 
    DspA_SetFloodSelectParams ( this->mapDisplays[this->mnLastClickedArea], 
				&settings );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetFloodSelectParams: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetSegBrushInfo ( ClientData  ipClientData, 
			       Tcl_Interp* ipInterp,
			       int         argc,
			       char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  DspA_tSegBrushSettings settings;

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 6 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a color and 3d */
  settings.mnPaintValue = (int) atoi( argv[1] );
  settings.mb3D         = (int) atoi( argv[2] );
  settings.mSrc         = (tkm_tVolumeTarget) atoi( argv[3] );
  settings.mnFuzzy      = (int) atoi( argv[4] );
  settings.mnDistance   = (int) atoi( argv[5] );

  /* call on the last clicked display. */
  eDispResult = 
    DspA_SetSegBrushInfo ( this->mapDisplays[this->mnLastClickedArea], 
			    &settings );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetSegBrushInfo: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSelectCurrentSegLabel ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* pass on to the last clicked display. */
  eDispResult = DspA_SelectCurrentSegLabel
    ( this->mapDisplays[this->mnLastClickedArea] );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSelectCurrentSegLabel: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclGraphCurrentSegLabelAvg ( ClientData  iClientData, 
         Tcl_Interp* ipInterp,
         int         argc,
         char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* pass on to the last clicked display. */
  eDispResult = DspA_GraphCurrentSegLabelAvg
    ( this->mapDisplays[this->mnLastClickedArea] );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclGraphCurrentSegLabelAvg: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}


int MWin_TclSetSurfaceDistanceAtCursor ( ClientData  iClientData, 
           Tcl_Interp* ipInterp,
           int         argc,
           char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* pass on to the last clicked display. */
  eDispResult = DspA_SetSurfaceDistanceAtCursor
    ( this->mapDisplays[this->mnLastClickedArea] );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetSurfaceDistanceAtCursor: %s\n",
	      eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSmartCutAtCursor ( ClientData  iClientData, 
			       Tcl_Interp* ipInterp,
			       int         argc,
			       char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) iClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* pass on to the last clicked display. */
  eDispResult = 
    DspA_SmartCutAtCursor ( this->mapDisplays[this->mnLastClickedArea] );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSmartCutAtCursor: %s\n",
	      eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}


int MWin_TclRedrawAll ( ClientData  ipClientData, 
      Tcl_Interp* ipInterp,
      int         argc,
      char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  char              sError[256]  = "";       

  /* grab us from the client data ptr */
  this = (tkmMeditWindowRef) ipClientData;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* if not accepting commands yet, return. */
  if( !this->mbAcceptingTclCommands )
    goto cleanup;

  /* redraw all. */
  eResult = MWin_RedrawAll( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclRedrawAll: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, MWin_GetErrorString(eResult), TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}
