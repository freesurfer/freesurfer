#include "tkmMeditWindow.h"
#include "tkmDisplayArea.h"

char MWin_ksaErrorStrings [MWin_knNumErrorCodes][256] = {

  "No error",
  "Allocation failed.",
  "Allocation of display area failed.",
  "Error accessing GLUT window.",
  "Invalid ptr to object (was probably NULL).",
  "Invalid signature found while verifying object.",
  "Invalid display index.",
  "Invalid display area.",
  "Invalid display configuration",
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

  /* set the size */
  this->mnWidth  = inWidth;
  this->mnHeight = inHeight;

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
  }

  /* default last clicked area. */
  this->mnLastClickedArea = 0;

  /* set the default configuration and position the displays */
  eResult =
    MWin_SetDisplayConfiguration ( this, MWin_tDisplayConfiguration_1x1 );
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
    DebugPrint "Error %d in MWin_New: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_Delete: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetWindowTitle ( tkmMeditWindowRef this,
        char*             isTitle ) {

  MWin_tErr  eResult  = MWin_tErr_NoErr;
  xGWin_tErr eGWin    = xGWin_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* set the gl window title */
  eGWin = xGWin_SetWindowTitle( this->mpGLutWindow, isTitle );
  if( eGWin != xGWin_tErr_NoErr ) {
    eResult = MWin_tErr_ErrorAccessingWindow;
    goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetWindowTitle: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

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
    DebugPrint "Error %d in MWin_SetWindowTitle: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetDisplayConfiguration ( tkmMeditWindowRef          this,
           MWin_tDisplayConfiguration iConfig ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;

  /* verify us. */
  eResult = MWin_Verify ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* save the configuration */
  this->mConfiguration = iConfig;

  /* position the subpanes */
  eResult = MWin_PositionDisplays_ ( this );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetDisplayConfiguration: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_PositionDisplays_ ( tkmMeditWindowRef this ) {

  MWin_tErr eResult     = MWin_tErr_NoErr;
  DspA_tErr eDispResult = DspA_tErr_NoErr;
  int       nDisplay    = 0;
  xPoint2n  location    = {0, 0};
  tkm_tOrientation orientation = tkm_tOrientation_None;

  switch ( this->mConfiguration ) {

  case MWin_tDisplayConfiguration_1x1:

    /* put the first one in the upper left corner */
    location.mnX = 0;
    location.mnY = 0;
    eDispResult = DspA_SetPosition ( this->mapDisplays[0],
             location, this->mnWidth, this->mnHeight );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }

    /* set the other ones to 0 width and height */
    for ( nDisplay = 1; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      eDispResult = DspA_SetPosition ( this->mapDisplays[nDisplay],
               location, 0, 0 );
      if ( DspA_tErr_NoErr != eDispResult ) {
  eResult = MWin_tErr_ErrorAccessingDisplay;
  goto error;
      }

    }

    break;

  case MWin_tDisplayConfiguration_2x2:

    /* put them in the four corners, with width and height of half of the
       window width and height */
    for ( nDisplay = 0; nDisplay < MWin_knMaxNumAreas; nDisplay++ ) {
      
      switch ( nDisplay ) {
      case 0:
  location.mnX = 0;
  location.mnY = 0;
  orientation = tkm_tOrientation_Coronal;
  break;
      case 1:
  location.mnX = this->mnWidth/2;
  location.mnY = 0;
  orientation = tkm_tOrientation_Sagittal;
  break;
      case 2:
  location.mnX = 0;
  location.mnY = this->mnHeight/2;
  orientation = tkm_tOrientation_Horizontal;
  break;
      case 3:
  location.mnX = this->mnWidth/2;
  location.mnY = this->mnHeight/2;
  orientation = tkm_tOrientation_Coronal;
  break;
      }
      
      eDispResult = 
  DspA_SetPosition ( this->mapDisplays[nDisplay],
         location, this->mnWidth / 2, this->mnHeight / 2);
      if ( DspA_tErr_NoErr != eDispResult ) {
  eResult = MWin_tErr_ErrorAccessingDisplay;
  goto error;
      }

      /* also set orientations */
      eDispResult = 
  DspA_SetOrientation ( this->mapDisplays[nDisplay], orientation );
      if ( DspA_tErr_NoErr != eDispResult ) {
  eResult = MWin_tErr_ErrorAccessingDisplay;
  goto error;
      }
      
    }
    break;

  default:
    eResult = MWin_tErr_InvalidDisplayConfiguration;
    goto error;
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_PositionDisplays_: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetVolume ( tkmMeditWindowRef this,
         int               inDispIndex,
         tVolumeRef        ipVolume,
         int               inSize ) {

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

    eDispResult = DspA_SetVolume ( this->mapDisplays[nDispIndex],
           ipVolume, inSize );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetAuxVolume ( tkmMeditWindowRef this,
            int               inDispIndex,
            tVolumeRef        ipVolume,
            int               inSize ) {

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
              ipVolume, inSize );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetAuxVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetParcellationVolume ( tkmMeditWindowRef this,
               int               inDispIndex,
               tVolumeRef        ipVolume,
               int               inSize ) {

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

    eDispResult = DspA_SetParcellationVolume ( this->mapDisplays[nDispIndex],
                 ipVolume, inSize );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetParcellationVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}
MWin_tErr MWin_SetSurface ( tkmMeditWindowRef this, 
          int               inDispIndex,
          MRI_SURFACE*      ipSurface ) {

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
            ipSurface );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetSurface: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_SetOverlayVolume: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetControlPointsSpace ( tkmMeditWindowRef this,
               int               inDispIndex,
               VoxelSpaceRef     ipVoxels ) {

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
    DebugPrint "Error %d in MWin_SetControlPointsSpace: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetControlPointsSelectionList ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 VoxelListRef      ipVoxels ) {

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

  /* set the selected control points */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = 
      DspA_SetControlPointsSelectionList ( this->mapDisplays[nDispIndex],
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
    DebugPrint "Error %d in MWin_SetControlPointsSelectionList: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


MWin_tErr MWin_SetSelectionSpace ( tkmMeditWindowRef this, 
           int               inDispIndex,
           VoxelSpaceRef     ipVoxels ) {
  
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
    DebugPrint "Error %d in MWin_SetSelectionSpace: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_SetLinkedCursorFlag: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_ToggleLinkedCursorFlag: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetCursor ( tkmMeditWindowRef this, 
         int               inDispIndex,
         VoxelRef          ipCursor ) {

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
    DebugPrint "Error %d in MWin_SetCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetOrientation ( tkmMeditWindowRef this, 
        int               inDispIndex,
        tkm_tOrientation  iOrientation ) {

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

  /* set the orientation */
  for ( nDispIndex = nDispIndexMin; 
  nDispIndex < nDispIndexMax; 
  nDispIndex++ ) {

    eDispResult = DspA_SetOrientation ( this->mapDisplays[nDispIndex],
          iOrientation );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_SetOrientation: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_SetZoomCenter ( tkmMeditWindowRef this, 
             int               inDispIndex,
             VoxelRef          ipCenter ) {

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
    DebugPrint "Error %d in MWin_SetZoomCenter: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_SetZoomCenterToCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_HiliteSurfaceVertex ( tkmMeditWindowRef this,
             int               inDispIndex,
             tkm_tSurfaceType  inSurface,
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
    DebugPrint "Error %d in MWin_HiliteSurfaceVertex: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_SetDisplayFlag: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_GetCursor ( tkmMeditWindowRef this,
         VoxelRef          opCursor ) {
  
  MWin_tErr eResult      = MWin_tErr_NoErr;
  DspA_tErr eDispResult  = DspA_tErr_NoErr;
  VoxelRef  pCursor      = NULL;

  /* new cursor. */
  Voxel_New( &pCursor );

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
  Voxel_Copy( opCursor, pCursor );

  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_GetCursor: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }
  
 cleanup:
  
  /* delete cursor */
  Voxel_Delete( &pCursor );
  
  return eResult;
}

MWin_tErr MWin_GetOrientation ( tkmMeditWindowRef this,
        tkm_tOrientation*   oOrientation ) {

  MWin_tErr        eResult      = MWin_tErr_NoErr;
  DspA_tErr        eDispResult  = DspA_tErr_NoErr;
  tkm_tOrientation orientation  = tkm_tOrientation_None;

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
    DebugPrint "Error %d in MWin_GetOrientation: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }
  
 cleanup:
  
  return eResult;
}

MWin_tErr MWin_CursorChanged  ( tkmMeditWindowRef this,
        tkmDisplayAreaRef ipDisplay,
        VoxelRef          ipCursor ) {

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

    /* if this is the display that called us... */
    if( ipDisplay != (this->mapDisplays[nDisplay]) ) {

      /* set the cursor. */
      eDispResult = DspA_SetCursor ( this->mapDisplays[nDisplay], ipCursor );
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
    DebugPrint "Error %d in MWin_CursorChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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

    /* if this is the display that called us... */
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
    DebugPrint "Error %d in MWin_DisplayFlagChanged: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
      DebugPrint "MWin_HandleEvent: Got\n" EndDebugPrint;
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

    /* flip the y. i hate this. */
    nFlippedY = this->mnHeight - ipEvent->mWhere.mnY;
    
    /* these events have valid mouse points. only send the events to areas */
    /* contain the point. */
    /* go thru the display areas... */
    for ( nDispIndex = 0; 
    nDispIndex < MWin_knMaxNumAreas; 
    nDispIndex++ ) {
      
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

  /* if this was a mouse event, see if this was not the last clicked
     area... */
  if( nDispIndex != this->mnLastClickedArea ) {

    /* focus on this display. */
    eDispResult = DspA_Focus( this->mapDisplays[nDispIndex] );
    if ( DspA_tErr_NoErr != eDispResult ) {
      eResult = MWin_tErr_ErrorAccessingDisplay;
      goto error;
    }
  }

  /* save new focused display. */
  this->mnLastClickedArea = nDispIndex;
  
  /* pass the event along */
  eDispResult = DspA_HandleEvent ( this->mapDisplays[nDispIndex],
           ipEvent );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
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
    DebugPrint "Error %d in MWin_Event: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_Redraw: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

MWin_tErr MWin_RedrawAll ( tkmMeditWindowRef this ) {

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

  /* now draw ourselves. */
  MWin_Redraw( this );

  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_RedrawAll: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_HandleDraw_: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    
  goto cleanup;

  goto error;
 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_PlaceToolWindow_: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
  Tcl_CreateCommand ( ipInterp, "SetBrush",
          MWin_TclSetBrush,
          (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "SetBrushThreshold",
          MWin_TclSetBrushThreshold,
          (ClientData) this, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( ipInterp, "RedrawAll",
          MWin_TclRedrawAll,
          (ClientData) this, (Tcl_CmdDeleteProc*) NULL );

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in MWin_RegisterTclCommands: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
    DebugPrint "Error %d in MWin_GetOrientation: %s\n",
      eResult, MWin_GetErrorString(eResult) EndDebugPrint;
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
  VoxelRef          pCursor      = NULL;

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
    
    sprintf ( sError, "Error %d in MWin_TclSetCursor: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;
    
    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }
  
  eTclResult = TCL_ERROR;
  
 cleanup:

  /* delete the voxel */
  Voxel_Delete( &pCursor );

  return eTclResult;
}

int MWin_TclSetCursor ( ClientData  iClientData, 
      Tcl_Interp* ipInterp,
      int         argc,
      char*       argv[] ) {
  
  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  VoxelRef          pCursor      = NULL;

  /* new the voxel */
  Voxel_New( &pCursor );

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
  if ( argc < 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }
  
  /* parse the args and set the cursor. */
  Voxel_Set( pCursor, atoi( argv[1] ), atoi( argv[2] ), atoi( argv[3] ) );

  /* set the cursor of the last clicked display. */
  eDispResult = DspA_SetCursor ( this->mapDisplays[this->mnLastClickedArea],
         pCursor );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    
    sprintf ( sError, "Error %d in MWin_TclSetCursor: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;
    
    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }
  
  eTclResult = TCL_ERROR;
  
 cleanup:

  /* delete the voxel */
  Voxel_Delete( &pCursor );

  return eTclResult;
}

int MWin_TclSetOrientation ( ClientData  ipClientData, 
           Tcl_Interp* ipInterp,
           int         argc,
           char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  tkm_tOrientation  orientation  = tkm_tOrientation_None;

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

  /* parse the args and get an orientation. */
  orientation = (tkm_tOrientation) atoi( argv[1] );

  /* set the zoom level of the last clicked display. */
  eDispResult = 
    DspA_SetOrientation ( this->mapDisplays[this->mnLastClickedArea],
        orientation );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetOrientation: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
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
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
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

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a zoom level. */
  nZoomLevel = atoi( argv[1] );

  /* set the zoom level of the last clicked display. */
  eDispResult = 
    DspA_SetZoomLevel ( this->mapDisplays[this->mnLastClickedArea],
      nZoomLevel );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetZoomLevel: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
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
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  VoxelRef          pCenter      = NULL;

  /* new the voxel */
  Voxel_New( &pCenter );

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
  if ( argc < 4 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }
  
  /* parse the args and set the Center. */
  Voxel_Set( pCenter, atoi( argv[1] ), atoi( argv[2] ), atoi( argv[3] ) );

  /* set the Center of the last clicked display. */
  eDispResult = DspA_SetZoomCenter( this->mapDisplays[this->mnLastClickedArea],
            pCenter );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  
  goto cleanup;
  
 error:
  
  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {
    
    sprintf ( sError, "Error %d in MWin_TclSetZoomCenter: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;
    
    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }
  
  eTclResult = TCL_ERROR;
  
 cleanup:

  /* delete the voxel */
  Voxel_Delete( &pCenter );

  return eTclResult;
}


int MWin_TclSetDisplayConfig ( ClientData  ipClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] ) {

  tkmMeditWindowRef          this          = NULL;
  int                        eTclResult    = TCL_OK;
  MWin_tErr                  eResult       = MWin_tErr_NoErr;
  char                       sError[256]   = "";       
  MWin_tDisplayConfiguration configuration = MWin_tDisplayConfiguration_None;

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

  /* parse the args and get an configuration. */
  configuration = (MWin_tDisplayConfiguration) atoi( argv[1] );

  /* set our configuration */
  eResult = MWin_SetDisplayConfiguration( this, configuration );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetDisplayConfig: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
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
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
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

  /* verify the last clicked display area index. */
  eResult = MWin_VerifyDisplayIndex ( this, this->mnLastClickedArea );
  if ( MWin_tErr_NoErr != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 3 ) {
    eResult = MWin_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse the args and get a flag and a value. */
  flag   = (DspA_tDisplayFlag) atoi( argv[1] );
  bValue = (tBoolean) atoi( argv[2] );

  /* set the flag of the last clicked display. */
  eDispResult = 
    DspA_SetDisplayFlag ( this->mapDisplays[this->mnLastClickedArea],
        flag, bValue );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetDisplayFlag: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
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

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetBrush ( ClientData  ipClientData, 
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
    DspA_SetBrush ( this->mapDisplays[this->mnLastClickedArea],
        nRadius, shape, b3D );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetBrush: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}

int MWin_TclSetBrushThreshold ( ClientData  ipClientData, 
        Tcl_Interp* ipInterp,
        int         argc,
        char*       argv[] ) {

  tkmMeditWindowRef this         = NULL;
  int               eTclResult   = TCL_OK;
  MWin_tErr         eResult      = MWin_tErr_NoErr;
  DspA_tErr         eDispResult  = DspA_tErr_NoErr;
  char              sError[256]  = "";       
  int               nLow         = 0;
  int               nHigh        = 0;
  int               nNewValue    = 0;

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

  /* parse the args and get a low, high, and new value */
  nLow      = (int) atoi( argv[1] );
  nHigh     = (int) atoi( argv[2] );
  nNewValue = (int) atoi( argv[3] );

  /* set the brush of the last clicked display. */
  eDispResult = 
    DspA_SetBrushThreshold ( this->mapDisplays[this->mnLastClickedArea],
           nLow, nHigh, nNewValue );
  if ( DspA_tErr_NoErr != eDispResult ) {
    eResult = MWin_tErr_ErrorAccessingDisplay;
    goto error;
  }
  goto cleanup;

 error:

  /* print error message */
  if ( MWin_tErr_NoErr != eResult ) {

    sprintf ( sError, "Error %d in MWin_TclSetBrushThreshold: %s\n",
        eResult, MWin_GetErrorString(eResult) );

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
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

    DebugPrint sError EndDebugPrint;

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

 cleanup:

  return eTclResult;
}
