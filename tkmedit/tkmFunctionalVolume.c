/**
 * @file  tkmFunctionalVolume.c
 * @brief Manages functional volume overlay and timecourse plotting
 *
 * Provides an interface for functional volumes as an overlay or
 * timecourse plot. The same volume can be both, or you can have
 * different volumes for each purpose. Handles loading via the
 * mriFunctionaDataAccess code. Keeps track of the current time point
 * and condition, and sends the Tcl commands to interface with
 * tkm_functional.tcl to draw the graph.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/01 01:41:22 $
 *    $Revision: 1.66 $
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


#include "tkmFunctionalVolume.h"
#include <stdlib.h>
#include <math.h>
#include "mri_transform.h"
#include "xUtilities.h"
#include "mri.h"
#include "utils.h"
#include "error.h"

#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)

/* about the script */
#define ksEnvVariable_UseLocalDirectoryForScript "DONT_USE_LOCAL_TKMFUNCTIONAL_TCL"
#define ksFileName_InterfaceScript "tkm_functional.tcl"
#define ksDir_LibraryPath          getenv ( "FREESURFER_HOME" )
#define ksDir_LibrarySubPath       "/tktools/"

#define knLengthOfGraphDataItem              18 // for "100.1 1000.12345 "
#define knLengthOfGraphDataHeader            20 // for header + cond + {}
#define knMaxCommandLength                   50

/* error strings */
char *FunV_ksaErrorString [FunV_tErr_knNumErrorCodes] = {
  "No error.",
  "Internal allocation failed.",
  "Internal deletion failed.",
  "Error parsing the given path and stem. Make sure it begins with a backslash and ends with the stem, i.e. \"/path/to/data/stem\".",
  "Couldn't load volume.",
  "Error accessing internal volume structure.",
  "Error accessing internal list structure.",
  "Error converting the second to a time point.",
  "Couldn't allocate overlay cache.",
  "Couldn't open file.",
  "Overlay data not loaded.",
  "Time course data not loaded.",
  "Trying to init graph window when it is already inited.",
  "Trying to draw graph when window is not inited.",
  "Couldn't find tcl/tk script file for graph.",
  "Couldn't load tcl/tk script file for graph.",
  "Invlid pointer to function volume (was probably NULL).",
  "Invalid signature in function volume (memory probably trashed).",
  "Invalid parameter.",
  "Invalid time point.",
  "Invalid condition.",
  "Invalid anatomical  voxel, data doesn't exist in functional space.",
  "Invalid threshold, min must be less than mid.",
  "Invalid display flag.",
  "Wrong number of arguments.",
  "Error accessing anatomical volume.",
  "Couldn't load brain volume to use ask mask.",
  "Invalid error code."
};

char *FunV_ksaTclCommand [FunV_knNumTclCommands] = {
  "Overlay_DoConfigDlog",
  "Overlay_UpdateNumTimePoints",
  "Overlay_UpdateNumConditions",
  "Overlay_UpdateDisplayFlag",
  "Overlay_UpdateDataName",
  "Overlay_UpdateTimePoint",
  "Overlay_UpdateCondition",
  "Overlay_UpdateThreshold",
  "Overlay_UpdateRange",
  "Overlay_ShowOffsetOptions",
  "TimeCourse_DoConfigDlog",
  "TimeCourse_BeginDrawingGraph",
  "TimeCourse_EndDrawingGraph",
  "TimeCourse_DrawGraph",
  "TimeCourse_ClearData",
  "TimeCourse_UpdateNumConditions",
  "TimeCourse_UpdateNumTimePoints",
  "TimeCourse_UpdateTimePoint",
  "TimeCourse_UpdateNumPreStimPoints",
  "TimeCourse_UpdateTimeResolution",
  "TimeCourse_UpdateDisplayFlag",
  "TimeCourse_UpdateDataName",
  "TimeCourse_UpdateLocationName",
  "TimeCourse_UpdateGraphData",
  "TimeCourse_UpdateErrorData",
  "TimeCourse_ShowOffsetOptions"
};

FunV_tErr FunV_New ( tkmFunctionalVolumeRef* oppVolume,
                     void(*ipOverlayChangedFunction)(void),
                     void(*ipSendTkmeditCmdFunction)(tkm_tTclCommand,char*),
                     char*(*ipSendTclCommandFunction)(char*) ) {

  FunV_tErr              eResult = FunV_tErr_NoError;
  tkmFunctionalVolumeRef this    = NULL;
  xList_tErr             eList   = xList_tErr_NoErr;
  int                    nFlag   = 0;

  /* allocate us */
  this = (tkmFunctionalVolumeRef) malloc ( sizeof(tkmFunctionalVolume) );
  if ( NULL == this ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* set signature */
  this->mSignature = FunV_kSignature;

  /* set volumes to null */
  this->mpOverlayVolume           = NULL;
  this->mpTimeCourseVolume        = NULL;
  this->mpOverlayOffsetVolume     = NULL;
  this->mpTimeCourseOffsetVolume  = NULL;
  this->mOverlayCache            = NULL;

  /* allocate voxel list */
  this->mpSelectedVoxels = NULL;
  eList = xList_New( &(this->mpSelectedVoxels) );
  if ( xList_tErr_NoErr != eList ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* set default state values */
  this->mbUseOverlayCache      = FALSE;
  this->mnCachedTimePoint      = -1;
  this->mnCachedCondition      = -1;
  this->mnTimePoint            = 0;
  this->mnCondition            = 0;
  this->mThresholdMin          = 0;
  this->mThresholdMid          = 1;
  this->mThresholdSlope        = 1;
  this->mbGraphInited          = FALSE;
  this->mbRegistrationEnabled  = FALSE;

  /* set flags to false */
  for ( nFlag = 0; nFlag < FunV_knNumDisplayFlags; nFlag++ )
    this->mabDisplayFlags[nFlag] = FALSE;

  /* default flag values */
  this->mabDisplayFlags[FunV_tDisplayFlag_Ol_Opaque] = TRUE;

  /* set functions */
  this->mpOverlayChangedFunction    = ipOverlayChangedFunction;
  this->mpSendTkmeditTclCmdFunction = ipSendTkmeditCmdFunction;
  this->mpSendTclCommandFunction    = ipSendTclCommandFunction;

  /* send a message telling the interface to hide the functional value */
  eResult = FunV_SendTkmeditTclCommand_( this,
                                         tkm_tTclCommand_ShowFuncValue, 
                                         "0" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );

  /* disable func overlay display options */
  eResult = 
    FunV_SendTkmeditTclCommand_( this,
                                 tkm_tTclCommand_ShowFuncOverlayOptions, 
                                 "0" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );

  /* disable time course display options */
  eResult = 
    FunV_SendTkmeditTclCommand_( this,
                                 tkm_tTclCommand_ShowFuncTimeCourseOptions, 
                                 "0" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadTimeCourse", __LINE__, eResult );

  /* return the volume */
  *oppVolume = this;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_New: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

  /* if we allocated us, delete us. */
  if ( NULL != this->mpSelectedVoxels )
    xList_Delete( &(this->mpSelectedVoxels) );
  if ( NULL != this )
    free( this );

cleanup:

  return eResult;
}

FunV_tErr FunV_Delete ( tkmFunctionalVolumeRef* ioppVolume ) {

  FunV_tErr              eResult = FunV_tErr_NoError;
  tkmFunctionalVolumeRef this    = NULL;
  xList_tErr             eList   = xList_tErr_NoErr;
  FunD_tErr              eVolume = FunD_tErr_NoError;

  /* get us */
  this = *ioppVolume;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* delete the volumes that exist. */
  if ( NULL != this->mpOverlayVolume ) {
    eVolume = FunD_Delete( &(this->mpOverlayVolume) );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_DeletionFailed;
      goto error;
    }
  }
  if ( NULL != this->mpTimeCourseVolume ) {
    eVolume = FunD_Delete( &(this->mpTimeCourseVolume) );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_DeletionFailed;
      goto error;
    }
  }
  if ( NULL != this->mpOverlayOffsetVolume ) {
    eVolume = FunD_Delete( &(this->mpOverlayOffsetVolume) );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_DeletionFailed;
      goto error;
    }
  }
  if ( NULL != this->mpTimeCourseOffsetVolume ) {
    eVolume = FunD_Delete( &(this->mpTimeCourseOffsetVolume) );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_DeletionFailed;
      goto error;
    }
  }

  /* clear the list */
  FunV_BeginSelectionRange( this );
  FunV_EndSelectionRange( this );

  /* delete the list */
  eList = xList_Delete( &(this->mpSelectedVoxels) );
  if ( xList_tErr_NoErr != eList ) {
    eResult = FunV_tErr_DeletionFailed;
    goto error;
  }

  /* trash signature */
  this->mSignature = 0x1;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_Delete: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_LoadOverlay ( tkmFunctionalVolumeRef this,
                             char*                  isFileName,
                             char*                  isOffsetFileName,
                             FunV_tRegistrationType iRegistrationType,
                             char*                  isRegistration,
                             mriVolumeRef           iAnatomicalVolume ) {

  FunV_tErr              eResult           = FunV_tErr_NoError;
  int                    nNumConditions    = 0;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* attempt to load the volume */
  eResult = FunV_LoadFunctionalVolume_( this, &(this->mpOverlayVolume),
                                        isFileName,
                                        NULL,
                                        iRegistrationType, isRegistration,
                                        iAnatomicalVolume, TRUE );
  if ( FunV_tErr_NoError != eResult ) {
    goto error;
  }

  /* if we have more than one condition, set condition to 1. */
  FunD_GetNumConditions( this->mpOverlayVolume, &nNumConditions );
  if ( nNumConditions > 1 ) {
    this->mnCondition = 1;
  }

  /* if we have an offset path... */
  if ( NULL != isOffsetFileName ) {

    /* attempt to load offset volume, if there is one */
    eResult = FunV_LoadFunctionalVolume_( this,
                                          &(this->mpOverlayOffsetVolume),
                                          isOffsetFileName,
                                          NULL,
                                          iRegistrationType, isRegistration,
                                          iAnatomicalVolume, FALSE );
    if ( FunV_tErr_NoError != eResult ) {
      /* no offset, that's fine. */
      this->mpOverlayOffsetVolume = NULL;
    }
  }

  /* initialize that overlay cache */
  if ( this->mbUseOverlayCache ) {
    eResult = FunV_InitOverlayCache_( this );
    if ( FunV_tErr_NoError != eResult ) {
      FunV_Signal( "FunV_LoadFunctionalVolume_: loading cache",
                   __LINE__, eResult );
      eResult = FunV_tErr_NoError;
    }
  }

  /* send a message telling the interface to show the functional value */
  eResult = FunV_SendTkmeditTclCommand_( this,
                                         tkm_tTclCommand_ShowFuncValue, 
                                         "1" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );

  /* show func overlay display options */
  eResult = 
    FunV_SendTkmeditTclCommand_( this,
                                 tkm_tTclCommand_ShowFuncOverlayOptions, 
                                 "1" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );

  /* if we have offset data, display offset options as well. */
  if ( NULL != this->mpOverlayOffsetVolume ) {
    eResult = 
      FunV_SendTclCommand_( this,
                            FunV_tTclCommand_Ol_ShowOffsetOptions, 
                            "1" );
    if ( FunV_tErr_NoError != eResult )
      FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );
  }

  /* send all our view state to tcl */
  eResult = FunV_SendViewStateToTcl( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* notify interested parties that the overlay is changed. */
  FunV_OverlayChanged_( this );

  goto cleanup;

error:

  /* send a message telling the interface to hide the functional value */
  eResult = FunV_SendTkmeditTclCommand_( this,
                                         tkm_tTclCommand_ShowFuncValue, "0" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );

  /* hide func overlay display options */
  eResult = 
    FunV_SendTkmeditTclCommand_( this,
                                 tkm_tTclCommand_ShowFuncOverlayOptions, "0" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );

  /* hide offset options */
  if ( NULL != this->mpOverlayOffsetVolume ) {
    eResult = 
      FunV_SendTclCommand_( this,
                            FunV_tTclCommand_Ol_ShowOffsetOptions, "0" );
    if ( FunV_tErr_NoError != eResult )
      FunV_Signal ( "FunV_LoadOverlay", __LINE__, eResult );
  }

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_LoadOverlay: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_LoadTimeCourse ( tkmFunctionalVolumeRef this,
                                char*                  isFileName,
                                char*                  isOffsetFileName,
                                FunV_tRegistrationType iRegistrationType,
                                char*                  isRegistration,
                                mriVolumeRef           iAnatomicalVolume ) {

  FunV_tErr              eResult           = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* attempt to load the volume */
  eResult = FunV_LoadFunctionalVolume_( this,
                                        &(this->mpTimeCourseVolume),
                                        isFileName,
                                        NULL,
                                        iRegistrationType, isRegistration,
                                        iAnatomicalVolume, TRUE );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if we have an offset path... */
  if ( NULL != isOffsetFileName ) {

    /* attempt to load offset volume, if there is one */
    eResult = FunV_LoadFunctionalVolume_( this,
                                          &(this->mpTimeCourseOffsetVolume),
                                          isOffsetFileName,
                                          NULL,
                                          iRegistrationType, isRegistration,
                                          iAnatomicalVolume, FALSE );
    if ( FunV_tErr_NoError != eResult ) {
      /* no offset, that's fine. */
      this->mpTimeCourseOffsetVolume = NULL;
    }
  }

  /* show the graph window */
  eResult = FunV_SetDisplayFlag( this, FunV_tDisplayFlag_TC_GraphWindowOpen,
                                 TRUE );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadTimeCourse", __LINE__, eResult );

  /* show time course display options */
  eResult = 
    FunV_SendTkmeditTclCommand_( this,
                                 tkm_tTclCommand_ShowFuncTimeCourseOptions, 
                                 "1" );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal ( "FunV_LoadTimeCourse", __LINE__, eResult );

  /* if we have offset data, display offset options as well. */
  if ( NULL != this->mpTimeCourseOffsetVolume ) {
    eResult = 
      FunV_SendTclCommand_( this,
                            FunV_tTclCommand_TC_ShowOffsetOptions, "1" );
    if ( FunV_tErr_NoError != eResult )
      FunV_Signal ( "FunV_LoadTimeCourse", __LINE__, eResult );
  }

  /* send all our view state to tcl */
  eResult = FunV_SendViewStateToTcl( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_LoadTimeCourse: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_LoadFunctionalVolume_ ( tkmFunctionalVolumeRef this,
                                       mriFunctionalDataRef*  ioppVolume,
                                       char*                  isFileName,
                                       char*                  isHeaderStem,
                                       FunV_tRegistrationType iRegistrationType,
                                       char*                  isRegistration,
                                       mriVolumeRef          iAnatomicalVolume,
                                       tBoolean              ibReportErrors ) {

  FunV_tErr            eResult            = FunV_tErr_NoError;
  FunD_tRegistrationType regType           = FunD_tRegistration_None;
  FunD_tErr            eVolume            = FunD_tErr_NoError;
  mriFunctionalDataRef pVolume            = NULL;
  char                 sError[tkm_knErrStringLen] = "";

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* Convert the registration type. */
  switch ( iRegistrationType ) {
  case FunV_tRegistration_File:
    regType = FunD_tRegistration_File;
    break;
  case FunV_tRegistration_Find:
    regType = FunD_tRegistration_Find;
    break;
  case FunV_tRegistration_Identity:
    regType = FunD_tRegistration_Identity;
    break;
  default:
    eResult = FunV_tErr_InvalidParameter;
    goto error;
  }

  /* delete volume if already allocated */
  if ( NULL != *ioppVolume ) {
    eVolume = FunD_Delete( ioppVolume );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_DeletionFailed;
      goto error;
    }
  }

  /* load the volume */
  eVolume = FunD_New( &pVolume,
                      isFileName,
                      regType,
                      isRegistration,
                      -1,  /* Don't try to be a scalar volume */
                      iAnatomicalVolume );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorLoadingVolume;
    goto error;
  }

  /* return the volume */
  *ioppVolume = pVolume;

  goto cleanup;

error:

  /* put up an error dlog if the volume didn't load */
  if ( FunD_tErr_NoError != eVolume
       && ibReportErrors ) {
    xUtil_snprintf( sError, sizeof(sError),
                    "Loading functional overlay %s.", isFileName );
    tkm_DisplayError( sError, FunD_GetErrorString( eVolume ),
                      "Tkmedit couldn't read the functional volume you "
                      "specified. This could be because the volume "
                      "wasn't found or was unreadable, or because a "
                      "valid header type couldn't be find, or a "
                      "registration file couldn't be found or opened." );
  }

  /* print error message */
  if ( FunV_tErr_NoError != eResult
       && ibReportErrors ) {
    DebugPrint( ("Error %d in FunV_LoadFunctionalVolume_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  /* if no we don't want to print errors, don't return them */
  if ( !ibReportErrors )
    eResult = FunV_tErr_NoError;

  return eResult;
}

FunV_tErr FunV_SetConversionMethod ( tkmFunctionalVolumeRef this,
                                     FunD_tConversionMethod iMethod ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  FunD_tErr eVolume = FunD_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* save this as the new default method */
  this->mDefaultConvMethod = iMethod;

  /* set the method in our volumes */
  if ( NULL != this->mpOverlayVolume ) {
    eVolume = FunD_SetConversionMethod( this->mpOverlayVolume, iMethod );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_ErrorAccessingInternalVolume;
      goto error;
    }
  }
  if ( NULL != this->mpOverlayOffsetVolume ) {
    eVolume = FunD_SetConversionMethod( this->mpOverlayOffsetVolume, iMethod );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_ErrorAccessingInternalVolume;
      goto error;
    }
  }
  if ( NULL != this->mpTimeCourseVolume ) {
    eVolume = FunD_SetConversionMethod( this->mpTimeCourseVolume, iMethod );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_ErrorAccessingInternalVolume;
      goto error;
    }
  }
  if ( NULL != this->mpTimeCourseOffsetVolume ) {
    eVolume = FunD_SetConversionMethod( this->mpTimeCourseOffsetVolume,
                                        iMethod );
    if ( FunD_tErr_NoError != eVolume ) {
      eResult = FunV_tErr_ErrorAccessingInternalVolume;
      goto error;
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetConversionMethod: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_InitOverlayCache_ ( tkmFunctionalVolumeRef this ) {

  FunV_tErr             eResult    = FunV_tErr_NoError;
  int                   nCacheSize = 0;
  xVoxelRef              pVoxel     = NULL;
  xVoxelRef              pMin     = NULL;
  xVoxelRef              pMax     = NULL;
  int                   nX         = 0;
  int                   nY         = 0;
  int                   nZ         = 0;
  FunV_tFunctionalValue funcValue  = 0;
  int                   nIndex     = 0;
  xDbg_tDebuggingState  debugState;

  xVoxl_New( &pVoxel );
  xVoxl_New( &pMin );
  xVoxl_New( &pMax );

  /* if we're off, return */
  if ( FALSE == this->mbUseOverlayCache ) {
    goto cleanup;
  }

  /* if we already have this time point and condition cached, exit */
  if ( this->mnCachedCondition == this->mnCondition
       && this->mnCachedTimePoint == this->mnTimePoint ) {
    goto cleanup;
  }

  /* if the cache volume already exists, kill it. */
  if ( NULL != this->mOverlayCache ) {
    free( this->mOverlayCache );
    this->mOverlayCache = NULL;
  }

  /* get our func bounds in antomical space */
  FunD_GetBoundsInClientSpace( this->mpOverlayVolume, pMin, pMax );

  /* calc dimensions */
  this->manCacheDimensions[0] = xVoxl_GetX(pMax) - xVoxl_GetX(pMin);
  this->manCacheDimensions[1] = xVoxl_GetY(pMax) - xVoxl_GetY(pMin);
  this->manCacheDimensions[2] = xVoxl_GetZ(pMax) - xVoxl_GetZ(pMin);

  /* save offsets */
  this->manCacheOffsets[0] = xVoxl_GetX(pMin);
  this->manCacheOffsets[1] = xVoxl_GetY(pMin);
  this->manCacheOffsets[2] = xVoxl_GetZ(pMin);

  /* allocate the cache. */
  nCacheSize = 
    this->manCacheDimensions[0] *
    this->manCacheDimensions[1] * 
    this->manCacheDimensions[2] * sizeof(float);
  this->mOverlayCache = (float*) malloc( nCacheSize );
  if ( NULL == this->mOverlayCache ) {
    eResult = FunV_tErr_ErrorAllocatingOverlayCache;
    goto error;
  }
  bzero( (this->mOverlayCache), nCacheSize );


  /* disable the cache so we actually calcuate the values instead of getting
     them out of our freshly bzeroed cache. */
  FunV_UseOverlayCache( this, FALSE );

  GetDebuggingState( &debugState );
  DisableDebuggingOutput;

  /* get our initial index */
  nIndex = 0;

  /* for every voxel from one corner to the other... */
  for ( nZ = 0; nZ < this->manCacheDimensions[2]; nZ++ ) {

    OutputPrint "\rBuilding cache... %.2f%%",
    ( (float)nZ / (float)this->manCacheDimensions[2] * 100.0 ) EndOutputPrint;

    for ( nY = 0; nY < this->manCacheDimensions[1]; nY++ ) {
      for ( nX = 0; nX < this->manCacheDimensions[0]; nX++ ) {

        /* set our voxel with the offsets */
        xVoxl_Set( pVoxel, nX + this->manCacheOffsets[0],
                   nY + this->manCacheOffsets[1],
                   nZ + this->manCacheOffsets[2] );

        /* get the functional value at this voxel */
        FunV_GetValueAtMRIIdx( this, pVoxel, TRUE, &funcValue );

        /* set the cache value */
        this->mOverlayCache[nIndex] = funcValue;

        nIndex ++;
      }
    }
  }

  /* mark what time poitn and condition we have cached */
  this->mnCachedTimePoint = this->mnTimePoint;
  this->mnCachedCondition = this->mnCondition;

  /* reenable cache. note that this will call the init function again,
     which is bad design - i admit - but since we have just saved the currently
     cached tp and cond, it will exit right away. */
  FunV_UseOverlayCache( this, TRUE );

  SetDebuggingState( debugState );

  OutputPrint " done!\n" EndOutputPrint;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_InitOverlayCache_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &pVoxel );
  xVoxl_Delete( &pMin );
  xVoxl_Delete( &pMax );

  return eResult;
}

FunV_tErr FunV_SetOverlayCacheValue_ ( tkmFunctionalVolumeRef this,
                                       xVoxelRef               ipVoxel,
                                       FunV_tFunctionalValue  iValue ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  int       nIndex  = 0;

  nIndex =
    ((xVoxl_GetZ(ipVoxel) - this->manCacheOffsets[2]) *
     (this->manCacheDimensions[1] * this->manCacheDimensions[0])) +
    ((xVoxl_GetY(ipVoxel) - this->manCacheOffsets[1]) *
     this->manCacheDimensions[0]) +
    (xVoxl_GetX(ipVoxel) - this->manCacheOffsets[0]) ;

  this->mOverlayCache[ nIndex ] = iValue;

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetOverlayCacheValue_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;

}

FunV_tErr FunV_GetOverlayCacheValue_ ( tkmFunctionalVolumeRef this,
                                       xVoxelRef               ipVoxel,
                                       FunV_tFunctionalValue* oValue ) {


  FunV_tErr eResult = FunV_tErr_NoError;
  int       nIndex  = 0;

  /* check bounds */
  if ( xVoxl_GetX(ipVoxel) < this->manCacheOffsets[0]
       || xVoxl_GetX(ipVoxel) > (this->manCacheOffsets[0] + 
                                 this->manCacheDimensions[0])
       || xVoxl_GetY(ipVoxel) < this->manCacheOffsets[1]
       || xVoxl_GetY(ipVoxel) > (this->manCacheOffsets[1] + 
                                 this->manCacheDimensions[1])
       || xVoxl_GetZ(ipVoxel) < this->manCacheOffsets[2]
       || xVoxl_GetZ(ipVoxel) > (this->manCacheOffsets[2] + 
                                 this->manCacheDimensions[2]) ) {
    *oValue = 0;
    eResult = FunV_tErr_InvalidMRIIdx;
    goto cleanup;
  }

  nIndex =
    ((xVoxl_GetZ(ipVoxel) - this->manCacheOffsets[2]) *
     (this->manCacheDimensions[1] * this->manCacheDimensions[0])) +
    ((xVoxl_GetY(ipVoxel) - this->manCacheOffsets[1]) *
     this->manCacheDimensions[0]) +
    (xVoxl_GetX(ipVoxel) - this->manCacheOffsets[0]) ;

  *oValue = this->mOverlayCache[ nIndex ];

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_GetOverlayCacheValue_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


FunV_tErr FunV_UseOverlayCache ( tkmFunctionalVolumeRef this,
                                 tBoolean               ibUseCache ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* set flag */
  this->mbUseOverlayCache = ibUseCache;

  /* if we're on, init the cache. */
  if ( this->mbUseOverlayCache ) {

    eResult = FunV_InitOverlayCache_( this );
    if ( FunV_tErr_NoError != eResult ) {

      /* turn cache off if it didn't work. */
      this->mbUseOverlayCache = FALSE;
      goto error;
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_UseOverlayCache: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;

}

FunV_tErr FunV_IsOverlayPresent    ( tkmFunctionalVolumeRef this,
                                     tBoolean*              obIsLoaded ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  tBoolean  bLoaded = FALSE;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* check if loaded */
  if ( NULL != this->mpOverlayVolume ) {
    bLoaded = TRUE;
  }

  /* return status */
  *obIsLoaded = bLoaded;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_IsOverlayPresent: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SmoothOverlayData ( tkmFunctionalVolumeRef this,
                                   float                  ifSigma ) {


  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* smooth if loaded */
  if ( NULL != this->mpOverlayVolume ) {
    FunD_Smooth( this->mpOverlayVolume, 0, 0, ifSigma );
  } else {
    eResult = FunV_tErr_OverlayNotLoaded;
    goto error;
  }

  /* reinit the cache and call for a redraw. */
  FunV_InitOverlayCache_( this );
  FunV_OverlayChanged_( this );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SmoothOverlayData: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;

}

FunV_tErr FunV_IsTimeCoursePresent ( tkmFunctionalVolumeRef this,
                                     tBoolean*              obIsLoaded ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  tBoolean  bLoaded = FALSE;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* check if loaded */
  if ( NULL != this->mpTimeCourseVolume ) {
    bLoaded = TRUE;
  }

  /* return status */
  *obIsLoaded = bLoaded;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_IsTimeCoursePresent: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


FunV_tErr FunV_SetTimeResolution ( tkmFunctionalVolumeRef this,
                                   float                  inTimeResolution ) {

  FunV_tErr         eResult            = FunV_tErr_NoError;
  FunD_tErr         eVolume            = FunD_tErr_NoError;
  float             nTimeRes           = 0;
  char              sTclArguments[tkm_knTclCmdLen] = "";
  char              sError[tkm_knErrStringLen] = "";

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* set value in volume */
  eVolume = FunD_SetTimeResolution( this->mpTimeCourseVolume,
                                    inTimeResolution );
  if ( FunD_tErr_NoError != eVolume
       && FunD_tErr_InvalidTimeResolution != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* get the new time res. */
  eVolume = FunD_GetTimeResolution( this->mpTimeCourseVolume,
                                    &nTimeRes );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* if they weren't equal, display an error message. */
  if ( nTimeRes != inTimeResolution ) {
    xUtil_snprintf( sError, sizeof(sError),
                    "Changing time potin to %d.", inTimeResolution );
    tkm_DisplayError( sError, "Invalid time point.",
                      "Tkmedit couldn't set the time point to the value "
                      "you specified. Please make sure it is in range, "
                      "from 0 to the number of time points minus one." );
  }

  /* send the new value to tcl */
  sprintf( sTclArguments, "%f", nTimeRes );
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateTimeResolution,
                        sTclArguments );

  /* udpate the graph */
  eResult = FunV_DrawGraph( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetTimeResolution: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SetNumPreStimPoints ( tkmFunctionalVolumeRef this,
                                     int inNumPreStimPoints ) {

  FunV_tErr         eResult            = FunV_tErr_NoError;
  FunD_tErr         eVolume            = FunD_tErr_NoError;
  int               nNumPreStimPoints  = 0;
  int               nNumPoints         = 0;
  char              sTclArguments[tkm_knTclCmdLen] = "";
  char              sError[tkm_knErrStringLen] = "";
  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* set value in volume */
  eVolume = FunD_SetNumPreStimTimePoints( this->mpTimeCourseVolume,
                                          inNumPreStimPoints );
  if ( FunD_tErr_NoError != eVolume
       && FunD_tErr_InvalidNumPreStimTimePoints != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* get the new num pts. */
  eVolume = FunD_GetNumPreStimTimePoints( this->mpTimeCourseVolume,
                                          &nNumPreStimPoints );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* if they weren't equal, display an error message. */
  if ( nNumPreStimPoints != inNumPreStimPoints ) {

    /* get number of time points */
    FunD_GetNumTimePoints( this->mpTimeCourseVolume, &nNumPoints );

    xUtil_snprintf( sError, sizeof(sError),
                    "Setting number of pre-stim time points to %d.",
                    inNumPreStimPoints );
    tkm_DisplayError
      ( sError, "Invalid number of points.",
        "Tkmedit couldn't set the number of pre-stim time points "
        "to the value you specified. Please make sure it is "
        "greater than 0 and less than the total number of time "
        "points." );
  }

  /* send the new value to tcl */
  sprintf( sTclArguments, "%d", nNumPreStimPoints );
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateNumPreStimPoints,
                        sTclArguments );

  /* udpate the graph */
  eResult = FunV_DrawGraph( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetNumPreStimPoints: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SetTimeSecond ( tkmFunctionalVolumeRef this,
                               int                    inSecond ) {

  FunV_tErr         eResult     = FunV_tErr_NoError;
  FunD_tErr  eVolume     = FunD_tErr_NoError;
  int               nTimePoint  = 0;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* don't do anything if we don't have overlay data. */
  if ( NULL == this->mpOverlayVolume )
    goto cleanup;

  /* convert to time point. */
  eVolume = FunD_ConvertSecondToTimePoint( this->mpOverlayVolume,
            inSecond, &nTimePoint );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorConvertingSecondToTimePoint;
    goto error;
  }

  /* set time point */
  eResult = FunV_SetTimePoint( this, nTimePoint );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetTimeSecond: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SetTimePoint ( tkmFunctionalVolumeRef this,
                              int                    inTimePoint ) {

  FunV_tErr         eResult            = FunV_tErr_NoError;
  FunD_tErr         eVolume            = FunD_tErr_NoError;
  float             fSecond            = 0;
  int               nNumTimePoints     = 0;
  char              sTclArguments[tkm_knTclCmdLen] = "";
  char              sError[tkm_knErrStringLen] = "";

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* don't do anything if we don't have overlay data. */
  if ( NULL == this->mpOverlayVolume )
    goto cleanup;

  if (inTimePoint >= this->mpOverlayVolume->mpData->nframes) inTimePoint = 0;
  if (inTimePoint < 0) inTimePoint = this->mpOverlayVolume->mpData->nframes-1;

  /* if it's valid, set it */
  eVolume = FunD_VerifyTimePoint( this->mpOverlayVolume, inTimePoint );
  if ( FunD_tErr_NoError == eVolume ) {

    /* set it */
    this->mnTimePoint = inTimePoint;

  } else {

    /* get number of time points. */
    FunD_GetNumTimePoints( this->mpOverlayVolume, &nNumTimePoints );

    /* display an error message */
    xUtil_snprintf( sError, sizeof(sError),
                    "Setting time point to %d.", inTimePoint );
    tkm_DisplayError( sError, "Invalid time point.",
                      "Tkmedit couldn't set the time point to the value you "
                      "specfied. Please make sure it is greater than 0 and "
                      "less than the total number of time points." );

    goto error;
  }


  /* Send time point update to the Overlay. */
  sprintf( sTclArguments, "%d", this->mnTimePoint );
  FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateTimePoint,
                        sTclArguments );

  /* If we have a time course, send time point update to the time
     course (we also convert the TP to a second and send that along
     too). */
  if( NULL != this->mpTimeCourseVolume ) {
    
    eVolume = FunD_ConvertTimePointToSecond( this->mpTimeCourseVolume,
					     this->mnTimePoint, &fSecond );
    sprintf( sTclArguments, "%d %f", this->mnTimePoint, fSecond );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateTimePoint,
			  sTclArguments );
  }

  /* reinit the cache */
  eResult = FunV_InitOverlayCache_( this );
  if ( FunV_tErr_NoError != eResult ) {

    /* turn cache off */
    FunV_UseOverlayCache( this, FALSE );

    /* clear error */
    eResult = FunV_tErr_NoError;
  }

  /* udpate the overlay */
  eResult = FunV_OverlayChanged_( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetTimePoint: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SetCondition ( tkmFunctionalVolumeRef this,
                              int                    inCondition ) {

  FunV_tErr         eResult            = FunV_tErr_NoError;
  FunD_tErr         eVolume            = FunD_tErr_NoError;
  int               nNumConditions     = 0;
  char              sTclArguments[tkm_knTclCmdLen] = "";
  char              sError[tkm_knErrStringLen] = "";

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* don't do anything if we don't have overlay data. */
  if ( NULL == this->mpOverlayVolume )
    goto cleanup;

  /* if valid.. */
  eVolume = FunD_VerifyCondition( this->mpOverlayVolume, inCondition );
  if ( FunD_tErr_NoError == eVolume ) {

    /* set it */
    this->mnCondition = inCondition;

  } else {

    /* get number of conditions */
    FunD_GetNumConditions( this->mpOverlayVolume, &nNumConditions );

    /* display an error message */
    xUtil_snprintf( sError, sizeof(sError),
                    "Setting condition to %d.", inCondition );
    tkm_DisplayError( sError, "Invalid condition.",
                      "Tkmedit couldn't set the time point to the value you "
                      "specfied. Please make sure it is greater than 0 and "
                      "less than the total number of conditions." );

    goto error;
  }

  /* send the value to tcl */
  sprintf( sTclArguments, "%d", this->mnCondition );
  FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateCondition,
                        sTclArguments );

  /* reinit the cache */
  eResult = FunV_InitOverlayCache_( this );
  if ( FunV_tErr_NoError != eResult ) {

    /* turn cache off */
    FunV_UseOverlayCache( this, FALSE );

    /* clear error */
    eResult = FunV_tErr_NoError;
  }

  /* reinit the cache */
  eResult = FunV_InitOverlayCache_( this );
  if ( FunV_tErr_NoError != eResult ) {

    /* turn cache off */
    FunV_UseOverlayCache( this, FALSE );

    /* clear error */
    eResult = FunV_tErr_NoError;
  }

  /* udpate the overlay */
  eResult = FunV_OverlayChanged_( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetCondition: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SetThreshold ( tkmFunctionalVolumeRef this,
                              FunV_tFunctionalValue  iMin,
                              FunV_tFunctionalValue  iMid,
                              FunV_tFunctionalValue  iSlope ) {

  FunV_tErr         eResult            = FunV_tErr_NoError;
  char              sTclArguments[256] = "";

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if they are the same, return */
  if ( this->mThresholdMin      == iMin
       && this->mThresholdMid   == iMid
       && this->mThresholdSlope == iSlope ) {
    goto cleanup;
  }

  /* if valid.. */
  if ( iMin <= iMid
       && iSlope > 0 ) {

    /* set everything */
    this->mThresholdMin   = iMin;
    this->mThresholdMid   = iMid;
    this->mThresholdSlope = iSlope;

  } else {

    /* display an error message */
    tkm_DisplayError( "Setting threshold values", "Invalid values.",
                      "Tkmedit couldn't set the threshold to the values you "
                      "specfied. Please make sure the minimum is less than "
                      "or equal to the midpoint, and the slope is greater "
                      "than 0." );
  }

  /* send the new value to tcl */
  sprintf( sTclArguments, "%f %f %f",
           (float)(this->mThresholdMin),
           (float)(this->mThresholdMid),
           (float)(this->mThresholdSlope) );
  FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateThreshold,
                        sTclArguments );

  /* udpate the overlay */
  eResult = FunV_OverlayChanged_( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetThreshold: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SetThresholdUsingFDR ( tkmFunctionalVolumeRef this,
                                      float                  iRate,
                                      tBoolean               ibMaskToBrain ) {

  FunV_tErr  eResult   = FunV_tErr_NoError;
  FunD_tErr  eFunD     = FunD_tErr_NoError;
  int        nSignFlag = 0;
  float      newMin    = 0;
  char       fnBrainVol[1024] = "";
  char       fnMaskStem[1024] = ""; 
  char       *pTmpStr = NULL;
  MRI*       pMaskVol  = NULL;

  DebugEnterFunction( ("FunV_SetThresholdUsingFDR( this=%p, iRate=%f, "
                       "ibMaskToBrain=%d)", this, iRate, ibMaskToBrain) );

  /* verify us */
  DebugNote( ("Verifying this") );
  eResult = FunV_Verify( this );
  DebugAssertThrow( (FunV_tErr_NoError == eResult) );

  /* verify the rate */
  DebugAssertThrowX( (iRate > 0 && iRate < 1),
                     eResult, tkm_tErr_InvalidParameter );

  /* Load up the brain volume and pass that in as the mask. We can
     pass this in 'client' space (our space) as mriFunctionalData will
     take care of converting it properly. */
  if ( ibMaskToBrain ) {
    DebugNote( ("Making filename for mask volume") );
    pTmpStr = getenv("TKM_FDR_MASK");
    if(pTmpStr == NULL) sprintf(fnMaskStem,"%s","brain.mgz");
    else                sprintf(fnMaskStem,"%s",pTmpStr);
    tkm_MakeFileName( fnMaskStem, tkm_tFileName_Volume,
                      fnBrainVol, sizeof(fnBrainVol) );
    printf("FDR Mask is %s\n",fnBrainVol);
    DebugNote( ("Reading brain volume from %s", fnBrainVol) );
    pMaskVol = MRIread( fnBrainVol );
    if(pMaskVol == NULL) printf("ERROR: could not load %s\n",fnBrainVol);
      
    DebugAssertThrowX( (NULL != pMaskVol),
                       eResult, FunV_tErr_CouldntLoadBrainMask );
  }

  /* Determine a sign. */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_TruncateNegative] ) {
    nSignFlag = 1;
  }
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_TruncatePositive] ) {
    nSignFlag = -1;
  }
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_ReversePhase] ) {
    nSignFlag = -nSignFlag;
  }

  /* Calc the FDR. */
  DebugNote( ("Calcing FDR") );
  eFunD = FunD_CalcFDRThreshold( this->mpOverlayVolume,
                                 this->mnCondition, this->mnTimePoint,
                                 nSignFlag, iRate, pMaskVol,
                                 &newMin );
  DebugAssertThrowX( (FunD_tErr_NoError == eFunD),
                     eResult, FunV_tErr_ErrorAccessingInternalVolume );

  /* Set threshold based on return value. */
  DebugNote( ("Setting threshold") );
  eResult = FunV_SetThreshold( this, newMin, newMin + 1.5, 0.66 );
  DebugAssertThrow( FunV_tErr_NoError == eResult );

  DebugCatch;
  DebugCatchError( eResult, FunV_tErr_NoError, FunV_GetErrorString );
  EndDebugCatch;

  if ( NULL != pMaskVol ) {
    DebugNote( ("Freeing mask vol") );
    MRIfree( &pMaskVol );
  }

  DebugExitFunction;

  return eResult;
}

FunV_tErr FunV_SetDisplayFlag ( tkmFunctionalVolumeRef this,
                                FunV_tDisplayFlag      iFlag,
                                tBoolean               iNewValue ) {

  FunV_tErr         eResult            = FunV_tErr_NoError;
  char              sTclArguments[256] = "";
  tBoolean          bNewValue          = FALSE;
  tBoolean          bUpdateOverlay     = FALSE;
  tBoolean          bUpdateGraph       = FALSE;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the flag */
  if ( iFlag <= FunV_tDisplayFlag_None
       || iFlag >= FunV_knNumDisplayFlags ) {
    eResult = FunV_tErr_InvalidDisplayFlag;
    goto error;
  }

  /* get the new value */
  bNewValue = iNewValue;

  /* act on the flag */
  switch ( iFlag ) {

  case FunV_tDisplayFlag_TC_GraphWindowOpen:

    /* if no time course data, set to false. */
    if ( NULL == this->mpTimeCourseVolume ) {
      bNewValue = FALSE;
    }

    /* if showing, need to update graph */
    if ( TRUE == bNewValue ) {
      bUpdateGraph = TRUE;
    }

    break;

  case FunV_tDisplayFlag_Ol_TruncateNegative:
  case FunV_tDisplayFlag_Ol_TruncatePositive:
  case FunV_tDisplayFlag_Ol_ReversePhase:
  case FunV_tDisplayFlag_Ol_IgnoreThreshold:
  case FunV_tDisplayFlag_Ol_Grayscale:
  case FunV_tDisplayFlag_Ol_Opaque:

    /* if no overlay data, set to false */
    if ( NULL == this->mpOverlayVolume ) {
      bNewValue = FALSE;

    } else {

      /* if values are different, need to update overlay */
      if ( this->mabDisplayFlags[iFlag] != bNewValue ) {
        bUpdateOverlay = TRUE;
      }
    }

    break;

  case FunV_tDisplayFlag_Ol_OffsetValues:

    /* if no offset data, set to false. */
    if ( NULL == this->mpOverlayOffsetVolume ) {
      bNewValue = FALSE;

    } else {

      /* if values are different, need to update overlay */
      if ( this->mabDisplayFlags[iFlag] != bNewValue ) {
        bUpdateOverlay = TRUE;
      }
    }

    break;

  case FunV_tDisplayFlag_TC_OffsetValues:
  case FunV_tDisplayFlag_TC_PreStimOffset:

    /* if no offset data, set to false. */
    if ( NULL == this->mpTimeCourseOffsetVolume ) {
      bNewValue = FALSE;

    } else {

      /* if values are different, need to update graph */
      if ( this->mabDisplayFlags[iFlag] != bNewValue ) {
        bUpdateGraph = TRUE;
      }
    }

    break;

  default:

    eResult = FunV_tErr_InvalidDisplayFlag;
    goto error;
  }

  /* set the value */
  this->mabDisplayFlags[iFlag] = bNewValue;

  /* send the tcl update */
  sprintf( sTclArguments, "%d %d", (int)iFlag, (int)bNewValue );

  /* is it an overlay flag or a time course flag? */
  if ( iFlag >= FunV_knFirstOverlayDisplayFlag
       && iFlag <= FunV_knLastOverlayDisplayFlag ) {
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateDisplayFlag,
                          sTclArguments );
  } else {
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateDisplayFlag,
                          sTclArguments );
  }
  /* udpate the overlay */
  if ( bUpdateOverlay ) {
    eResult = FunV_OverlayChanged_( this );
    if ( FunV_tErr_NoError != eResult )
      goto error;
  }

  /* udpate the graph */
  if ( bUpdateGraph ) {
    eResult = FunV_DrawGraph( this );
    if ( FunV_tErr_NoError != eResult )
      goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetDisplayFlag: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_EnableRegistration  ( tkmFunctionalVolumeRef this,
                                     tBoolean               iNewValue ) {

  FunV_tErr         eResult        = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* set the new value */
  this->mbRegistrationEnabled = iNewValue;

  /* enable or disable the menu items */
  eResult = 
    FunV_SendTkmeditTclCommand_
    ( this,
      tkm_tTclCommand_ShowOverlayRegistrationOptions,
      this->mbRegistrationEnabled ? "1" : "0" );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_EnableRegistration: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}


FunV_tErr FunV_ChangeTimePointBy ( tkmFunctionalVolumeRef this,
                                   int                    inDelta ) {

  FunV_tErr         eResult        = FunV_tErr_NoError;
  int               nTimePoint     = 0;
  int               nNumTimePoints = 0;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* don't do anything if we don't have overlay data. */
  if ( NULL == this->mpOverlayVolume )
    goto cleanup;

  /* get the time point */
  nTimePoint = this->mnTimePoint;

  /* get the number of time points */
  FunD_GetNumTimePoints( this->mpOverlayVolume, &nNumTimePoints );
  if ( nNumTimePoints == 0 )
    goto cleanup;

  /* add the delta */
  nTimePoint += inDelta;

  /* if out of range, bounce it back */
  while ( nTimePoint < 0 ) {
    nTimePoint += nNumTimePoints;
  }
  while ( nTimePoint >= nNumTimePoints ) {
    nTimePoint -= nNumTimePoints;
  }

  /* set time point */
  eResult = FunV_SetTimePoint( this, nTimePoint );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_ChangeTimePointBy: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_MRIIdxClicked ( tkmFunctionalVolumeRef this,
                               xVoxelRef              iMRIIdx ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  if ( NULL == iMRIIdx ) {
    eResult = FunV_tErr_InvalidParameter;
    goto error;
  }

  /* if we have no data, return */
  if ( NULL == this->mpOverlayVolume
       && NULL == this->mpTimeCourseVolume ) {
    goto cleanup;
  }

  /* set our selection range to just this voxel */
  eResult = FunV_BeginSelectionRange( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  eResult = FunV_AddMRIIdxToSelectionRange( this, iMRIIdx );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  eResult = FunV_EndSelectionRange( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_MRIIdxClicked: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_ApplyTransformToOverlay ( tkmFunctionalVolumeRef this,
    MATRIX*                iTransform ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  FunD_tErr eVolume = FunD_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  if ( NULL == iTransform ) {
    eResult = FunV_tErr_InvalidParameter;
    goto error;
  }

  /* if we have no data, return */
  if ( NULL == this->mpOverlayVolume ) {
    goto cleanup;
  }

  /* apply the transform */
  eVolume = FunD_ApplyTransformToRegistration( this->mpOverlayVolume,
            iTransform );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_ApplyTransformToOverlay: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_TranslateOverlayRegistration ( tkmFunctionalVolumeRef this,
    float                  ifDistance,
    tAxis                  iAxis ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  FunD_tErr eVolume = FunD_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if we have no data, return */
  if ( NULL == this->mpOverlayVolume ) {
    goto cleanup;
  }

  /* translate the overlay registration */
  eVolume = FunD_TranslateRegistration( this->mpOverlayVolume,
                                        ifDistance, iAxis );
  if ( FunV_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_TranslateOverlayRegistration: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_RotateOverlayRegistration ( tkmFunctionalVolumeRef this,
    float                  ifDegrees,
    tAxis                  iAxis,
    xVoxelRef         iCenterMRIIdx ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  FunD_tErr eVolume = FunD_tErr_NoError;
  xVoxel    centerFuncRAS;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if we have no data, return */
  if ( NULL == this->mpOverlayVolume ) {
    goto cleanup;
  }

  /* convert ana idx center to func idx */
  FunV_ConvertMRIIdxToFuncRAS( this, iCenterMRIIdx, &centerFuncRAS );

  /* rotate the overlay registration */
  eVolume = FunD_RotateRegistration( this->mpOverlayVolume,
                                     ifDegrees, iAxis, &centerFuncRAS );

  if ( FunV_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_RotateOverlayRegistration: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_ScaleOverlayRegistration ( tkmFunctionalVolumeRef this,
    float                  ifFactor,
    tAxis                  iAxis ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  FunD_tErr eVolume = FunD_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if we have no data, return */
  if ( NULL == this->mpOverlayVolume ) {
    goto cleanup;
  }

  /* scale the overlay registration */
  eVolume = FunD_ScaleRegistration( this->mpOverlayVolume,
                                    ifFactor, iAxis );
  if ( FunV_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_ScaleOverlayRegistration: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_GetValueAtMRIIdx ( tkmFunctionalVolumeRef this,
                                  xVoxelRef              ipVoxel,
                                  tBoolean               iSampled,
                                  FunV_tFunctionalValue* opValue ) {

  FunV_tErr        eResult = FunV_tErr_NoError;
  FunD_tErr        eVolume = FunD_tErr_NoError;
  float            fValue  = 0;
  float            fOffset = 0;
  xVoxelRef        RASVox  = NULL;

  xVoxl_New( &RASVox );

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* make sure we have overlay data */
  if ( NULL == this->mpOverlayVolume ) {
    eResult = FunV_tErr_OverlayNotLoaded;
    goto error;
  }

  /* if we have a cache and we're using it... */
  if ( NULL != this->mOverlayCache
       && TRUE == this->mbUseOverlayCache) {

    /* get the cached value */
    eResult = FunV_GetOverlayCacheValue_( this, ipVoxel, opValue );
    goto cleanup;
  }

  /* get the data */
  if ( iSampled ) {
    eVolume = FunD_GetSampledData( this->mpOverlayVolume, ipVoxel,
                                   this->mnCondition, this->mnTimePoint,
                                   &fValue );
  } else {
    eVolume = FunD_GetData( this->mpOverlayVolume, ipVoxel,
                            this->mnCondition, this->mnTimePoint,
                            &fValue );
  }
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_InvalidMRIIdx;
    *opValue = 0;
    goto error;
  }

  /* if we are displaying offsets and we have offset data... */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_OffsetValues]
       && NULL != this->mpOverlayOffsetVolume ) {

    /* get the offset at this value. only one plane in offset data. */
    eVolume = FunD_GetData( this->mpOverlayOffsetVolume, ipVoxel,
                            0, 0, &fOffset );
    if ( FunD_tErr_NoError == eVolume ) {

      /* divide the functional value by the offset and mult by 100 to
      get a percent */
      fValue = (fValue / fOffset) * 100.0;
    } else {
      DebugPrint( ("FunV_GetValueAtMRIIdx: %d,%d,%d: %s\n",
                   xVoxl_ExpandInt(ipVoxel), FunD_GetErrorString(eVolume) ) );
    }
  }

  /* return it */
  *opValue = (FunV_tFunctionalValue)fValue;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_GetValueAtMRIIdx: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  xVoxl_Delete( &RASVox );

  return eResult;
}

FunV_tErr FunV_InitColorCache_ ( tkmFunctionalVolumeRef this ) {

  FunV_tErr eResult         = FunV_tErr_NoError;
  int       nCacheEntry     = 0;
  float     fStep           = 0;
  float     fValue          = 0;
  float     fThreshDistance = 0;

  /* if we already have this one cached, return. */
  if ( this->mCachedThresholdMin == this->mThresholdMin
       && this->mCachedThresholdSlope == this->mThresholdSlope
       && this->mCachedThresholdMid == this->mThresholdMid ) {
    goto cleanup;
  }

  /* precalc some values */
  fThreshDistance = (1.0 / (float)(this->mThresholdSlope)) -
                    this->mThresholdMin;
  fStep = fThreshDistance / (float)FunV_knColorCacheSize;

  /* from the min to the max... */
  for ( nCacheEntry = 0; nCacheEntry < FunV_knColorCacheSize; nCacheEntry++ ) {

    /* get the value fort his step. */
    fValue = (float)nCacheEntry * fStep;

    /* cache the color */
    /*
      FunV_CalcColorValue( this, (FunV_tFunctionalValue)fValue,
      &(this->mafColorCache[nCacheEntry][0]),
      &(this->mafColorCache[nCacheEntry][2]),
      &(this->mafColorCache[nCacheEntry][3]) );
    */
  }

  /* saved cached values. */
  this->mCachedThresholdMin   = this->mThresholdMin;
  this->mCachedThresholdMid   = this->mThresholdMid;
  this->mCachedThresholdSlope = this->mThresholdSlope;

  goto cleanup;

  goto error;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_InitColorCache_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_GetColorForValue ( tkmFunctionalVolumeRef this,
                                  FunV_tFunctionalValue  iValue,
                                  xColor3fRef            iBaseColor,
                                  xColor3fRef            oColor ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  float     f       = 0; /* functional value */
  float     r       = 0; /* final color */
  float     g       = 0;
  float     b       = 0;
  float     br      = 0; /* base color */
  float     bg      = 0;
  float     bb      = 0;
  float     or      = 0; /* offset color (base color * scale) */
  float     og      = 0;
  float     ob      = 0;
  float     min     = 0; /* threshold */
  float     mid     = 0;
  float     max     = 0;
  float     tmp     = 0;
  int       bOpaque = this->mabDisplayFlags[FunV_tDisplayFlag_Ol_Opaque];

  /* set up values */
  f = (float)iValue;
  r = g = b = 0;
  br = iBaseColor->mfRed;
  bg = iBaseColor->mfGreen;
  bb = iBaseColor->mfBlue;

  /* if we're ignoreing the threshold,use the max and min values, else
     use the threshold */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_IgnoreThreshold] ) {
    FunD_GetValueRange( this->mpOverlayVolume, &min, &max );
    mid = (max-min) / 2.0;
  } else {
    min = (float)(this->mThresholdMin);
    if( bOpaque ) {
      mid = 1.001*min;
      max = 1 / (float)(this->mThresholdSlope) + min;
    }
    else {
      mid = (float)(this->mThresholdMid);
      max = 0.5 / (float)(this->mThresholdSlope) + mid;
    }
  }

  /* if we're truncating values, modify approriatly */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_TruncateNegative]
       && f < 0 ) {
    f = 0;
  }

  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_TruncatePositive]
       && f > 0 ) {
    f = 0;
  }

  /* if reversing, do so */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_ReversePhase] ) {
    f = -f;
  }

  /* at this pt, if the abs value of f is below the min, don't do
     any more processing. */
  if ( fabs(f) < min ) {
    oColor->mfRed   = br;
    oColor->mfGreen = bg;
    oColor->mfBlue  = bb;
    goto cleanup;
  }

  /* this puts values between min and mid on a nice scale */
  if ( fabs(f) > min && fabs(f) < mid ) {
    tmp = fabs(f);
    tmp = (1.0/(mid-min)) * (tmp-min)*(tmp-min) + min;
    f = (f<0) ? -tmp : tmp;
  }

  /* calc the color */
  if ( f >= 0 ) {

    /* the offset is a portion of the color that is 'blended' into
       the functional color. the rest is a standard interpolated
       color scale. */
    or = br * ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
    og = bg * ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
    ob = bb * ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
    r = or + ((f<min) ? 0.0 : (f<mid) ? (f-min)/(mid-min) : 1.0);
    g = og + ((f<mid) ? 0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
    b = ob;
    
  } else {
    f = -f;

    or = br * ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
    og = bg * ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
    ob = bb * ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
    b = ob + ((f<min) ? 0.0 : (f<mid) ? (f-min)/(mid-min) : 1.0);
    g = og + ((f<mid) ? 0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
    r = or;
  }

  /* cap values at 1 just in case */
  if ( r > 1.0 ) r = 1;
  if ( g > 1.0 ) g = 1;
  if ( b > 1.0 ) b = 1;

  /* if we're in grayscale, calc a grayscale value */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_Grayscale] ) {
    r = (r + g + b) / 3.0;
    g = b = r;
  }

  /* return values. */
  oColor->mfRed   = r;
  oColor->mfGreen = g;
  oColor->mfBlue  = b;

  goto cleanup;

cleanup:

  return eResult;
}

FunV_tErr FunV_GetAvgFunctionalValue ( tkmFunctionalVolumeRef this,
                                       FunV_tFunctionalValue* oValue,
                                       xVoxelRef              oFuncIdx,
                                       xVoxelRef              oFuncRAS,
                                       tBoolean*              obIsSelection ) {

  FunV_tErr             eResult            = FunV_tErr_NoError;
  FunD_tErr             eVolume            = FunD_tErr_NoError;
  xList_tErr            eList              = xList_tErr_NoErr;
  xVoxelRef             pVoxel             = NULL;
  FunV_tFunctionalValue value              = 0;
  FunV_tFunctionalValue sum                = 0;
  int                   nNumValues         = 0;
  FunV_tFunctionalValue average            = 0;
  xVoxel                funcIdx;
  xVoxel                funcRAS;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* make sure we have overlay data.. */
  if ( NULL == this->mpOverlayVolume ) {
    eResult =  FunV_tErr_OverlayNotLoaded;
    goto error;
  }

  /* for every voxel in our selection... */
  sum = 0;
  xList_ResetPosition( this->mpSelectedVoxels );
  eList = xList_tErr_NoErr;
  while ( xList_tErr_NoErr == eList ) {

    /* try to get a voxel */
    void* pvoid = (void*) &pVoxel;
    eList = xList_GetNextItemFromPosition( this->mpSelectedVoxels,
                                           (void**)pvoid );

    /* if we got one */
    if ( NULL != pVoxel ) {

      /* get the value to create an average. it's okay if this call
      fails, that just means one of the selected voxels are out
      of bounds. */
      eResult = FunV_GetValueAtMRIIdx( this, pVoxel, FALSE, &value );
      if ( FunV_tErr_NoError == eResult ) {

        /* add value to sum and inc count */
        sum += value;
        nNumValues++;
      } else {
        /* don't use this value, but clear the error flag. */
        eResult = FunV_tErr_NoError;
      }
    }
  }

  /* if we got something, divide sum by count */
  if ( nNumValues > 0 ) {
    average = sum / (FunV_tFunctionalValue)nNumValues;
  } else {
    average = 0;
  }

  /* return the value */
  *oValue = average;

  /* if they want the voxel location... */
  if ( NULL != oFuncIdx ||
       NULL != oFuncRAS ) {

    /* get the first voxel in the selection */
    void* pvoid = (void*) &pVoxel;
    eList = xList_GetFirstItem( this->mpSelectedVoxels, (void**)pvoid );
    if ( NULL != pVoxel ) {

      /* convert to functional index */
      FunD_ConvertClientToFuncIdx_( this->mpOverlayVolume, pVoxel, &funcIdx );

      /* if valid, return it. */
      eVolume = FunD_VerifyFuncIdx_( this->mpOverlayVolume, &funcIdx );
      if ( FunD_tErr_NoError == eVolume ) {

        if ( NULL != oFuncIdx ) {
          xVoxl_Copy( oFuncIdx, &funcIdx );
        }

        /* if they want the ras too, convert and return it. */
        if ( NULL != oFuncRAS ) {
          FunD_ConvertClientToFuncRAS_( this->mpOverlayVolume,
                                        pVoxel, &funcRAS );
          xVoxl_Copy( oFuncRAS, &funcRAS );
        }

      } else {

        /* not valid, return -1s for voxels they wanted. */
        if ( NULL != oFuncIdx ) {
          xVoxl_Set( oFuncIdx, -1, -1, -1 );
        }
        if ( NULL != oFuncRAS ) {
          xVoxl_Set( oFuncRAS, -1, -1, -1 );
        }
      }
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_GetAvgFunctionalValue: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_ConvertMRIIdxToFuncIdx ( tkmFunctionalVolumeRef this,
                                        xVoxelRef              iMRIIdx,
                                        xVoxelRef              oFuncIdx ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* make sure we have overlay data */
  if ( NULL == this->mpOverlayVolume ) {
    eResult = FunV_tErr_OverlayNotLoaded;
    goto error;
  }

  /* do the conversion */
  FunD_ConvertClientToFuncIdx_( this->mpOverlayVolume, iMRIIdx, oFuncIdx );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_ConvertMRIIdxToFuncIdx: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_ConvertMRIIdxToFuncRAS ( tkmFunctionalVolumeRef this,
                                        xVoxelRef              iMRIIdx,
                                        xVoxelRef              oFuncRAS ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* make sure we have overlay data */
  if ( NULL == this->mpOverlayVolume ) {
    eResult = FunV_tErr_OverlayNotLoaded;
    goto error;
  }

  /* do the conversion */
  FunD_ConvertClientToFuncRAS_( this->mpOverlayVolume, iMRIIdx, oFuncRAS );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_ConvertMRIIdxToFuncRAS: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_InitGraphWindow ( tkmFunctionalVolumeRef this,
                                 Tcl_Interp*            pInterp ) {

  FunV_tErr eResult               = FunV_tErr_NoError;
  tBoolean  bLookInLocalDirectory = FALSE;
  char      sFileName[256]        = "";
  FILE*     pScript               = NULL;
  int       eTcl                  = TCL_OK;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* only do this is we haven't already */
  if ( this->mbGraphInited ) {
    eResult = FunV_tErr_GraphWindowAlreadyInited;
    goto cleanup;
  }

  /* find if we can look in the local directory */
  bLookInLocalDirectory =
    ! (getenv( ksEnvVariable_UseLocalDirectoryForScript ));

  /* if we can look in the local directory... */
  if ( bLookInLocalDirectory ) {

    /* build the file name and try to open */
    sprintf( sFileName, "%s", ksFileName_InterfaceScript );
    pScript = fopen( sFileName, "r" );
  }

  /* if no file, build file name from default path and try to open. */
  if ( NULL == pScript ) {

    sprintf( sFileName, "%s%s%s", ksDir_LibraryPath, ksDir_LibrarySubPath,
             ksFileName_InterfaceScript );
    pScript = fopen( sFileName, "r" );
  }

  /* if still no file, fail. */
  if ( NULL == pScript ) {
    eResult = FunV_tErr_ErrorFindingScriptTclFile;
    goto error;
  }

  /* attempt to source the file. */
  eTcl = Tcl_EvalFile( pInterp, sFileName );
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("FunV_InitGraphWindow: error parsing %s file at line %d, result was:\n%s\n", sFileName, pInterp->errorLine, pInterp->result ) );
    eResult = FunV_tErr_ErrorParsingScriptTclFile;
    goto error;
  }

  /* set inited flag */
  this->mbGraphInited = TRUE;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_InitGraphWindow: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_BeginSelectionRange ( tkmFunctionalVolumeRef this ) {

  FunV_tErr  eResult = FunV_tErr_NoError;
  xList_tErr eList   = xList_tErr_NoErr;
  xVoxelRef   pVoxel  = NULL;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* travese the list and delete the voxels */
  xList_ResetPosition( this->mpSelectedVoxels );
  eList = xList_tErr_NoErr;
  while ( xList_tErr_NoErr == eList ) {

    /* try to get a voxel */
    void* pvoid = (void*) &pVoxel;
    eList = xList_GetNextItemFromPosition( this->mpSelectedVoxels,
                                           (void**)pvoid );

    /* if we got one */
    if ( NULL != pVoxel ) {

      /* delete it */
      xVoxl_Delete( &pVoxel );
    }
  }

  /* clear the selection list */
  eList = xList_Clear( this->mpSelectedVoxels );
  if ( xList_tErr_NoErr != eList ) {
    eResult = FunV_tErr_ErrorAccessingSelectionList;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_BeginSelectionRange: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_AddMRIIdxToSelectionRange ( tkmFunctionalVolumeRef this,
    xVoxelRef              ipVoxel ) {

  FunV_tErr  eResult = FunV_tErr_NoError;
  xVoxelRef   pVoxel  = NULL;
  xList_tErr eList   = xList_tErr_NoErr;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* create a voxel and copy into it */
  xVoxl_New( &pVoxel );
  xVoxl_Copy( pVoxel, ipVoxel );

  /* add voxel to list */
  eList = xList_InsertItem( this->mpSelectedVoxels, pVoxel );
  if ( xList_tErr_NoErr != eList ) {
    eResult = FunV_tErr_ErrorAccessingSelectionList;
    goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_AddMRIIdxToSelectionRange: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_EndSelectionRange ( tkmFunctionalVolumeRef this ) {

  FunV_tErr  eResult            = FunV_tErr_NoError;
  xList_tErr eList              = xList_tErr_NoErr;
  int        nNumValues         = 0;
  char       sTclArguments[256] = "";
  xVoxelRef  pVoxel             = NULL;
  xVoxel     funcRAS;

  DebugEnterFunction( ("FunV_EndSelectionRange( this=%p )", this) );

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  if ( NULL != this->mpTimeCourseVolume ) {

    /* find out how many selected voxels we have and set the location
       name in the time course window accordingly. */
    xList_GetCount( this->mpSelectedVoxels, &nNumValues );
    if ( 0 == nNumValues ) {
      /* Keep last update */
    } else if ( 1 == nNumValues ) {
      void* pvoid = (void*) &pVoxel;
      eList = xList_GetFirstItem( this->mpSelectedVoxels, (void**)pvoid );
      if ( NULL != pVoxel ) {
        FunD_ConvertClientToFuncRAS_( this->mpTimeCourseVolume,
                                      pVoxel, &funcRAS );
        sprintf( sTclArguments, "\"%.2f, %.2f, %.2f\"",
                 xVoxl_ExpandFloat( &funcRAS ) );
      }
    } else {
      strcpy( sTclArguments, "Selection" );
    }
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateLocationName,
                          sTclArguments );


    /* redraw the graph if we have time course data */
    eResult = FunV_DrawGraph( this );
    if ( FunV_tErr_NoError != eResult )
      goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_EndSelectionRange: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  DebugExitFunction;

  return eResult;
}
double round(double);

// DNG's attempt to circumvent the normal tkmedit plotter
// Plan:
//   1. Create a PlotInfo struct and reader 
//   2. Add plot info file/struct to tkmFunctionalVolumeRef struct
//   3. Hijack FunD_FindAndParseStemHeader_() in utils/mriFunctionalDataAccess.c
//   4. Hijack FunV_DrawGraph
//   5. Something with FunV_tErr FunV_InitGraphWindow?
//   6. Or maybe FunV_SendGraphData_?
FunV_tErr DNGFunV_DrawGraph(tkmFunctionalVolumeRef this ) ;
FunV_tErr DNGFunV_DrawGraph(tkmFunctionalVolumeRef this ) 
{
  static float afValues[1000];
  static char sTclArguments[10000];
  int n = 0, m = 0;

  xVoxelRef pVoxel  = NULL;
  FunV_tErr eResult = FunV_tErr_NoError;

  /* if graph window isn't open, bail */
  if( !(this->mbGraphInited) ) {
    eResult = FunV_tErr_GraphWindowNotInited;
    return(eResult);
  }

  this->nRawPlot = CountItemsInString(getenv("USEDNGDRAW"));
  for(n=0; n < this->nRawPlot; n++){
    this->RawPlotFiles[n] = GetNthItemFromString(getenv("USEDNGDRAW"),n);
    if(this->RawPlotFiles[n]){
      if(this->RawPlotVols[n] == NULL){
	printf("Reading %s\n",this->RawPlotFiles[n]);
	this->RawPlotVols[n] = MRIread(this->RawPlotFiles[n]);
	if(this->RawPlotVols[n] == NULL) exit(1);
      }
    }
  }

  void* pvoid = (void*) &pVoxel;
  xList_GetNextItemFromPosition( this->mpSelectedVoxels, (void**)pvoid );
  printf("%f %f %f\n",pVoxel->mfX,pVoxel->mfY,pVoxel->mfZ);

  /* send the command for starting to draw the graph */
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_BeginDrawingGraph, "" );

  for(n=0; n < this->nRawPlot; n++){
    memset(sTclArguments,'\0',strlen(sTclArguments));
    sprintf( sTclArguments, "%d {", n+1);
    for(m=0; m < this->RawPlotVols[n]->nframes; m++){ 
      afValues[m] = MRIgetVoxVal(this->RawPlotVols[n],
	 round(pVoxel->mfX),round(pVoxel->mfY),round(pVoxel->mfZ),m);
      sprintf(sTclArguments, "%s %1.1f %2.5f", 
	      sTclArguments, m*this->RawPlotVols[n]->tr/1000.0, afValues[m] );
    }
    sprintf( sTclArguments, "%s}", sTclArguments );
    FunV_SendTclCommand_(this, FunV_tTclCommand_TC_UpdateGraphData,sTclArguments );
    //FunV_SendGraphData_( this, n, 10+n, afValues );
    //FunV_SendGraphData_( this, nCondition, nNumTimePoints, afValues );
  }

  /* finish drawing */
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_EndDrawingGraph, "" );

  return(0);
}

FunV_tErr FunV_DrawGraph ( tkmFunctionalVolumeRef this ) {

  FunV_tErr             eResult         = FunV_tErr_NoError;
  FunD_tErr             eVolume         = FunD_tErr_NoError;
  tBoolean              bPresent        = FALSE;
  int                   nNumTimePoints  = 0;
  float*                afValues        = NULL;
  float*                afDeviations    = NULL;
  int                   nNumConditions  = 0;
  int                   nCondition      = 0;
  int                   nNumGoodValues  = 0;

  if(getenv("USEDNGDRAW")){
    printf("Using dng\n");
    eResult = DNGFunV_DrawGraph(this);
    return(eResult);
  }

  DebugEnterFunction( ("FunV_DrawGraph( this=%p )", this ) );

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )  goto error;

  /* if graph window isn't open, bail */
  if ( !(this->mbGraphInited) ) {
    eResult = FunV_tErr_GraphWindowNotInited;
    goto error;
  }

  /* find the number of time points. */
  eVolume = FunD_GetNumTimePoints( this->mpTimeCourseVolume,
                                   &nNumTimePoints );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* if there's only one, there's nothing to graph, so bail. */
  if ( nNumTimePoints <= 1 ) {
    goto cleanup;
  }

  /* alllocate our arrays to hold values. this will be the array we use
     to get all the values at a specific voxel. */
  afValues = (float*) malloc( nNumTimePoints *  sizeof(float) );
  if ( NULL == afValues ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* this will hold the deviations at a certain voxel */
  afDeviations = (float*) malloc( nNumTimePoints * sizeof(float) );
  if ( NULL == afDeviations ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* get the number of conditions */
  eVolume =
    FunD_GetNumConditions( this->mpTimeCourseVolume, &nNumConditions );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* send the command for starting to draw the graph */
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_BeginDrawingGraph, "" );

  /* for each condition.. */
  for ( nCondition = 0; nCondition < nNumConditions; nCondition++ ) {

    eResult = FunV_CalcTimeCourseAverages_( this, nCondition,
                                            &nNumGoodValues,
                                            afValues, afDeviations );
    if ( FunV_tErr_NoError != eResult )
      goto error;

    /* if we don't have any values at this point, our whole selections
       is out of range. */
    if ( 0 == nNumGoodValues ) {
      goto cleanup;
    }

    /* send the values to the graph */
    /* This is what actually does the plotting (dng), eg, the values
       can be changed with:
       for(nth = 0; nth < nNumTimePoints; nth++) afValues[nth] = nth;
    */
    FunV_SendGraphData_( this, nCondition, nNumTimePoints, afValues );

    /* if there is error data present.. */
    FunD_IsErrorDataPresent( this->mpTimeCourseVolume, &bPresent );
    if ( bPresent ) {

      /* send the deviations to the graph */
      FunV_SendGraphErrorBars_( this, nCondition, nNumTimePoints,
                                afDeviations );
    }
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_DrawGraph: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  /* finish drawing */
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_EndDrawingGraph, "" );

  /* dispose of float arrays */
  if ( NULL != afValues )
    free( afValues );
  if ( NULL != afDeviations )
    free( afDeviations );

  DebugExitFunction;

  return eResult;
}

FunV_tErr FunV_SendGraphData_ ( tkmFunctionalVolumeRef this,
                                int                    inCondition,
                                int                    inNumValues,
                                float*                 iafValues ) {

  FunV_tErr        eResult       = FunV_tErr_NoError;
  char*            sTclArguments = NULL;
  int              nTimePoint    = 0;
  float            fTimeSecond   = 0;
  FunD_tErr        eVolume       = FunD_tErr_NoError;

  /* allocate the argument string.  */
  sTclArguments = (char*) malloc( sizeof(char) *
                                  ((inNumValues * knLengthOfGraphDataItem) +
                                   knLengthOfGraphDataHeader) );
  if ( NULL == sTclArguments ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* write the condition number and first brace */
  sprintf( sTclArguments, "%d {", inCondition );

  /* for each time point... */
  for ( nTimePoint = 0; nTimePoint < inNumValues; nTimePoint++ ) {

    /* convert to a second. */
    eVolume = FunD_ConvertTimePointToSecond( this->mpTimeCourseVolume,
              nTimePoint, &fTimeSecond );

    /* write the second and value to the arg list */
    sprintf( sTclArguments, "%s %1.1f %2.5f", sTclArguments,
             fTimeSecond, iafValues[nTimePoint] );
  }

  /* write the last brace */
  sprintf( sTclArguments, "%s}", sTclArguments );

  /* send it to tcl */
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateGraphData,
                        sTclArguments );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SendGraphData_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  free( sTclArguments );

  return eResult;
}

FunV_tErr FunV_SendGraphErrorBars_ ( tkmFunctionalVolumeRef this,
                                     int                    inCondition,
                                     int                    inNumValues,
                                     float*                 iafBars ) {

  FunV_tErr eResult       = FunV_tErr_NoError;
  char*     sTclArguments = NULL;
  int       nTimePoint    = 0;

  /* allocate the argument string. */
  sTclArguments = (char*) malloc( sizeof(char) *
                                  ((inNumValues * knLengthOfGraphDataItem) +
                                   knLengthOfGraphDataHeader) );
  if ( NULL == sTclArguments ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* write the condition number and first brace */
  sprintf( sTclArguments, "%d {", inCondition );

  /* for each time point... */
  for ( nTimePoint = 0; nTimePoint < inNumValues; nTimePoint++ ) {

    /* write the deviation value */
    sprintf( sTclArguments, "%s %2.5f", sTclArguments, iafBars[nTimePoint] );
  }

  /* write the last brace */
  sprintf( sTclArguments, "%s}", sTclArguments );

  /* send it to tcl */
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateErrorData,
                        sTclArguments );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SendGraphErrorBars_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  free( sTclArguments );

  return eResult;
}

FunV_tErr FunV_CalcTimeCourseAverages_ ( tkmFunctionalVolumeRef this,
    int                    inCondition,
    int*                   onNumValues,
    float*                 oafValues,
    float*                 oafDeviations) {

  FunV_tErr             eResult           = FunV_tErr_NoError;
  FunD_tErr             eVolume           = FunD_tErr_NoError;
  int                   nNumValues        = 0;
  int                   nNumTimePoints    = 0;
  float*                afSums            = NULL;
  float                 fOffset           = 0;
  float                 fOffsetSum        = 0;
  xList_tErr            eList             = xList_tErr_NoErr;
  xVoxelRef             pVoxel            = NULL;
  int                   nValue            = 0;
  int                   nNumPreStimPoints = 0;
  float                 fPreStimSum       = 0;
  float                 fPreStimOffset    = 0;
  tBoolean              bPresent          = FALSE;

  DebugEnterFunction( ("FunV_CalcTimeCourseAverages_( this=%p, "
                       "inCondition=%d, onNumValues=%p, oafValues=%p, "
                       "oafDeviations=%p )", this, inCondition, onNumValues,
                       oafValues, oafDeviations) );

  /* find the number of time points. */
  eVolume = FunD_GetNumTimePoints( this->mpTimeCourseVolume,
                                   &nNumTimePoints );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* allocate sums */
  afSums = (float*) calloc( nNumTimePoints, sizeof(float) );
  if ( NULL == afSums ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* no good values yet */
  nNumValues = 0;

  /* get the values at all time points. add them
     in the sums array. keep track of how many we have. */
  xList_ResetPosition( this->mpSelectedVoxels );
  eList = xList_tErr_NoErr;
  while ( xList_tErr_NoErr == eList ) {

    /* try to get a voxel */
    void* pvoid = (void*) &pVoxel;
    eList = xList_GetNextItemFromPosition( this->mpSelectedVoxels,
                                           (void**)pvoid );

    /* if we got one */
    if ( NULL != pVoxel ) {

      /* get all values at this voxel */
      eVolume =
        FunD_GetDataForAllTimePoints( this->mpTimeCourseVolume,
                                      pVoxel, inCondition,
                                      oafValues );

      /* if it wasn't out of bounds... */
      if ( FunD_tErr_NoError == eVolume ) {

        /* if we are displaying offsets and we have offset data... */
        if ( this->mabDisplayFlags[FunV_tDisplayFlag_TC_OffsetValues]
             && NULL != this->mpTimeCourseOffsetVolume ) {

          /* get the offset at this value. only one plane in offset data. */
          eVolume =
            FunD_GetData( this->mpTimeCourseOffsetVolume,
                          pVoxel, 0, 0, &fOffset );
          if ( FunD_tErr_NoError == eVolume ) {

            /* divide all functional values by the offset and mult by 100 to
               get a percent */
            for ( nValue = 0; nValue < nNumTimePoints; nValue++ ) {
              oafValues[nValue] = (oafValues[nValue] / fOffset) * 100.0;
            }
          }
        }

        /* add all values to our sums array */
        for ( nValue = 0; nValue < nNumTimePoints; nValue++ ) {
          afSums[nValue] += oafValues[nValue];
        }


        /* if we have error data, we'll need to divide our error bar
           values by the average offset. */
        if ( this->mabDisplayFlags[FunV_tDisplayFlag_TC_OffsetValues]
             && NULL != this->mpOverlayOffsetVolume ) {

          /* get the offset at this value and add it to the sum. */
          eVolume = FunD_GetData( this->mpTimeCourseOffsetVolume,
                                  pVoxel, 0, 0, &fOffset );
          if ( FunD_tErr_NoError == eVolume ) {
            fOffsetSum += fOffset;
          }
        }

        /* inc our count */
        nNumValues++;
      }
    }
  }

  /* if we don't have any values at this point, our whole selections
     is out of range. */
  if ( 0 == nNumValues ) {

    /* return that fact */
    *onNumValues = nNumValues;
    goto cleanup;
  }

  /* divide everything by the number of values to find the average */
  for ( nValue = 0; nValue < nNumTimePoints; nValue++ ) {
    oafValues[nValue] = afSums[nValue] / (float)nNumValues;
  }

  /* if we have offset values, divide the offset sum by the number of
     values. */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_TC_OffsetValues]
       && NULL != this->mpOverlayOffsetVolume ) {
    fOffset = fOffsetSum / (float) nNumValues;
  }


  /* if we are subtracting the prestim average */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_TC_PreStimOffset] ) {

    /* get the number of prestim points */
    FunD_GetNumPreStimTimePoints( this->mpTimeCourseVolume,
                                  &nNumPreStimPoints );

    /* calc the offset by getting the sum of the values before the stim
       and dividing by the number of prestim points. */
    fPreStimSum = 0;
    for ( nValue = 0; nValue < nNumPreStimPoints; nValue++ ) {
      fPreStimSum += oafValues[nValue];
    }
    fPreStimOffset = fPreStimSum / nNumPreStimPoints;

    /* now subtract the prestim offset from all values */
    for ( nValue = 0; nValue < nNumTimePoints; nValue++ ) {
      oafValues[nValue] -= fPreStimOffset;
    }
  }


  /* if there is error data present.. */
  FunD_IsErrorDataPresent( this->mpTimeCourseVolume, &bPresent );
  if ( bPresent ) {

    /* Go through the voxel list again. */
    xList_ResetPosition( this->mpSelectedVoxels );
    eList = xList_tErr_NoErr;
    while ( xList_tErr_NoErr == eList ) {

      /* try to get a voxel */
      void* pvoid = (void*) &pVoxel;
      eList = xList_GetNextItemFromPosition( this->mpSelectedVoxels,
                                             (void**)pvoid );

      /* if we got one */
      if ( NULL != pVoxel ) {

        /* get the deviations at all time points */
        eVolume = FunD_GetDeviationForAllTimePoints( this->mpTimeCourseVolume,
                  pVoxel, inCondition,
                  oafDeviations );
        if ( FunD_tErr_NoError != eVolume ) {
          eResult = FunV_tErr_ErrorAccessingInternalVolume;
          goto error;
        }

        /* if we have offset values... */
        if ( this->mabDisplayFlags[FunV_tDisplayFlag_TC_OffsetValues]
             && NULL != this->mpTimeCourseOffsetVolume ) {

          /* divide all deviations by the offset and mult by 100 to
             get a percent */
          for ( nValue = 0; nValue < nNumTimePoints; nValue++ ) {
            oafDeviations[nValue] = (oafDeviations[nValue] / fOffset) * 100.0;
          }
        }
      }
    }
  } else {

    /* fill deviations with 0s */
    for ( nValue = 0; nValue < nNumTimePoints; nValue++ ) {
      oafDeviations[nValue] = 0;
    }
  }


  /* return the number of values */
  *onNumValues = nNumValues;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_CalcTimeCourseAverages_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  /* dispose of float arrays */
  if ( NULL != afSums )
    free( afSums );

  DebugExitFunction;

  return eResult;
}

FunV_tErr FunV_SetLocationString ( tkmFunctionalVolumeRef this,
                                   char*                  isLabel ) {

  FunV_tErr eResult                        = FunV_tErr_NoError;
  char      sTclArguments[tkm_knTclCmdLen] = "";

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* check the label */
  if ( NULL == isLabel ) {
    eResult = FunV_tErr_InvalidParameter;
    goto error;
  }

  /* check the graph */
  if ( !this->mbGraphInited )
    goto cleanup;

  /* set the label */
  xUtil_strncpy( sTclArguments, isLabel, sizeof(sTclArguments) );
  FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateLocationName,
                        sTclArguments );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SetLocationString_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_PrintSelectionRangeToFile ( tkmFunctionalVolumeRef this,
    char*                 isFileName ) {

  FunV_tErr             eResult            = FunV_tErr_NoError;
  FILE*      file = NULL;
  FunD_tErr             eVolume         = FunD_tErr_NoError;
  int                   nNumTimePoints  = 0;
  float*                afValues        = NULL;
  float*                afDeviations    = NULL;
  int                   nNumConditions  = 0;
  int                   nCondition      = 0;
  int                   nNumGoodValues  = 0;
  int                   nTimePoint = 0;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* attempt to open file */
  file = fopen( isFileName, "w" );
  if ( NULL == file ) {
    eResult = FunV_tErr_ErrorOpeningFile;
    goto error;
  }

  /* print a header */
  fprintf( file, "Voxel\t# voxels\tvalue\tstd dev\n" );

  /* find the number of time points. */
  eVolume = FunD_GetNumTimePoints( this->mpTimeCourseVolume,
                                   &nNumTimePoints );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* if there's only one, there's nothing to graph, so bail. */
  if ( nNumTimePoints <= 1 ) {
    goto cleanup;
  }

  /* alllocate our arrays to hold values. this will be the array we use
     to get all the values at a specific voxel. */
  afValues = (float*) malloc( nNumTimePoints *  sizeof(float) );
  if ( NULL == afValues ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* this will hold the deviations at a certain voxel */
  afDeviations = (float*) malloc( nNumTimePoints * sizeof(float) );
  if ( NULL == afDeviations ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* get the number of conditions */
  eVolume =
    FunD_GetNumConditions( this->mpTimeCourseVolume, &nNumConditions );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* for each condition.. */
  for ( nCondition = 0; nCondition < nNumConditions; nCondition++ ) {

    /* print condition header */
    fprintf( file, "Condition %d\n", nCondition );

    eResult = FunV_CalcTimeCourseAverages_( this, nCondition,
                                            &nNumGoodValues,
                                            afValues, afDeviations );
    if ( FunV_tErr_NoError != eResult )
      goto error;

    /* if we don't have any values at this point, our whole selections
       is out of range. */
    if ( 0 == nNumGoodValues ) {
      goto cleanup;
    }

    /* for each time point.. */
    for ( nTimePoint = 0; nTimePoint < nNumTimePoints; nTimePoint++ ) {

      /* print a line */
      fprintf( file, "%d\t%d\t%f\t%f\n",
               nTimePoint, nNumGoodValues, afValues[nTimePoint],
               afDeviations[nTimePoint] );
    }

    /* condition footer */
    fprintf( file, "End condition %d\n\n", nCondition );
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_PrintSelectionRangeToFile: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  if ( NULL != file )
    fclose( file );

  /* dispose of float arrays */
  if ( NULL != afValues )
    free( afValues );
  if ( NULL != afDeviations )
    free( afDeviations );

  return eResult;

}

FunV_tErr FunV_FloodSelect ( tkmFunctionalVolumeRef this,
                             xVoxelRef              iSeedMRIIdx,
                             tkm_tVolumeType        iVolume,
                             int                    inDistance,
                             FunV_tFindStatsComp    iCompare ) {

  FunV_tErr                     eResult      = FunV_tErr_NoError;
  Volm_tFloodParams             params;
  FunV_tFloodSelectCallbackData callbackData;
  mriVolumeRef                  sourceVolume = NULL;
  Volm_tErr                     eVolume      = Volm_tErr_NoErr;
  FunV_tFunctionalValue         funcValue    = 0;

  DebugEnterFunction( ("FunV_FloodSelect( this=%p, iSeedMRIIdx=%p, "
                       "iVolume=%d, inDistance=%d, iCompare=%d )",
                       this, iSeedMRIIdx, iVolume, inDistance,(int)iCompare) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iSeedMRIIdx),
                     eResult, tkm_tErr_InvalidParameter );

  DebugNote( ("Verifying this") );
  eResult = FunV_Verify( this );
  DebugAssertThrow( (FunV_tErr_NoError == eResult) );

  xVoxl_Copy( &params.mSourceIdx, iSeedMRIIdx );
  params.mfSourceValue           = 0;
  params.mfFuzziness             = 0;
  params.mfMaxDistance           = inDistance;
  params.mb3D                    = TRUE;

  /* Set the callback function data. Tell it to use the callback data
     we just initialized. */
  params.mpFunction     = FunV_FloodSelectCallback;
  params.mpFunctionData = (void*)&callbackData;

  /* Also pass in a comparator function because we want to do the
     value comparing ourselves, and use the the same callback data. */
  params.mComparatorFunc         = FunV_CompareMRIAndFuncValues;
  params.mComparatorFuncData     = (void*)&callbackData;

  /* Get the value at the starting voxel */
  eResult = FunV_GetValueAtMRIIdx( this, iSeedMRIIdx, FALSE, &funcValue );
  DebugAssertThrow( (FunV_tErr_NoError == eResult) );

  /* Initialize the callback data. */
  callbackData.mThis        = this;
  callbackData.mStartValue  = funcValue;
  callbackData.mCompareType = iCompare;
  callbackData.mnCount      = 0;

  /* Get the source functional volume. */
  tkm_GetAnatomicalVolume( iVolume, &sourceVolume );
  DebugAssertThrowX( (NULL != sourceVolume),
                     eResult, FunV_tErr_ErrorAccessingAnatomicalVolume );

  /* Start listening for a cancel. */
  xUtil_StartListeningForUserCancel();

  /* Do it! */
  eVolume = Volm_Flood( sourceVolume, &params );

  /* If we selected more than 1000 voxels, we printed a message and
     started printing update dots. Now close off the message. */
  if ( callbackData.mnCount > 1000 ) {
    printf( "done. %d voxels selected. \n", callbackData.mnCount );
  }

  /* Stop listening for the cancel. */
  xUtil_StopListeningForUserCancel();

  /* update the overlay */
  eResult = FunV_OverlayChanged_( this );
  DebugAssertThrow( FunV_tErr_NoError == eResult );

  DebugCatch;
  DebugCatchError( eResult, FunV_tErr_NoError, FunV_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

tBoolean FunV_CompareMRIAndFuncValues ( xVoxelRef  iMRIIdx,
                                        float      iValue,
                                        void*      ipData ) {

  FunV_tErr                      eResult      = FunV_tErr_NoError;
  tkmFunctionalVolumeRef         this         = NULL;
  FunV_tFunctionalValue          funcValue    = 0;
  tBoolean                       ebGood       = FALSE;
  FunV_tFloodSelectCallbackData* callbackData = NULL;

  callbackData = (FunV_tFloodSelectCallbackData*)ipData;
  this = callbackData->mThis;

  /* Get the functional value at this voxel. */
  DisableDebuggingOutput;
  eResult = FunV_GetValueAtMRIIdx( this, iMRIIdx, FALSE, &funcValue );
  DebugAssertThrow( (FunV_tErr_NoError == eResult) );
  EnableDebuggingOutput;

  /* Do the comparison and return the result. */
  switch ( callbackData->mCompareType ) {
  case FunV_tFindStatsComp_GTEoLTE:
    if ( callbackData->mStartValue > 0 ) {
      ebGood = (funcValue >= callbackData->mStartValue);
    } else {
      ebGood = (funcValue <= callbackData->mStartValue);
    }
    break;
  case FunV_tFindStatsComp_EQ:
    ebGood = (funcValue == callbackData->mStartValue);
    break;
  case FunV_tFindStatsComp_GTEThresholdMin:
    if ( this->mThresholdSlope > 0 ) {
      if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_ReversePhase] ) {
        ebGood = (funcValue <= -this->mThresholdMin);
      } else {
        ebGood = (funcValue >= this->mThresholdMin);
      }
    } else {
      ebGood = (funcValue <= this->mThresholdMin);
    }
    break;
  default:
    ebGood = FALSE;
  }

  DebugCatch;
  DebugCatchError( eResult, FunV_tErr_NoError, FunV_GetErrorString );
  EndDebugCatch;

  return ebGood;
}

Volm_tVisitCommand FunV_FloodSelectCallback ( xVoxelRef iMRIIdx,
    float     iValue,
    void*     iData ) {

  FunV_tFloodSelectCallbackData* callbackData;

  callbackData = (FunV_tFloodSelectCallbackData*)iData;

  /* Incremenet our count. If it's over 1000, print a message saying
      the user can cancel and start printing update dots. */
  callbackData->mnCount++;
  if ( callbackData->mnCount == 1000 ) {
    printf( "Selecting (press ctrl-c to cancel) " );
  }
  if ( callbackData->mnCount > 1000 &&
       callbackData->mnCount % 100 == 0 ) {
    printf( "." );
    fflush( stdout );
  }

  /* We only get here if our comparator, defined above, said this
     voxel was good, so just add it to the selection. */
  tkm_SelectVoxel( iMRIIdx );

  /* Check the user cancel. If they canceled, stop. */
  if ( xUtil_DidUserCancel() ) {
    return Volm_tVisitComm_Stop;
  }

  return Volm_tVisitComm_Continue;
}

int FunV_TclOlSaveRegistration ( ClientData iClientData,
                                 Tcl_Interp *ipInterp,
                                 int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* if we have overlay data... */
  if ( this->mpOverlayVolume ) {

    /* save the regsitration */
    FunD_SaveRegistration( this->mpOverlayVolume );
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {

    if ( IsDebugging ) {
      sprintf ( sError, "Error %d in FunV_TclOlSaveRegistration: %s\n",
                eResult, FunV_GetErrorString(eResult) );
      DebugPrint( (sError ) );
    } else {
      sprintf( sError, "Error %d while saving registration: %s\n",
               eResult, FunV_GetErrorString(eResult) );
    }

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

cleanup:

  return eTclResult;
}

int FunV_TclOlRestoreRegistration ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* if we have overlay data... */
  if ( this->mpOverlayVolume ) {

    /* restore the regsitration */
    FunD_RestoreRegistration( this->mpOverlayVolume );

    /* udpate the overlay */
    eResult = FunV_OverlayChanged_( this );
    if ( FunV_tErr_NoError != eResult )
      goto error;

  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {

    if ( IsDebugging ) {
      sprintf ( sError, "Error %d in FunV_TclOlRestoreRegistration: %s\n",
                eResult, FunV_GetErrorString(eResult) );
      DebugPrint( (sError ) );
    } else {
      sprintf( sError, "Error %d while restoring registration: %s\n",
               eResult, FunV_GetErrorString(eResult) );
    }

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

cleanup:

  return eTclResult;
}

int FunV_TclOlSetRegistrationToIdentity ( ClientData iClientData,
    Tcl_Interp *ipInterp,
    int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 1 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* if we have overlay data... */
  if ( this->mpOverlayVolume ) {

    /* set the regsitration to identity matrix */
    FunD_SetRegistrationToIdentity( this->mpOverlayVolume );

    /* udpate the overlay */
    eResult = FunV_OverlayChanged_( this );
    if ( FunV_tErr_NoError != eResult )
      goto error;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {

    if ( IsDebugging ) {
      sprintf ( sError, "Error %d in FunV_TclOlRestoreRegistration: %s\n",
                eResult, FunV_GetErrorString(eResult) );
      DebugPrint( (sError ) );
    } else {
      sprintf( sError, "Error %d while restoring registration: %s\n",
               eResult, FunV_GetErrorString(eResult) );
    }

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );
  }

  eTclResult = TCL_ERROR;

cleanup:

  return eTclResult;
}

int FunV_TclOlSetTimePoint ( ClientData iClientData,
                             Tcl_Interp *ipInterp, int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    nTimePoint   = 0;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  nTimePoint = rint( atof( argv[1] ) );
  eResult = FunV_SetTimePoint( this, nTimePoint );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclOlSetTimePoint: %s\n",
              eResult, FunV_GetErrorString(eResult) );
    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclOlSetCondition ( ClientData iClientData,
                             Tcl_Interp *ipInterp, int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    nCondition   = 0;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  nCondition = atoi( argv[1] );
  eResult = FunV_SetCondition( this, nCondition );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclOlSetCondition: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclOlSetDisplayFlag ( ClientData iClientData,
                               Tcl_Interp *ipInterp, int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  FunV_tDisplayFlag      flag         = FunV_tDisplayFlag_None;
  tBoolean               bValue       = FALSE;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 3 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  flag   = (FunV_tDisplayFlag) atoi( argv[1] );
  bValue = (tBoolean) atoi( argv[2] );

  /* set the value */
  eResult = FunV_SetDisplayFlag( this, flag, bValue );
  if ( FunV_tErr_NoError != eResult )
    goto error;


  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclOlSetDisplayFlag: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclOlSetThreshold ( ClientData iClientData,
                             Tcl_Interp *ipInterp, int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  FunV_tFunctionalValue  min          = 0;
  FunV_tFunctionalValue  mid          = 0;
  FunV_tFunctionalValue  slope        = 0;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 4 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  min   = (FunV_tFunctionalValue) atof( argv[1] );
  mid   = (FunV_tFunctionalValue) atof( argv[2] );
  slope = (FunV_tFunctionalValue) atof( argv[3] );

  /* set value */
  eResult = FunV_SetThreshold( this, min, mid, slope );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclOlSetThreshold: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclOlSetThresholdUsingFDR ( ClientData iClientData,
                                     Tcl_Interp *ipInterp,
                                     int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  FunV_tFunctionalValue  rate         = 0;
  tBoolean               bMask        = FALSE;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc != 3 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  rate   = (FunV_tFunctionalValue) atof( argv[1] );
  bMask  = (tBoolean) atoi( argv[2] );

  /* call functio */
  eResult = FunV_SetThresholdUsingFDR( this, rate, bMask );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclOlSetThresholdUsingFDR: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclOlSetSampleType ( ClientData iClientData,
                              Tcl_Interp *ipInterp, int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  FunD_tSampleType       sampleType   = FunD_tSampleType_Nearest;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  sampleType = (FunD_tSampleType) atoi( argv[1] );

  /* set value */
  eResult = FunD_SetSampleType( this->mpOverlayVolume, sampleType );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclOlSetSampleType: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclTCSetNumPreStimPoints ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    nNumPoints   = 0;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  nNumPoints = atoi( argv[1] );
  eResult = FunV_SetNumPreStimPoints( this, nNumPoints );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclTCSetNumPreStimPoints: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclTCSetTimeResolution ( ClientData iClientData,
                                  Tcl_Interp *ipInterp,
                                  int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  float                  fTimeRes     = 0;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 1 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  fTimeRes = atof( argv[1] );
  eResult = FunV_SetTimeResolution( this, fTimeRes );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclTCSetTimeResolution: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclTCSetDisplayFlag ( ClientData iClientData,
                               Tcl_Interp *ipInterp, int argc, char *argv[] ) {

  FunV_tErr              eResult      = FunV_tErr_NoError;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* treat this just like an overlay flag */
  eTclResult = FunV_TclOlSetDisplayFlag( iClientData, ipInterp, argc, argv );

  goto cleanup;

  goto error;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclTCSetDisplayFlag: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclTCPrintSelectionRangeToFile ( ClientData iClientData,
    Tcl_Interp *ipInterp,
    int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  char                   sFileName[256] = "";
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";

  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 2 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  strcpy( sFileName, argv[1] );
  eResult = FunV_PrintSelectionRangeToFile( this, sFileName );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  goto cleanup;

  goto error;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclTCPrintSelectionRangeToFile: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}

int FunV_TclTCPrintTimeCourseData ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] ) {

  tkmFunctionalVolumeRef this         = NULL;
  FunV_tErr              eResult      = FunV_tErr_NoError;
  FunD_tErr              eVolume         = FunD_tErr_NoError;
  xVoxel                 idx;
  xVoxel                 client;
  int                    nCondition      = 0;
  int                    eTclResult   = TCL_OK;
  char                   sError[256]  = "";
  int                    nNumTimePoints  = 0;
  int                    nTimePoint       = 0;
  float*                 afValues        = NULL;
  float*                 afDeviations    = NULL;


  /* grab us from the client data ptr */
  this = (tkmFunctionalVolumeRef) iClientData;

  /* verify us. */
  eResult = FunV_Verify ( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* verify the number of arguments. */
  if ( argc < 5 ) {
    eResult = FunV_tErr_WrongNumberArgs;
    goto error;
  }

  /* parse args */
  nCondition = atoi(argv[1]);
  xVoxl_Set( &idx, atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));


  /* we got func idx, but we have to convert to our client space to
     pass to the function. */
  FunD_ConvertFuncIdxToClient_( this->mpTimeCourseVolume, &idx, &client );

  /* find the number of time points. */
  eVolume = FunD_GetNumTimePoints( this->mpTimeCourseVolume,
                                   &nNumTimePoints );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  /* alllocate our arrays to hold values. this will be the array we use
     to get all the values at a specific voxel. */
  afValues = (float*) malloc( nNumTimePoints *  sizeof(float) );
  if ( NULL == afValues ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  /* this will hold the deviations at a certain voxel */
  afDeviations = (float*) malloc( nNumTimePoints * sizeof(float) );
  if ( NULL == afDeviations ) {
    eResult = FunV_tErr_AllocationFailed;
    goto error;
  }

  eVolume = FunD_GetDataForAllTimePoints( this->mpTimeCourseVolume,
                                          &client, nCondition, afValues );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  eVolume = FunD_GetDeviationForAllTimePoints( this->mpTimeCourseVolume,
            &client, nCondition, afDeviations );
  if ( FunD_tErr_NoError != eVolume ) {
    eResult = FunV_tErr_ErrorAccessingInternalVolume;
    goto error;
  }

  printf( "condition %d, %d,%d,%d (%.2f,%.2f,%.2f) \n tp, val, dev\n",
          nCondition, xVoxl_ExpandInt(&idx), xVoxl_ExpandFloat(&client) );

  /* for each time point.. */
  for ( nTimePoint = 0; nTimePoint < nNumTimePoints; nTimePoint++ ) {

    /* print a line */
    printf( "%d\t%f\t%f\n",
            nTimePoint, afValues[nTimePoint], afDeviations[nTimePoint] );
  }

  goto cleanup;

  goto error;

error:

  /* print error message if debugging is enabled. otherwise, the user got
     their error message already. */
  if ( IsDebugging && FunV_tErr_NoError != eResult ) {

    sprintf ( sError, "Error %d in FunV_TclTCPrintTimeCourseData: %s\n",
              eResult, FunV_GetErrorString(eResult) );

    DebugPrint( (sError ) );

    /* set tcl result, volatile so tcl will make a copy of it. */
    Tcl_SetResult( ipInterp, sError, TCL_VOLATILE );

    eTclResult = TCL_ERROR;
  }

cleanup:

  return eTclResult;
}


FunV_tErr FunV_SendViewStateToTcl ( tkmFunctionalVolumeRef this ) {

  FunV_tErr  eResult            = FunV_tErr_NoError;
  xList_tErr eList             = xList_tErr_NoErr;
  int        nValue             = 0;
  float      fValue             = 0;
  int        nFlag              = 0;
  char       sValue[256]        = "";
  char       sTclArguments[256] = "";
  float      fMin               = 0;
  float      fMax               = 0;
  xVoxelRef  pVoxel             = NULL;
  xVoxel     funcIdx;
  xVoxel     funcRAS;

  funcRAS.mfX = 0;
  funcRAS.mfY = 0;
  funcRAS.mfZ = 0;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if we have overlay data.. */
  if ( NULL != this->mpOverlayVolume ) {

    /* send the number of time points */
    FunD_GetNumTimePoints( this->mpOverlayVolume, &nValue );
    sprintf( sTclArguments, "%d", nValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateNumTimePoints,
                          sTclArguments );

    /* send the number of conditions */
    FunD_GetNumConditions( this->mpOverlayVolume, &nValue );
    sprintf( sTclArguments, "%d", nValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateNumConditions,
                          sTclArguments );

    /* send all the display flags */
    for ( nFlag = FunV_knFirstOverlayDisplayFlag;
          nFlag <= FunV_knLastOverlayDisplayFlag;
          nFlag ++ ) {
      sprintf( sTclArguments, "%d %d", nFlag, this->mabDisplayFlags[nFlag] );
      FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateDisplayFlag,
                            sTclArguments );
    }

    /* send the name and stem as the data name */
    FunD_GetSubjectName( this->mpOverlayVolume, sValue );
    sprintf( sTclArguments, "\"%s\"", sValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateDataName,
                          sTclArguments );

    /* send the current time point */
    sprintf( sTclArguments, "%d", this->mnTimePoint );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateTimePoint,
                          sTclArguments );

    /* the the current condition */
    sprintf( sTclArguments, "%d", this->mnCondition );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateCondition,
                          sTclArguments );

    /* send the threshold data */
    sprintf( sTclArguments, "%f %f %f",
             (float)(this->mThresholdMin),
             (float)(this->mThresholdMid),
             (float)(this->mThresholdSlope) );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateThreshold,
                          sTclArguments );

    /* send the volume range */
    FunD_GetValueRange( this->mpOverlayVolume, &fMin, &fMax );
    sprintf( sTclArguments, "%f %f", fMin, fMax );
    FunV_SendTclCommand_( this, FunV_tTclCommand_Ol_UpdateRange,
                          sTclArguments );
  }

  /* if we have time course data... */
  if ( NULL != this->mpTimeCourseVolume ) {

    /* Send the number of time points. */
    FunD_GetNumTimePoints( this->mpTimeCourseVolume, &nValue );
    sprintf( sTclArguments, "%d", nValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateNumTimePoints,
			  sTclArguments );

    /* send the number of conditions */
    FunD_GetNumConditions( this->mpTimeCourseVolume, &nValue );
    sprintf( sTclArguments, "%d", nValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateNumConditions,
                          sTclArguments );

    /* send the number of pre stim points */
    FunD_GetNumPreStimTimePoints( this->mpTimeCourseVolume, &nValue );
    sprintf( sTclArguments, "%d", nValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateNumPreStimPoints,
                          sTclArguments );

    /* send the time resolution */
    FunD_GetTimeResolution( this->mpTimeCourseVolume, &fValue );
    sprintf( sTclArguments, "%f", fValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateTimeResolution,
                          sTclArguments );

    /* Send time point update to the time course (we also convert the
       TP to a second and send that along too). */
    FunD_ConvertTimePointToSecond( this->mpTimeCourseVolume,
				   this->mnTimePoint, &fValue );
    sprintf( sTclArguments, "%d %f", this->mnTimePoint, fValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateTimePoint,
			  sTclArguments );

    /* send all the display flags */
    for ( nFlag = FunV_knFirstTimeCourseDisplayFlag;
          nFlag <= FunV_knLastTimeCourseDisplayFlag;
          nFlag ++ ) {
      sprintf( sTclArguments, "%d %d", nFlag, this->mabDisplayFlags[nFlag] );
      FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateDisplayFlag,
                            sTclArguments );
    }

    /* send the name and stem as the data name */
    FunD_GetSubjectName( this->mpTimeCourseVolume, sValue );
    sprintf( sTclArguments, "\"%s\"", sValue );
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateDataName,
                          sTclArguments );

    /* find out how many selected voxels we have and set the location
       name accordingly. */
    xList_GetCount( this->mpSelectedVoxels, &nValue );
    if ( nValue <= 1 ) {
      void* pvoid = (void*) &pVoxel;
      eList = xList_GetFirstItem( this->mpSelectedVoxels, (void**)pvoid );
      if ( NULL != pVoxel ) {
        FunD_ConvertClientToFuncRAS_( this->mpTimeCourseVolume,
                                      pVoxel, &funcIdx );
        sprintf( sTclArguments, "\"%.2f, %.2f, %.2f\"",
                 xVoxl_ExpandFloat( &funcRAS ) );
      }
    } else {
      strcpy( sTclArguments, "Selection" );
    }
    FunV_SendTclCommand_( this, FunV_tTclCommand_TC_UpdateLocationName,
                          sTclArguments );
  }


  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SendViewStateToTcl: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_RegisterTclCommands ( tkmFunctionalVolumeRef this,
                                     Tcl_Interp*            pInterp ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* register all our commands. */
  Tcl_CreateCommand( pInterp, "Overlay_SaveRegistration",
                     (Tcl_CmdProc*) FunV_TclOlSaveRegistration,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_RestoreRegistration",
                     (Tcl_CmdProc*) FunV_TclOlRestoreRegistration,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetRegistrationToIdentity",
                     (Tcl_CmdProc*) FunV_TclOlSetRegistrationToIdentity,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetTimePoint",
                     (Tcl_CmdProc*) FunV_TclOlSetTimePoint,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetCondition",
                     (Tcl_CmdProc*) FunV_TclOlSetCondition,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetDisplayFlag",
                     (Tcl_CmdProc*) FunV_TclOlSetDisplayFlag,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetThreshold",
                     (Tcl_CmdProc*) FunV_TclOlSetThreshold,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetThresholdUsingFDR",
                     (Tcl_CmdProc*) FunV_TclOlSetThresholdUsingFDR,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "Overlay_SetVolumeSampleType",
                     (Tcl_CmdProc*) FunV_TclOlSetSampleType,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "TimeCourse_SetNumPreStimPoints",
                     (Tcl_CmdProc*) FunV_TclTCSetNumPreStimPoints,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "TimeCourse_SetTimeResolution",
                     (Tcl_CmdProc*) FunV_TclTCSetTimeResolution,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "TimeCourse_SetDisplayFlag",
                     (Tcl_CmdProc*) FunV_TclTCSetDisplayFlag,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "TimeCourse_PrintSelectionRangeToFile",
                     (Tcl_CmdProc*) FunV_TclTCPrintSelectionRangeToFile,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );
  Tcl_CreateCommand( pInterp, "TimeCourse_PrintTimeCourseData",
                     (Tcl_CmdProc*) FunV_TclTCPrintTimeCourseData,
                     (ClientData) this, (Tcl_CmdDeleteProc *) NULL );

  /* send view state for good measure */
  eResult = FunV_SendViewStateToTcl( this );
  if ( FunV_tErr_NoError != eResult )
    FunV_Signal( "FunV_RegisterTclCommands", __LINE__, eResult );

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_RegisterTclCommands: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_GetThresholdMax ( tkmFunctionalVolumeRef this,
                                 FunV_tFunctionalValue* oValue ) {

  FunV_tErr eResult = FunV_tErr_NoError;
  float     min     = 0;
  float     max     = 0;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* if we're ignoring the threshold, get the max value of the volume,
     else use the threshold value */
  if ( this->mabDisplayFlags[FunV_tDisplayFlag_Ol_IgnoreThreshold] ) {
    FunD_GetValueRange( this->mpOverlayVolume, &min, &max );
    *oValue = max;
  } else {
    *oValue = 0.5 / this->mThresholdSlope + this->mThresholdMid;
  }

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_GetThreshold ( tkmFunctionalVolumeRef this,
                              FunV_tFunctionalValue* oMin,
                              FunV_tFunctionalValue* oMid,
                              FunV_tFunctionalValue* oSlope ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* verify us */
  eResult = FunV_Verify( this );
  if ( FunV_tErr_NoError != eResult )
    goto error;

  /* return the values */
  *oMin   = this->mThresholdMin;
  *oMid   = this->mThresholdMid;
  *oSlope = this->mThresholdSlope;

  goto cleanup;

error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_GetThreshold: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}



FunV_tErr FunV_OverlayChanged_ ( tkmFunctionalVolumeRef this ) {

  FunV_tErr              eResult            = FunV_tErr_NoError;
  xList_tErr             eList              = xList_tErr_NoErr;
  xVoxelRef               pVoxel             = NULL;
  FunV_tFunctionalValue  value              = 0;
  FunV_tErr              eValueResult       = FunV_tErr_NoError;
  char                   sTclArguments[256] = "";

  /* try to get the first selected voxel */
  void* pvoid = (void*) &pVoxel;
  eList = xList_GetFirstItem ( this->mpSelectedVoxels, (void**)pvoid );
  if ( xList_tErr_NoErr == eList ) {

    /* try to get a value for it */
    eValueResult = FunV_GetValueAtMRIIdx( this, pVoxel, FALSE, &value );
    if ( FunV_tErr_NoError == eValueResult ) {

      /* if we got one, send it to the tcl window */
      sprintf( sTclArguments, "cursor %f", (float)value );
      tkm_SendTclCommand( tkm_tTclCommand_UpdateFunctionalValue,
                          sTclArguments );
    }
  }

  /* if we have the function... */
  if ( NULL != (this->mpOverlayVolume) ) {

    /* call it */
    (this->mpOverlayChangedFunction)();
  }

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_OverlayChanged_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SendTkmeditTclCommand_ ( tkmFunctionalVolumeRef this,
                                        tkm_tTclCommand        iCommand,
                                        char*                  isArguments ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* if we have the function... */
  if ( NULL != (this->mpSendTkmeditTclCmdFunction) ) {

    /* call it */
    (this->mpSendTkmeditTclCmdFunction)( iCommand, isArguments );
  }

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SendTkmeditTclCommand_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  return eResult;
}

FunV_tErr FunV_SendTclCommand_ ( tkmFunctionalVolumeRef this,
                                 FunV_tTclCommand       iCommand,
                                 char*                  isArguments ) {

  FunV_tErr eResult   = FunV_tErr_NoError;
  char*     sCommand  = NULL;

  /* if we have the function... */
  if ( NULL != (this->mpSendTclCommandFunction) ) {

    /* allocate enough string. */
    sCommand = (char*) malloc (sizeof(char) *
                               (strlen(isArguments) + knMaxCommandLength) );

    /* build the command */
    sprintf( sCommand, "%s %s",
             FunV_ksaTclCommand[iCommand], isArguments );

    /* call it */
    (this->mpSendTclCommandFunction)( sCommand );
  }

  goto cleanup;

  goto error;
error:

  /* print error message */
  if ( FunV_tErr_NoError != eResult ) {
    DebugPrint( ("Error %d in FunV_SendTclCommand_: %s\n",
                 eResult, FunV_GetErrorString(eResult) ) );
  }

cleanup:

  free( sCommand );

  return eResult;
}

FunV_tErr FunV_Verify ( tkmFunctionalVolumeRef this ) {

  FunV_tErr eResult = FunV_tErr_NoError;

  /* check for null ptr */
  if ( NULL == this ) {
    eResult = FunV_tErr_InvalidPointer;
    goto cleanup;
  }

  /* check signature */
  if ( FunV_kSignature != this->mSignature ) {
    eResult = FunV_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;
}

void FunV_Signal ( char* isFuncName, int inLineNum, FunV_tErr ieCode ) {

  DebugPrint( ("Signal in %s, line %d: %d, %s\n",
               isFuncName, inLineNum, ieCode, FunV_GetErrorString(ieCode) ) );
}

char* FunV_GetErrorString ( FunV_tErr ieCode ) {

  FunV_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= FunV_tErr_knNumErrorCodes ) {
    eCode = FunV_tErr_InvalidErrorCode;
  }

  return FunV_ksaErrorString [eCode];
}
