#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mriHeadPointList.h"

char HPtL_ksaErrorStrings [HPtL_knNumErrorCodes][256] = {

  "No error.",
  "Invalid object.",
  "Invalid parameter.",
  "Invalid signature.",
  "Allocation failed.",
  "Error opening or finding head point file.",
  "Error parsing head point file.",
  "Error opening or finding transform file.",
  "Error parsing transform file.",
  "Error creating transform object.",
  "Error accessing client transform object.",
  "Couldn't access transform.",
  "Last point.",
  "Invalid error code."
};

HPtL_tErr HPtL_New ( mriHeadPointListRef* oList,
         char*                isListName,
         char*                isTransformName,
         mriTransformRef      iClientTransform ) {

  mriHeadPointListRef this    = NULL;
  HPtL_tErr           eResult = HPtL_tErr_NoErr;

  /* check for valid params */
  if( NULL == isListName ) {
    eResult = HPtL_tErr_InvalidParameter;
    goto error;
  }

  /* allocate us */
  this = (mriHeadPointListRef) malloc( sizeof( mriHeadPointList ) );
  if( NULL == this ) {
    eResult = HPtL_tErr_AllocationFailed;
    goto error;
  }

  /* set signature */
  this->mSignature = HPtL_kSignature;

  /* set stuff to initial values */
  strcpy( this->msPointFile,     "" );
  this->maPoints        = NULL;
  this->mnNumPoints     = 0;
  this->mTransform      = NULL;
  this->mnCurPoint      = 0;

  /* transform name could be null */
  if( NULL != isTransformName )
    strcpy( this->msTransformFile, "" );

  /* read in the head list file */
  eResult = HPtL_ReadHeadListFile_( this,
            isListName );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* read in the transform file. makes a transform object for us with
     the clients a->ras and the transform file as the b->ras */
  eResult = HPtL_CreateTransform_( this,
           isTransformName,
           iClientTransform );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* convert all our points to client space */
  eResult = HPtL_ConvertListToClientSpace_( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* set current point to the beginning */
  this->mnCurPoint = 0;

  /* return us */
  *oList = this;

  goto cleanup;

 error:
  
  if( NULL != this )
    free( this );
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_New: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

HPtL_tErr HPtL_Delete ( mriHeadPointListRef* iopList ) {

  mriHeadPointListRef this    = NULL;
  HPtL_tErr           eResult = HPtL_tErr_NoErr;

  /* check the param. */
  if( NULL == iopList ) {
    eResult = HPtL_tErr_InvalidParameter;
    goto error;
  }
  
  /* get and verify us. */
  this = *iopList;
  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* delete our point storage */
  if( NULL != this->maPoints )
     free( this->maPoints );

  /* mangle sig */
  this->mSignature = 0x1;

  /* delete us. */
  free( this );

  /* return nil */
  *iopList = NULL;

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_Delete: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;

}
        
HPtL_tErr HPtL_ReadHeadListFile_ ( mriHeadPointListRef this,
           char*               isListName ) {

  HPtL_tErr        eResult          = HPtL_tErr_NoErr;
  FILE*            file             = NULL;
  char             sLine[1024]      = "";
  int              nNumPoints       = 0;
  HPtL_tHeadPoint* aStorage         = NULL;
  int              nPoint           = 0;
  tBoolean         bGood            = FALSE;
  char             sPointLabel[256] = "";
  int              nPointIndex      = 0;
  float            fPointX          = 0;
  float            fPointY          = 0;
  float            fPointZ          = 0;

  /* try to open the file */
  file = fopen( isListName, "r" );
  if( NULL == file ) {
    eResult = HPtL_tErr_ErrorOpeningHeadPointFile;
    goto error;
  }

  /* scan it and count the number of points */
  nNumPoints = 0;
  while( !feof( file ) ) {
    fgets( sLine, 1024, file );
    nNumPoints++;
  }

  /* allocate storage */
  aStorage = (HPtL_tHeadPoint*) calloc( nNumPoints,
          sizeof( HPtL_tHeadPoint ) );
  if( NULL == aStorage ) {
    eResult = HPtL_tErr_AllocationFailed;
    goto error;
  }

  /* for each point */
  rewind( file );
  for( nPoint = 0; nPoint < nNumPoints; nPoint++ ) {

    /* read and parse line */
    fgets( sLine, 1024, file );
    bGood = sscanf( sLine, "%s %d %f %f %f",
        sPointLabel, &nPointIndex, &fPointX, &fPointY, &fPointZ );
    if( !bGood ) {
      eResult = HPtL_tErr_ErrorParsingHeadPointFile;
      fprintf( stdout, "ERROR: Error parsing head point file %s on line %d.\n",
         isListName, nPoint );
      goto error;
    }

    /* set the data */
    strcpy( aStorage[nPoint].msLabel, sPointLabel );
    aStorage[nPoint].mnIndex = nPointIndex;
    xVoxl_SetFloat( &aStorage[nPoint].mPoint,
        fPointX, fPointY, fPointZ );
  }

  /* save the file name, number of points, and storage */
  if( NULL != this->maPoints )
    free( this->maPoints );
  this->maPoints = aStorage;
  this->mnNumPoints = nNumPoints;
  strcpy( this->msPointFile, isListName );

  goto cleanup;

 error:

  if( NULL != aStorage ) 
    free( aStorage );

  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_ReadHeadListFile_: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  if( NULL != file )
    fclose( file );

  return eResult;
}

HPtL_tErr HPtL_CreateTransform_ ( mriHeadPointListRef this,
          char*               isTransformName,
          mriTransformRef     iClientTransform ) {

  HPtL_tErr       eResult       = HPtL_tErr_NoErr;
  char            sTransformName[256] = "";
  char*           psSuffix      = NULL;
  Trns_tErr       eTransform    = Trns_tErr_NoErr;
  mriTransformRef transform     = NULL;
  FILE*           file          = NULL;
  MATRIX*         mTmp          = NULL;  
  int             nRow          = 0;
  int             nCol          = 0;
  float           fValue        = 0;
  MATRIX*         mClientAtoRAS = NULL;
  tBoolean        bGood         = FALSE;

  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* create a new transform */
  eTransform = Trns_New( &transform );
  if( Trns_tErr_NoErr != eTransform ) {
    eResult = HPtL_tErr_AllocationFailed;
    goto error;
  }

  /* if we don't have a transform name, grab the head points file name
     and replace the suffix with .trans */
  if( NULL == isTransformName ) {
    strcpy( sTransformName, this->msPointFile );
    psSuffix = strstr( sTransformName, ".hpts" );
    if( NULL != psSuffix ) {
      strcpy( psSuffix, ".trans" );
    } else {
      eResult = HPtL_tErr_InvalidParameter;
      goto error;
    }
  } else {
    strcpy( sTransformName, isTransformName );
  }

  /* try to open the transform file */
  file = fopen( sTransformName, "r" );
  if( NULL == file ) {
    eResult = HPtL_tErr_ErrorOpeningTransformFile;
    goto error;
  }

  /* create a new matrix */
  mTmp  = MatrixAlloc( 4, 4, MATRIX_REAL );
  MatrixClear( mTmp );
  
  /* read in the transform file */
  for( nRow = 1; nRow <= 4; nRow++ ) {
    for( nCol = 1; nCol <= 4; nCol++ ) {

      bGood = fscanf( file, "%f", &fValue );
      if( !bGood ) {
  eResult = HPtL_tErr_ErrorParsingTransformFile;
  goto error;
      }

      *MATRIX_RELT(mTmp,nRow,nCol) = fValue;
    }
  }

  /* copy the matrix into b->ras */
  eTransform = Trns_CopyBtoRAS( transform, mTmp );
  if( Trns_tErr_NoErr != eTransform ) {
    eResult = HPtL_tErr_ErrorCreatingTransform;
    goto error;
  }

  /* if there is a client matrix.. */
  if( NULL != iClientTransform ) {

    /* try to get their a->ras */
    eTransform = Trns_GetAtoRAS( iClientTransform, &mClientAtoRAS );
    if( Trns_tErr_NoErr != eTransform ) {
      eResult = HPtL_tErr_ErrorAccessingClientTransform;
      goto error;
    }
    
    /* if we got one... */
    if( NULL != mClientAtoRAS ) {

      /* copy in their a->ras as our a->ras */
      eTransform = Trns_CopyAtoRAS( transform, mClientAtoRAS );
      if( Trns_tErr_NoErr != eTransform ) {
  eResult = HPtL_tErr_ErrorCreatingTransform;
  goto error;
      }
    }
  }

  /* use an identity matrix as our a->b */
  MatrixIdentity( 4, mTmp );
  eTransform = Trns_CopyARAStoBRAS( transform, mTmp );
  if( Trns_tErr_NoErr != eTransform ) {
    eResult = HPtL_tErr_ErrorCreatingTransform;
    goto error;
  }

  /* set our matrix and transform name */
  if( NULL != this->mTransform ) 
    Trns_Delete( &this->mTransform );
  this->mTransform = transform;
  strcpy( this->msTransformFile, sTransformName );

  goto cleanup;

 error:
  
  if( NULL != transform )
    Trns_Delete( &transform );

  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_CreateTransform_: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  if( NULL != file )
    fclose( file );

  if( NULL != mTmp )
    MatrixFree( &mTmp );

  return eResult;
}

HPtL_tErr HPtL_ConvertListToClientSpace_ ( mriHeadPointListRef this ) {

  HPtL_tErr          eResult = HPtL_tErr_NoErr;
  HPtL_tHeadPointRef pHeadPt = NULL;

  /* reset the iterator */
  eResult = HPtL_ResetIterator( this, HPtL_tIterationPlane_All, 0, 0 );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  while( (eResult = HPtL_NextPoint( this, &pHeadPt ))
   == HPtL_tErr_NoErr ) {
  
    /* first scale the point up */
    xVoxl_SetFloat( &(pHeadPt->mClientPoint),
        xVoxl_GetFloatX( &(pHeadPt->mPoint) ) * 1000.0,
        xVoxl_GetFloatY( &(pHeadPt->mPoint) ) * 1000.0,
        xVoxl_GetFloatZ( &(pHeadPt->mPoint) ) * 1000.0 );

    /* run the transform */
    Trns_ConvertBtoA( this->mTransform, 
          &(pHeadPt->mClientPoint),
          &(pHeadPt->mClientPoint) );
  }

  if( eResult != HPtL_tErr_LastPoint )
    goto error;

  /* clear error */
  eResult = HPtL_tErr_NoErr;

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_ConvertListToClientSpace_: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:

  return eResult;
}

HPtL_tErr HPtL_WriteTransform ( mriHeadPointListRef this,
        char*               isDest ) {

  HPtL_tErr eResult         = HPtL_tErr_NoErr;
  char      sTransform[256] = "";
  MATRIX*   mTmp            = NULL;
  MATRIX*   mARAStoBRAS     = NULL;
  MATRIX*   mBRAStoARAS     = NULL;
  MATRIX*   mBtoRAS         = NULL;
  FILE*     file            = NULL;
  int       nRow            = 0;
  int       nCol            = 0;
  float     fValue          = 0;

  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* if null arg, use saved transform name */
  if( NULL == isDest ) {
    strcpy( sTransform, this->msTransformFile );
  } else {
    strcpy( sTransform, isDest );
  }

  /* get our matrices */
  Trns_GetARAStoBRAS( this->mTransform, &mARAStoBRAS );
  Trns_GetBtoRAS( this->mTransform, &mBtoRAS );

  /* we actually want BRAS->ARAS to invert it */
  mBRAStoARAS = MatrixInverse( mARAStoBRAS, NULL );

  /* compose the BRAStoARAS and BtoRAS matrices */
  mTmp = MatrixMultiply( mBRAStoARAS, mBtoRAS, NULL );

  /* attempt to open the file */
  file = fopen( sTransform, "w" );
  if( NULL == file ) {
    eResult = HPtL_tErr_ErrorOpeningTransformFile;
    goto error;
  }

  /* write the transform file */
  for( nRow = 1; nRow <= 4; nRow++ ) {
    for( nCol = 1; nCol <= 4; nCol++ ) {

      fValue = *MATRIX_RELT(mTmp,nRow,nCol);
      fprintf( file, "%f ", fValue );
    }
    fprintf( file, "\n" );
  }

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_WriteTransform: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  if( NULL != mTmp )
    MatrixFree( &mTmp );
  if( NULL != mBRAStoARAS )
    MatrixFree( &mBRAStoARAS );

  if( NULL != file )
    fclose( file );

  return eResult;
}

HPtL_tErr HPtL_WriteHeadPointFile ( mriHeadPointListRef this,
            char*               isDest ) {

  HPtL_tErr eResult         = HPtL_tErr_NoErr;
  char      sFilename[256]  = "";
  FILE*     file            = NULL;
  int       nHeadPt         = 0;

  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* if null arg, use saved name */
  if( NULL == isDest ) {
    strcpy( sFilename, this->msTransformFile );
  } else {
    strcpy( sFilename, isDest );
  }
  /* attempt to open the file */
  file = fopen( sFilename, "w" );
  if( NULL == file ) {
    eResult = HPtL_tErr_ErrorOpeningHeadPointFile;
    goto error;
  }

  /* write the list */
  for( nHeadPt = 0; nHeadPt < this->mnNumPoints; nHeadPt++ ) {

    //    fprintf( file, "%s %d % \n" );
  }

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_WriteHeadPointFile: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  if( NULL != file )
    fclose( file );

  return eResult;
}

HPtL_tErr HPtL_ResetIterator ( mriHeadPointListRef  this,
             HPtL_tIterationPlane iPlane,
             float                ifPlaneNumber,
             float                ifPlaneRange ) {

  HPtL_tErr eResult = HPtL_tErr_NoErr;

  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* save iteration desires */
  this->mIterPlane        = iPlane;
  this->mfIterPlaneNumber = ifPlaneNumber;
  this->mfIterPlaneRange  = ifPlaneRange;

  /* set cur point to -1 */
  this->mnCurPoint = -1;

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_ResetIterator: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  return eResult;
}

HPtL_tErr HPtL_NextPoint ( mriHeadPointListRef this,
         HPtL_tHeadPointRef* opPoint ) {

  HPtL_tErr          eResult    = HPtL_tErr_NoErr;
  HPtL_tHeadPointRef point      = NULL;
  tBoolean           bGoodPoint = FALSE;
  
  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* while we haven't found a good point */
  bGoodPoint = FALSE;
  while( !bGoodPoint
   && this->mnCurPoint < this->mnNumPoints ) {

    /* get the next point */
    this->mnCurPoint++;
    point = &(this->maPoints[this->mnCurPoint]);
    
    /* check it */
    switch( this->mIterPlane ) {
    case HPtL_tIterationPlane_X:
      if( abs( xVoxl_GetFloatX( &(point->mClientPoint) ) - 
         this->mfIterPlaneNumber ) <= this->mfIterPlaneRange )
  bGoodPoint = TRUE;
      break;
    case HPtL_tIterationPlane_Y:
      if( abs( xVoxl_GetFloatY( &(point->mClientPoint) ) - 
         this->mfIterPlaneNumber ) <= this->mfIterPlaneRange )
  bGoodPoint = TRUE;
      break;
    case HPtL_tIterationPlane_Z:
      if( abs( xVoxl_GetFloatZ( &(point->mClientPoint) ) - 
         this->mfIterPlaneNumber ) <= this->mfIterPlaneRange )
  bGoodPoint = TRUE;
      break;
    case HPtL_tIterationPlane_All:
      bGoodPoint = TRUE;
      break;
    default:
      break;
    }
  }

  if( this->mnCurPoint >= this->mnNumPoints ) {
    eResult = HPtL_tErr_LastPoint;
    goto cleanup;
  }

  /* return the point */
  *opPoint = point;

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_NextPoint: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  return eResult;
}
HPtL_tErr HPtL_FindNearestPoint ( mriHeadPointListRef  this,
          HPtL_tIterationPlane iPlane,
          float                ifPlaneRange,
          xVoxelRef            iWhere,
          HPtL_tHeadPointRef*  opPoint ) {

  HPtL_tErr          eResult    = HPtL_tErr_NoErr;
  float              fPlane     = 0;
  HPtL_tHeadPointRef pPoint     = NULL;
  float              fDistance  = 0;
  HPtL_tHeadPointRef pNearestPoint = NULL;
  float              fNearestDistance = 0;
  
  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* get the plane we're looking on */
  switch( iPlane ) {
  case HPtL_tIterationPlane_X:
    fPlane = xVoxl_GetFloatX( iWhere );
    break;
  case HPtL_tIterationPlane_Y:
    fPlane = xVoxl_GetFloatY( iWhere );
    break;
  case HPtL_tIterationPlane_Z:
    fPlane = xVoxl_GetFloatZ( iWhere );
    break;
  case HPtL_tIterationPlane_All:
    /* no range necessary */
    break;
  default:
    eResult = HPtL_tErr_InvalidParameter;
    goto error;
  }

  /* set iterator */
  eResult = HPtL_ResetIterator( this, iPlane, fPlane, ifPlaneRange );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* for every point... */
  fNearestDistance = 999;
  while( (eResult = HPtL_NextPoint( this, &pPoint ))
   == HPtL_tErr_NoErr ) {

    /* calc the distance in client space don't bother square rooting
       because we don't care about the real distance */
    fDistance = pow( xVoxl_GetFloatX(&(pPoint->mClientPoint)) -
         xVoxl_GetFloatX(iWhere), 2 ) +
      pow( xVoxl_GetFloatY(&(pPoint->mClientPoint)) -
     xVoxl_GetFloatY(iWhere), 2 ) +
      pow( xVoxl_GetFloatZ(&(pPoint->mClientPoint)) -
     xVoxl_GetFloatZ(iWhere), 2 );
    
    /* if less than min so far, save it */
    if( fDistance < fNearestDistance ) {
      fNearestDistance = fDistance;
      pNearestPoint    = pPoint;
    }
  }

  if( eResult != HPtL_tErr_LastPoint )
    goto error;

  /* clear error */
  eResult = HPtL_tErr_NoErr;

  /* return found point */
  *opPoint = pNearestPoint;

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_FindNearestPoint: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  return eResult;
}

HPtL_tErr HPtL_RestoreTransform ( mriHeadPointListRef this ) {

  HPtL_tErr  eResult    = HPtL_tErr_NoErr;
  Trns_tErr  eTransform = Trns_tErr_NoErr;
  MATRIX*    mTransform = NULL;

  mTransform = MatrixIdentity( 4, NULL );

  /* set the transform again */
  eTransform = Trns_CopyARAStoBRAS( this->mTransform, mTransform );

  /* reconvert all the points */
  eResult = HPtL_ConvertListToClientSpace_( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:
  
  if( Trns_tErr_NoErr != eTransform )
    eResult = HPtL_tErr_ErrorAccessingTransform;

  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_RestoreTransform: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  if( NULL != mTransform )
    MatrixFree( &mTransform );

  return eResult;
}

HPtL_tErr HPtL_ApplyTransform ( mriHeadPointListRef this,
        MATRIX*             iTransform ) {

  HPtL_tErr  eResult    = HPtL_tErr_NoErr;
  Trns_tErr  eTransform = Trns_tErr_NoErr;

  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  if( NULL == iTransform ) {
    eResult = HPtL_tErr_InvalidParameter;
    goto error;
  }

  /* apply the transform */
  eTransform = Trns_ApplyTransform( this->mTransform, iTransform );
  if( Trns_tErr_NoErr != eTransform )
    goto error;

  /* reconvert all the points */
  eResult = HPtL_ConvertListToClientSpace_( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:
  
  if( Trns_tErr_NoErr != eTransform ) 
    eResult = HPtL_tErr_ErrorAccessingTransform;

  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_ApplyTransform: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  return eResult;
}
  
HPtL_tErr HPtL_AlignPointToClientVoxel ( mriHeadPointListRef this,
           HPtL_tHeadPointRef  iPoint,
           xVoxelRef           iClientVox ) {

  HPtL_tErr  eResult      = HPtL_tErr_NoErr;
  xVoxel     destPt;
  xVoxel     srcPt;
  float      fDeltaX      = 0;
  float      fDeltaY      = 0;
  float      fDeltaZ      = 0;
  MATRIX*    mTranslation = NULL;

  if( NULL == iPoint
      || NULL == iClientVox ) {
    eResult = HPtL_tErr_InvalidParameter;
    goto error;
  }
  
  eResult = HPtL_Verify( this );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  /* convert the client vox and our head pt to ras */
  Trns_ConvertBtoRAS( this->mTransform, &(iPoint->mPoint), &srcPt );
  Trns_ConvertAtoRAS( this->mTransform, iClientVox, &destPt );

  /* get the x,y, and z difference */
  fDeltaX = xVoxl_GetFloatX( &destPt ) - xVoxl_GetFloatX( &srcPt );
  fDeltaY = xVoxl_GetFloatY( &destPt ) - xVoxl_GetFloatY( &srcPt );
  fDeltaZ = xVoxl_GetFloatZ( &destPt ) - xVoxl_GetFloatZ( &srcPt );

  DebugPrint "ana dest: %d %d %d local src %.2f %.2f %.2f\nras dest %.2f %.2f %.2f ras src %.2f %.2f %.2f delta %.2f %.2f %.2f\n",
    xVoxl_ExpandInt( iClientVox ), xVoxl_ExpandFloat( &(iPoint->mPoint) ), 
    xVoxl_ExpandFloat( &destPt ), xVoxl_ExpandFloat( &srcPt ),
    fDeltaX, fDeltaY, fDeltaZ EndDebugPrint;
    
  /* create a trans matrix */
  mTranslation = MatrixIdentity( 4, NULL );
  *MATRIX_RELT(mTranslation,1,4) = fDeltaX;
  *MATRIX_RELT(mTranslation,2,4) = fDeltaY;
  *MATRIX_RELT(mTranslation,3,4) = fDeltaZ;

  DebugPrint "Align:\n" EndDebugPrint;
  MatrixPrint( stderr, mTranslation );

  /* apply it */
  eResult = HPtL_ApplyTransform( this, mTranslation );
  if( HPtL_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:
  
  if( HPtL_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in HPtL_AlignPointToClientVoxel: %s\n",
      eResult, HPtL_GetErrorString( eResult ) EndDebugPrint;
  }
  
 cleanup:
  
  if( NULL != mTranslation )
    MatrixFree( &mTranslation );

  return eResult;
}

HPtL_tErr HPtL_Verify ( mriHeadPointListRef this ) {

  HPtL_tErr eResult = HPtL_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this ) {
    eResult = HPtL_tErr_InvalidObject;
    goto cleanup;
  }
  
  /* check signature */
  if ( HPtL_kSignature != this->mSignature ) {
    eResult = HPtL_tErr_InvalidSignature;
    goto cleanup;
  }

 cleanup:

  return eResult;

}

char* HPtL_GetErrorString ( HPtL_tErr ieCode ) {

  HPtL_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= HPtL_knNumErrorCodes ) {
    eCode = HPtL_tErr_InvalidErrorCode;
  }

  return HPtL_ksaErrorStrings [eCode];
}
