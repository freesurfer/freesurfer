#include "mriSurface.h"
#include "error.h"

char Surf_ksaErrorStrings [Surf_knNumErrorCodes][256] = {

  "No error.",
  "Invalid pointer to object.",
  "Invalid paramter.",
  "Invalid signature.",
  "Memory allocation failed.",
  "Error loading surface (MRISread).",
  "Error loading vertex set.",
  "Error accessing surface.",
  "Last face.",
  "Last vertex.",
  "Invalid error code."
};

char Surf_ksaVertexSets [Surf_knNumVertexSets][256] = {
  "main", "original", "pial"
};

static xVoxel sTmpVertex;

Surf_tErr Surf_New ( mriSurfaceRef*  opSurface,
         char*           isFileName,
         mriTransformRef iTransform ) {

  mriSurfaceRef this    = NULL;
  Surf_tErr     eResult = Surf_tErr_NoErr;

  /* allocate us */
  this = (mriSurfaceRef) malloc( sizeof( mriSurface ) );
  if( NULL == this ) {
    eResult = Surf_tErr_AllocationFailed;
    goto error;
  }
  
  /* set the signature */
  this->mSignature = Surf_kSignature;

  /* read in the main surface */
  this->mSurface = MRISread( isFileName );
  if( NULL == (this->mSurface) ) {
    eResult = Surf_tErr_ErrorLoadingSurface;
    goto error;
  }

  /* set longest face to nothing */
  this->mfLongestEdge = 0;

  /* save the transformation */
  this->mTransform = iTransform;

  /* set the loaded flags */
  this->mabVertexSetLoaded[ Surf_tVertexSet_Main ]     = TRUE; 
  this->mabVertexSetLoaded[ Surf_tVertexSet_Original ] = FALSE;
  this->mabVertexSetLoaded[ Surf_tVertexSet_Pial ]     = FALSE;

  /* init iterators */
  this->mnCurFace   = 0;
  this->mnCurVertex = 0;

  /* allocate our faces */
  this->maFaces = (Surf_tFaceRef) calloc( this->mSurface->nfaces,
            sizeof( Surf_tFace ) );
  if( NULL == this->maFaces ) {
    eResult = Surf_tErr_AllocationFailed;
    goto error;
  }

  /* convert our surface */
  Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Main );

  /* return us. */
  *opSurface = this;

  goto cleanup;
  
 error:

  if( NULL != this->mSurface )
    MRISfree( &this->mSurface );

  if( NULL != this )
    free( this );

  if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_New: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Surf_tErr Surf_Delete ( mriSurfaceRef* iopSurface ) {

  mriSurfaceRef this    = NULL;
  Surf_tErr     eResult = Surf_tErr_NoErr;

  if( NULL == iopSurface ) {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }

  /* get and verify us */
  this = *iopSurface;
  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* free surface */
  MRISfree( &(this->mSurface) );
  
  /* mangle signature */
  this->mSignature = 0x1;

  /* free us */
  free( this );

  /* return null */
  *iopSurface = NULL;
    
  goto cleanup;

 error:

  if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_Delete: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


Surf_tErr Surf_LoadVertexSet ( mriSurfaceRef   this,
             char*           isName,
             Surf_tVertexSet iSet ) {

  Surf_tErr eResult = Surf_tErr_NoErr;
  int       eMRIS   = NO_ERROR;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;


  /* save current vertex set to tmp */
  MRISsaveVertexPositions( this->mSurface, TMP_VERTICES );

  /* read in the vertices */
  eMRIS = MRISreadVertexPositions( this->mSurface, isName );
  if( eMRIS != NO_ERROR ) {
    
    /* restore the vertices */
    MRISrestoreVertexPositions( this->mSurface, TMP_VERTICES );
    eResult = Surf_tErr_ErrorLoadingVertexSet;
    goto error;
  }
  
  /* save them to the correct position and restore the vertices from
     temp if necessary. */
  switch( iSet ) {
  case Surf_tVertexSet_Main:
    /* do nothing here, leave the vertices we just read in in the main
       position */
    break;
  case Surf_tVertexSet_Original:
    MRISsaveVertexPositions( this->mSurface, ORIG_VERTICES );
    MRISrestoreVertexPositions( this->mSurface, TMP_VERTICES );
    break;
  case Surf_tVertexSet_Pial:
    MRISsaveVertexPositions( this->mSurface, CANONICAL_VERTICES );
    MRISrestoreVertexPositions( this->mSurface, TMP_VERTICES );
    break;
  default:
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }
  
  /* conert this set */
  Surf_ConvertSurfaceToClientSpace_( this, iSet );

  /* set load status flag */
  this->mabVertexSetLoaded[ iSet ] = TRUE;

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_LoadVertexSet: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
   }

 cleanup:

  return eResult;
}

Surf_tErr Surf_IsVertexSetLoaded ( mriSurfaceRef   this,
           Surf_tVertexSet iSet,
           tBoolean*       obIsLoaded ) {

  Surf_tErr     eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* return flag */
  *obIsLoaded = this->mabVertexSetLoaded[ iSet ];

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_IsVertexSetLoaded: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


Surf_tErr Surf_ConvertSurfaceToClientSpace_ ( mriSurfaceRef   this,
                Surf_tVertexSet iSet ) {

  Surf_tErr     eResult            = Surf_tErr_NoErr;
  int           nFaceIdx           = 0;
  face_type*    face               = NULL;
  Surf_tFaceRef localFace          = NULL;
  int           nVertexIdx         = 0;
  int           nNeighborVertexIdx = 0;
  vertex_type*  vertex             = NULL;
  vertex_type*  neighborVertex     = NULL;
  xVoxelRef     localVox           = NULL;
  float         fDistance          = 0;

  for( nFaceIdx = 0; nFaceIdx < this->mSurface->nfaces; nFaceIdx++ ) {

    face = &(this->mSurface->faces[nFaceIdx]);
    localFace = &(this->maFaces[nFaceIdx]);

    for( nVertexIdx = 0; nVertexIdx < VERTICES_PER_FACE; nVertexIdx++ ) {
      
      vertex = &(this->mSurface->vertices[face->v[nVertexIdx]]);
      localVox = &(localFace->mVoxel[(int)iSet][nVertexIdx]);

      Surf_ConvertVertexToVoxel( vertex, iSet, this->mTransform, localVox );

      /* get the next vertex */
      nNeighborVertexIdx = nVertexIdx==0 ? VERTICES_PER_FACE-1 : nVertexIdx-1;
      neighborVertex =&(this->mSurface->vertices[face->v[nNeighborVertexIdx]]);

      /* calc the distance */
      fDistance = sqrt( 
 pow( (Surf_GetVertexValue( vertex,         iSet, Surf_tOrientation_X ) -
       Surf_GetVertexValue( neighborVertex, iSet, Surf_tOrientation_X )), 2 ) +
 pow( (Surf_GetVertexValue( vertex,         iSet, Surf_tOrientation_Y ) -
       Surf_GetVertexValue( neighborVertex, iSet, Surf_tOrientation_Y )), 2 ) +
 pow( (Surf_GetVertexValue( vertex,         iSet, Surf_tOrientation_Z ) -
       Surf_GetVertexValue( neighborVertex, iSet, Surf_tOrientation_Z )), 2 ) );
      
      /* if longer than the longest edge, set it */
      if( fDistance > this->mfLongestEdge )
  this->mfLongestEdge = fDistance;
    }

    if( !(nFaceIdx % 1000) ) {
      fprintf( stdout, "\rConverting %s surface: %.2f%% done", 
         Surf_ksaVertexSets[iSet],
         ((float)nFaceIdx / (float)this->mSurface->nfaces) * 100.0 );
      }
  }
  
  fprintf( stdout, "\rConverting %s surface: 100%% done.       \n", 
     Surf_ksaVertexSets[iSet] );

  goto cleanup;

  goto error;
 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_ConvertSurfaceToClientSpace_: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


Surf_tErr Surf_SetIteratorPosition ( mriSurfaceRef    this,
             xVoxelRef        plane ) {

  Surf_tErr     eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* see what orientation we're in and save the plane */
  if( xVoxl_GetFloatX( plane ) != 0 ) {
    this->mIterOrientation = Surf_tOrientation_X;
    this->mfIterPlane      = xVoxl_GetFloatX( plane );
  } else if( xVoxl_GetFloatY( plane ) != 0 ) {
    this->mIterOrientation = Surf_tOrientation_Y;
    this->mfIterPlane      = xVoxl_GetFloatY( plane );
  } else if( xVoxl_GetFloatZ( plane ) != 0 ) {
    this->mIterOrientation = Surf_tOrientation_Z;
    this->mfIterPlane      = xVoxl_GetFloatZ( plane );
  }

  /* right now just start on first face */
  this->mnCurFace   = -1;
  this->mnCurVertex = -1;

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_SetIteratorPosition: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Surf_tErr Surf_GetNextAndNeighborVertex ( mriSurfaceRef    this,
             Surf_tVertexSet iSet, 
             xVoxelRef       oNextVoxel,
             xVoxelRef       oNeighborVoxel ) {

  Surf_tErr     eResult         = Surf_tErr_NoErr;
  tBoolean      bFaceFound      = FALSE;
  Surf_tFaceRef face            = NULL;
  float         fVertexPlane    = 0;
  int           nNeighborVertex = 0;
  xVoxelRef     vertex          = NULL;
  xVoxelRef     neighborVertex  = NULL;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* if vertex is -1, search for next face. */
  if( -1 == this->mnCurVertex ) {

    this->mnCurFace++;
    if( this->mnCurFace >= this->mSurface->nfaces ) {
      eResult = Surf_tErr_LastFace;
      goto cleanup;
    }

    bFaceFound = FALSE;
    while( this->mnCurFace < this->mSurface->nfaces &&
     !bFaceFound) {
      
      /* get the face and the first vertex on the face. */
      face = &(this->maFaces[this->mnCurFace]);
      vertex = &(face->mVoxel[iSet][0]);

      /* get the plane it's on relative to our iteration orientation */
      switch( this->mIterOrientation ) {
      case Surf_tOrientation_X:
  fVertexPlane = xVoxl_GetFloatX( vertex );
  break;
      case Surf_tOrientation_Y:
  fVertexPlane = xVoxl_GetFloatY( vertex );
  break;
      case Surf_tOrientation_Z:
  fVertexPlane = xVoxl_GetFloatZ( vertex );
  break;
      default:
  eResult = Surf_tErr_InvalidObject;
  goto error;
      }

      /* get the distance from the vertex to the iteration plane. if it's 
   less than some amount, use this face. */
      if( fabs( fVertexPlane - this->mfIterPlane ) <= this->mfLongestEdge ) {
  bFaceFound = TRUE;
      } else {
  this->mnCurFace++;
      }
    }

    /* if not found, we're at the last face. */
    if( this->mnCurFace >= this->mSurface->nfaces ) {
      eResult = Surf_tErr_LastFace;
      goto cleanup;
    }

    /*
      DebugPrint "Looking at face %d, %.2f,%.2f,%.2f\n",
      this->mnCurFace, xVoxl_ExpandFloat( vertex ) EndDebugPrint;
    */

    /* now we have a face to check. set the vertex to 0. */
    this->mnCurVertex = 0;
  }

  /* get the current face. */
  face = &(this->maFaces[this->mnCurFace]);

  /* get the index of the last vertex. */
  nNeighborVertex = this->mnCurVertex - 1;
  if( nNeighborVertex < 0 )
    nNeighborVertex = VERTICES_PER_FACE - 1;

  /* get this vertex and the last one */
  vertex         = &(face->mVoxel[iSet][this->mnCurVertex]);
  neighborVertex = &(face->mVoxel[iSet][nNeighborVertex]);

  /* copy them out */
  xVoxl_Copy( oNextVoxel,     vertex );
  xVoxl_Copy( oNeighborVoxel, neighborVertex );

  /* inc current vertex. if out of bounds, set to -1. will look for new
     face next time. return Surf_tErr_LastVertex to let the client know
     that this is the last vertex in the face. */
  this->mnCurVertex++;
  if( this->mnCurVertex >= VERTICES_PER_FACE ) {
    this->mnCurVertex = -1;
    eResult = Surf_tErr_LastVertex;
  }

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Surf_GetNextAndNeighborVertex: %s\n",
      eResult, Surf_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

float Surf_GetVertexValue ( vertex_type*      iVertex,
          Surf_tVertexSet   iSet,
          Surf_tOrientation iOrientation ) {

  switch( iSet ) {
  case Surf_tVertexSet_Main:
    switch( iOrientation ) {
    case Surf_tOrientation_X:
      return iVertex->x;
      break;
    case Surf_tOrientation_Y:
      return iVertex->y;
      break;
    case Surf_tOrientation_Z:
      return iVertex->z;
      break;
    default:
      return -1;
    }
    break;
  case Surf_tVertexSet_Original:
    switch( iOrientation ) {
    case Surf_tOrientation_X:
      return iVertex->origx;
      break;
    case Surf_tOrientation_Y:
      return iVertex->origy;
      break;
    case Surf_tOrientation_Z:
      return iVertex->origz;
      break;
    default:
      return -1;
    }
    break;
  case Surf_tVertexSet_Pial:
    switch( iOrientation ) {
    case Surf_tOrientation_X:
      return iVertex->cx;
      break;
    case Surf_tOrientation_Y:
      return iVertex->cy;
      break;
    case Surf_tOrientation_Z:
      return iVertex->cz;
      break;
    default:
      return -1;
    }
    break;
  default:
    return -1;
  }
}

void Surf_ConvertVertexToVoxel ( vertex_type*    iVertex,
         Surf_tVertexSet iSet,
         mriTransformRef iTransform,
         xVoxelRef       oVoxel ) {

  /* if we don't have a transform, just copy vertex into voxel */
  if( NULL == iTransform ) {
    xVoxl_SetFloat( oVoxel, 
        Surf_GetVertexValue( iVertex, iSet, Surf_tOrientation_X ),
        Surf_GetVertexValue( iVertex, iSet, Surf_tOrientation_Y ),
        Surf_GetVertexValue( iVertex, iSet, Surf_tOrientation_Z ));
  } else {

    /* stuff the vertex into a voxel */
    xVoxl_SetFloat( &sTmpVertex, 
        Surf_GetVertexValue( iVertex, iSet, Surf_tOrientation_X ),
        Surf_GetVertexValue( iVertex, iSet, Surf_tOrientation_Y ),
        Surf_GetVertexValue( iVertex, iSet, Surf_tOrientation_Z ));

    /* transform voxel */
    Trns_ConvertBtoA( iTransform, &sTmpVertex, oVoxel );
  }
 
}

Surf_tErr Surf_Verify ( mriSurfaceRef this ) {

  Surf_tErr eResult = Surf_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this ) {
    eResult = Surf_tErr_InvalidObject;
    goto cleanup;
  }
  
  /* check signature */
  if ( Surf_kSignature != this->mSignature ) {
    eResult = Surf_tErr_InvalidSignature;
    goto cleanup;
  }

 cleanup:

  return eResult;

}

char* Surf_GetErrorString ( Surf_tErr ieCode ) {

  Surf_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= Surf_knNumErrorCodes ) {
    eCode = Surf_tErr_InvalidErrorCode;
  }

  return Surf_ksaErrorStrings [eCode];
}
