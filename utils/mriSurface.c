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

  /* save the file name. */
  this->msFileName = strdup( isFileName );

  /* return us. */
  *opSurface = this;

  goto cleanup;
  
 error:

  if( NULL != this->mSurface )
    MRISfree( &this->mSurface );

  if( NULL != this )
    free( this );

  if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_New: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
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
  
  /* free file name */
  free( this->msFileName );

  /* mangle signature */
  this->mSignature = 0x1;

  /* free us */
  free( this );

  /* return null */
  *iopSurface = NULL;
    
  goto cleanup;

 error:

  if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_Delete: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
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
    DebugPrint( ("Error %d in Surf_LoadVertexSet: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
   }

 cleanup:

  return eResult;
}

Surf_tErr Surf_WriteValues ( mriSurfaceRef this,
           char*         isFileName ) {
  
  Surf_tErr eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* Write the val field */
  MRISwriteValues( this->mSurface, isFileName );
  printf( "Surface values written to %s.\n", isFileName );
  
  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_WriteValueSet: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
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
    DebugPrint( ("Error %d in Surf_IsVertexSetLoaded: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
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
 pow( (Surf_GetVertexCoord( vertex,         iSet, Surf_tOrientation_X ) -
       Surf_GetVertexCoord( neighborVertex, iSet, Surf_tOrientation_X )), 2 ) +
 pow( (Surf_GetVertexCoord( vertex,         iSet, Surf_tOrientation_Y ) -
       Surf_GetVertexCoord( neighborVertex, iSet, Surf_tOrientation_Y )), 2 ) +
 pow( (Surf_GetVertexCoord( vertex,         iSet, Surf_tOrientation_Z ) -
       Surf_GetVertexCoord( neighborVertex, iSet, Surf_tOrientation_Z )), 2 ) );
      
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
    DebugPrint( ("Error %d in Surf_ConvertSurfaceToClientSpace_: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}


Surf_tErr Surf_SetIteratorPosition ( mriSurfaceRef    this,
             xVoxelRef        iPlane ) {

  Surf_tErr     eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* see what orientation we're in and save the plane */
  if( xVoxl_GetFloatX( iPlane ) != 0 ) {
    this->mIterOrientation  = Surf_tOrientation_X;
    this->mfIterPlane       = xVoxl_GetFloatX( iPlane );
  } else if( xVoxl_GetFloatY( iPlane ) != 0 ) {
    this->mIterOrientation  = Surf_tOrientation_Y;
    this->mfIterPlane       = xVoxl_GetFloatY( iPlane );
  } else if( xVoxl_GetFloatZ( iPlane ) != 0 ) {
    this->mIterOrientation  = Surf_tOrientation_Z;
    this->mfIterPlane       = xVoxl_GetFloatZ( iPlane );
  }

  /* right now just start on first face */
  this->mnCurFace   = -1;
  this->mnCurVertex = -1;

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_SetIteratorPosition: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
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
    DebugPrint( ("Error %d in Surf_GetNextAndNeighborVertex: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}


Surf_tErr Surf_GetNthVertex ( mriSurfaceRef   this,
            Surf_tVertexSet iSet,
            int             inIndex,
            xVoxelRef       oVoxel,
            char*           osDescription ) {

  Surf_tErr     eResult = Surf_tErr_NoErr;
  vertex_type*  vertex  = NULL;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* verify the index */
  if( inIndex < 0 ||
      inIndex >= this->mSurface->nvertices ) {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }

  /* grab the voxel */
  vertex = &(this->mSurface->vertices[inIndex]);
  if( NULL == vertex ) {
    eResult = Surf_tErr_ErrorAccesssingSurface;
    goto error;
  }

  /* get the voxel in client space */
  Surf_ConvertVertexToVoxel( vertex, iSet,
           this->mTransform, oVoxel );
  
  /* make a string of info if they want it */
  if( NULL != osDescription ) {
    sprintf( osDescription, "RAS Coords: %.2f %.2f %.2f",
       Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_X ), 
       Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Y ), 
       Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Z ) );
  }

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_GetNthVertex: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;

}

Surf_tErr Surf_GetClosestVertexVoxel ( mriSurfaceRef   this,
               Surf_tVertexSet iSet,
               xVoxelRef       iClientVoxel,
               xVoxelRef       oClientVoxel,
               char*           osDescription ) {

  Surf_tErr     eResult       = Surf_tErr_NoErr;
  xVoxel        surfaceVoxel;
  vertex_type*  vertex        = NULL;
  int           nIndex        = 0;
  float         fDistance     = 0;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* since we want to return results with the vertex index, it is better
     to search for the closest vertex in the surface space, not in client
     space. so convert the target point to surface space */
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel, this->mTransform,
           &surfaceVoxel );

  /* get the closest vertex. */
  Surf_GetClosestVertex( this, iSet, &surfaceVoxel, &vertex, 
       &nIndex, &fDistance );

  /* convert it to client space */
  Surf_ConvertVertexToVoxel( vertex, iSet, 
           this->mTransform, oClientVoxel );

  /* make a string of info */
  if( NULL != osDescription ) {
    sprintf( osDescription, "Index: %d Distance: %.2f RAS Coords: %.2f %.2f %.2f", nIndex, fDistance, Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_X ), Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Y ), Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Z ) );
  }

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_GetClosestVertex: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}

float Surf_GetVertexCoord ( vertex_type*      iVertex,
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

Surf_tErr Surf_GetSurfaceSetName ( Surf_tVertexSet iSet,
           char*           osName ) {

  Surf_tErr eResult = Surf_tErr_NoErr;
  
  if( iSet < 0 || iSet >= Surf_knNumVertexSets ||
      NULL == osName  ) {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }
  
  strcpy( osName, Surf_ksaVertexSets[iSet] );
  
  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_GetSurfaceSetName: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}

Surf_tErr Surf_SetVertexValue ( mriSurfaceRef   this,
        Surf_tVertexSet iVertexSet,
        Surf_tValueSet  iValueSet,
        xVoxelRef       iClientVoxel,
        float           iValue ) {
  
  Surf_tErr eResult = Surf_tErr_NoErr;
  xVoxel    surfaceVoxel;
  VERTEX*   vertex = NULL;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* Get the closest vertex to this voxel. First get the voxel in
     surface space, then find the closest vertex. */
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel, this->mTransform,
           &surfaceVoxel );
  Surf_GetClosestVertex( this, iVertexSet, &surfaceVoxel, &vertex,
       NULL, NULL );
  if( NULL == vertex ) {
    eResult = Surf_tErr_ErrorAccesssingSurface;
    goto error;
  }
  
  /* Now set the value in the vertex. */
  switch( iValueSet ) {
  case Surf_tValueSet_Val:
    vertex->val = iValue;
    printf( "vertex %f,%f,%f val set to %f\n", 
      vertex->x, vertex->y, vertex->z, iValue );
    break;
  default:
    return Surf_tErr_InvalidParameter;
    break;
  }

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_SetVertexValue: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}

Surf_tErr Surf_GetVertexValue ( mriSurfaceRef   this,
        Surf_tVertexSet iVertexSet,
        Surf_tValueSet  iValueSet,
        xVoxelRef       iClientVoxel,
        float*          opValue ) {

  Surf_tErr eResult = Surf_tErr_NoErr;
  xVoxel    surfaceVoxel;
  VERTEX*   vertex = NULL;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* Get the closest vertex to this voxel. First get the voxel in
     surface space, then find the closest vertex. */
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel, this->mTransform,
           &surfaceVoxel );
  Surf_GetClosestVertex( this, iVertexSet, &surfaceVoxel, &vertex,
       NULL, NULL );
  if( NULL == vertex ) {
    eResult = Surf_tErr_ErrorAccesssingSurface;
    goto error;
  }
  
  /* Now get the value in the vertex. */
  switch( iValueSet ) {
  case Surf_tValueSet_Val:
    *opValue = vertex->val;
    break;
  default:
    return Surf_tErr_InvalidParameter;
    break;
  }

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_GetVertexValue: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}

Surf_tErr Surf_GetDistance ( mriSurfaceRef this,
           xVoxelRef     iClientVoxel1,
           xVoxelRef     iClientVoxel2,
           float*        ofDistance ) {

  
  Surf_tErr eResult = Surf_tErr_NoErr;
  xVoxel    surfaceVoxel1;
  xVoxel    surfaceVoxel2;
  float     fDistanceX  = 0;
  float     fDistanceY  = 0;
  float     fDistanceZ  = 0;
  float     fDistance;

  eResult = Surf_Verify( this );
  if( Surf_tErr_NoErr != eResult ) 
    goto error;

  /* Get both voxels in surface space, then calc the distance. */
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel1, this->mTransform,
           &surfaceVoxel1 );
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel2, this->mTransform,
           &surfaceVoxel2 );

  fDistanceX = xVoxl_GetFloatX(&surfaceVoxel1) - 
    xVoxl_GetFloatX(&surfaceVoxel2);
  fDistanceY = xVoxl_GetFloatY(&surfaceVoxel1) - 
    xVoxl_GetFloatY(&surfaceVoxel2);
  fDistanceZ = xVoxl_GetFloatZ(&surfaceVoxel1) - 
    xVoxl_GetFloatZ(&surfaceVoxel2);

  fDistance = sqrt( fDistanceX * fDistanceX +
        fDistanceY * fDistanceY +
        fDistanceZ * fDistanceZ );

  /* return the distance. */
  *ofDistance = fDistance;

  goto cleanup;

 error:

   if( Surf_tErr_NoErr != eResult ) {
    DebugPrint( ("Error %d in Surf_GetVertexValue: %s\n",
      eResult, Surf_GetErrorString( eResult ) ) );
  }

 cleanup:

  return eResult;
}


void Surf_GetClosestVertex ( mriSurfaceRef   this,
           Surf_tVertexSet iSet,
           xVoxelRef       iSurfaceVoxel,
           vertex_type**   opVertex,
           int*            onIndex,
           float*          ofDistance) {

  vertex_type*  currentVertex   = NULL;
  int           nVertex         = 0;
  int           nBestVertex     = -1;
  float         dx              = 0;
  float         dy              = 0;
  float         dz              = 0;
  float         fDistance       = 0;
  float         fLowestDistance = 0;

  /* start high */
  fLowestDistance = 256*256*256;

  /* for every vertex... */
  for( nVertex = 0; nVertex < this->mSurface->nvertices; nVertex++ ) {
    
    /* get the vertex */
    currentVertex = &(this->mSurface->vertices[nVertex]);

    /* calc distance */
    dx = xVoxl_GetFloatX( iSurfaceVoxel ) - 
      Surf_GetVertexCoord( currentVertex, iSet, Surf_tOrientation_X );
    dy = xVoxl_GetFloatY( iSurfaceVoxel ) - 
      Surf_GetVertexCoord( currentVertex, iSet, Surf_tOrientation_Y );
    dz = xVoxl_GetFloatZ( iSurfaceVoxel ) - 
      Surf_GetVertexCoord( currentVertex, iSet, Surf_tOrientation_Z );
    fDistance = dx*dx + dy*dy + dz*dz;

    /* save the lowest */
    if( fDistance < fLowestDistance ) {
      fLowestDistance = fDistance;
      nBestVertex = nVertex;
    }
  }

  if( -1 == nBestVertex ) {
    goto error;
  }
    
  /* get the best vertex and convert it to client space */
  *opVertex = &(this->mSurface->vertices[nBestVertex]);
  
    
  /* return the index if they want it. */
  if( NULL != onIndex ) {
    *onIndex = nBestVertex;
  }

  /* find and return the real distance if they want it. */
  if( NULL != ofDistance ) {
    *ofDistance = sqrt( fDistance );
  }    

  goto cleanup;

 error:

   *opVertex = NULL;

 cleanup:

  return;
}


void Surf_ConvertVertexToVoxel ( vertex_type*    iVertex,
         Surf_tVertexSet iSet,
         mriTransformRef iTransform,
         xVoxelRef       oVoxel ) {

  /* if we don't have a transform, just copy vertex into voxel */
  if( NULL == iTransform ) {
    xVoxl_SetFloat( oVoxel, 
        Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_X ),
        Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Y ),
        Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Z ));
  } else {

    /* stuff the vertex into a voxel */
    xVoxl_SetFloat( &sTmpVertex, 
        Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_X ),
        Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Y ),
        Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Z ));

    /* transform voxel */
    Trns_ConvertBtoA( iTransform, &sTmpVertex, oVoxel );
  }
 
}

void Surf_ConvertVoxelToSurfaceSpace ( xVoxelRef       iVoxel,
               mriTransformRef iTransform,
               xVoxelRef       oSurfVox ) {

  /* if we don't have a transform, just copy vertex into voxel */
  if( NULL == iTransform ) {
    xVoxl_Copy( oSurfVox, iVoxel );
  } else {

    /* transform voxel */
    Trns_ConvertAtoB( iTransform, iVoxel, oSurfVox );
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
