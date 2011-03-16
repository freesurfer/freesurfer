/**
 * @file  mriSurface.c
 * @brief surface utils
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/03/16 17:31:48 $
 *    $Revision: 1.39 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include "mriSurface.h"
#include "error.h"

char *Surf_ksaErrorStrings [Surf_knNumErrorCodes] =
{

  "No error.",
  "Invalid pointer to object.",
  "Invalid paramter.",
  "Invalid signature.",
  "Memory allocation failed.",
  "Error loading surface (MRISread).",
  "Error loading vertex set.",
  "Error loading annotation.",
  "Error accessing surface.",
  "Last face.",
  "Last vertex.",
  "Invalid error code."
};

char *Surf_ksaVertexSets [Surf_knNumVertexSets] =
{
  "main", "original", "pial"
};

static xVoxel sTmpVertex;

Surf_tErr Surf_New ( mriSurfaceRef*  opSurface,
                     char*           isFileName )
{

  mriSurfaceRef   this     = NULL;
  Surf_tErr       eResult  = Surf_tErr_NoErr;
  mriTransformRef transform = NULL;
  MATRIX*         identity = NULL;

  /* allocate us */
  this = (mriSurfaceRef) malloc( sizeof( mriSurface ) );
  if ( NULL == this )
  {
    eResult = Surf_tErr_AllocationFailed;
    goto error;
  }

  /* set the signature */
  this->mSignature = Surf_kSignature;

  /* read in the main surface */
  this->mSurface = MRISread( isFileName );
  if ( NULL == (this->mSurface) )
  {
    eResult = Surf_tErr_ErrorLoadingSurface;
    goto error;
  }

  /* set longest face to nothing */
  this->mfLongestEdge = 0;

  /* set the loaded flags */
  this->mabVertexSetLoaded[ Surf_tVertexSet_Main ]     = TRUE;
  this->mabVertexSetLoaded[ Surf_tVertexSet_Original ] = FALSE;
  this->mabVertexSetLoaded[ Surf_tVertexSet_Pial ]     = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Main ]     = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Original ] = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Pial ]     = FALSE;

  /* init iterators */
  this->mnCurFace   = 0;
  this->mnCurVertex = 0;

  /* allocate our faces */
  this->maFaces = (Surf_tFaceRef) calloc( this->mSurface->nfaces,
                                          sizeof( Surf_tFace ) );
  if ( NULL == this->maFaces )
  {
    eResult = Surf_tErr_AllocationFailed;
    goto error;
  }

  /* save the file name. */
  this->msFileName = strdup( isFileName );

  /* Start with identity transform. */
  identity = MatrixIdentity( 4, NULL );
  Trns_New( &transform );
  Trns_CopyAtoRAS( transform, identity );
  Trns_CopyBtoRAS( transform, identity );
  Trns_CopyARAStoBRAS( transform, identity );
  Surf_SetTransform( this, transform );

  /* return us. */
  *opSurface = this;

  goto cleanup;

error:

  if ( NULL != this->mSurface )
  {
    MRISfree( &this->mSurface );
  }

  if ( NULL != this )
  {
    free( this );
  }

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_New: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  if ( NULL != transform )
  {
    Trns_Delete( &transform );
  }

  return eResult;
}

Surf_tErr Surf_Delete ( mriSurfaceRef* iopSurface )
{

  mriSurfaceRef this    = NULL;
  Surf_tErr     eResult = Surf_tErr_NoErr;

  if ( NULL == iopSurface )
  {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }

  /* get and verify us */
  this = *iopSurface;
  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* free surface */
  MRISfree( &(this->mSurface) );

  /* free faces */
  free( this->maFaces );

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

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_Delete: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_SetTransform ( mriSurfaceRef this,
                              mriTransformRef iTransform )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* copy the transformation */
  Trns_DeepClone( iTransform, &this->mTransform );

  /* Clear our converted flags. */
  this->mabVertexSetConverted[ Surf_tVertexSet_Main ]     = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Original ] = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Pial ]     = FALSE;

  /* convert our surface */
  Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Main );
  if ( this->mabVertexSetLoaded[Surf_tVertexSet_Pial] )
  {
    Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Pial );
  }
  if ( this->mabVertexSetLoaded[Surf_tVertexSet_Original] )
  {
    Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Original );
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_SetTransform: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_LoadVertexSet ( mriSurfaceRef   this,
                               char*           isName,
                               Surf_tVertexSet iSet )
{

  Surf_tErr eResult = Surf_tErr_NoErr;
  int       eMRIS   = NO_ERROR;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* save current vertex set to tmp */
  MRISsaveVertexPositions( this->mSurface, TMP_VERTICES );

  /* read in the vertices */
  eMRIS = MRISreadVertexPositions( this->mSurface, isName );
  if ( eMRIS != NO_ERROR )
  {

    /* restore the vertices */
    MRISrestoreVertexPositions( this->mSurface, TMP_VERTICES );
    eResult = Surf_tErr_ErrorLoadingVertexSet;
    goto error;
  }

  /* save them to the correct position and restore the vertices from
     temp if necessary. */
  switch ( iSet )
  {
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
  this->mabVertexSetConverted[iSet] = FALSE;
  Surf_ConvertSurfaceToClientSpace_( this, iSet );

  /* set load status flag */
  this->mabVertexSetLoaded[ iSet ] = TRUE;

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_LoadVertexSet: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_WriteValues ( mriSurfaceRef this,
                             char*         isFileName )
{

  Surf_tErr eResult = Surf_tErr_NoErr;
  int     vno;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  for ( vno = 0; vno < this->mSurface->nvertices; vno++ )
  {
    if ( this->mSurface->vertices[vno].val != 0 )
    {
      fprintf( stderr, "vertex %d = %.2f\n", vno,
               this->mSurface->vertices[vno].val );
    }
  }

  /* Write the val field */
  MRISwriteValues( this->mSurface, isFileName );
  printf( "Surface values written to %s.\n", isFileName );

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_WriteValueSet: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_IsVertexSetLoaded ( mriSurfaceRef   this,
                                   Surf_tVertexSet iSet,
                                   tBoolean*       obIsLoaded )
{

  Surf_tErr     eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return flag */
  *obIsLoaded = this->mabVertexSetLoaded[ iSet ];

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_IsVertexSetLoaded: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_LoadAnnotation ( mriSurfaceRef   this,
                                char*           isFileName )
{

  Surf_tErr eResult = Surf_tErr_NoErr;
  int       eMRIS   = NO_ERROR;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* read in the annotation */
  eMRIS = MRISreadAnnotation( this->mSurface, isFileName );
  if ( eMRIS != NO_ERROR )
  {

    /* restore the vertices */
    eResult = Surf_tErr_ErrorLoadingAnnotation;
    goto error;
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_LoadAnnotation: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_IsInternalColorTablePresent ( mriSurfaceRef this,
    tBoolean*     obPresent )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Return whether or not we have a ct. */
  *obPresent = (this->mSurface->ct != NULL);

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_IsInternalColorTablePresent: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_NewColorTableFromInternal   ( mriSurfaceRef   this,
    COLOR_TABLE**   opTable )
{

  Surf_tErr    eResult    = Surf_tErr_NoErr;
  COLOR_TABLE* colorTable = NULL;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Make a COPY our ct, if we have one. */
  if ( this->mSurface->ct )
  {

    colorTable = CTABdeepCopy( this->mSurface->ct );
    if ( NULL == colorTable )
    {
      goto error;
    }

    *opTable = colorTable;
  }


  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_NewColorTableFromInternal: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }
  if ( NULL == colorTable )
  {
    DebugPrint( ("Error in CTABdeepCopy: Couldn't copy color table\n") );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_ConvertSurfaceToClientSpace_ ( mriSurfaceRef   this,
    Surf_tVertexSet iSet )
{

  Surf_tErr     eResult            = Surf_tErr_NoErr;
  int           nFaceIdx           = 0;
  face_type*    face               = NULL;
  Surf_tFaceRef localFace          = NULL;
  int           nVertexIdx         = 0;
  int           nNeighborVertexIdx = 0;
  vertex_type*  vertex             = NULL;
  vertex_type*  neighborVertex     = NULL;
  xVoxelRef     localVox           = NULL;
  float         dx                 = 0;
  float         dy                 = 0;
  float         dz                 = 0;
  float         fDistance          = 0;

  if ( this->mabVertexSetConverted[iSet] )
  {
    goto cleanup;
  }

  for ( nFaceIdx = 0; nFaceIdx < this->mSurface->nfaces; nFaceIdx++ )
  {

    face = &(this->mSurface->faces[nFaceIdx]);
    localFace = &(this->maFaces[nFaceIdx]);

    for ( nVertexIdx = 0; nVertexIdx < VERTICES_PER_FACE; nVertexIdx++ )
    {

      vertex = &(this->mSurface->vertices[face->v[nVertexIdx]]);
      localVox = &(localFace->mVoxel[(int)iSet][nVertexIdx]);

      /* copy the vertex into the voxel. also set the index in the face. */
      Surf_ConvertVertexToVoxel( vertex, iSet, this->mTransform, localVox );
      localFace->mnVertexIndex[(int)iSet][nVertexIdx] = face->v[nVertexIdx];

      /* get the next vertex */
      nNeighborVertexIdx = nVertexIdx==0 ? VERTICES_PER_FACE-1 : nVertexIdx-1;
      neighborVertex =&(this->mSurface->vertices[face->v[nNeighborVertexIdx]]);

      /* calc the distance */
      dx = (Surf_GetVertexCoord( vertex,         iSet, Surf_tOrientation_X ) -
            Surf_GetVertexCoord( neighborVertex, iSet, Surf_tOrientation_X ));
      dy = (Surf_GetVertexCoord( vertex,         iSet, Surf_tOrientation_Y ) -
            Surf_GetVertexCoord( neighborVertex, iSet, Surf_tOrientation_Y ));
      dz = (Surf_GetVertexCoord( vertex,         iSet, Surf_tOrientation_Z ) -
            Surf_GetVertexCoord( neighborVertex, iSet, Surf_tOrientation_Z ));
      fDistance = sqrt( dx*dx + dy*dy + dz*dz );

      /* if longer than the longest edge, set it */
      if ( fDistance > this->mfLongestEdge )
      {
        this->mfLongestEdge = fDistance;
      }
    }

    if ( !(nFaceIdx % 1000) )
    {
      fprintf( stdout, "\rConverting %s surface: %.2f%% done",
               Surf_ksaVertexSets[iSet],
               ((float)nFaceIdx / (float)this->mSurface->nfaces) * 100.0 );
    }
  }

  fprintf( stdout, "\rConverting %s surface: 100%% done.       \n",
           Surf_ksaVertexSets[iSet] );

  this->mabVertexSetConverted[iSet] = TRUE;

  goto cleanup;

  goto error;
error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_ConvertSurfaceToClientSpace_: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_SetIteratorPosition ( mriSurfaceRef    this,
                                     xVoxelRef        iPlane )
{

  Surf_tErr     eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* see what orientation we're in and save the plane */
  if ( xVoxl_GetFloatX( iPlane ) != 0 )
  {
    this->mIterOrientation  = Surf_tOrientation_X;
    this->mfIterPlane       = xVoxl_GetFloatX( iPlane );
  }
  else if ( xVoxl_GetFloatY( iPlane ) != 0 )
  {
    this->mIterOrientation  = Surf_tOrientation_Y;
    this->mfIterPlane       = xVoxl_GetFloatY( iPlane );
  }
  else if ( xVoxl_GetFloatZ( iPlane ) != 0 )
  {
    this->mIterOrientation  = Surf_tOrientation_Z;
    this->mfIterPlane       = xVoxl_GetFloatZ( iPlane );
  }

  /* right now just start on first face */
  this->mnCurFace   = -1;
  this->mnCurVertex = -1;

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_SetIteratorPosition: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_GetNextAndNeighborVertex ( mriSurfaceRef    this,
    Surf_tVertexSet iSet,
    xVoxelRef       oNextVoxel,
    int*            onNextIndex,
    xVoxelRef       oNeighborVoxel,
    int*            onNeighborIndex )
{

  Surf_tErr     eResult         = Surf_tErr_NoErr;
  tBoolean      bFaceFound      = FALSE;
  Surf_tFaceRef face            = NULL;
  float         fVertexPlane    = 0;
  int           nNeighborVertex = 0;
  xVoxelRef     vertex          = NULL;
  xVoxelRef     neighborVertex  = NULL;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if vertex is -1, search for next face. */
  if ( -1 == this->mnCurVertex )
  {

    this->mnCurFace++;
    if ( this->mnCurFace >= this->mSurface->nfaces )
    {
      eResult = Surf_tErr_LastFace;
      goto cleanup;
    }

    bFaceFound = FALSE;
    while ( this->mnCurFace < this->mSurface->nfaces &&
            !bFaceFound)
    {

      /* get the face and the first vertex on the face. */
      face = &(this->maFaces[this->mnCurFace]);
      vertex = &(face->mVoxel[iSet][0]);

      /* get the plane it's on relative to our iteration orientation */
      switch ( this->mIterOrientation )
      {
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
      if ( fabs( fVertexPlane - this->mfIterPlane ) <= this->mfLongestEdge )
      {
        bFaceFound = TRUE;
      }
      else
      {
        this->mnCurFace++;
      }
    }

    /* if not found, we're at the last face. */
    if ( this->mnCurFace >= this->mSurface->nfaces )
    {
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
  if ( nNeighborVertex < 0 )
  {
    nNeighborVertex = VERTICES_PER_FACE - 1;
  }

  /* get this vertex and the last one */
  vertex         = &(face->mVoxel[iSet][this->mnCurVertex]);
  neighborVertex = &(face->mVoxel[iSet][nNeighborVertex]);

  /* copy them out */
  xVoxl_Copy( oNextVoxel,     vertex );
  xVoxl_Copy( oNeighborVoxel, neighborVertex );

  /* if they want the indices, copy those too. */
  if ( NULL != onNextIndex )
  {
    *onNextIndex = face->mnVertexIndex[iSet][this->mnCurVertex];
  }
  if ( NULL != onNeighborIndex )
  {
    *onNeighborIndex = face->mnVertexIndex[iSet][nNeighborVertex];
  }

  /* inc current vertex. if out of bounds, set to -1. will look for new
     face next time. return Surf_tErr_LastVertex to let the client know
     that this is the last vertex in the face. */
  this->mnCurVertex++;
  if ( this->mnCurVertex >= VERTICES_PER_FACE )
  {
    this->mnCurVertex = -1;
    eResult = Surf_tErr_LastVertex;
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetNextAndNeighborVertex: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_CopyGeometryInformation ( mriSurfaceRef this,
    VOL_GEOM*     ioVolumeGeometry )
{

  Surf_tErr eResult   = Surf_tErr_NoErr;

  DebugEnterFunction( ("Surf_CopyGeometryInformation( this=%p, "
                       "ioVolumeGeometry=%p )", this, ioVolumeGeometry) );

  DebugNote( ("Verifying volume") );
  eResult = Surf_Verify( this );
  DebugAssertThrow( (eResult == Surf_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != ioVolumeGeometry),
                     eResult, Surf_tErr_InvalidParameter );

  /* Copy the information out. */
  memmove( ioVolumeGeometry, &(this->mSurface->vg), sizeof(VOL_GEOM) );

  DebugCatch;
  DebugCatchError( eResult, Surf_tErr_NoErr, Surf_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Surf_tErr Surf_TransformToVolumeGeometry ( mriSurfaceRef this,
    VOL_GEOM*     iVolumeGeometry )
{

  Surf_tErr eResult  = Surf_tErr_NoErr;
  int       eMRIS    = ERROR_NONE;
  MRI*      mri      = NULL;

  DebugEnterFunction( ("Surf_TransformToVolumeGeometry( this=%p, "
                       "iVolumeGeometry=%p )", this, iVolumeGeometry) );

  DebugNote( ("Verifying volume") );
  eResult = Surf_Verify( this );
  DebugAssertThrow( (eResult == Surf_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iVolumeGeometry),
                     eResult, Surf_tErr_InvalidParameter );

  /* Make a fake MRI from the volume geometry we got. */
  DebugNote( ("Allocting fake header.") );
  mri = MRIallocHeader( iVolumeGeometry->width,
                        iVolumeGeometry->height,
                        iVolumeGeometry->depth, MRI_VOLUME_TYPE_UNKNOWN, 1);
  DebugAssertThrowX( (NULL != mri),
                     eResult, Surf_tErr_AllocationFailed );

  /* Copy geometry information to the fake MRI. */
  DebugNote( ("Putting geometry in fake header with useVolGeomToMRI.") );
  useVolGeomToMRI( iVolumeGeometry, mri );

  /* Do the transform. */
  DebugNote( ("Running MRISsurf2surf") );
  eMRIS = MRISsurf2surf( this->mSurface, mri, NULL );
  DebugAssertThrowX( (ERROR_NONE == eMRIS),
                     eResult, Surf_tErr_ErrorAccesssingSurface );

  /* Since that changes the vertex positions, reconvert our client
     space cache. Clear our converted flags first. */
  this->mabVertexSetConverted[ Surf_tVertexSet_Main ]     = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Original ] = FALSE;
  this->mabVertexSetConverted[ Surf_tVertexSet_Pial ]     = FALSE;

  Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Main );
  if ( this->mabVertexSetLoaded[Surf_tVertexSet_Pial] )
  {
    Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Pial );
  }
  if ( this->mabVertexSetLoaded[Surf_tVertexSet_Original] )
  {
    Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Original );
  }

  DebugCatch;
  DebugCatchError( eResult, Surf_tErr_NoErr, Surf_GetErrorString );
  EndDebugCatch;

  if ( NULL != mri )
  {
    MRIfree( &mri );
  }

  DebugExitFunction;

  return eResult;
}


Surf_tErr Surf_GetNthVertex ( mriSurfaceRef   this,
                              Surf_tVertexSet iSet,
                              int             inIndex,
                              xVoxelRef       oVoxel,
                              char*           osDescription )
{

  Surf_tErr     eResult = Surf_tErr_NoErr;
  vertex_type*  vertex  = NULL;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* verify the index */
  if ( inIndex < 0 ||
       inIndex >= this->mSurface->nvertices )
  {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }

  /* grab the voxel */
  vertex = &(this->mSurface->vertices[inIndex]);
  if ( NULL == vertex )
  {
    eResult = Surf_tErr_ErrorAccesssingSurface;
    goto error;
  }

  /* get the voxel in client space */
  Surf_ConvertVertexToVoxel( vertex, iSet,
                             this->mTransform, oVoxel );

  /* make a string of info if they want it */
  if ( NULL != osDescription )
  {
    sprintf( osDescription, "RAS Coords: %.2f %.2f %.2f",
             Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_X ),
             Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Y ),
             Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Z ) );
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
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
                                       char*           osDescription )
{

  Surf_tErr     eResult       = Surf_tErr_NoErr;
  xVoxel        surfaceVoxel;
  vertex_type*  vertex        = NULL;
  int           nIndex        = 0;
  float         fDistance     = 0;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

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
  if ( NULL != osDescription )
  {
    sprintf( osDescription,
             "Index: %d Distance: %.2f RAS Coords: %.2f %.2f %.2f",
             nIndex, fDistance,
             Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_X ),
             Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Y ),
             Surf_GetVertexCoord( vertex, iSet, Surf_tOrientation_Z ) );
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetClosestVertex: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

float Surf_GetVertexCoord ( vertex_type*      iVertex,
                            Surf_tVertexSet   iSet,
                            Surf_tOrientation iOrientation )
{

  switch ( iSet )
  {
  case Surf_tVertexSet_Main:
    switch ( iOrientation )
    {
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
    switch ( iOrientation )
    {
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
    switch ( iOrientation )
    {
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
  return -1;
}

Surf_tErr Surf_SetVertexValue ( mriSurfaceRef   this,
                                Surf_tVertexSet iVertexSet,
                                Surf_tValueSet  iValueSet,
                                xVoxelRef       iClientVoxel,
                                float           iValue )
{

  Surf_tErr eResult = Surf_tErr_NoErr;
  xVoxel    surfaceVoxel;
  VERTEX*   vertex = NULL;
  int       nVno   = 0;
  float     fDistance = 0;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Get the closest vertex to this voxel. First get the voxel in
     surface space, then find the closest vertex. */
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel, this->mTransform,
                                   &surfaceVoxel );
  Surf_GetClosestVertex( this, iVertexSet, &surfaceVoxel, &vertex,
                         &nVno, &fDistance );
  if ( NULL == vertex )
  {
    eResult = Surf_tErr_ErrorAccesssingSurface;
    goto error;
  }

  /* Now set the value in the vertex. */
  switch ( iValueSet )
  {
  case Surf_tValueSet_Val:
    vertex->val = iValue;
    fprintf( stderr, "vertex %d (%.2f,%.2f,%.2f, d=%.2f) set to %.2f\n",
             nVno, vertex->x, vertex->y, vertex->z, fDistance, iValue );
    break;
  default:
    return Surf_tErr_InvalidParameter;
    break;
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
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
                                float*          opValue )
{

  Surf_tErr eResult = Surf_tErr_NoErr;
  xVoxel    surfaceVoxel;
  VERTEX*   vertex = NULL;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Get the closest vertex to this voxel. First get the voxel in
     surface space, then find the closest vertex. */
  Surf_ConvertVoxelToSurfaceSpace( iClientVoxel, this->mTransform,
                                   &surfaceVoxel );
  Surf_GetClosestVertex( this, iVertexSet, &surfaceVoxel, &vertex,
                         NULL, NULL );
  if ( NULL == vertex )
  {
    eResult = Surf_tErr_ErrorAccesssingSurface;
    goto error;
  }

  /* Now get the value in the vertex. */
  switch ( iValueSet )
  {
  case Surf_tValueSet_Val:
    *opValue = vertex->val;
    break;
  default:
    return Surf_tErr_InvalidParameter;
    break;
  }

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetVertexValue: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_GetVertexAnnotationByIndex ( mriSurfaceRef   this,
    int             iVertexIndex,
    int*            oAnnotation )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Verify the index. */
  if ( iVertexIndex < 0 || iVertexIndex >= this->mSurface->nvertices )
  {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }

  /* Return the annotation value. */
  *oAnnotation = this->mSurface->vertices[iVertexIndex].annotation;

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetVertexAnnotationByIndex: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_GetDistance ( mriSurfaceRef this,
                             xVoxelRef     iClientVoxel1,
                             xVoxelRef     iClientVoxel2,
                             float*        ofDistance )
{


  Surf_tErr eResult = Surf_tErr_NoErr;
  xVoxel    surfaceVoxel1;
  xVoxel    surfaceVoxel2;
  float     fDistanceX  = 0;
  float     fDistanceY  = 0;
  float     fDistanceZ  = 0;
  float     fDistance;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

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

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetVertexValue: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_GetMRIS ( mriSurfaceRef   this,
                         MRIS**          opSurface )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* pass back the surface. */
  *opSurface = this->mSurface;

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetMRIS: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Surf_tErr Surf_GetSurfaceSetName ( Surf_tVertexSet iSet,
                                   char*           osName )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  if ( iSet < 0 || iSet >= Surf_knNumVertexSets ||
       NULL == osName  )
  {
    eResult = Surf_tErr_InvalidParameter;
    goto error;
  }

  strcpy( osName, Surf_ksaVertexSets[iSet] );

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_GetSurfaceSetName: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Surf_tErr Surf_UsesRealRAS ( mriSurfaceRef this,
                             tBoolean*     obUseRealRAS )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Return the useRealRAS flag. */
  *obUseRealRAS = this->mSurface->useRealRAS;

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_UsesRealRAS: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }
cleanup:

  return eResult;
}

Surf_tErr Surf_AverageVertexPositions ( mriSurfaceRef this,
                                        int           inNumAverages )
{

  Surf_tErr eResult = Surf_tErr_NoErr;
  int       eMRIS   = NO_ERROR;

  eResult = Surf_Verify( this );
  if ( Surf_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* Call the MRIS function. */
  eMRIS = MRISaverageVertexPositions( this->mSurface, inNumAverages );
  if ( NO_ERROR != eMRIS )
  {
    goto error;
  }

  /* Since that changes the vertex positions, reconvert our client
     space cache. */
  Surf_ConvertSurfaceToClientSpace_( this, Surf_tVertexSet_Main );

  goto cleanup;

error:

  if ( Surf_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Surf_AverageVertexPositions: %s\n",
                 eResult, Surf_GetErrorString( eResult ) ) );
  }

  if ( NO_ERROR != eMRIS )
  {
    DebugPrint( ("Error %d (from MRIS) in Surf_AverageVertexPositions\n",
                 eMRIS) );
  }

cleanup:

  return eResult;
}

void Surf_GetClosestVertex ( mriSurfaceRef   this,
                             Surf_tVertexSet iSet,
                             xVoxelRef       iSurfaceVoxel,
                             vertex_type**   opVertex,
                             int*            onIndex,
                             float*          ofDistance)
{

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
  for ( nVertex = 0; nVertex < this->mSurface->nvertices; nVertex++ )
  {

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
    if ( fDistance < fLowestDistance )
    {
      fLowestDistance = fDistance;
      nBestVertex = nVertex;
    }
  }

  if ( -1 == nBestVertex )
  {
    goto error;
  }

  /* get the best vertex and convert it to client space */
  *opVertex = &(this->mSurface->vertices[nBestVertex]);


  /* return the index if they want it. */
  if ( NULL != onIndex )
  {
    *onIndex = nBestVertex;
  }

  /* find and return the real distance if they want it. */
  if ( NULL != ofDistance )
  {
    *ofDistance = sqrt( fLowestDistance );
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
                                 xVoxelRef       oVoxel )
{

  /* if we don't have a transform, just copy vertex into voxel */
  if ( NULL == iTransform )
  {
    oVoxel->mfX = ((iSet == Surf_tVertexSet_Main) ? iVertex->x :
                   (iSet == Surf_tVertexSet_Original) ? iVertex->origx :
                   iVertex->cx );
    oVoxel->mfY = ((iSet == Surf_tVertexSet_Main) ? iVertex->y :
                   (iSet == Surf_tVertexSet_Original) ? iVertex->origy :
                   iVertex->cy );
    oVoxel->mfZ = ((iSet == Surf_tVertexSet_Main) ? iVertex->z :
                   (iSet == Surf_tVertexSet_Original) ? iVertex->origz :
                   iVertex->cz );
#if 0
    xVoxl_SetFloat( oVoxel,
                    Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_X ),
                    Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Y ),
                    Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Z ));
#endif
  }
  else
  {

    /* stuff the vertex into a voxel */
#if 0
    xVoxl_SetFloat( &sTmpVertex,
                    Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_X ),
                    Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Y ),
                    Surf_GetVertexCoord( iVertex, iSet, Surf_tOrientation_Z ));
#endif

    sTmpVertex.mfX = ((iSet == Surf_tVertexSet_Main) ? iVertex->x :
                      (iSet == Surf_tVertexSet_Original) ? iVertex->origx :
                      iVertex->cx );
    sTmpVertex.mfY = ((iSet == Surf_tVertexSet_Main) ? iVertex->y :
                      (iSet == Surf_tVertexSet_Original) ? iVertex->origy :
                      iVertex->cy );
    sTmpVertex.mfZ = ((iSet == Surf_tVertexSet_Main) ? iVertex->z :
                      (iSet == Surf_tVertexSet_Original) ? iVertex->origz :
                      iVertex->cz );

    /* transform voxel */
    //    Trns_ConvertBtoA(iTransform, &sTmpVertex, oVoxel);
    Trns_ConvertBRAStoB(iTransform, &sTmpVertex, oVoxel);
#if 0
    printf("vertex: (%.2f, %.2f, %.2f), localVox: (%.2f, %.2f, %.2f)\n",
           sTmpVertex.mfX, sTmpVertex.mfY, sTmpVertex.mfZ,
           oVoxel->mfX, oVoxel->mfY, oVoxel->mfZ);
#endif

  }

}

void Surf_ConvertVoxelToSurfaceSpace ( xVoxelRef       iVoxel,
                                       mriTransformRef iTransform,
                                       xVoxelRef       oSurfVox )
{

#if 0
  static MATRIX * BtoRAS = NULL;
  static MATRIX * tmp1   = NULL;
  static MATRIX * tmp2   = NULL;
#endif

  /* if we don't have a transform, just copy vertex into voxel */
  if ( NULL == iTransform )
  {
    xVoxl_Copy( oSurfVox, iVoxel );
  }
  else
  {

    /* transform voxel */
    /* RKT: In the opposite of this function, CovertVertexToVoxel, we
       use CovertBRAStoB, so we'll use ConvertBtoRAS here. Even though
       at some point we should probably fix this. */
    //    Trns_ConvertAtoB( iTransform, iVoxel, oSurfVox );
    Trns_ConvertBtoRAS( iTransform, iVoxel, oSurfVox );


#if 0
    /* RKT: This doesn't work because tkmedit.c messes with the BtoRAS
       of the conversion matrix, so we'll calc it here manually. */
    /* RKT: Somehow I got it working with needing to calc a new
       BtoRAS. */
    if ( NULL == BtoRAS )
    {
      BtoRAS = MatrixInverse( iTransform->mRAStoB, NULL );
      tmp1 = MatrixAlloc( 4, 1, MATRIX_REAL );
      tmp2 = MatrixAlloc( 4, 1, MATRIX_REAL );

      DebugPrint(("iTransform->mBtoRAS\n"));
      MatrixPrint(stderr,iTransform->mBtoRAS);
      DebugPrint(("calc'd BtoRAS\n"));
      MatrixPrint(stderr,BtoRAS);
    }

    *MATRIX_RELT(tmp1,1,1) = xVoxl_GetFloatX( iVoxel );
    *MATRIX_RELT(tmp1,2,1) = xVoxl_GetFloatY( iVoxel );
    *MATRIX_RELT(tmp1,3,1) = xVoxl_GetFloatZ( iVoxel );
    *MATRIX_RELT(tmp1,4,1) = 1.0;

    MatrixMultiply( BtoRAS, tmp1, tmp2 );

    xVoxl_SetFloat( oSurfVox,
                    *MATRIX_RELT(tmp2,1,1),
                    *MATRIX_RELT(tmp2,2,1),
                    *MATRIX_RELT(tmp2,3,1) );
#endif
  }

}

Surf_tErr Surf_Verify ( mriSurfaceRef this )
{

  Surf_tErr eResult = Surf_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this )
  {
    eResult = Surf_tErr_InvalidObject;
    goto cleanup;
  }

  /* check signature */
  if ( Surf_kSignature != this->mSignature )
  {
    eResult = Surf_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;

}

char* Surf_GetErrorString ( Surf_tErr ieCode )
{

  Surf_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= Surf_knNumErrorCodes )
  {
    eCode = Surf_tErr_InvalidErrorCode;
  }

  return Surf_ksaErrorStrings [eCode];
}
