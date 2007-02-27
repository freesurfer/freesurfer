/**
 * @file  vtkFSVolumeSource.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/02/27 22:37:29 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdexcept>
#include "vtkFSVolumeSource.h"
#include "vtkObjectFactory.h"
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"

using namespace std;

vtkStandardNewMacro( vtkFSVolumeSource );
vtkCxxRevisionMacro( vtkFSVolumeSource, "$Revision: 1.4 $" );

vtkFSVolumeSource::vtkFSVolumeSource () :
    mMRI( NULL ),
    mbBoundsCacheDirty( true ) {
  mImageData = vtkImageData::New();
}

void
vtkFSVolumeSource::MRIRead ( char const* ifn ) {

  char* fn = strdup( ifn );
  mMRI = ::MRIread( fn );
  free( fn );

  if ( mMRI == NULL ) {
    throw runtime_error( "MRIread failed" );
  }

  this->CopyMRIToImage();
  this->CopyMatricesFromMRI();
}

void
vtkFSVolumeSource::MRIWrite ( char const* ifn ) {

  char* fn = strdup( ifn );
  int err = ::MRIwrite( mMRI, fn );
  free( fn );

  if ( err != 0 ) {
    throw runtime_error( "MRIwrite failed" );
  }
}

void
vtkFSVolumeSource::MRIWrite () {

  int err = ::MRIwrite( mMRI, mMRI->fname );
  if ( err != 0 ) {
    throw runtime_error( "MRIwrite failed" );
  }
}

float
vtkFSVolumeSource::GetRASCenterX () {

  if ( mMRI )
    return mMRI->c_r;
  else
    return -1;
}

float
vtkFSVolumeSource::GetRASCenterY () {

  if ( mMRI )
    return mMRI->c_a;
  else
    return -1;
}

float
vtkFSVolumeSource::GetRASCenterZ () {

  if ( mMRI )
    return mMRI->c_s;
  else
    return -1;
}

void
vtkFSVolumeSource::SetVoxelToRASMatrix ( vtkMatrix4x4& iMatrix ) {

  MATRIX* m = MatrixIdentity( 4, NULL );
  for ( int r = 0; r < 4; r++ )
    for ( int c = 0; c < 4; c++ )
      *MATRIX_RELT((m),r+1,c+1) = (Real) iMatrix[r][c];

  MRIsetVoxelToRasXform( mMRI, m );

  MatrixFree( &m );

  this->CopyMRIToImage();
  this->CopyMatricesFromMRI();
}


int
vtkFSVolumeSource::ConvertIndexToRAS ( float iIdxX, float iIdxY, float iIdxZ,
                                       float& oRASX, float& oRASY, float& oRASZ  ) {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  Real ix, iy, iz, wx, wy, wz;
  int r;

  ix = iIdxX;
  iy = iIdxY;
  iz = iIdxZ;
  r = ::MRIvoxelToWorld( mMRI, ix, iy, iz, &wx, &wy, &wz );
  oRASX = wx;
  oRASY = wy;
  oRASZ = wz;

  return r;
}

int
vtkFSVolumeSource::ConvertRASToIndex( float iRASX, float iRASY, float iRASZ,
                                      float& oIdxX, float& oIdxY, float& oIdxZ ) {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  Real ix, iy, iz, wx, wy, wz;
  int r;

  wx = iRASX;
  wy = iRASY;
  wz = iRASZ;
  r = ::MRIworldToVoxel( mMRI, wx, wy, wz, &ix, &iy, &iz );
  oIdxX = ix;
  oIdxY = iy;
  oIdxZ = iz;

  return r;
}

int
vtkFSVolumeSource::ConvertRASToIndex( float iRASX, float iRASY, float iRASZ,
                                      int& oIdxX, int& oIdxY, int& oIdxZ ) {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  Real wx, wy, wz;
  int ix, iy, iz;
  int r;

  wx = iRASX;
  wy = iRASY;
  wz = iRASZ;
  r = ::MRIworldToVoxelIndex( mMRI, wx, wy, wz, &ix, &iy, &iz );
  oIdxX = ix;
  oIdxY = iy;
  oIdxZ = iz;

  return r;
}

void
vtkFSVolumeSource::Execute () {

  vtkStructuredPoints* output = this->GetOutput();
  if ( output == NULL ) {
    return;
  }

  // Set relevant information from the ImageData in the output.
  output->SetSpacing( mImageData->GetSpacing() );
  output->SetOrigin( mImageData->GetOrigin() );
  output->SetDimensions( mImageData->GetDimensions() );

  // Just pass the output a pointer to our scalars.
  output->GetPointData()->SetScalars
  ( mImageData->GetPointData()->GetScalars() );
}



void
vtkFSVolumeSource::ExecuteInformation () {

  vtkStructuredPoints* output = this->GetOutput();
  if ( output == NULL ) {
    return;
  }

  // Set relevant information from the ImageData in the output.
  output->SetWholeExtent( mImageData->GetWholeExtent() );
  output->SetScalarType( mImageData->GetScalarType() );
  output->SetNumberOfScalarComponents
  ( mImageData->GetNumberOfScalarComponents() );
  output->SetSpacing( mImageData->GetSpacing() );
  output->SetOrigin( mImageData->GetOrigin() );
}


float
vtkFSVolumeSource::GetValueAtIndex (float iIdxX, float iIdxY, float iIdxZ ) {


  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  int aDimensions[3];
  mImageData->GetDimensions( aDimensions );
  if ( iIdxX < 0 || iIdxX >= aDimensions[0] ||
       iIdxY < 0 || iIdxY >= aDimensions[1] ||
       iIdxZ < 0 || iIdxZ >= aDimensions[2] ) {
    return -1;
  }

  float oValue;
  switch ( mMRI->type ) {
  case MRI_UCHAR:
    oValue = MRIvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ );
    break;
  case MRI_INT:
    oValue = MRIIvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ );
    break;
  case MRI_LONG:
    oValue = MRILvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ );
    break;
  case MRI_FLOAT:
    oValue = MRIFvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ );
    break;
  case MRI_SHORT:
    oValue = MRISvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ );
    break;
  default:
    oValue = 0;
    break ;
  }

  return oValue;
}

void
vtkFSVolumeSource::SetValueAtIndex ( float iIdxX, float iIdxY, float iIdxZ,
                                     float iValue ) {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return;
  }

  // Make sure the index is in our bounds.
  int aDimensions[3];
  mImageData->GetDimensions( aDimensions );
  if ( iIdxX < 0 || iIdxX >= aDimensions[0] ||
       iIdxY < 0 || iIdxY >= aDimensions[1] ||
       iIdxZ < 0 || iIdxZ >= aDimensions[2] ) {
    return;
  }

  // Get our scalars.
  vtkDataArray *scalars = mImageData->GetPointData()->GetScalars();
  if ( !scalars ) {
    vtkErrorMacro( << "Couldn't get scalars" );
  }

  // Calculate a flat index from our coordinates.
  int index = (int)iIdxX +
              ((int)iIdxY * aDimensions[0]) +
              ((int)iIdxZ * aDimensions[0] * aDimensions[1]);

  // Set the value in the scalars array.
  scalars->SetComponent( index, 0, iValue );

  // We are now modified.
  this->Modified();

  // Set the value in the MRI.
  switch ( mMRI->type ) {
  case MRI_UCHAR:
    MRIvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ ) = (unsigned char)iValue;
    break;
  case MRI_INT:
    MRIIvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ ) = (int)iValue;
    break;
  case MRI_LONG:
    MRILvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ ) = (long)iValue;
    break;
  case MRI_FLOAT:
    MRIFvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ ) = iValue;
    break;
  case MRI_SHORT:
    MRISvox( mMRI, (int)iIdxX, (int)iIdxY, (int)iIdxZ ) = (short)iValue;
    break;
  default:
    break ;
  }
}

float
vtkFSVolumeSource::GetMinValue () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  float min, max;
  MRIvalRange( mMRI, &min, &max );
  return min;
}

float
vtkFSVolumeSource::GetMaxValue () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  float min, max;
  MRIvalRange( mMRI, &min, &max );
  return max;
}

float
vtkFSVolumeSource::GetPixelSizeX () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  return mMRI->xsize;
}

float
vtkFSVolumeSource::GetPixelSizeY () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  return mMRI->ysize;
}

float
vtkFSVolumeSource::GetPixelSizeZ () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  return mMRI->zsize;
}

float
vtkFSVolumeSource::GetXDimension () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  return mMRI->width;
}

float
vtkFSVolumeSource::GetYDimension () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  return mMRI->height;
}

float
vtkFSVolumeSource::GetZDimension () {

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return -1;
  }

  return mMRI->depth;
}

int*
vtkFSVolumeSource::GetDimensions () {

  if ( NULL != mImageData )
    return mImageData->GetDimensions();
  else
    return NULL;
}

void
vtkFSVolumeSource::CopyMRIToImage () {

  vtkDataArray *scalars = NULL;
  vtkUnsignedCharArray  *ucharScalars = NULL;
  vtkIntArray           *intScalars = NULL;
  vtkShortArray         *shortScalars = NULL;
  vtkLongArray          *longScalars = NULL;
  vtkFloatArray         *floatScalars = NULL;
  void* pixels;
  int cValues;
  int nZ, nY;
  int cRead;
  int zElement=0;

  if ( mMRI == NULL ) {
    vtkErrorMacro( << "No MRI is present." );
    return;
  }

  int zX = mMRI->width;
  int zY = mMRI->height;
  int zZ = mMRI->depth;
  int zFrames = mMRI->nframes;

  // This object's output space is in voxel coordinates.
  mImageData->SetDimensions( zX, zY, zZ );
  mImageData->SetSpacing( 1, 1, 1 );
  mImageData->SetOrigin( -zX/2, -zY/2, -zZ/2 );
  mImageData->SetWholeExtent( 0, zX-1, 0, zY-1, 0, zZ-1 );
  mImageData->SetNumberOfScalarComponents( zFrames );

  // create the scalars for all of the images. set the element size
  // for the data we will read.
  vtkDebugMacro (<< "Creating vtkUnsignedCharArray");
  switch ( mMRI->type ) {
  case MRI_UCHAR:
    mImageData->SetScalarTypeToUnsignedChar();
    ucharScalars = vtkUnsignedCharArray::New();
    scalars = (vtkDataArray*) ucharScalars;
    zElement = sizeof( unsigned char );
    break;
  case MRI_INT:
    mImageData->SetScalarTypeToInt();
    intScalars = vtkIntArray::New();
    scalars = (vtkDataArray*) intScalars;
    zElement = sizeof( int );
    break;
  case MRI_LONG:
    mImageData->SetScalarTypeToLong();
    longScalars = vtkLongArray::New();
    scalars = (vtkDataArray*) longScalars;
    zElement = sizeof( long );
    break;
  case MRI_FLOAT:
    mImageData->SetScalarTypeToFloat();
    floatScalars = vtkFloatArray::New();
    scalars = (vtkDataArray*) floatScalars;
    zElement = sizeof( float );
    break;
  case MRI_SHORT:
    mImageData->SetScalarTypeToShort();
    shortScalars = vtkShortArray::New();
    scalars = (vtkDataArray*) shortScalars;
    zElement = sizeof( short );
    break;
  default:
    break ;
  }
  if ( NULL == scalars ) {
    vtkErrorMacro(<< "Couldn't allocate scalars array.");
    return;
  }

  cValues = zX * zY * zZ;
  scalars->Allocate( cValues );

  // Copy the slice data into the scalars.
  vtkDebugMacro (<< "Copying " << cValues << " values into the scalars array");
  float* tuple = (float*) malloc( sizeof(float) * zFrames );
  cRead = 0;
  int nTuple = 0;
  pixels = NULL;
  for ( nZ = 0; nZ < zZ; nZ++ ) {
    for ( nY = 0; nY < zY; nY++ ) {
#if 1
      for ( int nX = 0; nX < zY; nX++ ) {
	for ( int nFrame = 0; nFrame < zFrames; nFrame++ ) {
	  
	  switch ( mMRI->type ) {
	  case MRI_UCHAR:
	    scalars->InsertComponent( nTuple, nFrame, 
				      MRIseq_vox( mMRI, nX, nY, nZ, nFrame ) );
	    break;
	  case MRI_INT:
	    scalars->InsertComponent( nTuple, nFrame, 
				      MRIIseq_vox( mMRI, nX, nY, nZ, nFrame ) );
	    break;
	  case MRI_LONG:
	    scalars->InsertComponent( nTuple, nFrame, 
				      MRILseq_vox( mMRI, nX, nY, nZ, nFrame ) );
	    break;
	  case MRI_FLOAT:
	    scalars->InsertComponent( nTuple, nFrame, 
				      MRIFseq_vox( mMRI, nX, nY, nZ, nFrame ) );
	    break;
	  case MRI_SHORT:
	    scalars->InsertComponent( nTuple, nFrame, 
				      MRISseq_vox( mMRI, nX, nY, nZ, nFrame ) );
	    break;
	  default:
	    break;
	  }
	  
	  nTuple++;
	}
      }
      
#else
      // Not sure how this will work when reading in frames.
      vtkDebugMacro (<< "Getting a write pointer for " <<
                     cValues << " values");
      switch ( mMRI->type ) {
      case MRI_UCHAR:
        pixels = (void*) ucharScalars->WritePointer( cRead, zX * zElement );
        break;
      case MRI_INT:
        pixels = (void*) intScalars->WritePointer( cRead, zX * zElement );
        break;
      case MRI_LONG:
        pixels = (void*) longScalars->WritePointer( cRead, zX * zElement );
        break;
      case MRI_FLOAT:
        pixels = (void*) floatScalars->WritePointer( cRead, zX * zElement );
        break;
      case MRI_SHORT:
        pixels = (void*) shortScalars->WritePointer( cRead, zX * zElement );
        break;
      default:
        break ;
      }
      if ( NULL == pixels ) {
        vtkErrorMacro (<< "Couldn't get a write pointer");
        return;
      }
      memcpy( pixels, mMRI->slices[nZ][nY], zX* zElement );
      cRead += zX;
#endif
    }
  }
  free( tuple );
  
  // Assign the scalars array to the image.
  mImageData->GetPointData()->SetScalars( scalars );

  scalars->Delete();

  this->Modified();

  // Need to recalc our bounds at some point.
  mbBoundsCacheDirty = true;
}

void
vtkFSVolumeSource::GetRASBounds ( float oRASBounds[6] ) {

  if ( NULL == mMRI ) {

    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
                                      oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;

    return;
  }

  // If we need to rebuild the cache...
  if ( mbBoundsCacheDirty ) {

    mRASBounds[0] = mRASBounds[2] = mRASBounds[4] = 999999;
    mRASBounds[1] = mRASBounds[3] = mRASBounds[5] = -999999;

    // For each corner, convert to RAS, and get the bounds.
    float cornerFactor[3];
    float RAS[3];
    for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ ) {
      for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ ) {
        for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ ) {

          this->ConvertIndexToRAS( cornerFactor[0]*mMRI->width,
                                   cornerFactor[1]*mMRI->height,
                                   cornerFactor[2]*mMRI->depth,
                                   RAS[0], RAS[1], RAS[2] );

          if ( RAS[0] < mRASBounds[0] ) mRASBounds[0] = RAS[0];
          if ( RAS[0] > mRASBounds[1] ) mRASBounds[1] = RAS[0];
          if ( RAS[1] < mRASBounds[2] ) mRASBounds[2] = RAS[1];
          if ( RAS[1] > mRASBounds[3] ) mRASBounds[3] = RAS[1];
          if ( RAS[2] < mRASBounds[4] ) mRASBounds[4] = RAS[2];
          if ( RAS[2] > mRASBounds[5] ) mRASBounds[5] = RAS[2];
        }
      }
    }

    mbBoundsCacheDirty = false;
  }

  oRASBounds[0] = mRASBounds[0];
  oRASBounds[1] = mRASBounds[1];
  oRASBounds[2] = mRASBounds[2];
  oRASBounds[3] = mRASBounds[3];
  oRASBounds[4] = mRASBounds[4];
  oRASBounds[5] = mRASBounds[5];
}

void
vtkFSVolumeSource::CopyMatricesFromMRI () {

  if ( NULL == mMRI )
    return;

  MATRIX* m = extract_i_to_r( mMRI );
  for ( int i = 0; i < 16; i++ ) {
    mVoxelToRASMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &m );

  m = extract_r_to_i( mMRI );
  for ( int i = 0; i < 16; i++ ) {
    mRASToVoxelMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &m );
}

float
vtkFSVolumeSource::GetPreferredValueIncrement () {

  if ( !mMRI )
    return 1;

  // First look at the range. If it's an int or char type, increment
  // is 1.
  switch ( mMRI->type ) {

  case MRI_UCHAR:
  case MRI_INT:
  case MRI_LONG:
  case MRI_SHORT:
    return 1;

  default:

    // Look at the range. If it's > 100, inc is 1, 10-100, inc is .1,
    // 1-10, inc is .01, etc.
    float range = this->GetMaxValue() - this->GetMinValue();
    float inc = 1;
    if ( range >= 1000000 )        {
      inc = 1000;
    } else if ( range >=  100000 )        {
      inc =  100;
    } else if ( range >=   10000 )        {
      inc =   10;
    } else if ( range >=    1000 )        {
      inc =    1;
    } else if ( range >=      10 )        {
      inc =    0.1;
    } else if ( range >=       1 )        {
      inc =    0.01;
    } else if ( range >=       0.1 )      {
      inc =    0.001;
    } else if ( range >=       0.01 )     {
      inc =    0.0001;
    } else if ( range >=       0.001 )    {
      inc =    0.00001;
    } else if ( range >=       0.0001 )   {
      inc =    0.000001;
    } else if ( range >=       0.00001 )  {
      inc =    0.0000001;
    } else if ( range >=       0.000001 ) {
      inc =    0.00000001;
    } else                          {
      inc =    0.000000001;
    }

    return inc;
  }
}
