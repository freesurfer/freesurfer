/**
 * @file  FSVolume.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:27:25 $
 *    $Revision: 1.17 $
 *
 * Copyright (C) 2008-2009,
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

#include <wx/wx.h>
#include <wx/ffile.h>
#include "FSVolume.h"
#include "MyUtils.h"
#include <stdexcept>
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkImageReslice.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageChangeInformation.h"
#include "vtkMath.h"

extern "C"
{
#include "registerio.h"
#include "transform.h"
}

#define DEFAULT_RESAMPLE_ALGORITHM SAMPLE_NEAREST

using namespace std;

FSVolume::FSVolume( FSVolume* ref ) :
    m_MRI( NULL ),
    m_MRITarget( NULL ),
    m_MRIRef( NULL ),
    m_MRIOrigTarget( NULL ),
    m_matReg( NULL ),
    m_volumeRef( ref ),
    m_fMinValue( 0 ),
    m_fMaxValue( 1 ),
    m_bResampleToRAS( true ),
    m_bBoundsCacheDirty( true )
{
  m_imageData = NULL;
  if ( ref )
  {
    SetMRI( m_MRIRef, ref->m_MRI );
    SetMRI( m_MRIOrigTarget, ref->m_MRIOrigTarget );

    m_bResampleToRAS = ref->m_bResampleToRAS;

    if ( !m_bResampleToRAS )
      SetMRI( m_MRITarget, ref->m_MRITarget );
  }
}

FSVolume::~FSVolume()
{
  if ( m_MRI )
    ::MRIfree( &m_MRI );

  if ( m_MRITarget )
    ::MRIfree( &m_MRITarget );

  if ( m_MRIRef )
    ::MRIfree( &m_MRIRef );

  if ( m_MRIOrigTarget )
    ::MRIfree( &m_MRIOrigTarget );

  if ( m_matReg )
    ::MatrixFree( &m_matReg );
}

bool FSVolume::MRIRead( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_MRI )
    ::MRIfree( &m_MRI );

  event.SetInt( event.GetInt() + 1 );
  wxPostEvent( wnd, event );

  char* fn = strdup( filename );
  m_MRI = ::MRIread( fn );
  free( fn );

  if ( m_MRI == NULL )
  {
    cerr << "MRIread failed" << endl;
    return false;
  }

  // read registration matrix
  if ( reg_filename && !LoadRegistrationMatrix( reg_filename ) )
  {
    cerr << "Read registration failed" << endl;
    return false;
  }

  MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );

  this->CopyMatricesFromMRI();
  this->MapMRIToImage( wnd, event );

  return true;
}

/*
bool FSVolume::LoadRegistrationMatrix( const char* filename )
{
 if ( m_matReg )
  ::MatrixFree( &m_matReg );
 m_matReg = NULL;

 if ( MyUtils::HasExtension( filename, "xfm" ) )  // MNI style
 {
  if ( regio_read_mincxfm( (char*)filename, &m_matReg, NULL ) != 0 )
   return false;

  m_nRegType = REG_MNI;
 }
 else if ( MyUtils::HasExtension( filename, "lta" ) ) // LTA style
 {
  TRANSFORM* FSXform = TransformRead( (char*)filename );
  if ( FSXform == NULL )
   return false;
  LTA* lta = (LTA*) FSXform->xform;
  if ( lta->type != LINEAR_RAS_TO_RAS )
  {
   cout << "INFO: LTA input is not RAS to RAS...converting..." << endl;
   lta = LTAchangeType( lta, LINEAR_RAS_TO_RAS );
  }
  if ( lta->type != LINEAR_RAS_TO_RAS )
  {
   cerr << "ERROR: LTA input is not RAS to RAS" << endl;
   TransformFree( &FSXform );
   return false;
  }
     // Assume RAS2RAS and uses vox2ras from input volumes:
     // Note: This ignores the volume geometry in the LTA file.
  m_matReg = MatrixCopy( lta->xforms[0].m_L, NULL );
  TransformFree( &FSXform );

  m_nRegType = REG_LTA;
 }
 else  // tkregister style
 {
  if ( MyUtils::HasExtension( filename, "mat" ) )
  {
   wxFFile file( filename );
   wxString strg;
   if ( !file.ReadAll( &strg ) )
    return false;
   strg.Replace( "\n", " " );
   wxArrayString ar = MyUtils::SplitString( strg, " " );
   if ( ar.Count() < 16 )
    return false;

   MATRIX* m = MatrixAlloc( 4, 4, MATRIX_REAL );
   double val;
   for ( int i = 0; i < 16; i++ )
   {
    ar[i].ToDouble( &val );
    *MATRIX_RELT(m, (i/4)+1, (i%4)+1) = val;
   }
   m_matReg = MRIfsl2TkReg( m_MRIRef, m_MRI, m );
  // MatrixFree( &m );
  }
  else
  {
   char* subject;
   float inplaneres, betplaneres, intensity;
   int float2int;
   if ( regio_read_register( (char*)filename, &subject, &inplaneres, &betplaneres, &intensity, &m_matReg, &float2int ) != 0 )
    return false;

   free( subject );
  }
  m_nRegType = REG_TKREGISTER;
 }

 return true;
}
*/

// read in registration file and convert it to tkreg style
bool FSVolume::LoadRegistrationMatrix( const char* filename )
{
  if ( m_matReg )
    ::MatrixFree( &m_matReg );
  m_matReg = NULL;

  if ( MyUtils::HasExtension( filename, "xfm" ) )  // MNI style
  {
    MATRIX* m = NULL;
    if ( regio_read_mincxfm( (char*)filename, &m, NULL ) != 0 )
      return false;

    m_matReg = MRItkRegMtx( m_MRIRef, m_MRI, m );
    MatrixFree( &m );
  }
  else if ( MyUtils::HasExtension( filename, "mat" ) )  // fsl style
  {
    wxFFile file( filename );
    wxString strg;
    if ( !file.ReadAll( &strg ) )
      return false;
    strg.Replace( "\n", " " );
    wxArrayString ar = MyUtils::SplitString( strg, " " );
    if ( ar.Count() < 16 )
      return false;

    MATRIX* m = MatrixAlloc( 4, 4, MATRIX_REAL );
    double val;
    for ( int i = 0; i < 16; i++ )
    {
      ar[i].ToDouble( &val );
      *MATRIX_RELT(m, (i/4)+1, (i%4)+1) = val;
    }
    m_matReg = MRIfsl2TkReg( m_MRIRef, m_MRI, m );
    MatrixFree( &m );
  }
  else if ( MyUtils::HasExtension( filename, "dat" ) )  // tkregister style
  {
    char* subject = NULL;
    float inplaneres, betplaneres, intensity;
    int float2int;
    if ( regio_read_register( (char*)filename, &subject, &inplaneres, &betplaneres, &intensity, &m_matReg, &float2int ) != 0 )
      return false;

    free( subject );
  }
  else  // LTA style & all possible other styles
  {
    TRANSFORM* FSXform = TransformRead( (char*)filename );
    if ( FSXform == NULL )
      return false;
    LTA* lta = (LTA*) FSXform->xform;
    if ( lta->type != LINEAR_RAS_TO_RAS )
    {
      cout << "INFO: LTA input is not RAS to RAS...converting..." << endl;
      lta = LTAchangeType( lta, LINEAR_RAS_TO_RAS );
    }
    if ( lta->type != LINEAR_RAS_TO_RAS )
    {
      cerr << "ERROR: LTA input is not RAS to RAS" << endl;
      TransformFree( &FSXform );
      return false;
    }
    // Assume RAS2RAS and uses vox2ras from input volumes:
    // Note: This ignores the volume geometry in the LTA file.
    m_matReg = MRItkRegMtx( m_MRIRef, m_MRI, lta->xforms[0].m_L );
    TransformFree( &FSXform );
  }


  return true;
}

void FSVolume::Create( FSVolume* src_vol, bool bCopyVoxelData )
{
  if ( m_MRI )
    ::MRIfree( &m_MRI );
  if ( m_matReg )
    ::MatrixFree( &m_matReg );

  SetMRI( m_MRIRef, src_vol->GetMRI() );
  memcpy( m_RASToVoxelMatrix, src_vol->m_RASToVoxelMatrix, 16 * sizeof( double ) );
  memcpy( m_VoxelToRASMatrix, src_vol->m_VoxelToRASMatrix, 16 * sizeof( double ) );
  memcpy( m_MRIToImageMatrix, src_vol->m_MRIToImageMatrix, 16 * sizeof( double ) );
  memcpy( m_VoxelToVoxelMatrix, src_vol->m_VoxelToVoxelMatrix, 16 * sizeof( double ) );
  memcpy( m_RASToRASMatrix, src_vol->m_RASToRASMatrix, 16 * sizeof( double ) );
  memcpy( m_RASToTkRegMatrix, src_vol->m_RASToTkRegMatrix, 16 * sizeof( double ) );

  m_matReg = MatrixCopy( src_vol->m_matReg, NULL );

  m_fMinValue = src_vol->m_fMinValue;
  m_fMaxValue = src_vol->m_fMaxValue;
  m_bResampleToRAS = src_vol->m_bResampleToRAS;

// SetOriginalOrigin( src_vol->m_fOriginalOrigin );

  if ( bCopyVoxelData )
  {
    m_MRI = MRIcopy( src_vol->m_MRI, NULL );
  }
  else
  {
    m_MRI = MRIallocSequence( src_vol->m_MRI->width, src_vol->m_MRI->height, src_vol->m_MRI->depth,
                              src_vol->m_MRI->type, 1 );
  }

  if ( NULL == m_MRI )
  {
    cerr << "Couldn't allocate new mri." << endl;
    return;
  }

  SetMRITarget( src_vol->m_MRITarget );

  if ( src_vol->m_MRIOrigTarget )
  {
    MRI* mri = src_vol->m_MRIOrigTarget;
    m_MRIOrigTarget = MRIallocHeader( mri->width, mri->height, mri->depth, mri->type );
    MRIcopyHeader( mri, m_MRIOrigTarget );
  }

  // Copy the header from the template into the new mri.
  if ( !bCopyVoxelData )
    MRIcopyHeader( src_vol->m_MRI, m_MRI );

// if ( !m_imageData.GetPointer() )
  if ( m_imageData == NULL )
  {
    m_imageData = vtkSmartPointer<vtkImageData>::New();
  }

  if ( bCopyVoxelData )
  {
    m_imageData->DeepCopy( src_vol->m_imageData );
  }
  else
  {
    m_imageData->SetNumberOfScalarComponents( 1 );
    switch ( m_MRI->type )
    {
    case MRI_UCHAR:
      m_imageData->SetScalarTypeToUnsignedChar();
      break;
    case MRI_INT:
      m_imageData->SetScalarTypeToInt();
      break;
    case MRI_LONG:
      m_imageData->SetScalarTypeToLong();
      break;
    case MRI_FLOAT:
      m_imageData->SetScalarTypeToFloat();
      break;
    case MRI_SHORT:
      m_imageData->SetScalarTypeToShort();
      ;
      break;
    default:
      break;
    }

    m_imageData->SetOrigin( src_vol->m_imageData->GetOrigin() );
    m_imageData->SetSpacing( src_vol->m_imageData->GetSpacing() );
    m_imageData->SetDimensions( src_vol->m_imageData->GetDimensions() );
    m_imageData->AllocateScalars();
    char* ptr = ( char* )m_imageData->GetScalarPointer();
    int* nDim = m_imageData->GetDimensions();
    // cout << nDim[0] << ", " << nDim[1] << ", " << nDim[2] << endl;
    memset( ptr, 0, m_imageData->GetScalarSize() * nDim[0] * nDim[1] * nDim[2] );
  }
}

void FSVolume::SetMRI( MRI*& mri_out, MRI* mri_in  )
{
  if ( mri_out )
    ::MRIfree( &mri_out );

  if ( mri_in )
  {
    mri_out = MRIallocHeader( mri_in->width, mri_in->height, mri_in->depth, mri_in->type );
    MRIcopyHeader( mri_in, mri_out );
  }
  else
  {
    mri_out = NULL;
  }
}

void FSVolume::SetMRITarget( MRI* mri )
{
  if ( mri )
  {
    if ( m_MRITarget )
      ::MRIfree( &m_MRITarget );

    m_MRITarget = MRIallocHeader( mri->width, mri->height, mri->depth, mri->type );
    MRIcopyHeader( mri, m_MRITarget );
  }
}

bool FSVolume::MRIWrite( const char* filename )
{
  char* fn = strdup( filename );
  int err = ::MRIwrite( m_MRI, fn );
  free( fn );

  if ( err != 0 )
  {
    cerr << "MRIwrite failed" << endl;
  }

  return err == 0;
}

bool FSVolume::MRIWrite()
{
  int err = ::MRIwrite( m_MRI, m_MRI->fname );
  if ( err != 0 )
  {
    cerr << "MRIwrite failed" << endl;
  }
  return err == 0;
}

/*
void FSVolume::UpdateMRIFromImage( vtkImageData* rasImage, wxWindow* wnd, wxCommandEvent& event )
{
 cout << "UpdateMRIFromImage begins" << endl;
 int nProgressStep = ( 99 - event.GetInt() ) / 5;
 vtkMatrix4x4* m = vtkMatrix4x4::New();
 m->DeepCopy( m_VoxelToVoxelMatrix );
 double p[4] = { 0, 0, 0, 1 };
 int* dim = rasImage->GetDimensions();
 int n[3];
 for ( int j = 0; j < m_MRI->height; j++ )
 {
  for ( int k = 0; k < m_MRI->depth; k++ )
  {
   for ( int i = 0; i < m_MRI->width; i++ )
   {
    p[0] = i;
    p[1] = j;
    p[2] = k;
    m->MultiplyPoint( p, p );
    for ( int nFrame = 0; nFrame < m_MRI->nframes; nFrame++ )
    {
     void* buf = (char*)&MRIseq_vox( m_MRI, 0, j, k, nFrame) + i * rasImage->GetScalarSize();
     n[0] = (int)(p[0]+0.5);
     n[1] = (int)(p[1]+0.5);
     n[2] = (int)(p[2]+0.5);
     if ( n[0] >= 0 && n[0] < dim[0] &&
       n[1] >= 0 && n[1] < dim[1] &&
       n[2] >= 0 && n[2] < dim[2] )
     {
      void* ptr = (char*)rasImage->GetScalarPointer( n[0], n[1], n[2] )
        + nFrame * rasImage->GetScalarSize();
      memcpy( buf, ptr, rasImage->GetScalarSize() );
     }
    }
   }
  }
  if ( m_MRI->height >= 5 && j%(m_MRI->height/5) == 0 )
  {
   event.SetInt( event.GetInt() + nProgressStep );
   wxPostEvent( wnd, event );
  }
 }
 m->Delete();
 cout << "UpdateMRIFromImage finished" << endl;

 MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );
}

*/

void FSVolume::UpdateMRIFromImage( vtkImageData* rasImage, wxWindow* wnd, wxCommandEvent& event )
{
  int nProgressStep = ( 30 - event.GetInt() ) / 5;

  MATRIX* vox2vox = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1) = m_VoxelToVoxelMatrix[i];
    ;
  }

  MRI* mri = MRIallocSequence( m_MRITarget->width, m_MRITarget->height, m_MRITarget->depth, m_MRITarget->type, m_MRI->nframes );
  MRIcopyHeader( m_MRITarget, mri );

// cout << "begin copy pixels" << endl;
  if ( mri->nframes > 1 )
  {
    for ( int j = 0; j < mri->height; j++ )
    {
      for ( int k = 0; k < mri->depth; k++ )
      {
        for ( int i = 0; i < mri->width; i++ )
        {
          for ( int nFrame = 0; nFrame < mri->nframes; nFrame++ )
          {
            void* buf = (char*)&MRIseq_vox( mri, i, j, k, nFrame);
            void* ptr = (char*)rasImage->GetScalarPointer( i, j, k )
                        + nFrame * rasImage->GetScalarSize();
            memcpy( buf, ptr, rasImage->GetScalarSize() );
          }
        }
      }
      if ( m_MRI->height >= 5 && j%(m_MRI->height/5) == 0 )
      {
        event.SetInt( event.GetInt() + nProgressStep );
        wxPostEvent( wnd, event );
      }
    }
  }
  else
  {
    for ( int k = 0; k < mri->depth; k++ )
    {
      void* ptr = rasImage->GetScalarPointer( 0, 0, k );
      BUFTYPE* buf = &MRIseq_vox( mri, 0, 0, k, 0);
      memcpy( buf, ptr, mri->bytes_per_slice );

      if ( mri->depth >= 5 && k%(mri->depth/5) == 0 )
      {
        event.SetInt( event.GetInt() + nProgressStep );
        wxPostEvent( wnd, event );
      }
    }
  }

// cout << "begin vol2vol" << endl;
  MRIvol2Vol( mri, m_MRI, vox2vox, DEFAULT_RESAMPLE_ALGORITHM, 0 );
// cout << "end vol2vol" << endl;

  MRIfree( &mri );
  MatrixFree( &vox2vox );
  MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );
}


int FSVolume::OriginalIndexToRAS( float iIdxX, float iIdxY, float iIdxZ,
                                  float& oRASX, float& oRASY, float& oRASZ )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present." << endl;
    return 1;
  }

  Real ix, iy, iz, wx, wy, wz;
  int r;

  ix = iIdxX;
  iy = iIdxY;
  iz = iIdxZ;
  r = ::MRIvoxelToWorld( m_MRI, ix, iy, iz, &wx, &wy, &wz );
  oRASX = wx;
  oRASY = wy;
  oRASZ = wz;

  return r;
}

int FSVolume::RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                                   float& oIdxX, float& oIdxY, float& oIdxZ )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present." << endl;
    return 1;
  }

  Real ix, iy, iz, wx, wy, wz;
  int r;

  wx = iRASX;
  wy = iRASY;
  wz = iRASZ;
  r = ::MRIworldToVoxel( m_MRI, wx, wy, wz, &ix, &iy, &iz );
  oIdxX = ix;
  oIdxY = iy;
  oIdxZ = iz;

  return r;
}

int FSVolume::RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                                   int& oIdxX, int& oIdxY, int& oIdxZ )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present." << endl;
    return 1;
  }

  Real ix, iy, iz;
  int r;

  r = ::MRIworldToVoxel( m_MRI, iRASX, iRASY, iRASZ, &ix, &iy, &iz );

  oIdxX = (int)( ix + 0.5 );
  oIdxY = (int)( iy + 0.5 );
  oIdxZ = (int)( iz + 0.5 );

  return r;
}


void FSVolume::MapMRIToImage( wxWindow* wnd, wxCommandEvent& event )
{
  event.SetInt( event.GetInt() + 1 );
  wxPostEvent( wnd, event );

  // first create target MRI
  float bounds[6];
  double voxelSize[3];
  this->GetPixelSize( voxelSize );
  int dim[3];

  MRI* rasMRI = NULL;
  MATRIX* m = MatrixZero( 4, 4, NULL );
  if ( m_matReg && m_MRIRef && m_volumeRef )
  {
    // if there is registration matrix, set target as the reference's target
    MRI* mri = m_volumeRef->m_MRITarget;
    rasMRI = MRIallocSequence( mri->width, mri->height, mri->depth, m_MRI->type, m_MRI->nframes );
    MRIcopyHeader( mri, rasMRI );
  }
  else if ( m_bResampleToRAS )
  {
    this->GetBounds( bounds );
    for ( int i = 0; i < 3; i++ )
      dim[i] = (int) ( ( bounds[i*2+1] - bounds[i*2] ) / voxelSize[i] + 0.5 );

    rasMRI = MRIallocSequence( dim[0], dim[1], dim[2], m_MRI->type, m_MRI->nframes );
    if ( rasMRI == NULL )
    {
      cerr << "Can not allocate memory for rasMRI" << endl;
      return;
    }
    MRIsetResolution( rasMRI, voxelSize[0], voxelSize[1], voxelSize[2] );

    *MATRIX_RELT( m, 1, 1 ) = voxelSize[0];
    *MATRIX_RELT( m, 2, 2 ) = voxelSize[1];
    *MATRIX_RELT( m, 3, 3 ) = voxelSize[2];
    *MATRIX_RELT( m, 4, 4 ) = 1;
    *MATRIX_RELT( m, 1, 4 ) = bounds[0];
    *MATRIX_RELT( m, 2, 4 ) = bounds[2];
    *MATRIX_RELT( m, 3, 4 ) = bounds[4];
    MRIsetVoxelToRasXform( rasMRI, m);
  }
  else
  {
    if ( !m_MRITarget )
    {
      double* rtv = this->GetVoxelToRASMatrix();
      int odim[3] = { m_MRI->width, m_MRI->height, m_MRI->depth };
      if ( fabs( rtv[0] ) > fabs( rtv[1] ) && fabs( rtv[0] ) > fabs( rtv[2] ) )
      {
        *MATRIX_RELT( m, 1, 1 ) = ( rtv[0] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( rtv[0] > 0 ? 0 : odim[0] - 1 );
        dim[0] = odim[0];
      }
      else if ( fabs( rtv[1] ) > fabs( rtv[2] ) )
      {
        *MATRIX_RELT( m, 1, 2 ) = ( rtv[1] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( rtv[1] > 0 ? 0 : odim[1] - 1 );
        dim[0] = odim[1];
      }
      else
      {
        *MATRIX_RELT( m, 1, 3 ) = (rtv[2] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( rtv[2] > 0 ? 0 : odim[2] - 1 );
        dim[0] = odim[2];
      }

      if ( fabs( rtv[4] ) > fabs( rtv[5] ) && fabs( rtv[4] ) > fabs( rtv[6] ) )
      {
        *MATRIX_RELT( m, 2, 1 ) = ( rtv[4] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( rtv[4] > 0 ? 0 : odim[0] - 1 );
        dim[1] = odim[0];
      }
      else if ( fabs( rtv[5] ) > fabs( rtv[6] ) )
      {
        *MATRIX_RELT( m, 2, 2 ) = ( rtv[5] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( rtv[5] > 0 ? 0 : odim[1] - 1 );
        dim[1] = odim[1];
      }
      else
      {
        *MATRIX_RELT( m, 2, 3 ) = ( rtv[6] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( rtv[6] > 0 ? 0 : odim[2] - 1 );
        dim[1] = odim[2];
      }

      if ( fabs( rtv[8] ) > fabs( rtv[9] ) && fabs( rtv[8] ) > fabs( rtv[10] ) )
      {
        *MATRIX_RELT( m, 3, 1 ) = ( rtv[8] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( rtv[8] > 0 ? 0 : odim[0] - 1 );
        dim[2] = odim[0];
      }
      else if ( fabs( rtv[9] ) > fabs( rtv[10] ) )
      {
        *MATRIX_RELT( m, 3, 2 ) = ( rtv[9] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( rtv[9] > 0 ? 0 : odim[1] - 1 );
        dim[2] = odim[1];
      }
      else

      {
        *MATRIX_RELT( m, 3, 3 ) = ( rtv[10] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( rtv[10] > 0 ? 0 : odim[2] - 1 );
        dim[2] = odim[2];
      }

      *MATRIX_RELT( m, 4, 4 ) = 1;

      rasMRI = MRIallocSequence( dim[0], dim[1], dim[2], m_MRI->type, m_MRI->nframes );
      if ( rasMRI == NULL )
      {
        cerr << "Can not allocate memory for rasMRI" << endl;
        return;
      }
      MRIsetResolution( rasMRI, voxelSize[0], voxelSize[1], voxelSize[2] );

      MATRIX* m1 = MRIgetVoxelToRasXform( m_MRI );
      MATRIX* m_inv = MatrixInverse( m, NULL );
      MATRIX* m2 = MatrixMultiply( m1, m_inv, NULL );

      MRIsetVoxelToRasXform( rasMRI, m2 );

      MatrixFree( &m1 );
      MatrixFree( &m2 );
      MatrixFree( &m_inv );
    }
    else
    {
      rasMRI = MRIallocSequence( m_MRITarget->width,
                                 m_MRITarget->height,
                                 m_MRITarget->depth,
                                 m_MRI->type,
                                 m_MRI->nframes );  // must use source data type and frames!
      if ( rasMRI == NULL )
      {
        cerr << "Can not allocate memory for rasMRI" << endl;
        return;
      }
      MRIcopyHeader( m_MRITarget, rasMRI );
      rasMRI->type = m_MRI->type;
      rasMRI->nframes = m_MRI->nframes;
    }
  }
  MatrixFree( &m );

  if ( m_matReg && !m_MRIRef )
  {
    cerr << "No target volume available! Can not use registration matrix." << endl;
  }

  if ( m_matReg && m_MRIRef )
  {
    // registration matrix is always converted to tkReg style now
    // if ( m_nRegType == REG_TKREGISTER )
    {
      MATRIX *vox2vox, *Tmov, *invTmov, *Ttarg;
      Tmov = MRIxfmCRS2XYZtkreg( m_MRI );
      invTmov = MatrixInverse( Tmov, NULL );
      Ttarg = MRIxfmCRS2XYZtkreg( m_MRIRef );
      // vox2vox = invTmov*R*Ttarg
      vox2vox = MatrixMultiply( invTmov, m_matReg, NULL );
      MatrixMultiply( vox2vox, Ttarg, vox2vox );

      MATRIX* t2r = MRIgetVoxelToVoxelXform( rasMRI, m_MRIRef );
      MatrixMultiply( vox2vox, t2r, t2r );

      MRIvol2Vol( m_MRI, rasMRI, t2r, DEFAULT_RESAMPLE_ALGORITHM, 0 );

      // copy vox2vox
      MatrixInverse( t2r, vox2vox );
      for ( int i = 0; i < 16; i++ )
      {
        m_VoxelToVoxelMatrix[i] = (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
      }

      MatrixFree( &vox2vox );
      MatrixFree( &Tmov );
      MatrixFree( &invTmov );
      MatrixFree( &Ttarg );
      MatrixFree( &t2r );
    }
    // Reg matrix is always converted to tkReg style at loading time,
    // but keep the following section of code for reference
    /*
    else if ( m_nRegType == REG_MNI || m_nRegType == REG_LTA )
    {
     MATRIX *vox2vox, *Tmov, *invTmov, *invReg, *Ttarg;
     Tmov = MRIgetVoxelToRasXform( m_MRI );
     Ttarg = MRIgetVoxelToRasXform( m_MRIRef );
     invTmov = MatrixInverse( Tmov, NULL );
     invReg = MatrixInverse( m_matReg, NULL );

     // vox2vox = invTmov*invR*Ttarg
     vox2vox = MatrixMultiply( invReg, Ttarg, NULL );
     MatrixMultiply( invTmov, vox2vox, vox2vox );

     MATRIX* t2r = MRIgetVoxelToVoxelXform( rasMRI, m_MRIRef );
     MatrixMultiply( vox2vox, t2r, t2r );

     MRIvol2Vol( m_MRI, rasMRI, t2r, DEFAULT_RESAMPLE_ALGORITHM, 0 );

     // copy vox2vox
     MatrixInverse( t2r, vox2vox );
     for ( int i = 0; i < 16; i++ )
     {
      m_VoxelToVoxelMatrix[i] = (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
     }

     MatrixFree( &vox2vox );
     MatrixFree( &Tmov );
     MatrixFree( &invTmov );
     MatrixFree( &invReg );
     MatrixFree( &Ttarg );
     MatrixFree( &t2r );
    }
    */
  }
  else
  {
    MRIvol2Vol( m_MRI, rasMRI, NULL, DEFAULT_RESAMPLE_ALGORITHM, 0 );
    MATRIX* vox2vox = MRIgetVoxelToVoxelXform( m_MRI, rasMRI );
    for ( int i = 0; i < 16; i++ )
    {
      m_VoxelToVoxelMatrix[i] = (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
    }
    MatrixFree( &vox2vox );
  }

// cout << "MRIvol2Vol finished" << endl;

  SetMRITarget( rasMRI );
  UpdateRASToRASMatrix();

  CreateImage( rasMRI, wnd, event );

  // copy mri pixel data to vtkImage we will use for display
  CopyMRIDataToImage( rasMRI, m_imageData, wnd, event );

  // Need to recalc our bounds at some point.
  m_bBoundsCacheDirty = true;

  ::MRIfree( &rasMRI );
}

void FSVolume::UpdateRASToRASMatrix()
{
  // compute ras2ras matrix
  MATRIX* r2r = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT( r2r, (i/4)+1, (i%4)+1 ) = m_VoxelToVoxelMatrix[i];
  }
  MATRIX* mMov = MRIgetVoxelToRasXform( m_MRI );
  MATRIX* mTarg = MRIgetVoxelToRasXform( m_MRITarget );
  MATRIX* mMov_inv = MatrixInverse( mMov, NULL );
  MatrixMultiply( r2r, mMov_inv, r2r );
  MatrixMultiply( mTarg, r2r, r2r );

  for ( int i = 0; i < 16; i++ )
  {
    m_RASToRASMatrix[i] = (double) *MATRIX_RELT((r2r),(i/4)+1,(i%4)+1);
  }

  MatrixFree( &r2r );
  MatrixFree( &mMov );
  MatrixFree( &mMov_inv );
  MatrixFree( &mTarg );
}

void FSVolume::CreateImage( MRI* rasMRI, wxWindow* wnd, wxCommandEvent& event )
{
  // first copy mri data to image
  vtkDataArray *scalars = NULL;
  vtkUnsignedCharArray  *ucharScalars = NULL;
  vtkIntArray           *intScalars = NULL;
  vtkShortArray         *shortScalars = NULL;
  vtkLongArray          *longScalars = NULL;
  vtkFloatArray         *floatScalars = NULL;
  int cValues;
  int zElement=0;

  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present." << endl;
    return;
  }

  int zX = rasMRI->width;
  int zY = rasMRI->height;
  int zZ = rasMRI->depth;
  int zFrames = rasMRI->nframes;

  m_imageData = vtkSmartPointer<vtkImageData>::New();
  vtkImageData* imageData = m_imageData;

  // This object's output space is in voxel coordinates.
  imageData->SetDimensions( zX, zY, zZ );

  double origin[3] = { 0, 0, 0 };
  if ( m_bResampleToRAS )
  {
    float bounds[6];
    GetBounds( bounds );
    origin[0] = bounds[0];
    origin[1] = bounds[2];
    origin[2] = bounds[4];
  }
  else if ( m_volumeRef )
  {
    m_volumeRef->m_imageData->GetOrigin( origin );
  }

// cout << "origin: " << origin[0] << "  " << origin[1] << "  " << origin[2] << endl
//   << "voxelsize: " << rasMRI->xsize << "  " << rasMRI->ysize << "  " << rasMRI->zsize << endl;
  // if map to RAS
  {
    imageData->SetSpacing( rasMRI->xsize, rasMRI->ysize, rasMRI->zsize );
    imageData->SetOrigin( origin[0], origin[1], origin[2] );
  }

  imageData->SetWholeExtent( 0, zX-1, 0, zY-1, 0, zZ-1 );
  imageData->SetNumberOfScalarComponents( zFrames );

  // create the scalars for all of the images. set the element size
  // for the data we will read.
  switch ( rasMRI->type )
  {
  case MRI_UCHAR:
    imageData->SetScalarTypeToUnsignedChar();
    ucharScalars = vtkUnsignedCharArray::New();
    scalars = (vtkDataArray*) ucharScalars;
    zElement = sizeof( unsigned char );
    break;
  case MRI_INT:
    imageData->SetScalarTypeToInt();
    intScalars = vtkIntArray::New();
    scalars = (vtkDataArray*) intScalars;
    zElement = sizeof( int );
    break;
  case MRI_LONG:
    imageData->SetScalarTypeToLong();
    longScalars = vtkLongArray::New();
    scalars = (vtkDataArray*) longScalars;
    zElement = sizeof( long );
    break;
  case MRI_FLOAT:
    imageData->SetScalarTypeToFloat();
    floatScalars = vtkFloatArray::New();
    scalars = (vtkDataArray*) floatScalars;
    zElement = sizeof( float );
    break;
  case MRI_SHORT:
    imageData->SetScalarTypeToShort();
    shortScalars = vtkShortArray::New();
    scalars = (vtkDataArray*) shortScalars;
    zElement = sizeof( short );
    break;
  default:
    break ;
  }

  if ( NULL == scalars )
  {
    cerr << "Couldn't allocate scalars array." << endl;
    return;
  }

  // change the number of components to store tuples
  if ( zFrames > 1 )
  {
    scalars->SetNumberOfComponents( zFrames );
  }

  cValues = zX * zY * zZ;
  scalars->Allocate( cValues );
  scalars->SetNumberOfTuples( zX*zY*zZ );

  // Assign the scalars array to the image.
  imageData->GetPointData()->SetScalars( scalars );
  scalars->Delete();

  if ( wnd )
  {
    event.SetInt( event.GetInt() + ( 100 - event.GetInt() ) / 10 );
    wxPostEvent( wnd, event );
  }
}


bool FSVolume::Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  if ( rotations.size() == 0 )
    return false;

  MRI* rasMRI;
  if ( rotations[0].Plane == -1 )   // restore
  {
    if ( !m_MRIOrigTarget )  // try to restore but no where to restore
      return false;

    rasMRI = MRIallocSequence( m_MRIOrigTarget->width,
                               m_MRIOrigTarget->height,
                               m_MRIOrigTarget->depth,
                               m_MRI->type,
                               m_MRI->nframes );
    MRIcopyHeader( m_MRIOrigTarget, rasMRI );
  }
  else    // rotate
  {
    MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );

    // calculate rotation matrix
    MATRIX* m_r = MatrixIdentity( 4, NULL ) ;
    for ( size_t i = 0; i < rotations.size(); i++ )
    {
      /* // convert RAS to voxel index
       Real vx, vy, vz;
       ::MRIworldToVoxel( m_MRITarget,
           rotations[i].Point[0],
           rotations[i].Point[1],
           rotations[i].Point[2],
           &vx, &vy, &vz );
       double pos[3] = { vx, vy, vz };
       MATRIX* m_tmp = GetRotationMatrix( rotations[i].Plane, rotations[i].Angle, pos );
       MATRIX* m_tmp1 = MatrixMultiply( m_tmp, m_r, NULL );
      */

      MATRIX* m_tmp = GetRotationMatrix( rotations[i].Plane, rotations[i].Angle, rotations[i].Point );
      MATRIX* m_tmp1 = MatrixMultiply( m_tmp, m_r, NULL );

      MatrixCopy( m_tmp1, m_r );
      MatrixFree( &m_tmp );
      MatrixFree( &m_tmp1 );

    }
    MATRIX* m_new = MatrixMultiply( m_r, m, NULL );

    // MRIsetVoxelToRasXform( m_MRITarget, m_new );

    rasMRI = MRIallocSequence( m_MRITarget->width,
                               m_MRITarget->height,
                               m_MRITarget->depth,
                               m_MRI->type,
                               m_MRI->nframes );
    // MRIcopyHeader( m_MRITarget, rasMRI );
    MRIsetResolution( rasMRI, m_MRITarget->xsize, m_MRITarget->ysize, m_MRITarget->zsize );
    MRIsetVoxelToRasXform( rasMRI, m_new );

    MatrixFree( &m );
    MatrixFree( &m_r );
    MatrixFree( &m_new );

    if ( !m_MRIOrigTarget )  // never rotated before, so save it as original
    {
      m_MRIOrigTarget = MRIallocHeader( m_MRITarget->width, m_MRITarget->height, m_MRITarget->depth, m_MRITarget->type );
      MRIcopyHeader( m_MRITarget, m_MRIOrigTarget );
    }
  }

  MATRIX* v2v = MRIgetVoxelToVoxelXform( m_MRITarget, rasMRI ); // old target to new target
  MATRIX* vm = MatrixAlloc( 4, 4, v2v->type );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT((vm),(i/4)+1,(i%4)+1) = m_VoxelToVoxelMatrix[i];
  }
  MATRIX* vox2vox = MatrixMultiply( v2v, vm, NULL );
  MatrixInverse( vox2vox, v2v );

  MRIvol2Vol( m_MRI, rasMRI, v2v, DEFAULT_RESAMPLE_ALGORITHM, 0 );

  // copy vox2vox
  for ( int i = 0; i < 16; i++ )
  {
    m_VoxelToVoxelMatrix[i] = (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &vox2vox );
  MatrixFree( &vm );
  MatrixFree( &v2v );

  SetMRITarget( rasMRI );
  UpdateRASToRASMatrix();

  // copy mri pixel data to vtkImage we will use for display
  CopyMRIDataToImage( rasMRI, m_imageData, wnd, event );

  // Need to recalc our bounds at some point.
  m_bBoundsCacheDirty = true;

  ::MRIfree( &rasMRI );

  m_bResampleToRAS = false;

  return true;
}


MATRIX* FSVolume::GetRotationMatrix( int nPlane, double angle, double* origin )
{
  // calculate rotation matrix
  MATRIX* m_r;
  switch ( nPlane )
  {
  case 0:
    m_r = MatrixAllocRotation( 4, angle/180 * vtkMath::Pi(), X_ROTATION );
    break;
  case 1:
    m_r = MatrixAllocRotation( 4, angle/180 * vtkMath::Pi(), Y_ROTATION );
    break;
  default:
    m_r = MatrixAllocRotation( 4, angle/180 * vtkMath::Pi(), Z_ROTATION );
    break;
  }

  MATRIX* m_t = MatrixZero( 4, 4, NULL );
  *MATRIX_RELT( m_t, 1, 1 ) = 1;
  *MATRIX_RELT( m_t, 2, 2 ) = 1;
  *MATRIX_RELT( m_t, 3, 3 ) = 1;
  *MATRIX_RELT( m_t, 4, 4 ) = 1;
  *MATRIX_RELT( m_t, 1, 4 ) = origin[0];
  *MATRIX_RELT( m_t, 2, 4 ) = origin[1];
  *MATRIX_RELT( m_t, 3, 4 ) = origin[2];

  MATRIX* m_t_inv = MatrixInverse( m_t, NULL );

  MATRIX* m1 = MatrixMultiply( m_r, m_t_inv, NULL );
  MATRIX* m = MatrixMultiply( m_t, m1, NULL );

  MatrixFree( &m1 );
  MatrixFree( &m_r );
  MatrixFree( &m_t );
  MatrixFree( &m_t_inv );

  return m;
}


void FSVolume::CopyMRIDataToImage( MRI* mri, vtkImageData* image, wxWindow* wnd, wxCommandEvent& event )
{
  // Copy the slice data into the scalars.
  int zX = mri->width;
  int zY = mri->height;
  int zZ = mri->depth;
  int zFrames = mri->nframes;

  int nTuple = 0;
  vtkDataArray *scalars = image->GetPointData()->GetScalars();
  int nProgressStep = ( 99-event.GetInt() ) / 5;
  for ( int nZ = 0; nZ < zZ; nZ++ )
  {
    for ( int nY = 0; nY < zY; nY++ )
    {
      for ( int nX = 0; nX < zX; nX++ )
      {
        for ( int nFrame = 0; nFrame < zFrames; nFrame++ )
        {
          switch ( mri->type )
          {
          case MRI_UCHAR:
            scalars->SetComponent( nTuple, nFrame,
                                   MRIseq_vox( mri, nX, nY, nZ, nFrame ) );
            break;
          case MRI_INT:
            scalars->SetComponent( nTuple, nFrame,
                                   MRIIseq_vox( mri, nX, nY, nZ, nFrame ) );
            break;
          case MRI_LONG:
            scalars->SetComponent( nTuple, nFrame,
                                   MRILseq_vox( mri, nX, nY, nZ, nFrame ) );
            break;
          case MRI_FLOAT:
            scalars->SetComponent( nTuple, nFrame,
                                   MRIFseq_vox( mri, nX, nY, nZ, nFrame ) );
            break;
          case MRI_SHORT:
            scalars->SetComponent( nTuple, nFrame,
                                   MRISseq_vox( mri, nX, nY, nZ, nFrame ) );
            break;
          default:
            break;
          }
        }
        nTuple++;
      }
    }

    if ( nZ%(max(1, zZ/5)) == 0 )
    {
      event.SetInt( event.GetInt() + nProgressStep );
      wxPostEvent( wnd, event );
    }
  }
}

/*

void FSVolume::CopyMRIDataToImage( MRI* mri, vtkImageData* image, wxWindow* wnd, wxCommandEvent& event )
{
  // Copy the slice data into the scalars.
// vtkDebugMacro (<< "Copying " << cValues << " values into the scalars array");
// float* tuple = (float*) malloc( sizeof(float) * zFrames );
 int zX = mri->width;
 int zY = mri->height;
 int zZ = mri->depth;
 int zFrames = mri->nframes;

 void* p = image->GetScalarPointer();
 memcpy( p, mri->chunk, image->GetScalarSize() * zX * zY * zZ * zFrames );
 event.SetInt( event.GetInt() + 10 );
 wxPostEvent( wnd, event );
}
*/
vtkImageData* FSVolume::GetImageOutput()
{
  return m_imageData;
}

void FSVolume::CopyMatricesFromMRI ()
{
  if ( NULL == m_MRI )
    return;

  MATRIX* m = MRIgetVoxelToRasXform( m_MRI );
  for ( int i = 0; i < 16; i++ )
  {
    m_VoxelToRASMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &m );

  m = MRIgetRasToVoxelXform( m_MRI );
  for ( int i = 0; i < 16; i++ )
  {
    m_RASToVoxelMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
  }

  MATRIX* tkreg = MRIxfmCRS2XYZtkreg( m_MRI );
  MATRIX* m1 = MatrixMultiply( tkreg, m, NULL );
  for ( int i = 0; i < 16; i++ )
  {
    m_RASToTkRegMatrix[i] = (double) *MATRIX_RELT((m1),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &tkreg );
  MatrixFree( &m );
  MatrixFree( &m1 );
}

void FSVolume::GetBounds ( float oRASBounds[6] )
{
  if ( NULL == m_MRI )
  {

    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
                                      oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;

    return;
  }

  if ( !m_bResampleToRAS )
  {
//  if ( m_imageData.GetPointer() )
    if ( m_imageData != NULL )
    {
      double* origin = m_imageData->GetOrigin();
      int* dim = m_imageData->GetDimensions();
      double* vs = m_imageData->GetSpacing();

      for ( int i = 0; i < 3; i++ )
      {
        oRASBounds[i*2] = origin[i];
        oRASBounds[i*2+1] = origin[i] + dim[i] * vs[i];
      }
    }
    else if ( m_MRITarget )
    {
      MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );
      oRASBounds[0] = *MATRIX_RELT( m, 1, 4 );
      oRASBounds[2] = *MATRIX_RELT( m, 2, 4 );
      oRASBounds[4] = *MATRIX_RELT( m, 3, 4 );
      oRASBounds[1] = m_MRITarget->width * m_MRITarget->xsize + oRASBounds[0];
      oRASBounds[3] = m_MRITarget->height * m_MRITarget->ysize + oRASBounds[2];
      oRASBounds[5] = m_MRITarget->depth * m_MRITarget->zsize + oRASBounds[4];
      MatrixFree( &m );
    }
    else
    {
      cout << "Did not find bounds." << endl;
    }
    return;
  }
  else if ( m_MRITarget )
  {
    MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );
    oRASBounds[0] = *MATRIX_RELT( m, 1, 4 );
    oRASBounds[2] = *MATRIX_RELT( m, 2, 4 );
    oRASBounds[4] = *MATRIX_RELT( m, 3, 4 );
    oRASBounds[1] = m_MRITarget->width * m_MRITarget->xsize + oRASBounds[0];
    oRASBounds[3] = m_MRITarget->height * m_MRITarget->ysize + oRASBounds[2];
    oRASBounds[5] = m_MRITarget->depth * m_MRITarget->zsize + oRASBounds[4];
    MatrixFree( &m );
    return;
  }
  else if ( m_bBoundsCacheDirty )
  {
    m_RASBounds[0] = m_RASBounds[2] = m_RASBounds[4] = 999999;
    m_RASBounds[1] = m_RASBounds[3] = m_RASBounds[5] = -999999;

    // For each corner, convert to RAS, and get the bounds.
    float cornerFactor[3];
    float RAS[3];
    for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
    {
      for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
      {
        for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
        {

          this->OriginalIndexToRAS( cornerFactor[0]*m_MRI->width,
                                    cornerFactor[1]*m_MRI->height,
                                    cornerFactor[2]*m_MRI->depth,
                                    RAS[0], RAS[1], RAS[2] );

          if ( RAS[0] < m_RASBounds[0] ) m_RASBounds[0] = RAS[0];
          if ( RAS[0] > m_RASBounds[1] ) m_RASBounds[1] = RAS[0];
          if ( RAS[1] < m_RASBounds[2] ) m_RASBounds[2] = RAS[1];
          if ( RAS[1] > m_RASBounds[3] ) m_RASBounds[3] = RAS[1];
          if ( RAS[2] < m_RASBounds[4] ) m_RASBounds[4] = RAS[2];
          if ( RAS[2] > m_RASBounds[5] ) m_RASBounds[5] = RAS[2];
        }
      }
    }

    m_bBoundsCacheDirty = false;
  }

  oRASBounds[0] = m_RASBounds[0];
  oRASBounds[1] = m_RASBounds[1];
  oRASBounds[2] = m_RASBounds[2];
  oRASBounds[3] = m_RASBounds[3];
  oRASBounds[4] = m_RASBounds[4];
  oRASBounds[5] = m_RASBounds[5];
}

void FSVolume::GetBounds( MRI* mri, float oRASBounds[6] )
{
  oRASBounds[0] = oRASBounds[2] = oRASBounds[4] = 1e10;
  oRASBounds[1] = oRASBounds[3] = oRASBounds[5] = -1e10;

  // For each corner, convert to RAS, and get the bounds.
  float cornerFactor[3];
  Real RAS[3];
  for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
  {
    for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
    {
      for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
      {
        ::MRIvoxelToWorld( mri,
                           cornerFactor[0]*mri->width,
                           cornerFactor[1]*mri->height,
                           cornerFactor[1]*mri->height,
                           &RAS[0], &RAS[1], &RAS[2] );

        if ( RAS[0] < oRASBounds[0] ) oRASBounds[0] = RAS[0];
        if ( RAS[0] > oRASBounds[1] ) oRASBounds[1] = RAS[0];
        if ( RAS[1] < oRASBounds[2] ) oRASBounds[2] = RAS[1];
        if ( RAS[1] > oRASBounds[3] ) oRASBounds[3] = RAS[1];
        if ( RAS[2] < oRASBounds[4] ) oRASBounds[4] = RAS[2];
        if ( RAS[2] > oRASBounds[5] ) oRASBounds[5] = RAS[2];
      }
    }
  }
}

int FSVolume::GetNumberOfFrames()
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present." << endl;
    return 0;
  }

  return m_MRI->nframes;
}

void FSVolume::GetPixelSize( double* pixelSize )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present." << endl;
    return;
  }


  if ( !m_bResampleToRAS )
  {
    if ( m_MRITarget )
    {
      pixelSize[0] = m_MRITarget->xsize;
      pixelSize[1] = m_MRITarget->ysize;
      pixelSize[2] = m_MRITarget->zsize;

      return;
    }
  }
  else if ( m_matReg && m_MRIRef )
  {
    MRI* mri = m_volumeRef->m_MRITarget;
    if ( mri )
    {
      pixelSize[0] = mri->xsize;
      pixelSize[1] = mri->ysize;
      pixelSize[2] = mri->zsize;

      return;
    }
  }

  double* m = GetVoxelToRASMatrix();

  if ( fabs( m[0] ) > fabs( m[1] ) && fabs( m[0] ) > fabs( m[2] ) )
    pixelSize[0] = m_MRI->xsize;
  else if ( fabs( m[1] ) > fabs( m[2] ) )
    pixelSize[0] = m_MRI->ysize;
  else
    pixelSize[0] = m_MRI->zsize;

  if ( fabs( m[4] ) > fabs( m[5] ) && fabs( m[4] ) > fabs( m[6] ) )
    pixelSize[1] = m_MRI->xsize;
  else if ( fabs( m[5] ) > fabs( m[6] ) )
    pixelSize[1] = m_MRI->ysize;
  else
    pixelSize[1] = m_MRI->zsize;

  if ( fabs( m[8] ) > fabs( m[9] ) && fabs( m[8] ) > fabs( m[10] ) )
    pixelSize[2] = m_MRI->xsize;
  else if ( fabs( m[9] ) > fabs( m[10] ) )
    pixelSize[2] = m_MRI->ysize;
  else
    pixelSize[2] = m_MRI->zsize;
}

void FSVolume::TargetToRAS( double x_in, double y_in, double z_in,
                            double& x_out, double& y_out, double& z_out )
{
  double pos[3] = { x_in, y_in, z_in };
  TargetToRAS( pos, pos );
  x_out = pos[0];
  y_out = pos[1];
  z_out = pos[2];
}

void FSVolume::TargetToRAS( const double* pos_in, double* pos_out )
{
  // if already resampled to standard RAS, no need to remap
  /* if ( m_bResampleToRAS )
    memcpy( pos_out, pos_in, sizeof( double ) * 3 );
   else
  */
  {
    double pos[4] = { 0 };
    double* vs = m_imageData->GetSpacing();
    double* origin = m_imageData->GetOrigin();
    for ( int i = 0; i < 3; i++ )
      pos[i] = ( pos_in[i] - origin[i] ) / vs[i];

    Real fpos[3];
    ::MRIvoxelToWorld( m_MRITarget, (float)pos[0], (float)pos[1], (float)pos[2], &fpos[0], &fpos[1], &fpos[2] );
//  cout << "out: " << fpos[0] << " " << fpos[1] << " " << fpos[2] << endl;
    for ( int i = 0; i < 3; i++ )
    {
      pos_out[i] = fpos[i];
    }
  }
}

void FSVolume::TargetToRAS( const float* pos_in, float* pos_out )
{
  double pos[3] = { 0 };
  TargetToRAS( pos_in[0], pos_in[1], pos_in[2], pos[0], pos[1], pos[2] );
  pos_out[0] = pos[0];
  pos_out[1] = pos[1];
  pos_out[2] = pos[2];
}

void FSVolume::RASToTarget( const double* pos_in, double* pos_out )
{
  /* if ( m_bResampleToRAS )
   {
    memcpy( pos_out, pos_in, sizeof( double ) * 3 );
   }
   else
  */
  {
    Real pos[4] = { 0 };
    ::MRIworldToVoxel( m_MRITarget, (float)pos_in[0], (float)pos_in[1], (float)pos_in[2],
                       &pos[0], &pos[1], &pos[2] );
    double* vs = m_imageData->GetSpacing();
    double* origin = m_imageData->GetOrigin();
    for ( int i = 0; i < 3; i++ )
    {
      pos_out[i] = vs[i] * pos[i] + origin[i];
    }
  }
}

void FSVolume::RASToTarget( const float* pos_in, float* pos_out )
{
  double p_in[3] = { pos_in[0], pos_in[1], pos_in[2] };
  double p_out[3];
  RASToTarget( p_in, p_out );
  pos_out[0] = p_out[0];
  pos_out[1] = p_out[1];
  pos_out[2] = p_out[2];
}

void FSVolume::RASToTargetIndex( const double* pos_in, int* index_out )
{
  double pos[3];
  RASToTarget( pos_in, pos );
  double* vs = m_imageData->GetSpacing();
  double* orig = m_imageData->GetOrigin();
  for ( int i = 0; i < 3; i++ )
  {
    index_out[i] = ( int )( ( pos[i] - orig[i] ) / vs[i] + 0.5 );
  }
}

void FSVolume::RASToNativeRAS( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  double m[16];
  vtkMatrix4x4::Invert( m_RASToRASMatrix, m );
  vtkMatrix4x4::MultiplyPoint( m, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::NativeRASToRAS( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  vtkMatrix4x4::MultiplyPoint( m_RASToRASMatrix, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::NativeRASToTkReg( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  vtkMatrix4x4::MultiplyPoint( m_RASToTkRegMatrix, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::TkRegToNativeRAS( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  double m[16];
  vtkMatrix4x4::Invert( m_RASToTkRegMatrix, m );
  vtkMatrix4x4::MultiplyPoint( m, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}
