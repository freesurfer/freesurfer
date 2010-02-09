/**
 * @file  FSVolume.h
 * @brief Base volume class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/02/09 20:24:51 $
 *    $Revision: 1.43 $
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
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"

extern "C"
{
#include "registerio.h"
#include "transform.h"
}

using namespace std;

FSVolume::FSVolume( FSVolume* ref ) :
    m_MRI( NULL ),
    m_MRITarget( NULL ),
    m_MRIRef( NULL ),
    m_MRIOrigTarget( NULL ),
    m_MRITemp( NULL ),
    m_matReg( NULL ),
    m_ctabEmbedded( NULL ),
    m_volumeRef( ref ),
    m_fMinValue( 0 ),
    m_fMaxValue( 1 ),
    m_bResampleToRAS( false ),
    m_bBoundsCacheDirty( true ),
    m_nInterpolationMethod( SAMPLE_NEAREST )
{
  m_imageData = NULL;
  if ( ref )
  {
    SetMRI( m_MRIRef, ref->m_MRI );
//    SetMRI( m_MRIOrigTarget, ref->m_MRIOrigTarget );

    m_bResampleToRAS = ref->m_bResampleToRAS;

//    if ( !m_bResampleToRAS )
//      SetMRI( m_MRITarget, ref->m_MRITarget );
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

  if ( m_MRITemp )
    ::MRIfree( &m_MRITemp );

  if ( m_matReg )
    ::MatrixFree( &m_matReg );
  
  if ( m_ctabEmbedded )
    ::CTABfree( &m_ctabEmbedded );
}

bool FSVolume::LoadMRI( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event )
{
  event.SetInt( event.GetInt() + 1 );
  wxPostEvent( wnd, event );
  
  // save old header to release later so there is no gap where m_MRI becomes NULL during re-loading process
  MRI* tempMRI = m_MRI;

  char* fn = strdup( filename );
  m_MRI = ::MRIread( fn );      // could be long process
  free( fn );

  if ( m_MRI == NULL )
  {
    cerr << "MRIread failed" << endl;
    if ( tempMRI )
      m_MRI = tempMRI;
    return false;
  }
  
  // if m_MRI successfully loaded, release old header.
  if ( tempMRI )
    ::MRIfree( &tempMRI );

  if ( m_MRI->ct != NULL )
    m_ctabEmbedded = CTABdeepCopy( m_MRI->ct );
  
  /*
  if ( m_MRI->AutoAlign != NULL )
  {
    MATRIX* M = m_MRI->AutoAlign;
    cout << M->rptr[1][1] << " " << M->rptr[1][2] << " " << M->rptr[1][3] << " " << M->rptr[1][4] << endl
        << M->rptr[2][1] << " " << M->rptr[2][2] << " " << M->rptr[2][3] << " " << M->rptr[2][4] << endl
        << M->rptr[3][1] << " " << M->rptr[3][2] << " " << M->rptr[3][3] << " " << M->rptr[3][4] << endl
        << M->rptr[4][1] << " " << M->rptr[4][2] << " " << M->rptr[4][3] << " " << M->rptr[4][4] << endl;
  }
  */
  
  if ( reg_filename && !m_MRIRef )
  {
    cerr << "Error: A target volume must be loaded first to apply registration matrix." << endl;
    return false;
  }
  
  // read registration matrix
  if ( reg_filename && !LoadRegistrationMatrix( reg_filename ) )
  {
    cerr << "Read registration failed" << endl;
    return false;
  }

  MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );

  return true;
}

bool FSVolume::MRIRead( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event )
{
  if ( LoadMRI( filename, reg_filename, wnd, event ) )
  {  
    
    this->CopyMatricesFromMRI();
    if ( !this->MapMRIToImage( wnd, event ) )
      return false;
  
    if ( m_volumeRef && m_volumeRef->m_MRIOrigTarget && !m_MRIOrigTarget )
    {
      m_MRIOrigTarget = CreateTargetMRI( m_MRI, m_volumeRef->m_MRIOrigTarget, false );
    } 
  
    // free MRI pixel space
    MRI* mri = MRIallocHeader( m_MRI->width, 
                              m_MRI->height, 
                              m_MRI->depth, 
                              m_MRI->type );
    MRIcopyHeader( m_MRI, mri );
    MRIfree( &m_MRI );
    m_MRI = mri;
  
    return true;
  }
  else
    return false;
}

bool FSVolume::Restore( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event )
{
  if ( LoadMRI( filename, reg_filename, wnd, event ) )
  {  
    // create m_MRITemp for save/rotate
    if ( m_MRITemp )
      MRIfree( &m_MRITemp );
    
    MRI* mri = MRIallocHeader( m_MRI->width, 
                               m_MRI->height, 
                               m_MRI->depth, 
                               m_MRI->type );
    MRIcopyHeader( m_MRI, mri );
    m_MRITemp = m_MRI;
    m_MRI = mri;
  
    return true;
  }
  else
  {
    cerr << "Restore failed." << endl;
    return false;
  }
}

// read in registration file and convert it to tkreg style
bool FSVolume::LoadRegistrationMatrix( const char* filename )
{
  if ( m_matReg )
    ::MatrixFree( &m_matReg );
  m_matReg = NULL;

  if ( MyUtils::HasExtension( wxString::FromAscii(filename),
			      _("xfm") ) )  // MNI style
  {
    MATRIX* m = NULL;
    if ( regio_read_mincxfm( (char*)filename, &m, NULL ) != 0 )
      return false;

    m_matReg = MRItkRegMtx( m_MRIRef, m_MRI, m );
    MatrixFree( &m );
  }
  else if ( MyUtils::HasExtension( wxString::FromAscii(filename), 
				   _("mat") ) )  // fsl style
  {
    wxFFile file( wxString::FromAscii(filename) );
    wxString strg;
    if ( !file.ReadAll( &strg ) )
      return false;
    strg.Replace( _("\n"), _(" ") );
    wxArrayString ar = MyUtils::SplitString( strg, _(" ") );
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
  else if ( MyUtils::HasExtension( wxString::FromAscii(filename), 
				   _("dat") ) )  // tkregister style
  {
    char* subject = NULL;
    float inplaneres, betplaneres, intensity;
    int float2int;
    if ( regio_read_register( (char*)filename, 
			      &subject, 
			      &inplaneres, 
			      &betplaneres, 
			      &intensity, 
			      &m_matReg, 
			      &float2int ) != 0 )
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

bool FSVolume::Create( FSVolume* src_vol, bool bCopyVoxelData, int data_type )
{
  if ( m_MRI )
    ::MRIfree( &m_MRI );
  if ( m_matReg )
    ::MatrixFree( &m_matReg );

  SetMRI( m_MRIRef, src_vol->GetMRI() );
  memcpy( m_RASToVoxelMatrix, 
	  src_vol->m_RASToVoxelMatrix, 
	  16 * sizeof( double ) );
  memcpy( m_VoxelToRASMatrix, 
	  src_vol->m_VoxelToRASMatrix, 
	  16 * sizeof( double ) );
  memcpy( m_MRIToImageMatrix, 
	  src_vol->m_MRIToImageMatrix, 
	  16 * sizeof( double ) );
  memcpy( m_VoxelToVoxelMatrix, 
	  src_vol->m_VoxelToVoxelMatrix, 
	  16 * sizeof( double ) );
  memcpy( m_RASToRASMatrix, 
	  src_vol->m_RASToRASMatrix, 
	  16 * sizeof( double ) );
  memcpy( m_RASToTkRegMatrix, 
	  src_vol->m_RASToTkRegMatrix, 
	  16 * sizeof( double ) );

  m_matReg = MatrixCopy( src_vol->m_matReg, NULL );

  m_fMinValue = src_vol->m_fMinValue;
  m_fMaxValue = src_vol->m_fMaxValue;
  m_bResampleToRAS = src_vol->m_bResampleToRAS;

// SetOriginalOrigin( src_vol->m_fOriginalOrigin );

  if ( data_type == -1 )
    data_type = src_vol->m_MRI->type;
  
  m_MRI = MRIallocSequence( src_vol->m_MRI->width, 
                            src_vol->m_MRI->height, 
                            src_vol->m_MRI->depth,
                            data_type, 1 );
  if ( NULL == m_MRI )
  {
    cerr << "Could not allocate new mri volume." << endl;
    return false;
  }  
  
  if ( bCopyVoxelData )
  {
    MRIcopy( src_vol->m_MRI, m_MRI );
  }
  
  SetMRITarget( src_vol->m_MRITarget );

  if ( src_vol->m_MRIOrigTarget )
  {
    MRI* mri = src_vol->m_MRIOrigTarget;
    m_MRIOrigTarget = MRIallocHeader( mri->width, 
				      mri->height, 
				      mri->depth, 
				      data_type );
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
    if ( data_type != src_vol->GetDataType() )
      cerr << "Can not copy voxel data with different data type." << endl;
    else
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
    if ( !bCopyVoxelData )
    {
      memset( ptr, 
	       0, 
	       m_imageData->GetScalarSize() * nDim[0] * nDim[1] * nDim[2] );
    }
  }
  return true;
}

void FSVolume::SetMRI( MRI*& mri_out, MRI* mri_in )
{
  if ( mri_out )
    ::MRIfree( &mri_out );

  if ( mri_in )
  {
    mri_out = MRIallocHeader( mri_in->width, 
			      mri_in->height, 
			      mri_in->depth, 
			      mri_in->type );
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

    // must use source data type!
    m_MRITarget = MRIallocHeader( mri->width, 
				  mri->height, 
				  mri->depth, 
          m_MRI->type );    
    MRIcopyHeader( mri, m_MRITarget );
  }
}

/*
bool FSVolume::MRIWrite( const char* filename, bool bSaveToOriginalSpace )
{
  char* fn = strdup( filename );
  int err;
  if ( bSaveToOriginalSpace )
    err = ::MRIwrite( m_MRI, fn );
  else
  {
    // save MRITarget temporarily
    MRI* mri = MRIallocHeader( m_MRITarget->width, 
                               m_MRITarget->height, 
                               m_MRITarget->depth, 
                               m_MRITarget->type );
    MRIcopyHeader( m_MRITarget, mri );
    
    // if rotated save data in original target space 
    if ( m_MRIOrigTarget) 
      MRIcopyHeader( m_MRIOrigTarget, m_MRITarget );
    
    err = ::MRIwrite( m_MRITarget, fn );
    // free pixel space
    MRIfreeFrames( m_MRITarget, 0 );
    
    // restore MRITarget saved at the beginning
    MRIcopyHeader( mri, m_MRITarget );
    MRIfree( &mri );
  }
  free( fn );

  if ( err != 0 )
  {
    cerr << "MRIwrite failed" << endl;
  }

  return err == 0;
}
*/

bool FSVolume::MRIWrite( const char* filename, bool bSaveToOriginalSpace )
{
  if ( !m_MRITemp )
  {
    cerr << "Volume not ready for save." << endl;
    return false;
  }
    
  LTA* lta = NULL;
  if ( !bSaveToOriginalSpace && m_MRIOrigTarget) 
  {
    // prepare lta to be saved later
    LINEAR_TRANSFORM *lt ;
    VOL_GEOM srcG, dstG;

    lta = LTAalloc(1, NULL) ;
    lt = &lta->xforms[0] ;
    lt->sigma = 1.0f ;
    lt->x0 = lt->y0 = lt->z0 = 0 ;
    lta->type = LINEAR_VOX_TO_VOX;
    lt->m_L = MRIgetVoxelToVoxelXform( m_MRI, m_MRITemp );

    MRIcopyHeader( m_MRIOrigTarget, m_MRITemp );
    
    getVolGeom(m_MRI,     &srcG);
    getVolGeom(m_MRITemp, &dstG);
    
    lta->xforms[0].src = srcG;
    lta->xforms[0].dst = dstG;
  }
  
  // check if file is writable
  FILE* fp = fopen( filename, "w" );
  if ( !fp )
  {
    cerr << "Failed to open file " << filename 
         << " for write. Please check if the directory exists and you have permission to write in that location." 
         << endl;
    MRIfree( &m_MRITemp );
    m_MRITemp = NULL;
    return false;
  }
  else
    fclose( fp );
  
  char* fn = strdup( filename );
  int err;
  err = ::MRIwrite( m_MRITemp, fn );
  free( fn );

  if ( err != 0 )
  {
    cerr << "MRIwrite failed" << endl;
  }
  else if ( lta )    // save lta file only if volume was saved successfully
  {
    string fname = filename;
    fname = fname + ".lta";
    fp = fopen(fname.c_str(),"w");
    if(!fp)
    {
      cerr << "ERROR: cannot open for writing: " << fname.c_str() << endl;
    }
    LTAprint(fp, lta);
    fclose( fp );
  }
  
  if ( lta )
    LTAfree(&lta);

  MRIfree( &m_MRITemp );
  m_MRITemp = NULL;
  
  return err == 0;
}

bool FSVolume::MRIWrite()
{
  if ( !m_MRITemp )
  {
    cerr << "Volume not ready for save." << endl;
    return false;
  }
  
  int err = ::MRIwrite( m_MRITemp, m_MRI->fname );
  if ( err != 0 )
  {
    cerr << "MRIwrite failed" << endl;
  }
  MRIfree( &m_MRITemp );
  m_MRITemp = NULL;
   
  return err == 0;
}

bool FSVolume::UpdateMRIFromImage( vtkImageData* rasImage, 
				   wxWindow* wnd, 
           wxCommandEvent& event,
           bool resampleToOriginal )
{
  int nProgressStep = ( 30 - event.GetInt() ) / 5;

  MATRIX* vox2vox = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1) = m_VoxelToVoxelMatrix[i];
  }

  // create a target volume
  MRI* mri = MRIallocSequence( m_MRITarget->width, 
			       m_MRITarget->height, 
			       m_MRITarget->depth, 
			       m_MRITarget->type, 
			       m_MRI->nframes );
  if ( mri == NULL )
  {
    cout << "Can not allocate mri volume for buffering." << endl;
    return false;
  }
  
  MRIcopyHeader( m_MRITarget, mri );

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
  
  // create m_MRItemp for writing or rotation
  if ( m_MRITemp )
    MRIfree( &m_MRITemp );
   
  if ( resampleToOriginal )
  {
    m_MRITemp = MRIallocSequence( m_MRI->width, 
                                m_MRI->height, 
                                m_MRI->depth, 
                                m_MRI->type,
                                m_MRI->nframes );
    if ( m_MRITemp == NULL )
    {
      cout << "Can not allocate mri volume for buffering." << endl;
      MRIfree( &mri );
      return false;
    }
    MRIcopyHeader( m_MRI, m_MRITemp );
    MRIvol2Vol( mri, m_MRITemp, vox2vox, m_nInterpolationMethod, 0 );
    MRIfree( &mri );
  }
  else
  {
    m_MRITemp = mri;
  }

  MatrixFree( &vox2vox );
  MRIvalRange( m_MRITemp, &m_fMinValue, &m_fMaxValue );
  return true;
}

double FSVolume::GetVoxelValue( int i, int j, int k, int frame )
{
  if ( i < 0 || i >= m_MRI->width ||
       j < 0 || j >= m_MRI->height ||
       k < 0 || k >= m_MRI->depth )
    return 0;
  else
    return MRIgetVoxVal( m_MRI, i, j, k, frame );
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

MRI* FSVolume::CreateTargetMRI( MRI* src, MRI* refTarget, bool bAllocatePixel )
{
  float cornerFactor[3];
  Real RAS[3], index[3];
  Real indexBounds[6];
  indexBounds[0] = indexBounds[2] = indexBounds[4] = 1e10;
  indexBounds[1] = indexBounds[3] = indexBounds[5] = -1e10;
  for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
  {
    for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
    {
      for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
      {
          ::MRIvoxelToWorld( src,
                             cornerFactor[0]*src->width-0.5,
                             cornerFactor[1]*src->height-0.5,
                             cornerFactor[2]*src->depth-0.5,
                             &RAS[0], &RAS[1], &RAS[2] );
          ::MRIworldToVoxel( refTarget,
                             RAS[0], RAS[1], RAS[2],
                             &index[0], &index[1], &index[2] );
          
          if ( index[0] < indexBounds[0] ) indexBounds[0] = index[0];
          if ( index[0] > indexBounds[1] ) indexBounds[1] = index[0];
          if ( index[1] < indexBounds[2] ) indexBounds[2] = index[1];
          if ( index[1] > indexBounds[3] ) indexBounds[3] = index[1];
          if ( index[2] < indexBounds[4] ) indexBounds[4] = index[2];
          if ( index[2] > indexBounds[5] ) indexBounds[5] = index[2];
      }
    }
  }
  
  // find out the voxelsize of the converted target volume
  MATRIX* mat = MRIgetVoxelToVoxelXform( src, refTarget );
  double m[16];
  for ( int i = 0; i < 16; i++ )
  {
    m[i] = (double) *MATRIX_RELT((mat),(i/4)+1,(i%4)+1);
  }
  double pixelSize[3] = { 1, 1, 1 };
  if ( fabs( m[0] ) > fabs( m[4] ) && fabs( m[0] ) > fabs( m[8] ) )
    pixelSize[0] = src->xsize;
  else if ( fabs( m[4] ) > fabs( m[8] ) )
    pixelSize[1] = src->xsize;
  else
    pixelSize[2] = src->xsize;

  if ( fabs( m[1] ) > fabs( m[5] ) && fabs( m[1] ) > fabs( m[9] ) )
    pixelSize[0] = src->ysize;
  else if ( fabs( m[5] ) > fabs( m[9] ) )
    pixelSize[1] = src->ysize;
  else
    pixelSize[2] = src->ysize;

  if ( fabs( m[2] ) > fabs( m[6] ) && fabs( m[2] ) > fabs( m[10] ) )
    pixelSize[0] = src->zsize;
  else if ( fabs( m[6] ) > fabs( m[10] ) )
    pixelSize[1] = src->zsize;
  else
    pixelSize[2] = src->zsize;
  
  // calculate dimension of the converted target volume
  int dim[3];
  dim[0] = (int)( ( indexBounds[1] - indexBounds[0] ) * refTarget->xsize / pixelSize[0] + 0.5 );
  dim[1] = (int)( ( indexBounds[3] - indexBounds[2] ) * refTarget->ysize / pixelSize[1] + 0.5 );
  dim[2] = (int)( ( indexBounds[5] - indexBounds[4] ) * refTarget->zsize / pixelSize[2] + 0.5 );
  
//  cout << "indexBounds: " << indexBounds[0] << " " << indexBounds[2] << " " << indexBounds[4] << endl;
  
  MRI* mri;
  if ( bAllocatePixel )
    mri = MRIallocSequence( dim[0], dim[1], dim[2], src->type, src->nframes );
  else
    mri = MRIallocHeader( dim[0], dim[1], dim[2], src->type );
  
  if ( mri == NULL )
    return NULL;
  
  MRIcopyHeader( refTarget, mri );
  MRIsetResolution( mri, pixelSize[0], pixelSize[1], pixelSize[2] );
  Real p0[3];
  ::MRIvoxelToWorld( refTarget, 
                     indexBounds[0]+0.5, 
                     indexBounds[2]+0.5, 
                     indexBounds[4]+0.5,
                     &p0[0], &p0[1], &p0[2] );
  MRIp0ToCRAS( mri, p0[0], p0[1], p0[2] );

  return mri;
}

bool FSVolume::MapMRIToImage( wxWindow* wnd, wxCommandEvent& event )
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
    rasMRI = MRIallocSequence( mri->width, 
			       mri->height, 
			       mri->depth, 
			       m_MRI->type, 
			       m_MRI->nframes );
    if ( rasMRI == NULL )
    {
      cerr << "Can not allocate memory for volume transformation" << endl;
      MatrixFree( &m );
      return false;
    }
    MRIcopyHeader( mri, rasMRI );
  }
  else if ( m_bResampleToRAS )
  {
    this->GetBounds( bounds );
    for ( int i = 0; i < 3; i++ )
      dim[i] = (int) ( ( bounds[i*2+1] - bounds[i*2] ) / voxelSize[i] + 0.5 );

    rasMRI = MRIallocSequence( dim[0], dim[1], dim[2], 
			       m_MRI->type, m_MRI->nframes );
    if ( rasMRI == NULL )
    {
      cerr << "Can not allocate memory for volume transformation" << endl;
      MatrixFree( &m );
      return false;
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
    if ( !m_volumeRef )
    {
      double* rtv = this->GetVoxelToRASMatrix();
      int odim[3] = { m_MRI->width, m_MRI->height, m_MRI->depth };
      if ( fabs( rtv[0] ) > fabs( rtv[4] ) && 
	   fabs( rtv[0] ) > fabs( rtv[8] ) )
      {
        *MATRIX_RELT( m, 1, 1 ) = ( rtv[0] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( rtv[0] > 0 ? 0 : odim[0] - 1 );
        dim[0] = odim[0];
      }
      else if ( fabs( rtv[4] ) > fabs( rtv[8] ) )
      {
        *MATRIX_RELT( m, 2, 1 ) = ( rtv[4] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( rtv[4] > 0 ? 0 : odim[0] - 1 );
        dim[1] = odim[0];
      }
      else
      {
        *MATRIX_RELT( m, 3, 1 ) = ( rtv[8] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( rtv[8] > 0 ? 0 : odim[0] - 1 );
        dim[2] = odim[0];
      }

      if ( fabs( rtv[1] ) > fabs( rtv[5] ) && 
	   fabs( rtv[1] ) > fabs( rtv[9] ) )
      {
        *MATRIX_RELT( m, 1, 2 ) = ( rtv[1] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( rtv[1] > 0 ? 0 : odim[1] - 1 );
        dim[0] = odim[1];
      }
      else if ( fabs( rtv[5] ) > fabs( rtv[9] ) )
      {
        *MATRIX_RELT( m, 2, 2 ) = ( rtv[5] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( rtv[5] > 0 ? 0 : odim[1] - 1 );
        dim[1] = odim[1];
      }
      else
      {
        *MATRIX_RELT( m, 3, 2 ) = ( rtv[9] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( rtv[9] > 0 ? 0 : odim[1] - 1 );
        dim[2] = odim[1];
      }

      if ( fabs( rtv[2] ) > fabs( rtv[6] ) && 
	   fabs( rtv[2] ) > fabs( rtv[10] ) )
      {
        *MATRIX_RELT( m, 1, 3 ) = ( rtv[2] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( rtv[2] > 0 ? 0 : odim[2] - 1 );
        dim[0] = odim[2];
      }
      else if ( fabs( rtv[6] ) > fabs( rtv[10] ) )
      {
        *MATRIX_RELT( m, 2, 3 ) = ( rtv[6] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( rtv[6] > 0 ? 0 : odim[2] - 1 );
        dim[1] = odim[2];
      }
      else
      {
        *MATRIX_RELT( m, 3, 3 ) = ( rtv[10] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( rtv[10] > 0 ? 0 : odim[2] - 1 );
        dim[2] = odim[2];
      }

      *MATRIX_RELT( m, 4, 4 ) = 1;
      rasMRI = MRIallocSequence( dim[0], dim[1], dim[2], 
				 m_MRI->type, m_MRI->nframes );
      if ( rasMRI == NULL )
      {
        cerr << "Can not allocate memory for volume transformation" << endl;
        MatrixFree( &m );
        return false;
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
      rasMRI = CreateTargetMRI( m_MRI, m_volumeRef->m_MRITarget );
      if ( rasMRI == NULL )
      {
        cerr << "Can not allocate memory for volume transformation" << endl;
        MatrixFree( &m );
        return false;
      }
    }
  }
  MatrixFree( &m );

  if ( m_matReg && !m_MRIRef )
  {
    cerr << "No target volume available! Cannot use registration matrix." 
	 << endl;
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

      MRIvol2Vol( m_MRI, rasMRI, t2r, m_nInterpolationMethod, 0 );

      // copy vox2vox
      MatrixInverse( t2r, vox2vox );
      for ( int i = 0; i < 16; i++ )
      {
        m_VoxelToVoxelMatrix[i] = 
	  (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
      }

      MatrixFree( &vox2vox );
      MatrixFree( &Tmov );
      MatrixFree( &invTmov );
      MatrixFree( &Ttarg );
      MatrixFree( &t2r );
    }
  }
  else
  {  
    MRIvol2Vol( m_MRI, rasMRI, NULL, m_nInterpolationMethod, 0 );
    MATRIX* vox2vox = MRIgetVoxelToVoxelXform( m_MRI, rasMRI );
    for ( int i = 0; i < 16; i++ )
    {
      m_VoxelToVoxelMatrix[i] = 
	(double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
    }
    MatrixFree( &vox2vox );
  }

  SetMRITarget( rasMRI );
  UpdateRASToRASMatrix();

  if ( !CreateImage( rasMRI, wnd, event ) )
    return false;

  // copy mri pixel data to vtkImage we will use for display
  CopyMRIDataToImage( rasMRI, m_imageData, wnd, event );

  // Need to recalc our bounds at some point.
  m_bBoundsCacheDirty = true;

  ::MRIfree( &rasMRI );
  
  return true;
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

bool FSVolume::CreateImage( MRI* rasMRI, wxWindow* wnd, wxCommandEvent& event )
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
    return false;
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
    Real ras[3], cindex[3];       
   
    ::MRIvoxelToWorld( rasMRI, 0., 0., 0., &ras[0], &ras[1], &ras[2] );
    ::MRIworldToVoxel( m_volumeRef->m_MRITarget, 
                       ras[0], ras[1], ras[2], 
                       &cindex[0], &cindex[1], &cindex[2] );
    
    m_volumeRef->GetImageOutput()->GetOrigin( origin );
    for ( int i = 0; i < 3; i++ )
    {
      if ( fabs(cindex[i]) < 1e-04 )
        cindex[i] = 0;
    } 
    
    origin[0] += cindex[0] * m_volumeRef->m_MRITarget->xsize;
    origin[1] += cindex[1] * m_volumeRef->m_MRITarget->ysize;
    origin[2] += cindex[2] * m_volumeRef->m_MRITarget->zsize;
  }

//  cout << "origin: " << origin[0] << "  " << origin[1] << "  " << origin[2] << endl; fflush(0);
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
    cerr << "Could not allocate scalars array." << endl;
    return false;
  }

  // change the number of components to store tuples
  if ( zFrames > 1 )
  {
    scalars->SetNumberOfComponents( zFrames );
  }

  cValues = zX * zY * zZ;
  if ( !scalars->Allocate( cValues ) )
  {
    cerr << "Could not allocate scalars array." << endl;
    return false;
  }
  scalars->SetNumberOfTuples( zX*zY*zZ );

  // Assign the scalars array to the image.
  imageData->GetPointData()->SetScalars( scalars );
  scalars->Delete();

  if ( wnd )
  {
    event.SetInt( event.GetInt() + ( 100 - event.GetInt() ) / 10 );
    wxPostEvent( wnd, event );
  }
  return true;
}


bool FSVolume::ResizeRotatedImage( MRI* rasMRI, MRI* refTarget, vtkImageData* refImageData,
                                   double* rasPoint,
                                   wxWindow* wnd, wxCommandEvent& event )
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
    return false;
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
  refImageData->GetOrigin( origin );
  
  /*
  float cornerFactor[3];
  Real RAS[3], index[3];
  Real indexBounds[6];
  indexBounds[0] = indexBounds[2] = indexBounds[4] = 1e10;
  indexBounds[1] = indexBounds[3] = indexBounds[5] = -1e10;
  for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
  {
    for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
    {
      for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
      {
          ::MRIvoxelToWorld( rasMRI,
                             cornerFactor[0]*rasMRI->width-0.5,
                             cornerFactor[1]*rasMRI->height-0.5,
                             cornerFactor[2]*rasMRI->depth-0.5,
                             &RAS[0], &RAS[1], &RAS[2] );
          ::MRIworldToVoxel( refTarget,
                             RAS[0], RAS[1], RAS[2],
                             &index[0], &index[1], &index[2] );
          
          if ( index[0] < indexBounds[0] ) indexBounds[0] = index[0];
          if ( index[0] > indexBounds[1] ) indexBounds[1] = index[0];
          if ( index[1] < indexBounds[2] ) indexBounds[2] = index[1];
          if ( index[1] > indexBounds[3] ) indexBounds[3] = index[1];
          if ( index[2] < indexBounds[4] ) indexBounds[4] = index[2];
          if ( index[2] > indexBounds[5] ) indexBounds[5] = index[2];
      }
    }
  } 
  
  origin[0] += indexBounds[0] * refTarget->xsize;;
  origin[1] += indexBounds[2] * refTarget->ysize;
  origin[2] += indexBounds[4] * refTarget->zsize;  
  */
  
  Real vox[3], tvox[3];
  ::MRIworldToVoxel( rasMRI,
                     rasPoint[0], rasPoint[1], rasPoint[2],
                     &vox[0], &vox[1], &vox[2] );
  ::MRIworldToVoxel( refTarget,
                     rasPoint[0], rasPoint[1], rasPoint[2],
                     &tvox[0], &tvox[1], &tvox[2] );
//  cout << "original: " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
  origin[0] -= ( vox[0] * rasMRI->xsize - tvox[0] * refTarget->xsize );
  origin[1] -= ( vox[1] * rasMRI->ysize - tvox[1] * refTarget->ysize );
  origin[2] -= ( vox[2] * rasMRI->zsize - tvox[2] * refTarget->zsize );
  
//  cout << "new: " << origin[0] << " " << origin[1] << " " << origin[2] << endl; fflush(0);
 
  imageData->SetSpacing( rasMRI->xsize, rasMRI->ysize, rasMRI->zsize );
  imageData->SetOrigin( origin[0], origin[1], origin[2] );

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
    cerr << "Could not allocate scalars array." << endl;
    return false;
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
  return true;
}

bool FSVolume::Rotate( std::vector<RotationElement>& rotations, 
		       wxWindow* wnd, 
		       wxCommandEvent& event,
           int nSampleMethod )
{
  if ( rotations.size() == 0 )
    return false;
  
  if ( !m_MRITemp )
  {
    cerr << "Volume not ready for rotation" << endl;
    return false;
  }

  MRI* rasMRI;
  if ( rotations[0].Plane == -1 )   // restore
  {
    if ( !m_MRIOrigTarget )  // try to restore but no where to restore
    {
      cout << "Can not restore because no original space was defined." << endl;
      return false;
    }

    rasMRI = MRIallocSequence( m_MRIOrigTarget->width,
                               m_MRIOrigTarget->height,
                               m_MRIOrigTarget->depth,
                               m_MRI->type,
                               m_MRI->nframes );
    if ( rasMRI == NULL )
    {
      cout << "Can not allocate memory for volume transformation." << endl;
      return false;
    }
    MRIcopyHeader( m_MRIOrigTarget, rasMRI );
  }
  else    // rotate
  {
    MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );

    // calculate rotation matrix
    MATRIX* m_r = MatrixIdentity( 4, NULL ) ;
    for ( size_t i = 0; i < rotations.size(); i++ )
    {
      MATRIX* m_tmp = GetRotationMatrix( rotations[i].Plane, 
					 rotations[i].Angle, 
           rotations[i].Point );
      MATRIX* m_tmp1 = MatrixMultiply( m_tmp, m_r, NULL );

      MatrixCopy( m_tmp1, m_r );
      MatrixFree( &m_tmp );
      MatrixFree( &m_tmp1 );
    }
    
    MATRIX* m_new = MatrixMultiply( m_r, m, NULL );

    rasMRI = MRIallocHeader( m_MRITarget->width,
                               m_MRITarget->height,
                               m_MRITarget->depth,
                               m_MRI->type );
    MRIsetResolution( rasMRI, 
		      m_MRITarget->xsize, 
		      m_MRITarget->ysize, 
		      m_MRITarget->zsize );
    MRIsetVoxelToRasXform( rasMRI, m_new );

    MatrixFree( &m );
    MatrixFree( &m_r );
    MatrixFree( &m_new );

    if ( !m_MRIOrigTarget )  // never rotated before, so save it as original
    {
      m_MRIOrigTarget = MRIallocHeader( m_MRITarget->width, 
					m_MRITarget->height, 
					m_MRITarget->depth, 
					m_MRITarget->type );
      MRIcopyHeader( m_MRITarget, m_MRIOrigTarget );
    }
  }

  /*
  MATRIX* v2v = 
    MRIgetVoxelToVoxelXform( m_MRITarget, rasMRI ); // old target to new target
  MATRIX* vm = MatrixAlloc( 4, 4, v2v->type );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT((vm),(i/4)+1,(i%4)+1) = m_VoxelToVoxelMatrix[i];
  }
  MATRIX* vox2vox = MatrixMultiply( v2v, vm, NULL );
  MatrixInverse( vox2vox, v2v );
  MatrixFree( &vm );
  MatrixFree( &v2v );
  */
  
  // change the field of view of the new target to full coverage of the image data
  if ( rotations[0].Plane != -1 )
  {
    MRI* mri = rasMRI;
    rasMRI = CreateTargetMRI( m_MRIOrigTarget, mri );
    ::MRIfree( &mri );
  }
    
  if ( nSampleMethod < 0 )
    nSampleMethod = rotations[0].SampleMethod;
  if ( rotations[0].Plane == -1 )
    MRIvol2Vol( m_MRITemp, rasMRI, NULL, m_nInterpolationMethod, 0 );
  else
    MRIvol2Vol( m_MRITemp, rasMRI, NULL, nSampleMethod, 0 );
  
  MRIfree( &m_MRITemp );
  m_MRITemp = NULL;

  // copy vox2vox
  MATRIX* vox2vox = MRIgetVoxelToVoxelXform( m_MRI, rasMRI );
  for ( int i = 0; i < 16; i++ )
  {
    m_VoxelToVoxelMatrix[i] = (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &vox2vox );
  
  bool bUpdateImage;
  if ( rotations[0].Plane == -1 )
    bUpdateImage = CreateImage( rasMRI, wnd, event );
  else
    // copy mri pixel data to vtkImage we will use for display
    bUpdateImage = ResizeRotatedImage( rasMRI, m_MRITarget, m_imageData, rotations[0].Point, wnd, event ); 
  
  if ( !bUpdateImage )
  {
    return false;
  }

  CopyMRIDataToImage( rasMRI, m_imageData, wnd, event );
  
  SetMRITarget( rasMRI );
  UpdateRASToRASMatrix();

  // Need to recalc our bounds at some point.
  m_bBoundsCacheDirty = true;

  ::MRIfree( &rasMRI );

  m_bResampleToRAS = false;

  return true;
}

/*
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
*/

MATRIX* FSVolume::GetRotationMatrix( int nPlane, double angle, double* origin )
{
  // calculate rotation matrix
  double p0[3] = { 0, 0, 0 }, p1[3] = { 0, 0, 0 };
  p1[nPlane] = 1;
  TargetToRAS( p0, p0 );
  TargetToRAS( p1, p1 );
  
  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->Translate( origin );
  tr->RotateWXYZ( angle, p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] );
  tr->Translate( -origin[0], -origin[1], -origin[2] );
  tr->Update();

  vtkSmartPointer<vtkMatrix4x4> m = vtkSmartPointer<vtkMatrix4x4>::New();
  tr->GetTranspose( m );
  
  MATRIX* m_r = MatrixZero( 4, 4, NULL );
  for ( int i = 0; i < 4; i++ )
    for ( int j = 0; j < 4; j++ )
      *MATRIX_RELT( m_r, i+1, j+1 ) = m->GetElement( i, j );
  
  return m_r;
}

void FSVolume::CopyMRIDataToImage( MRI* mri, 
				   vtkImageData* image, 
				   wxWindow* wnd, 
				   wxCommandEvent& event )
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

    if ( wnd && nZ%(max(1, zZ/5)) == 0 )
    {
      event.SetInt( event.GetInt() + nProgressStep );
      wxPostEvent( wnd, event );
    }
  }
}

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
                           cornerFactor[0]*(mri->width-1),
                           cornerFactor[1]*(mri->height-1),
                           cornerFactor[2]*(mri->depth-1),
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

  if ( fabs( m[0] ) > fabs( m[4] ) && fabs( m[0] ) > fabs( m[8] ) )
    pixelSize[0] = m_MRI->xsize;
  else if ( fabs( m[4] ) > fabs( m[8] ) )
    pixelSize[1] = m_MRI->xsize;
  else
    pixelSize[2] = m_MRI->xsize;

  if ( fabs( m[1] ) > fabs( m[5] ) && fabs( m[1] ) > fabs( m[9] ) )
    pixelSize[0] = m_MRI->ysize;
  else if ( fabs( m[5] ) > fabs( m[9] ) )
    pixelSize[1] = m_MRI->ysize;
  else
    pixelSize[2] = m_MRI->ysize;

  if ( fabs( m[2] ) > fabs( m[6] ) && fabs( m[2] ) > fabs( m[10] ) )
    pixelSize[0] = m_MRI->zsize;
  else if ( fabs( m[6] ) > fabs( m[10] ) )
    pixelSize[1] = m_MRI->zsize;
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
  {
    double pos[4] = { 0 };
    double* vs = m_imageData->GetSpacing();
    double* origin = m_imageData->GetOrigin();
    for ( int i = 0; i < 3; i++ )
      pos[i] = ( pos_in[i] - origin[i] ) / vs[i];

    Real fpos[3];
    ::MRIvoxelToWorld( m_MRITarget, 
		       (float)pos[0], (float)pos[1], (float)pos[2], 
		       &fpos[0], &fpos[1], &fpos[2] );
//  cout << "out: " << fpos[0] << " " << fpos[1] << " " << fpos[2] << endl;
    for ( int i = 0; i < 3; i++ )
    {
      pos_out[i] = fpos[i];
    }
  }
}

MATRIX* FSVolume::GetTargetToRASMatrix()
{
  double* vs = m_imageData->GetSpacing();
  double* origin = m_imageData->GetOrigin();
  MATRIX* m = MatrixZero( 4, 4, NULL );
  *MATRIX_RELT( m, 1, 1 ) = vs[0];
  *MATRIX_RELT( m, 2, 2 ) = vs[1];
  *MATRIX_RELT( m, 3, 3 ) = vs[2];
  *MATRIX_RELT( m, 4, 4 ) = 1;
  *MATRIX_RELT( m, 1, 4 ) = origin[0];
  *MATRIX_RELT( m, 2, 4 ) = origin[1];
  *MATRIX_RELT( m, 3, 4 ) = origin[2];
  MATRIX* m_invert = MatrixInverse( m, NULL );
  
  MATRIX* v2r = MRIgetVoxelToRasXform( m_MRITarget );
  MATRIX* t2r = MatrixMultiply( m_invert, v2r, NULL );
  MatrixFree( &m );
  MatrixFree( &v2r );
  MatrixFree( &m_invert );
  
  return t2r;
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
  Real pos[4] = { 0 };
  ::MRIworldToVoxel( m_MRITarget, 
	         (float)pos_in[0], 
		       (float)pos_in[1], 
		       (float)pos_in[2],
           &pos[0], 
		       &pos[1], 
		       &pos[2] );
  double* vs = m_imageData->GetSpacing();
  double* origin = m_imageData->GetOrigin();
  for ( int i = 0; i < 3; i++ )
  {
    pos_out[i] = vs[i] * pos[i] + origin[i];
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

void FSVolume::SetInterpolationMethod( int nMethod )
{
  m_nInterpolationMethod = nMethod;
}

int FSVolume::GetDataType()
{
  if ( m_MRI )
    return m_MRI->type;
  else
    return -1;
}
