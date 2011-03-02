/**
 * @file  mriio_nrrd_itk.cpp
 * @brief Provides Nrrd read support of diffusion data to Freesurfer
 *
 * Implements mriNrrdRead and mriNrrdWrite using ITK library.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:46 $
 *    $Revision: 1.4 $
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

#ifndef HAVE_ITK_LIBS
// if Freesurfer is configured *without* ITK libs (that is, if --with-itk-dir
// wss not specified on the configure line), then include these stubs:
extern "C" {
#include "error.h"
#include "mri.h"
MRI *mriNrrdReadDiffusion(char *fname, int read_volume)
{
  ErrorReturn
    (NULL,
     (ERROR_UNSUPPORTED,
      "mriNrrdReadDiffusion(): Nrrd input of diffusion data not supported!"));
  return NULL;
}
int mriNrrdWriteDiffusion(MRI *mri, char *fname)
{
  ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED, 
      "mriNrrdWriteDiffusion(): Nrrd diffusion data output not supported!"));
}
}
#else // ITK libs available, so implement IO...

#include <iostream>

#include <ctype.h>
#include <unistd.h>
#include <memory.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#include <dirent.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>


extern "C" {
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "minc_volume_io.h"
#include "region.h"
#include "machine.h"
#include "analyze.h"
#include "fio.h"
#include "mri_identify.h"
#include "signa.h"
#include "fio.h"
#include "matfile.h"
#include "math.h"
#include "matrix.h"
#include "diag.h"
#include "chklc.h"
}

#ifdef UCHAR
#undef UCHAR
#endif
#ifdef USHORT
#undef USHORT
#endif
#ifdef UINT
#undef UINT
#endif
#ifdef ULONG
#undef ULONG
#endif

#include "itkImage.h"
//#include "itkExceptionObject.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "itkDiffusionTensor3D.h"

using namespace std;

extern "C" {
  MRI *mriNrrdReadDiffusion(char *fname, int read_volume);
  int mriNrrdWriteDiffusion(MRI *mri, char *fname);
}

MRI *mriNrrdReadDiffusion(char *fname, int read_volume)
{
  typedef itk::DiffusionTensor3D<float> PixelType;
  typedef itk::Image<PixelType,3> DiffusionImageType;
  typedef itk::ImageFileReader<DiffusionImageType> ReaderType;

  MRI *mri = NULL;
  int type = MRI_UCHAR;

  // Read the file.
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName (fname);

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
  //  io->SetNrrdVectorType (nrrdKindList);
  io->SetFileType (itk::ImageIOBase::ASCII);

  reader->SetImageIO (io);
  reader->Update();

  // Get the component type.
  switch (io->GetComponentType()) 
  {
  case itk::ImageIOBase::UCHAR: type = MRI_UCHAR; break;
  case itk::ImageIOBase::USHORT:
  case itk::ImageIOBase::SHORT: type = MRI_SHORT; break;
  case itk::ImageIOBase::ULONG:
  case itk::ImageIOBase::LONG: type = MRI_LONG; break;
  case itk::ImageIOBase::DOUBLE:
  case itk::ImageIOBase::FLOAT: type = MRI_FLOAT; break;
  default:
    cerr << "Unsupported type: "
         << io->GetComponentTypeAsString(io->GetComponentType()) << endl;
    return NULL;
  }

  // Get a pointer to the image.
  DiffusionImageType::Pointer image = reader->GetOutput();

  // Get the bounds (the size of the largest possible region).
  const DiffusionImageType::RegionType& region =
    image->GetLargestPossibleRegion();
  const DiffusionImageType::SizeType size = 
    region.GetSize();
  cerr << "Size " << size[0] << " " 
       << size[1] << " " << size[2] << endl;

  // Get the pixel type and the number of frames.
  DiffusionImageType::PixelType value;
  cerr << "Size of value " << value.Size() << endl;

  // Make a mri with the correct number of frames.
  mri = MRIallocSequence( size[0], size[1], size[2], type, value.Size() );
  if( NULL == mri ) {
    cerr << "mriNrrdReadDiffusion(): Couldn't allocate MRI of size "
         << size[0] << " " << size[1] << " " << size[2] << endl;
    return NULL;
  }

  // Copy all the pixel data.
  DiffusionImageType::IndexType index;
  for( unsigned int z = 0; z < size[2]; z++ ) {
    index[2] = z;
    for( unsigned int y = 0; y < size[1]; y++ ) {
      index[1] = y;
      for( unsigned int x = 0; x < size[0]; x++ ) {
        index[0] = x;
        value = image->GetPixel( index );
        switch( type ) {
        case MRI_UCHAR:
          MRIseq_vox( mri, x, y, z, 0 ) = (BUFTYPE)value[0];
          MRIseq_vox( mri, x, y, z, 1 ) = (BUFTYPE)value[1];
          MRIseq_vox( mri, x, y, z, 2 ) = (BUFTYPE)value[2];
          break;
        case MRI_INT:
          MRIIseq_vox( mri, x, y, z, 0 ) = (int)value[0];
          MRIIseq_vox( mri, x, y, z, 1 ) = (int)value[1];
          MRIIseq_vox( mri, x, y, z, 2 ) = (int)value[2];
          break;
        case MRI_LONG:
          MRILseq_vox( mri, x, y, z, 0 ) = (long)value[0];
          MRILseq_vox( mri, x, y, z, 1 ) = (long)value[1];
          MRILseq_vox( mri, x, y, z, 2 ) = (long)value[2];
          break;
        case MRI_FLOAT:
          MRIFseq_vox( mri, x, y, z, 0 ) = value[0];
          MRIFseq_vox( mri, x, y, z, 1 ) = value[1];
          MRIFseq_vox( mri, x, y, z, 2 ) = value[2];
          break;
        case MRI_SHORT:
          MRISseq_vox( mri, x, y, z, 0 ) = (short)value[0];
          MRISseq_vox( mri, x, y, z, 1 ) = (short)value[1];
          MRISseq_vox( mri, x, y, z, 2 ) = (short)value[2];
          break;
        default:
          break;
        }
      }
    }
  }

  // Get and set the center.
  const DiffusionImageType::PointType origin =
    image->GetOrigin();
  cerr << "Origin " << origin[0] << " " 
       << origin[1] << " " << origin[2] << endl;
  mri->c_r = origin[0];
  mri->c_a = origin[1];
  mri->c_s = origin[2];

  // Get and set the spacing.
  const DiffusionImageType::SpacingType spacing =
    image->GetSpacing();
  cerr << "Spacing " << spacing[0] << " " 
       << spacing[1] << " " << spacing[2] << endl;
  mri->xsize = spacing[0];
  mri->ysize = spacing[1];
  mri->zsize = spacing[2];

  typedef itk::MetaDataDictionary DictionaryType;
  typedef itk::MetaDataObject<string> DictionaryStringType;

  DictionaryType& headers = io->GetMetaDataDictionary();

  DictionaryType::ConstIterator iHeader = headers.Begin();
  DictionaryType::ConstIterator end = headers.End();

  while( iHeader != end ) {
    itk::MetaDataObjectBase::Pointer entry = iHeader->second;

    DictionaryStringType::Pointer value =
      dynamic_cast<DictionaryStringType*>( entry.GetPointer() );

    if( value ) {
      string sValue = value->GetMetaDataObjectValue();
      string sKey = iHeader->first;
      cerr << "--- " << sKey << " = " << sValue << endl;
    }

    ++iHeader;
  }

  return mri;
}

int mriNrrdWriteDiffusion(MRI *mri, char *fname)
{
  //just give an error until write function is complete and tested
  ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED, 
      "mriNrrdWriteDiffusion(): "
      "Nrrd diffusion data output not yet supported"));
}

#endif // #ifndef HAVE_ITK_LIBS

