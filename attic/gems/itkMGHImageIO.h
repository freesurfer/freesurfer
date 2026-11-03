

#ifndef H_ITK_IMAGE_IO_H
#define H_ITK_IMAGE_IO_H

// STL includes

// ITK includes
#include "itkImageIOBase.h"
#include "itkIOCommon.h"
#if ITK_VERSION_MAJOR >= 5
#include "itkMacro.h"
#else
#include "itkExceptionObject.h"
#endif
#include "itkByteSwapper.h"
#include "itkMetaDataObject.h"
#include "itkMatrix.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>

#include <itk_zlib.h>

#define KVL_ORIENTATION_HACK 1 // This MGH reader/writer IO class seems to swap x and y-axis compared
                               // to NIFTI. Hard-coded a correction here; this really should be 
                               // investigated further.

// variables used in the IO
//
// because of the templated functions, need to declare them here
// try to avoid name scoping
namespace fs
{
const int MRI_UCHAR = 0;
const int MRI_INT   = 1;
const int MRI_FLOAT = 3;
const int MRI_SHORT = 4;
const int MRI_TENSOR= 6;

const int FS_DIMENSION_HEADER_SIZE = sizeof(int) * 7;
const int FS_RAS_HEADER_SIZE = (sizeof(float) * 15) + sizeof(short);
const int FS_UNUSED_HEADER_SIZE = 256 - FS_RAS_HEADER_SIZE;
const int FS_WHOLE_HEADER_SIZE =
  FS_RAS_HEADER_SIZE + FS_DIMENSION_HEADER_SIZE + FS_UNUSED_HEADER_SIZE;

};

namespace itk
{

  class OutputStreamWrapper;

  class ITK_EXPORT MGHImageIO : public ImageIOBase
  {
  public:
    typedef MGHImageIO              Self;
    typedef ImageIOBase             Superclass;
    typedef SmartPointer<Self>      Pointer;

    /** Method for creation through the object factory **/
    itkNewMacro(Self);
    /** RTTI (and related methods) **/
    itkTypeMacro(MGHImageIO, Superclass);

    
    /**--------------- Read the data----------------- **/
    virtual bool CanReadFile(const char* FileNameToRead);
    /* Set the spacing and dimension information for the set file name */
    virtual void ReadImageInformation();
    /* Read the data from the disk into provided memory buffer */
    virtual void Read(void* buffer);
    
    /**---------------Write the data------------------**/
    
    virtual bool CanWriteFile(const char* FileNameToWrite);
    /* Set the spacing and dimension information for the set file name */
    virtual void WriteImageInformation();
    /* Write the data to the disk from the provided memory buffer */
    virtual void Write(const void* buffer);

  protected:
    MGHImageIO();
    ~MGHImageIO();
    void PrintSelf(std::ostream& os, Indent indent) const;

    void ReadVolumeHeader(gzFile fp);

  private:
    // processes the actual data buffer
    void SwapBytesIfNecessary(void* buffer, 
            unsigned long numberOfPixels);
    // examines the direction cosines and creates encapsulation data
    //void MriDirCos();

    template<class Writer>
    void WriteHeader(Writer& writer);

    template<class Writer>
      void WriteData(Writer& writer, const void* buffer);

    void PermuteFrameValues(const void* buffer, char* tempmemory);

    unsigned int GetComponentSize() const;

  }; // end class declaration

  //-----------------------------------
  //
  // WriteHeader
  //
  //-----------------------------------
  template<class Writer>
  void
  MGHImageIO::WriteHeader(Writer& writer)
  {
    const int mghVersion = 1;
    
    // version
    writer.Write( mghVersion );

    // dimensions
    for(unsigned int ui=0; ui<3; ++ui)
      writer.Write( (int)m_Dimensions[ui] );

    // nframes
    writer.Write( (int)m_NumberOfComponents );

    // type
    switch( m_ComponentType )
      {
      case UCHAR:
        writer.Write( fs::MRI_UCHAR ); break;
      case USHORT:
        itkWarningMacro( << " Casting from USHORT to INT!!!" );
      case INT:
        writer.Write( fs::MRI_INT ); break;
      case DOUBLE:
        itkWarningMacro( << " Casting from DOUBLE to INT!!!" );
      case FLOAT:
        writer.Write( fs::MRI_FLOAT ); break;
      case SHORT:
        writer.Write( fs::MRI_SHORT ); break;

      // DJ -- added these cases to make the compiler shut up
      case UNKNOWNCOMPONENTTYPE:
      case CHAR:
      case UINT:
      case ULONG:
      case LONG:
        break;

      }

    // dof !?! -> default value = 1
    writer.Write( (int)1 );
    
    // write RAS and voxel size info
    writer.Write( (short)1 );
    
    // spacing
    for(unsigned int ui=0; ui<3; ++ui)
      writer.Write( (float)m_Spacing[ui] );

    //=================================================
    // get directions matrix
    std::vector<std::vector<double> > vvRas;
#if KVL_ORIENTATION_HACK
    for(unsigned int ui=0; ui<3; ++ui) 
      {
      std::vector< double > direction = GetDirection(ui);  
      direction[ 0 ] = -direction[ 0 ];
      direction[ 1 ] = -direction[ 1 ];
      vvRas.push_back( direction );
      }
#else    
    for(unsigned int ui=0; ui<3; ++ui) vvRas.push_back( GetDirection(ui) );
#endif    
    // transpose data before writing it
    std::vector<float> vBufRas;
    // transpose matrix
    for(unsigned int ui=0; ui<3; ++ui)
      for(unsigned int uj=0; uj<3; ++uj)
	vBufRas.push_back( (float)vvRas[ui][uj] );
    //==================================================

    for(std::vector<float>::const_iterator cit = vBufRas.begin();
	cit != vBufRas.end(); ++cit )
      writer.Write( *cit );

    // write c_r, c_a, c_s
    // defined as origin + DC x resolution x ( dim0/2 , dim1/2, dim2/2 )
    // 
    float crasBuf;
    for(unsigned int ui=0; ui<3; ++ui)
      {
	crasBuf = m_Origin[ui];
#if KVL_ORIENTATION_HACK
        if ( ui < 2 ) crasBuf = -crasBuf; 
#endif        
	for(unsigned int uj=0; uj<3; ++uj)
	  crasBuf += vvRas[uj][ui]*m_Spacing[uj]*(float)m_Dimensions[uj]/2.0f;
	writer.Write( crasBuf );
      } // next ui

    // fill the rest of the buffer with zeros
    char* buffer = new char[ fs::FS_UNUSED_HEADER_SIZE *sizeof(char) ];
    memset( buffer, 0, fs::FS_UNUSED_HEADER_SIZE *sizeof(char) );
    writer.WriteBuffer( buffer, fs::FS_UNUSED_HEADER_SIZE );
    delete[] buffer;
  }

  //--------------------------------------------
  //
  // WriteData
  //
  //--------------------------------------------
  template<class Writer>
    void
    MGHImageIO::WriteData(Writer& writer,
			  const void* buffer)
    {
      // swap bytes if necessary
      const unsigned long int numvalues = m_Dimensions[0]
	* m_Dimensions[1] * m_Dimensions[2] * m_NumberOfComponents;
      const unsigned long int numbytes = this->GetComponentSize() * numvalues;

      char* tempmemory = new char[numbytes];

      // re-arrange data in frames
      if ( m_NumberOfComponents > 1 )
	{
	  PermuteFrameValues(buffer, tempmemory);
	}
      else
	memcpy(tempmemory, buffer, numbytes);


      this->SwapBytesIfNecessary(tempmemory,numvalues);
    
      writer.WriteBuffer( tempmemory, numbytes );

      delete[] tempmemory;

      itk::MetaDataDictionary &thisDic=this->GetMetaDataDictionary();
    
      float fScanBuffer=0.0;

      if ( ExposeMetaData<float>(thisDic, "TR", fScanBuffer) )
	{
	  writer.Write(fScanBuffer);
	  if ( ExposeMetaData<float>(thisDic, "FlipAngle", fScanBuffer) )
	    {
	      writer.Write(fScanBuffer);
	      if ( ExposeMetaData<float>(thisDic, "TE", fScanBuffer) )
		{
		  writer.Write(fScanBuffer);
		  if ( ExposeMetaData<float>(thisDic, "TI", fScanBuffer) )
		    {
		      writer.Write(fScanBuffer);
		      if ( ExposeMetaData<float>(thisDic, "FoV", fScanBuffer) )
			writer.Write(fScanBuffer);
		    } // end TI
		} // end TE
	    } // end FlipAngle
	} // end TR

    }


} // end namespace itk


#endif // H_ITK_IMAGE_IO_H



