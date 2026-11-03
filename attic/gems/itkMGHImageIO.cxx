
#include <fstream>

#include "itkMGHImageIO.h"


#include <cmath>

#include <stdio.h>
#include <stdlib.h>

//-------------------------------
//
// Convert to BE
//
//-------------------------------

template<class T>
int TReadZ(gzFile iFile,
	   T&     out)
{
  T* pt = new T(0);
  int result;
  result = gzread(iFile,
		  pt,
		  sizeof(T));
  itk::ByteSwapper<T>::SwapFromSystemToBigEndian(pt);
  out = *pt;
  delete pt;
  return result;
}

//-------------------------------


static
std::string
GetExtension( const std::string& filename)
{
  const std::string::size_type pos = 
    filename.find_last_of(".");
  std::string extension(filename,
			pos+1,
			filename.length());
  return extension;
}

//============================================================

class STLWrapper
{
public:
  explicit STLWrapper(std::ofstream* p_ofs) 
    : m_ptrOfs(p_ofs), 
      m_own(false) {}
  explicit STLWrapper(const char* fname) 
    : m_ptrOfs(NULL),
      m_own(true)
  {
    m_ptrOfs = new std::ofstream(fname,
			      std::ios::out | std::ios::binary
			      );
  }
  virtual ~STLWrapper()
  {
    if ( m_own )
      delete m_ptrOfs;
  }
  
  template<class T>
  void Write(T value)
  {
    T* pt = new T(value);
    itk::ByteSwapper<T>::SwapFromSystemToBigEndian(pt);
    m_ptrOfs->write( (const char*)pt, sizeof(T) );
    delete pt;
  }
  virtual void WriteBuffer(char* buffer, 
			   unsigned int size)
  {
    m_ptrOfs->write(buffer, size *sizeof(char));
  }
private:
  std::ofstream* m_ptrOfs;
  bool m_own;
};

//========================

class GZWrapper
{
public:
  explicit GZWrapper(gzFile of) 
    : m_ofs(of), m_own(false) {}
  explicit GZWrapper(const char* fname) 
    : m_ofs(NULL), 
      m_own(true)
  {
    m_ofs = gzopen( fname, "wb" );
  }
  virtual ~GZWrapper() 
  {
    if ( m_own )
      gzclose(m_ofs);
  }

  template<class T>
  void Write(T value)
  {
    T* pt = new T(value);
    itk::ByteSwapper<T>::SwapFromSystemToBigEndian(pt);
    ::gzwrite(m_ofs, pt, sizeof(T));
    delete pt;
  }
  virtual void WriteBuffer(char* buffer, 
			   unsigned int size)
  {
    ::gzwrite( m_ofs, buffer, size *sizeof(char) );
  }
private:
  gzFile m_ofs;
  bool m_own;
};

//===========================================================================


//--------------------------------------
//
// MGHImageIO
//

namespace itk
{

  MGHImageIO::MGHImageIO() 
  {
    this->SetNumberOfDimensions(3);
    const unsigned int uzero = 0;
    m_Dimensions[0] = uzero;
    m_Dimensions[1] = uzero;
    m_Dimensions[2] = uzero;

    if ( ByteSwapper<int>::SystemIsBigEndian())
      m_ByteOrder = BigEndian;
    else
      m_ByteOrder = LittleEndian;

  }

  MGHImageIO::~MGHImageIO()
  {
  }
    

  bool
  MGHImageIO::CanReadFile(const char* FileNameToRead)
  {
    std::string filename(FileNameToRead);

    if ( filename == "" )
      {
	itkExceptionMacro(<<"A FileName must be specified.");
	return false;
      }

    
    // check if the correct extension is given by the user
    std::string extension = GetExtension(filename);
    if ( extension == std::string("mgh") ||
	 extension == std::string("mgz") )
      return true;
    
    if ( extension == std::string("gz") )
      {
	if ( filename.substr( filename.size() - 7 ) == std::string(".mgh.gz") )
	  return true;
      }
    
    return false;
  }
  
  void
  MGHImageIO::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    std::string strSep = ", ";

    os << indent 
       << "Data Dimensions: (" 
       << m_Dimensions[0] << strSep
       << m_Dimensions[1] << strSep
       << m_Dimensions[2] << ")\n"
       << indent
       << "Data Spacing: (" 
       << m_Spacing[0] << strSep
       << m_Spacing[1] << strSep
       << m_Spacing[2] << ")\n"
       << indent
       << "Scalar Type: " << m_ComponentType << std::endl
       << indent
       << "Number of Frames: " << m_NumberOfComponents << std::endl;
    
    os << indent << "RAS to IJK matrix: " << std::endl;
  }


  void
  MGHImageIO::ReadImageInformation()
  {
    gzFile     fp;
    
    fp = gzopen( m_FileName.c_str(), "rb");
    if ( !fp )
      {
	itkExceptionMacro(<<"Can't find/open file: " << m_FileName);
	return;
      }

    ReadVolumeHeader(fp);

    gzclose(fp);
  }

  void
  MGHImageIO::ReadVolumeHeader(gzFile fp)
  {
    int        version;
    int        bufInt; // buffer -> int type (most ITK types are unsigned)
    int        type;
    int        dof;
    short      RASgood;

    // check file reading

    if ( !fp )
      {
	itkExceptionMacro(<< "Can't find/open file: " << this->m_FileName);
	return;
      }
    TReadZ(fp, version);
    TReadZ(fp, bufInt);   m_Dimensions[0] = (unsigned int)bufInt;
    TReadZ(fp, bufInt);   m_Dimensions[1] = (unsigned int)bufInt;
    TReadZ(fp, bufInt);   m_Dimensions[2] = (unsigned int)bufInt;
    // next is nframes
    TReadZ(fp, bufInt);   m_NumberOfComponents = (unsigned int)bufInt;
    TReadZ(fp, type);
    TReadZ(fp, dof); // what does this do?

    // Convert type to an ITK type
    switch( type )
      {
      case fs::MRI_UCHAR: 
	m_ComponentType = UCHAR; break;
      case fs::MRI_INT:
	m_ComponentType = INT; break;
      case fs::MRI_FLOAT:
	m_ComponentType = FLOAT; break;
      case fs::MRI_SHORT:
	m_ComponentType = SHORT; break;
      case fs::MRI_TENSOR:
	m_ComponentType = FLOAT; m_NumberOfComponents = 9; break;
      default:
	itkExceptionMacro(<<" Unknown data type " << type << " using float by default.");
	m_ComponentType = FLOAT;
      }

    // Next short says whether RAS registration information is good. 
    // If so, read the voxel size and then the matrix
    TReadZ( fp, RASgood);
    float spacing;
    if ( RASgood )
      {
	for( int nSpacing = 0; nSpacing<3; ++nSpacing)
	  {
	    TReadZ( fp, spacing); // type is different
	    m_Spacing[nSpacing] = spacing;
	  }
	/*
	  From http://www.nmr.mgh.harvard.edu/~tosa/#coords:
	  To go from freesurfer voxel coordinates to RAS coordinates, they use:
	  translation:  t_r, t_a, t_s is defined using c_r, c_a, c_s centre voxel position in RAS
	  rotation: direction cosines x_(r,a,s), y_(r,a,s), z_(r,a,s)
	  voxel size for scale: s_x, s_y, s_z
    
	  [ x_r y_r z_r t_r][s_x  0   0  0]
	  [ x_a y_a z_a t_a][0   s_y  0  0]
	  [ x_s y_s z_s t_s][0    0  s_z 0]
	  [  0   0   0   1 ][0    0   0  1]
	  Voxel center is a column matrix, multipled from the right
	  [v_x]
	  [v_y]
	  [v_z]
	  [ 1 ]
    
	  In the MGH header, they hold:
	  x_r x_a x_s
	  y_r y_a y_s
	  z_r z_a z_s
	  c_r c_a c_s
	*/
	typedef itk::Matrix<double> MatrixType;
	MatrixType matrix;
  
	float fBuffer;
	float c[3];

	for(unsigned int uj=0; uj<3; ++uj)
	  {
	    for(unsigned int ui=0; ui<3; ++ui)
	      {
		TReadZ(fp, fBuffer);
		matrix[ui][uj] = fBuffer;
	      }
	  }

	for(unsigned int ui=0; ui<3; ++ui)
	  TReadZ(fp, c[ui]);

	for(unsigned int ui=0; ui<3; ++ui)
	  {
	    std::vector<double> vDir;
#if KVL_ORIENTATION_HACK
            vDir.push_back( -matrix[ 0 ][ui] );
            vDir.push_back( -matrix[ 1 ][ui] );
            vDir.push_back( matrix[ 2 ][ui] );
#else
            for(unsigned int uj=0; uj<3; ++uj)
	      vDir.push_back( matrix[uj][ui] );
#endif            
	    SetDirection( ui, vDir );
	  }

	float idxCenter[3];
	for(unsigned int ui=0; ui<3; ++ui)
	  idxCenter[ui] = static_cast<float>(m_Dimensions[ui])/2.0f;

	for( unsigned int ui=0; ui<3; ++ui)
	  {
	    m_Origin[ui] = c[ui];
	    for(unsigned int uj=0; uj<3; ++uj)
	      m_Origin[ui] -= matrix[ui][uj]*m_Spacing[uj]*idxCenter[uj];
	  }
#if KVL_ORIENTATION_HACK
            m_Origin[ 0 ] = -m_Origin[ 0 ];
            m_Origin[ 1 ] = -m_Origin[ 1 ];
#endif            
      }
   
    unsigned long numValues = m_Dimensions[0]*m_Dimensions[1]*m_Dimensions[2];
    gzseek(fp, fs::FS_WHOLE_HEADER_SIZE + 
	   ( m_NumberOfComponents * numValues * this->GetComponentSize() ),
	   SEEK_SET);

    float fBuf;
    // read TR, Flip, TE, FI, FOV
    if ( TReadZ( fp, fBuf) )
      {
	itk::MetaDataDictionary &thisDic = this->GetMetaDataDictionary();
	itk::EncapsulateMetaData<float>(thisDic,
					std::string("TR"),
					fBuf);

	// try to read flipAngle
	if ( TReadZ( fp, fBuf ) )
	  {
	    itk::EncapsulateMetaData<float>(thisDic,
					    std::string("FlipAngle"),
					    fBuf);
	    // TE
	    if ( TReadZ( fp, fBuf ) )
	      {
		itk::EncapsulateMetaData<float>(thisDic,
						std::string("TE"),
						fBuf);
		// TI
		if ( TReadZ(fp, fBuf) )
		  {
		    itk::EncapsulateMetaData<float>(thisDic,
						    std::string("TI"),
						    fBuf);
		    // FOV
		    if ( TReadZ(fp, fBuf) )
		      {
			itk::EncapsulateMetaData<float>(thisDic,
							std::string("FoV"),
							fBuf);
		      }
		  }
	      }
	  }
      } // end if

    //==================
    // read tags at the end of file
    
  }

  void
  MGHImageIO::Read(void* pData)
  {
    gzFile fp;
    fp = gzopen( m_FileName.c_str(), "rb");
    if ( !fp )
      {
	itkExceptionMacro(<<"Can't find/open file: " << m_FileName);
	return;
      }

    const unsigned long numPixels = m_Dimensions[0]*m_Dimensions[1]* m_Dimensions[2];

    const unsigned int componentSize( this->GetComponentSize() );
    
    // check that the offset is actually computed wrt. the beginning
    gzseek(fp, fs::FS_WHOLE_HEADER_SIZE, SEEK_SET );
    
    const unsigned int frameSize = numPixels*componentSize;

    if ( m_NumberOfComponents > 1  )
      {
	char* pBuffer = new char[ frameSize ];
  
	const unsigned int pixelSize = componentSize*m_NumberOfComponents;
  
	for(unsigned int frameIndex = 0;
	    frameIndex < m_NumberOfComponents;
	    ++frameIndex)
	  {
	    // read current frame
	    gzread( fp, pBuffer, frameSize ); 
	    // copy memory location in the final buffer
      
	    char* pSrc = (char*)pBuffer;
	    char* pDst = (char*)pData;
      
	    pDst += frameIndex * componentSize;
	    for( unsigned int ui=0;
		 ui < numPixels; 
		 ++ui, pSrc+=componentSize, pDst+=pixelSize)
	      {
		for(unsigned int byteCount = 0;
		    byteCount < componentSize; ++byteCount)
		  *(pDst+byteCount) = *(pSrc+byteCount);
	      } // next ui
	  } // next frameIndex

	// clear resources
	delete[] pBuffer;
      }
    else
      {
	gzread( fp, pData, frameSize);
      }
    
    gzclose(fp);
    
    SwapBytesIfNecessary( pData, numPixels*m_NumberOfComponents );

  } // end Read function

  void
  MGHImageIO::SwapBytesIfNecessary(void* buffer,
				   unsigned long numberOfPixels)
  {
    // NOTE: If machine order is little endian, and the data needs to be
    // swapped, the SwapFromBigEndianToSystem is equivalent to
    // SwapFromSystemToBigEndian.

    switch(m_ComponentType)
      {
      case UCHAR:
	ByteSwapper<unsigned char>::SwapRangeFromSystemToBigEndian((unsigned char*)buffer,
								   numberOfPixels);
	break;
      case SHORT:
	ByteSwapper<short>::SwapRangeFromSystemToBigEndian((short*)buffer,
							   numberOfPixels);
	break;
      case INT:
	ByteSwapper<int>::SwapRangeFromSystemToBigEndian((int*)buffer,
							 numberOfPixels);
	break;
      case FLOAT:
	ByteSwapper<float>::SwapRangeFromSystemToBigEndian((float*)buffer,
							   numberOfPixels);
	break;
      default:
	ExceptionObject exception(__FILE__,__LINE__);
	exception.SetDescription("Pixel Type Unknown");
	throw exception;
      } // end switch
  }


  bool
  MGHImageIO::CanWriteFile(const char* name)
  {
    std::string filename(name);

    if ( filename == "" )
      {
	itkExceptionMacro(<<"A FileName must be specified.");
	return false;
      }

    std::string extension = GetExtension(filename);
    if ( extension != std::string("mgh") &&
	 extension != std::string("mgz") )
      return false;

    return true;
  }

  void
  MGHImageIO::WriteImageInformation()
  {
    std::string extension = GetExtension(m_FileName);

    if ( extension == std::string("mgh") )
      {
	std::ofstream ofs( m_FileName.c_str(),
			   std::ios::out | std::ios::binary 
			   );
	STLWrapper writer(&ofs);
	this->WriteHeader( writer );
      }
    else
      {
	gzFile fp = gzopen(m_FileName.c_str(), "wb");
	if(!fp)
	  itkExceptionMacro(<<" Failed to open gzFile for writing");
	GZWrapper writer(fp);
	this->WriteHeader(writer);
	gzclose(fp);
      }
  }
  

  void
  MGHImageIO::Write(const void* buffer)
  {
    std::string extension = GetExtension(m_FileName);

    if ( extension == std::string("mgh") )
      {
	std::ofstream ofs( m_FileName.c_str(),
			   std::ios::out | std::ios::binary );
	if ( ofs.fail() )
	  {
	    itkExceptionMacro(<< "File cannot be opened for writing");
	  }
	STLWrapper writer(&ofs);
	this->WriteHeader(writer);
	this->WriteData(writer,buffer);
      }
    else
      {
	gzFile file_p = gzopen( m_FileName.c_str(), "wb");
	if (!file_p)
	  itkExceptionMacro(<<" Failed to open gzFile for writing");
	GZWrapper writer(file_p);
	this->WriteHeader(writer);
	this->WriteData(writer, buffer);
  
	gzclose(file_p);
      }

  }

  void
  MGHImageIO::PermuteFrameValues(const void* buffer,
				 char* tempmemory)
  {
    const unsigned int numPixels =  m_Dimensions[0]
      * m_Dimensions[1] * m_Dimensions[2];

// DJ -- these weren't used
//    const unsigned long int numvalues = numPixels * m_NumberOfComponents;
//    const unsigned long int numbytes = this->GetComponentSize() * numvalues;    

    const unsigned int valueSize( this->GetComponentSize() );

    const unsigned int frameSize = numPixels * valueSize;
    
    const char* pSrc = (const char*)buffer;
    char* pDst = (char*)tempmemory;
    
    for(unsigned int pixelIndex = 0;
	pixelIndex < numPixels; ++pixelIndex, pDst+=valueSize )
      {
	for(unsigned int componentIndex =0 ;
	    componentIndex < m_NumberOfComponents;
	    ++componentIndex, pSrc+=valueSize)
	  {
	    std::copy( pSrc, pSrc+ valueSize,
		       pDst + frameSize*componentIndex );
	  } // next component index
      } // next pixelIndex 
  }

  unsigned int
  MGHImageIO::GetComponentSize() const
  {
    unsigned int returnValue = 0;
    switch ( m_ComponentType )
      {
      case UCHAR: returnValue = sizeof(unsigned char); break;
      case SHORT: returnValue = sizeof(short); break;
      case INT:   returnValue = sizeof(int); break;
      case FLOAT: returnValue = sizeof(float); break;
      
      // DJ -- added this in to get the compiler to shut up
      case UNKNOWNCOMPONENTTYPE:
      case CHAR:
      case USHORT:
      case UINT:
      case ULONG:
      case LONG:
      case DOUBLE:
        break;
      }
    return returnValue;
  }

} // end NAMESPACE ITK

