#include "atlasmeshalphadrawercuda.hpp"

#include "atlasmeshalphadrawercudaimpl.hpp"

namespace kvl {
  namespace cuda {
    void AtlasMeshAlphaDrawerCUDA::SetRegions( const ImageType::RegionType&  region ) {
      	this->tSetRegions.Start();
	auto size = region.GetSize();
	
	Dimension<3, unsigned short> imageDims;
	imageDims[0] = size[2];
	imageDims[1] = size[1];
	imageDims[2] = size[0];

	this->d_Output.SetDimensions(imageDims);
	
	this->d_Output.SetMemory(0);
	
	// Set up the ITK image
	this->image = AtlasMeshAlphaDrawerCUDA::ImageType::New();
	this->image->SetRegions( region );
	this->image->Allocate();
	this->tSetRegions.Stop();
    }

    void AtlasMeshAlphaDrawerCUDA::Interpolate( const kvl::AtlasMesh* mesh ) {
      CudaTetrahedralMesh<double,unsigned long,float> ctm;
      this->tInterpolate.Start();

      this->tSendMesh.Start();
      ctm.Send(mesh);
      this->tSendMesh.Stop();

      this->tKernel.Start();
      RunAtlasMeshAlphaDrawerCUDA( this->d_Output, ctm, this->classNumber );
      this->tKernel.Stop();

      this->tInterpolate.Stop();
    }

    const AtlasMeshAlphaDrawerCUDA::ImageType* AtlasMeshAlphaDrawerCUDA::GetImage() const {
      std::vector<float> tmp;
      CudaImage<float,3,unsigned short>::DimensionType dims;
      
      this->tGetImage.Start();
      
      // Get the image data back
      this->tGetImageTransfer.Start();
      this->d_Output.Recv( tmp, dims );
      this->tGetImageTransfer.Stop();
      
      this->tGetImageUnpack.Start();
      for( unsigned short k=0; k<dims[0]; k++ ) {
	for( unsigned short j=0; j<dims[1]; j++ ) {
	  for( unsigned short i=0; i<dims[2]; i++ ) {
	    auto result = tmp.at(dims.GetLinearIndex(k,j,i));
	    ImageType::IndexType idx;
	    idx[0] = i;
	    idx[1] = j;
	    idx[2] = k;
	    
	    this->image->SetPixel(idx,result);
	  }
	}
      }
      this->tGetImageUnpack.Stop();
      
      this->tGetImage.Stop();
      
      return this->image;
    }
    
    void AtlasMeshAlphaDrawerCUDA::SetClassNumber( const int targetClass ) {
      this->classNumber = targetClass;
    }
  }
}
