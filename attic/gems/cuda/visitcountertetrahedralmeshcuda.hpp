#pragma once

#include "atlasmeshvisitcounter.hpp"
#include "cudatetrahedralmesh.hpp"

#include "visitcountertetrahedralmeshcudaimpl.hpp"

namespace kvl {
  namespace cuda {
    class VisitCounterTetrahedralMesh : public kvl::interfaces::AtlasMeshVisitCounter {
    public:
      virtual void SetRegions( const kvl::interfaces::AtlasMeshVisitCounter::ImageType::RegionType& region ) override {
	this->tSetRegions.Start();
	auto size = region.GetSize();
	
	Dimension<3, unsigned short> imageDims;
	imageDims[0] = size[2];
	imageDims[1] = size[1];
	imageDims[2] = size[0];

	this->d_Output.SetDimensions(imageDims);
	
	this->d_Output.SetMemory(0);
	
	// Set up the ITK image
	this->image = VisitCounterTetrahedralMesh::ImageType::New();
	this->image->SetRegions( region );
	this->image->Allocate();
	this->tSetRegions.Stop();
      }

      virtual void VisitCount( const kvl::AtlasMesh* mesh ) override {
	this->tVisitCount.Start();
	CudaTetrahedralMesh<double,unsigned long,float> ctm;

	this->tVisitCountPack.Start();
	ctm.Send(mesh);
	this->tVisitCountPack.Stop();

	this->tVisitCountKernel.Start();
	RunVisitCounterTetrahedralMeshCUDA( this->d_Output, ctm );
	this->tVisitCountKernel.Stop();

	this->tVisitCount.Stop();
      }

      virtual const VisitCounterTetrahedralMesh::ImageType* GetImage() const override {
	std::vector<int> tmp;
	CudaImage<int,3,unsigned short>::DimensionType dims;
	
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
      
      mutable kvl::Stopwatch tSetRegions, tVisitCount, tGetImage;

      mutable kvl::Stopwatch tVisitCountPack, tVisitCountTransfer, tVisitCountKernel;
      mutable kvl::Stopwatch tGetImageTransfer, tGetImageUnpack;

    private:
      CudaImage<int,3,unsigned short> d_Output;
      VisitCounterTetrahedralMesh::ImageType::Pointer image;
    };
  }
}
