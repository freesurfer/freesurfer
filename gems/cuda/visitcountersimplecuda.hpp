#pragma once

#include <vector>

#include "atlasmeshvisitcounter.hpp"

#include "dimensioncuda.hpp"
#include "cudaimage.hpp"
#include "stopwatch.hpp"

#include "visitcountersimplecudaimpl.hpp"

namespace kvl {
  namespace cuda {
    template<typename T,typename Internal>
    class VisitCounterSimple : public kvl::interfaces::AtlasMeshVisitCounter {
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
	this->image = VisitCounterSimple<T,Internal>::ImageType::New();
	this->image->SetRegions( region );
	this->image->Allocate();
	this->tSetRegions.Stop();
      }

      virtual void VisitCount( const kvl::AtlasMesh* mesh ) override {
	this->tVisitCount.Start();
	std::vector<AtlasMesh::CellIdentifier> tetrahedronIds;

	// Find the tetrahedra
	this->tVisitCountPack.Start();
	for( auto cellIt = mesh->GetCells()->Begin();
	     cellIt != mesh->GetCells()->End();
	     ++cellIt ) {
	  if( cellIt.Value()->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL ) {
	    tetrahedronIds.push_back( cellIt.Index() );
	  }
	}


	// Extract co-ordinates
	Dimension<3, unsigned long> tetArrDims;
	std::vector<T> tetrahedra;
	tetArrDims[0] = tetrahedronIds.size();
	tetArrDims[1] = nVertices;
	tetArrDims[2] = nDims;
	tetrahedra.resize(tetArrDims.ElementCount());

	for( int iTet=0; iTet<tetrahedronIds.size(); iTet++ ) {
	  AtlasMesh::CellAutoPointer cell;
	  mesh->GetCell( tetrahedronIds.at(iTet), cell );

	  auto pit = cell->PointIdsBegin();
	  unsigned long iVertex = 0;
	  for( auto pit = cell->PointIdsBegin(); pit != cell->PointIdsEnd(); ++pit ) {
	    AtlasMesh::PointType p;
	    mesh->GetPoint( *pit, &p );

	    for( unsigned long i=0; i<nDims; i++ ) {
	      size_t idx = tetArrDims.GetLinearIndex((unsigned long)iTet,iVertex,i);
	      tetrahedra.at(idx) = p[i];
	    }

	    iVertex++;
	  }
	}
	this->tVisitCountPack.Stop();

	CudaImage<T,3,size_t> d_tetrahedra;

	this->tVisitCountTransfer.Start();
	d_tetrahedra.Send(tetrahedra, tetArrDims);
	this->tVisitCountTransfer.Stop();

	this->tVisitCountKernel.Start();
	RunVisitCounterSimpleCUDA<T,Internal>( d_Output, d_tetrahedra );
	this->tVisitCountKernel.Stop();

	this->tVisitCount.Stop();
      };

      virtual const VisitCounterSimple<T,Internal>::ImageType* GetImage() const override {
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
      const int nDims = 3;
      const int nVertices = 4;

      CudaImage<int,3,unsigned short> d_Output;
      VisitCounterSimple<T,Internal>::ImageType::Pointer image;
    };
  }
}
