#pragma once

#include <vector>

#include "atlasmeshvisitcounter.hpp"

#include "dimensioncuda.hpp"
#include "cudaimage.hpp"

#include "visitcountersimplecudaimpl.hpp"

namespace kvl {
  namespace cuda {
    template<typename T>
    class VisitCounterSimple : public kvl::interfaces::AtlasMeshVisitCounter {
      virtual void SetRegions( const kvl::interfaces::AtlasMeshVisitCounter::ImageType::RegionType& region ) override {
	auto size = region.GetSize();

	Dimension<3, unsigned short> imageDims;
	imageDims[0] = size[2];
	imageDims[1] = size[1];
	imageDims[2] = size[0];

	this->d_Output.SetDimensions(imageDims);

	this->d_Output.SetMemory(0);
	
	// Set up the ITK image
	this->image = VisitCounterSimple<T>::ImageType::New();
	this->image->SetRegions( region );
	this->image->Allocate();
      }

      virtual void VisitCount( const kvl::AtlasMesh* mesh ) override {

	std::vector<AtlasMesh::CellIdentifier> tetrahedronIds;

	// Find the tetrahedra
	for( auto cellIt = mesh->GetCells()->Begin();
	     cellIt != mesh->GetCells()->End();
	     ++cellIt ) {
	  if( cellIt.Value()->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL ) {
	    tetrahedronIds.push_back( cellIt.Index() );
	  }
	}

	std::cout << "Found " << tetrahedronIds.size() << " tetrahedra" << std::endl;

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

	CudaImage<T,3,size_t> d_tetrahedra;

	d_tetrahedra.Send(tetrahedra, tetArrDims);

	RunVisitCounterSimpleCUDA( d_Output, d_tetrahedra );

	std::cout << __FUNCTION__ << ": Complete" << std::endl;
      };

      virtual const VisitCounterSimple<T>::ImageType* GetImage() const override {
	std::vector<int> tmp;
	CudaImage<int,3,unsigned short>::DimensionType dims;

	// Get the image data back
	this->d_Output.Recv( tmp, dims );

	for( unsigned short k=0; k<dims[0]; k++ ) {
	  for( unsigned short j=0; j<dims[1]; j++ ) {
	    for( unsigned short i=0; i<dims[2]; i++ ) {
	      int result = tmp.at(dims.GetLinearIndex(k,j,i));
	      ImageType::IndexType idx;
	      idx[0] = i;
	      idx[1] = j;
	      idx[2] = k;
	      
	      this->image->SetPixel(idx,result);
	    }
	  }
	}

	return this->image;
      }
      
    private:
      const int nDims = 3;
      const int nVertices = 4;

      CudaImage<int,3,unsigned short> d_Output;
      VisitCounterSimple<T>::ImageType::Pointer image;
    };
  }
}
