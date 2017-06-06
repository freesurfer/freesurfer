#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"

const unsigned int nDims = 3;
const unsigned int nVertices = 4;

template<typename ArgType,typename InvertType>
class SimpleSharedTetrahedron {
public:
  __device__
  SimpleSharedTetrahedron( ArgType tetrahedron[nVertices][nDims],
			   ArgType M[nDims][nDims] ) : tet(&tetrahedron[0][0]),
						       transf(&M[0][0]) {
    // In this class, we assume that the arguments passed in the constructor
    // are actually in shared memory
  } 
  
  __device__
  void LoadAndBoundingBox( const kvl::cuda::Image_GPU<ArgType,3,size_t> tetrahedra,
			   const size_t iTet,
			   unsigned short min[nDims],
			   unsigned short max[nDims] ) {
    // We presume that min and max are in shared memory
    if( (threadIdx.x < nDims) && (threadIdx.y==0) ) {
      for( unsigned int iVert=0; iVert<nVertices; iVert++ ) {
	this->tet[(iVert*nDims)+threadIdx.x] = tetrahedra(iTet,iVert,threadIdx.x);
      }
      
      // No need to sync since we're only using the first 3 threads
      ArgType minVal = this->tet[(0*nDims)+threadIdx.x];
      ArgType maxVal = this->tet[(0*nDims)+threadIdx.x];
      for( unsigned int iVert=1; iVert<nVertices; iVert++ ) {
	ArgType nxt = this->tet[(iVert*nDims)+threadIdx.x];
	if( nxt < minVal ) {
	  minVal = nxt;
	}
	if( nxt > maxVal ) {
	  maxVal = nxt;
	}
      }
      
      // Indices are always positive
      if( minVal < 0 ) {
	minVal = 0;
      }
      
      // Add one to the max to make the loops
      // simpler later
      // It also avoids some pathological cases of
      // planar tetrahedra with all integral vertices
      min[threadIdx.x] = floor(minVal);
      max[threadIdx.x] = ceil(maxVal)+1;
    }
    __syncthreads();
  }

  __device__
  void ComputeBarycentricTransform() {
    // Compute barycentric co-ordinate conversion matrix
    // This is taken from kvlTetrahedronInteriorConstIterator.hxx
    
    // Do single threaded
    if( (threadIdx.x==0) && (threadIdx.y==0) ) {
      // Start by computing locations relative to the first vertex
      // of the tetrahedron
      const InvertType a = this->tet[(1*nDims)+0] - this->tet[(0*nDims)+0];
      const InvertType b = this->tet[(2*nDims)+0] - this->tet[(0*nDims)+0];
      const InvertType c = this->tet[(3*nDims)+0] - this->tet[(0*nDims)+0];
      const InvertType d = this->tet[(1*nDims)+1] - this->tet[(0*nDims)+1];
      const InvertType e = this->tet[(2*nDims)+1] - this->tet[(0*nDims)+1];
      const InvertType f = this->tet[(3*nDims)+1] - this->tet[(0*nDims)+1];
      const InvertType g = this->tet[(1*nDims)+2] - this->tet[(0*nDims)+2];
      const InvertType h = this->tet[(2*nDims)+2] - this->tet[(0*nDims)+2];
      const InvertType i = this->tet[(3*nDims)+2] - this->tet[(0*nDims)+2];
    
      const InvertType A = ( e * i - f * h );
      const InvertType D = -( b * i - c * h );
      const InvertType G = ( b * f - c * e );
      const InvertType B = -(d * i - f * g );
      const InvertType E = ( a * i - c * g );
      const InvertType H = -( a * f - c * d );
      const InvertType C = ( d * h - e * g );
      const InvertType F = - (a * h - b * g );
      const InvertType I = ( a * e - b * d );
      
      const InvertType determinant = a * A + b * B + c * C;
    
      this->transf[(0*nDims) + 0] = A / determinant;
      this->transf[(1*nDims) + 0] = B / determinant;
      this->transf[(2*nDims) + 0] = C / determinant;
      this->transf[(0*nDims) + 1] = D / determinant;
      this->transf[(1*nDims) + 1] = E / determinant;
      this->transf[(2*nDims) + 1] = F / determinant;
      this->transf[(0*nDims) + 2] = G / determinant;
      this->transf[(1*nDims) + 2] = H / determinant;
      this->transf[(2*nDims) + 2] = I / determinant;
    }
    __syncthreads();
  }

  __device__
  void TransformToBarycentric( ArgType p[nDims], const ArgType x, const ArgType y, const ArgType z ) const {
    ArgType r[nDims];
    
    // Compute location relative to first vertex
    r[0] = x - this->tet[(0*nDims)+0];
    r[1] = y - this->tet[(0*nDims)+1];
    r[2] = z - this->tet[(0*nDims)+2];

    for( unsigned int i=0; i<nDims; i++ ) {
      p[i] = 0;
      for( unsigned int j=0; j<nDims; j++ ) {
	p[i] += this->transf[(i*nDims)+j] * r[j];
      }
    }
  }

  __device__
  bool PointInside( const ArgType x, const ArgType y, const ArgType z ) const {
    bool inside = true;
	  
    ArgType pTmp[nDims];
    
    this->TransformToBarycentric(pTmp, x, y, z);

    // Form the full set of barycentric co-ordinates
    ArgType p[nDims+1];
    for( unsigned int i=0; i< nDims; i++ ) {
      p[i+1] = pTmp[i];
    }
    p[0] = 1 - p[1] - p[2] - p[3];

    // Do the easiest cull
    if( (p[0] < 0) || (p[1] < 0) || (p[2] < 0) || (p[3] < 0) ) {
      inside = false;
    } else {
      // Have to handle case where one of the barycentric co-ordinates is zero

      // This comes down to seeing what would happen if the voxel were moved a little to the right, or up etc.
      ArgType nxtRowAdd[nVertices];
      nxtRowAdd[0] = -( this->transf[(0*nDims) + 0] + this->transf[(1*nDims) + 0] + this->transf[(2*nDims) + 0]);
      nxtRowAdd[1] = this->transf[(0*nDims) + 0];
      nxtRowAdd[2] = this->transf[(1*nDims) + 0];
      nxtRowAdd[3] = this->transf[(2*nDims) + 0];

      ArgType nxtColAdd[nVertices];
      nxtColAdd[0] = -( this->transf[(0*nDims) + 1] + this->transf[(1*nDims) + 1] + this->transf[(2*nDims) + 1]);
      nxtColAdd[1] = this->transf[(0*nDims) + 1];
      nxtColAdd[2] = this->transf[(1*nDims) + 1];
      nxtColAdd[3] = this->transf[(2*nDims) + 1];

      ArgType nxtSliceAdd[nVertices];
      nxtSliceAdd[0] = -( this->transf[(0*nDims) + 2] + this->transf[(1*nDims) + 2] + this->transf[(2*nDims) + 2]);
      nxtSliceAdd[1] = this->transf[(0*nDims) + 2];
      nxtSliceAdd[2] = this->transf[(1*nDims) + 2];
      nxtSliceAdd[3] = this->transf[(2*nDims) + 2];

      // Loop over the barycentric co-ords, checking each for zero
      for( unsigned int iVert=0; iVert<nVertices; iVert++ ) {
	if( inside && (p[iVert]==0) ) {
	  if( this->CheckBorder( nxtRowAdd[iVert], nxtColAdd[iVert], nxtSliceAdd[iVert] ) ) {
	    inside = false;
	  }
	}
      }
    }
    return inside;
  }

  __device__
  bool CheckBorder( ArgType a, ArgType b, ArgType c ) const {
    // The arguments a, b and c describe how the barycentric co-ordinate
    // would change in response to shifting the location slightly
    // along each co-ordinate axis
    if ( a < 0 ) {
      return true;  
    }
    
    if ( a == 0 ) {
      if ( b < 0 ) {
	return true;  
      }
      
      if ( b == 0 ) {
	if ( c < 0 ) {
	  return true;  
        }
      }
    }
    
    return false; 
  }

private:
  ArgType* tet;
  ArgType* transf;
};

template<typename T,typename Internal>
__global__
void SimpleVisitCounterKernel( kvl::cuda::Image_GPU<int,3,unsigned short> output,
			       const kvl::cuda::Image_GPU<T,3,size_t> tetrahedra ) {
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);
  
  // Check if this block has an assigned tetrahedron
  if( iTet >= tetrahedra.dims[0] ) {
    return;
  }

  // Load the tetrahedron and determine bounding box
  __shared__ T tetrahedron[nVertices][nDims];
  __shared__ unsigned short min[nDims], max[nDims];
  __shared__ T M[nDims][nDims];
  SimpleSharedTetrahedron<T,Internal> tet(tetrahedron, M);

  tet.LoadAndBoundingBox( tetrahedra, iTet, min, max );

  tet.ComputeBarycentricTransform();

  // Figure out how to cover the bounding box with the current thread block
  // We assume that each thread block is strictly 2D

  // Divide the bounding box into blocks equal to the blockDim
  for( unsigned short iyStart=min[1]; iyStart<max[1]; iyStart += blockDim.y ) {
    for( unsigned short ixStart=min[0]; ixStart<max[0]; ixStart += blockDim.x ) {
      const unsigned short ix = ixStart + threadIdx.x;
      const unsigned short iy = iyStart + threadIdx.y;

      // Could probably do this test a little better
      if( output.PointInRange(0,iy,ix) ) {

	for( unsigned short iz=min[2]; iz<max[2]; iz++ ) {
	  bool inside = tet.PointInside(ix,iy,iz);
	  
	  if( inside ) {
	    atomicAdd(&output(iz,iy,ix),1);
	  }
	}
      }
    }
  }
}


namespace kvl {
  namespace cuda {
    template<typename T,typename Internal>
    void SimpleVisitCounter( CudaImage<int,3,unsigned short>& d_output,
			     const CudaImage<T,3,size_t>& d_tetrahedra ) {
      const unsigned int nBlockx = 1024;

      const unsigned int nThreadsx = 32;
      const unsigned int nThreadsy = 32;
      const unsigned int nThreadsz = 1;

      dim3 grid, threads;

      const size_t nTetrahedra = d_tetrahedra.GetDimensions()[0];
      
      if( nTetrahedra > nBlockx ) {
	grid.x = nBlockx;
	grid.y = (nTetrahedra / grid.x)+1;
	if( (grid.y * grid.x) < nTetrahedra ) {
	  grid.y++;
	}
      } else {
	grid.x = nTetrahedra;
	grid.y = 1;
      }

      threads.x = nThreadsx;
      threads.y = nThreadsy;
      threads.z = nThreadsz;
      
      // Run the kernel
      auto err = cudaGetLastError();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
      SimpleVisitCounterKernel<T,Internal><<<grid,threads>>>( d_output.getArg(), d_tetrahedra.getArg() );
      err = cudaDeviceSynchronize();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
    }

    // -----------------------------------------------------------

    template<>
    void RunVisitCounterSimpleCUDA<float,float>( CudaImage<int,3,unsigned short>& d_output,
						 const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter<float,float>( d_output, d_tetrahedra );
    }

    template<>
    void RunVisitCounterSimpleCUDA<double,double>( CudaImage<int,3,unsigned short>& d_output,
						   const CudaImage<double,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter<double,double>( d_output, d_tetrahedra );
    }

    template<>
    void RunVisitCounterSimpleCUDA<float,double>( CudaImage<int,3,unsigned short>& d_output,
						  const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter<float,double>( d_output, d_tetrahedra );
    }
  }
}
