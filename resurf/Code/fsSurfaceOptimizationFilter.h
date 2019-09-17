#ifndef __fsSurfaceOptimizationFilter_h
#define __fsSurfaceOptimizationFilter_h

#include "itkMesh.h"
#include "itkMeshToMeshFilter.h"
#include "fsSurface.h"

using namespace itk;
namespace fs
{

	template <typename TSurfaceIn, typename  TSurfaceOut>
		class SurfaceOptimizationFilter : public itk::MeshToMeshFilter<TSurfaceIn, TSurfaceOut>
	{
		public:

			typedef SurfaceOptimizationFilter	                          Self;
			typedef SmartPointer<Self>       	                       Pointer;
			typedef itk::MeshToMeshFilter<TSurfaceIn, TSurfaceOut>       Superclass;
			typedef TSurfaceIn InSurfaceType;
			typedef typename InSurfaceType::Pointer InSurfacePointer;
			typedef TSurfaceOut OutSurfaceType;
			typedef typename OutSurfaceType::Pointer OutSurfacePointer;
			itkNewMacro(Self);
		protected:

			void GenerateData() override;	
		private:


	};


#include "fsSurfaceOptimizationFilter.txx"

}
#endif
