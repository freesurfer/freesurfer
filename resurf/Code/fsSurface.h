#ifndef __fsSurface_h
#define __fsSurface_h
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include "vtkSmartPointer.h"
extern "C" 
{
	#include "mrisurf.h"
}

using namespace itk;
namespace fs
{
struct Face
{
	int indexPoint[3];
	float area;
};
struct Edge 
{
	int indexPoint[2];
	float length;
};

template <typename TValueType, unsigned int  VDimension =3>
class Surface : public itk::Mesh<TValueType, VDimension>
{
public:
 
  typedef TValueType                                    ValueType;
  typedef Surface		                          Self;
  typedef SmartPointer<Self>                              Pointer;
  typedef itk::Mesh<TValueType>              Superclass;

  typedef typename  Superclass::PointsContainer PointsContainer;
  typedef typename Superclass::PointType PointType;
  typedef typename Superclass::PointIdentifier PointIdentifier;
  typedef typename Superclass::CellType CellType;
  typedef typename CellType::CellAutoPointer CellPointer;
	
  typedef itk::TriangleCell< CellType > TriangleType;
  //typedef typename TriangleType::Pointer TrianglePointer;
itkNewMacro(Self);
void Load(MRI_SURFACE *surf);
MRI_SURFACE* GetFSSurface(MRI_SURFACE *surf);

private:

std::vector<Face> faces;
std::vector<Edge> edges;
void AddFace(int idPoint1, int idPoint2, int idPoint3);
void AddEdge(int idPoint1, int idPoint2);

};


#include "fsSurface.txx"

}
#endif
