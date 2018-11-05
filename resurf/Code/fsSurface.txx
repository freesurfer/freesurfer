#ifndef __fsSurface_txx
#define __fsSurface_txx

template< typename TValueType ,unsigned int VDimension >
void fs::Surface<TValueType, VDimension>::Load(MRI_SURFACE *surf)
{
	//typename Surface::Pointer mesh = Self::New();
	typename PointsContainer::Pointer points = PointsContainer::New();
	points->Reserve( surf->nvertices );

	for( int i = 0; i < surf->nvertices; i++ )
	{
		PointType p;
		p[0]=surf->vertices[i].x;
		p[1]=surf->vertices[i].y;
		p[2]=surf->vertices[i].z;
		points->SetElement( i, p );
	}

	this->SetPoints( points );
	CellPointer cellpointer;
	//typedef typename CellType::CellAutoPointer CellPointer;

	for( int i = 0; i < surf->nfaces; i++ )
	{
		cellpointer.TakeOwnership( new TriangleType );
		for(int j=0;j<3;j++)
		{
			cellpointer->SetPointId( j, surf->faces[i].v[j]);
		}
		this->SetCell( i, cellpointer );
		
		this->AddFace(surf->faces[i].v[0],surf->faces[i].v[1],surf->faces[i].v[2]);
	}
	for( int i=0;i<surf->nedges; i++)
	{
		// int vtxno[4]; // vertex numbers of 2 ends + 2 opposites
		this->AddEdge(surf->edges[i].vtxno[0],surf->edges[i].vtxno[1]);
	}
}

template< typename TValueType ,unsigned int VDimension >
MRI_SURFACE* fs::Surface<TValueType, VDimension>::GetFSSurface(MRI_SURFACE *surf)
{
	//MRI_SURFACE *surf;
	for(int i=0; i< this->GetNumberOfPoints();i++)
	{
		typename Self::PointType p;
		if (this->GetPoint(i,&p))
		{
			std::cout << "point "<< p[0]+1 << std::endl;
			surf->vertices[i].x = p[0]+1;
			surf->vertices[i].y = p[1];
			surf->vertices[i].z = p[2];
		}
	}
	return surf;
}
template<typename TValueType, unsigned int VDimension>
void fs::Surface<TValueType, VDimension>::AddFace(int idPoint1, int idPoint2, int idPoint3)
{
	Face f;
	f.indexPoint[0]= idPoint1;
	f.indexPoint[1]= idPoint2;
	f.indexPoint[2]= idPoint3;
	faces.push_back(f);
}
template<typename TValueType, unsigned int VDimension>
void fs::Surface<TValueType, VDimension>::AddEdge(int idPoint1, int idPoint2)
{
	Edge e;
	e.indexPoint[0]= idPoint1;
	e.indexPoint[1]= idPoint2;
	e.length = this->GetPoint(idPoint1).EuclideanDistanceTo(this->GetPoint(idPoint2));
	edges.push_back(e);
}
#endif
