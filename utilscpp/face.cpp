#include "topology/face.h"

Face::Face(void)
{
	marked=0;
}

Face::~Face(void)
{
}

const Face &Face::operator=(const Face &face)
{
	for(int n = 0 ; n < 3 ; n++){
		v[n] = face.v[n];
		f[n] = face.f[n];
	}
	marked=face.marked;	
	nx=face.nx; ny=face.ny; nz=face.nz; 
	x=face.x; y=face.y; z=face.z; 

	return face;
}
