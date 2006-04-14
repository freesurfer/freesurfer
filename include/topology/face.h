#ifndef TOPOLOGY_FACE_H
#define TOPOLOGY_FACE_H

#ifdef __cplusplus

#include "globals.h"

class Face
{
public:
	int v[3]; // list of vertices
	int f[3]; // list of neighboring faces

	int marked;	//for computational purposes

	double nx,ny,nz; // the face normal
	double x,y,z; // the coordinates of the face

	//constructor/destructor
	Face(void);
	~Face(void);
	const Face& operator=(const Face &face);
};

#endif

#endif
