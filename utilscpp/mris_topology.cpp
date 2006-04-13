#include "mris_topology.h"
#include <iostream>
using namespace std;

#include "topology/loop.h"
#include "topology/segment.h"
#include "topology/surface.h"

extern "C" void MRIScutSurface(MRIS *surf_in){
	int i;
	i=0;

	Loop loop;
	Segment s;
	Surface surf;
	cout << "this is a test"<<endl;
}

