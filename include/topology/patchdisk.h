#ifndef TOPOLOGY_PATCHDISK_H
#define TOPOLOGY_PATCHDISK_H

#include "surface.h"

class PatchDisk
{
	void _Init();
public:
	Surface disk;
	Loop ring,init_ring;
	int *vtrans ;
	int *ftrans ;

	PatchDisk(int which_patch);
	PatchDisk(const string s):disk(s),init_ring(10){
		vtrans = new int[disk.nvertices];
		ftrans = new int[disk.nfaces];
		_Init();
	};
	~PatchDisk(void);

	void Init();
};

#endif
