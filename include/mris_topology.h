// include guard
#ifndef MRIS_TOPOLOGY_INCLUDED_H
#define MRIS_TOPOLOGY_INCLUDED_H

// The following is usable from C
#ifdef __cplusplus
extern "C" {
#endif

#include "mrisurf.h"

	//functions
	bool MRIScorrectDefect(MRIS *mris, int defect_number);
	int MRISgetEulerNumber(const MRIS *mris, const int *list_of_faces, int nfs);
	MRIP* MRIPextractFromMRIS(MRIS *mris, int defect_number);
	void MRISinitSurface(MRIS *mris);

#ifdef __cplusplus
}
#endif



// C++ portion starts here
#ifdef __cplusplus

extern "C" {
#include "mrisurf.h"
#include "error.h"
}
#include "topology/surface.h"

MRIP *MRIPalloc(int nvertices, int nfaces);
void MRIPfree(MRIP **mrip);
MRIP *MRIPclone(MRIP *src);

Surface *MRIStoSurface(MRIS *mris);
void SurfaceToMRIP(Surface &surface, MRIS &mris);
bool MRISaddMRIP(MRIS *mris_dst, MRIP *mrip);

#endif


#endif
