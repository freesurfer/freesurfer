#ifndef _TEST_CORRES_3D_H_
#define _TEST_CORRES_3D_H_
namespace hb {


/*! \file TestCorres3D.cc
    \brief The TestCorres3D module is used to test 3D correspondence field estimation using 
	synthetic data.  The module generates a synthetic volume and a synthetic 3D correspondence 
	field, then distorts the volume according to the field, and tries to estimate the 
	correspondences from the pair of volumes.
*/


// register commands, etc. defined in this module
void initTestCorres3D();


} // end namespace hb
#endif // _TEST_CORRES_3D_H_
