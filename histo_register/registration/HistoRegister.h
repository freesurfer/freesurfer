#ifndef _HISTO_REGISTER_H_
#define _HISTO_REGISTER_H_
namespace hb {


/*! \file HistoRegister.cc
	\brief The HistoRegister module registers histology with MR, assuming that the MR has 
	already been aligned in 3D using the blockface images.  This histology-to-MRI
	registration is performed using 2D optical flow (with the VarMotion algorithm).
*/


// register commands, etc. defined in this module
void initHistoRegister();


} // end namespace hb
#endif // _HISTO_REGISTER_H_
