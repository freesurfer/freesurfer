#ifndef _BLOCK_REGISTER_H_
#define _BLOCK_REGISTER_H_
namespace hb {


/*! \file BlockRegister.cc
	\brief The BlockRegister module registers a blockface volume with an MR volume.
	First the algorithm performs a linear registration and re-samples the MR 
	to match the blockface according to this registration.  Then the alogirhtm
	performs a non-linear registration (using the VarCorres3D algorithm) and
	re-samples the MR to match the blockface according to this registration. 
*/


// register commands, etc. defined in this module
void initBlockRegister();


} // end namespace hb
#endif // _BLOCK_REGISTER_H_
