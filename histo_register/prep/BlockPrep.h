#ifndef _BLOCK_PREP_H_
#define _BLOCK_PREP_H_
namespace hb {


/*! \file BlockPrep.cc
    \brief The BlockPrep module is used to prepare block-face images for subsequent segmentation.
	This includes segmenting the foreground from background, registering across microtome resets,
	and resizing/cropping to better match the MR data. 
*/


// register commands, etc. defined in this module
void initBlockPrep();


} // end namespace hb
#endif // _BLOCK_PREP_H_
