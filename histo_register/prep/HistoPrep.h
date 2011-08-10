#ifndef _HISTO_PREP_H_
#define _HISTO_PREP_H_
namespace hb {


/*! \file HistoPrep.cc
	\brief The HistoPrep module prepares histology images for registration.
	The code segments the foreground from background, then splits the tissue slices
	(assuming one or two slices per slide), then normalizes each slice.
*/


// register commands, etc. defined in this module
void initHistoPrep();


} // end namespace hb
#endif // _HISTO_PREP_H_
