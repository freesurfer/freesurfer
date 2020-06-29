/**
 * @brief Reduce the number of vertices and faces in a surface. 
 *
 * This tool reduces the number of triangles in a surface using the
 * the GNU Triangulated Surface Library documented at:
 *
 *            http://gts.sourceforge.net/reference/book1.html
 *
 * Please see the GTS documentation for details on the decimation algorithm.
 * mris_decimate will read in an existing surface and write out a new one
 * with less triangles.  The decimation level and other options can be provided
 * on the command-line.
 *
 */

/*
 * Original Author: Dan Ginsburg
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef MRIS_DECIMATE_H
#define MRIS_DECIMATE_H



#include "mrisurf.h"


///
//  Types
//

///
/// The follow structure contains options for the DECIMATION algorithm
///	provided by the GNU Trianglulated Surface Library (gts)
///
typedef struct
{
  /// A value between (0, 1.0) that controls how decimated the surface
  /// is.  A value of 1.0 means that the surface is not decimated at all
  /// and a value close to 0 means that it will be reduced to nearly no
  /// edges.  More specifically, this controls the number of edges that
  /// the decimated mesh will contain.
  float decimationLevel;
  
  // Set the desired number of faces and then compute the decimationLevel from that
  int desiredNumFaces;
  
  // Set the desired average face area and then compute the decimationLevel from that
  float desiredFaceArea;
  
  ///	The minimum angle between two neighboring triangles allowed by
  /// the decimation
  bool setMinimumAngle;
  float minimumAngle;
	
  // Run code to sort output vertices. This was needed when compiling
  // GTS with hashes instead of btrees to make the output
  // deterministic. Even then it did not always give the same output
  // when the input was ?h.orig.nofix, because the underlying
  // algorithm did not give the same vertices. The algorithm appears
  // to give deterministic output when using btrees, so the sorting
  // feature is off by default.
  int SortVertices;

} DECIMATION_OPTIONS;

///	Optional callback function 
typedef void (*DECIMATE_PROGRESS_FUNC)(float percentDone, const char *msg, void *userData);

///
/// \fn int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions,
///							DECIMATE_PROGRESS_FUNC decimateProgressFn)
/// \brief This function performs decimation on the input surface and outputs the new surface to a file.
/// \param mris Input loaded MRI_SURFACE to decimate
/// \param decimationOptions Options controlling the decimation arguments (see DECIMATION_OPTIONS)
/// \param decimateProgressFn If non-NULL, provides updates on decimation percentage complete and
///							  a status message that can, for example, be displayed in a GUI.
/// \param userData If decimateProgressFn is non-NULL, argument passed into decimateProgressFn
/// \return 0 on success, 1 on failure
///
int decimateSurface(MRI_SURFACE **mris, const DECIMATION_OPTIONS &decimationOptions, 
    				DECIMATE_PROGRESS_FUNC decimateProgressFn = NULL,
    				void *userData = NULL);

#endif // MRIS_DECIMATE_H

