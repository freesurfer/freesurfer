/**
 * @file  mris_decimate.h
 * @brief Reduce the number of vertices and faces in a surface. 
 *
 * This tool reduces the number of triangles in a surface using the
 * the vtkDecimatePro class documented at:
 *
 *            http://www.vtk.org/doc/release/5.2/html/a00310.htmls 
 *
 * Please see the VTK documentation for details on the decimation algorithm.
 * mris_decimate will read in an existing surface and write out a new one
 * with less triangles.  The decimation level and other options can be provided
 * on the command-line.
 *
 * By default, this tool is command-line driven.  If you wish to use a VTK-based
 * interactive visualization of the surface, uncomment the '#define BUILT_IN_VISUALIZER'
 * below.
 */

/*
 * Original Author: Dan Ginsburg
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/07/12 19:38:12 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2009,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef MRIS_DECIMATE_H
#define MRIS_DECIMATE_H

extern "C"
{
#include "mrisurf.h"
}

///
//  Types
//

///
/// The follow structure correspond to options for vtkDecimatePro and
/// are settable via the command-line of mris_decimate.  These options
/// all correspond to members of vtkDecimatePro as described here:
///
///   http://www.vtk.org/doc/release/5.2/html/a00310.htmls
///
typedef struct
{
    /// See vtkDecimatePro::SetTargetReduction()
    float reductionLevel;

    /// See vtkDecimatePro::SetPreserveTopology()
    bool   setPreserveTopology;
    int    preserveTopology;

    /// See vtkDecimatePro::SetSplitting()
    bool   setSplitting;
    int    splitting;

    /// See vtkDecimatePro::SetFeatureAngle()
    bool   setFeatureAngle;
    double featureAngle;

    /// See vtkDecimatePro::SetBoundaryVertexDeletion()
    bool   setBoundaryVertexDeletion;
    int    boundaryVertexDeletion;

    /// See vtkDecimatePro::SetDegree()   
    bool   setDegree;
    int    degree;

} DECIMATION_OPTIONS;


///
/// \fn int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions, char *outFilePath)
/// \brief This function performs decimation on the input surface and outputs the new surface to a file.
/// \param mris Input loaded MRI_SURFACE to decimate
/// \param decimationOptions Options controlling the decimation arguments (see DECIMATION_OPTIONS)
/// \param outFilePath Full path to output file to write decimated surface to
/// \return 0 on success, 1 on failure
///
int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions, char *outFilePath);

#endif // MRIS_DECIMATE_H

