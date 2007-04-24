/**
 * @file  vtkInflatePolyData.h
 * @brief Inflates (or deflates) a surface by moving verts along the normals
 *
 * This VTK filter class takes vtkPolyData as in put and outputs
 * vtkPolyData. It translates all vertices along the average normal
 * for all polygons that use that vertex by InflateFactor. A >0 value
 * will inflate the poly data, and a <0 value will deflate
 * it. Modified from vtkCleanPolyData.cxx.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/24 19:54:54 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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

#ifndef vtkInflatePolyData_h
#define vtkInflatePolyData_h

#include "vtkPolyDataAlgorithm.h"

class vtkInflatePolyData : public vtkPolyDataAlgorithm {

public:

  static vtkInflatePolyData *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeRevisionMacro(vtkInflatePolyData,vtkPolyDataAlgorithm);

  // Description:
  // The factor to use. >0 will inflate, and <0 will deflate.
  vtkSetClampMacro(InflateFactor,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(InflateFactor,double);

protected:
  vtkInflatePolyData();
 ~vtkInflatePolyData();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  double InflateFactor;

private:
  vtkInflatePolyData(const vtkInflatePolyData&);  // Not implemented.
  void operator=(const vtkInflatePolyData&);  // Not implemented.
};

#endif
