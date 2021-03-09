/**
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
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef vtkInflatePolyData_h
#define vtkInflatePolyData_h

#include "vtkPolyDataAlgorithm.h"

class vtkInflatePolyData : public vtkPolyDataAlgorithm {

public:

  static vtkInflatePolyData *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkInflatePolyData,vtkPolyDataAlgorithm);

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
