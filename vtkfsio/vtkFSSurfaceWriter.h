/**
 * @file  vtkFSSurfaceWriter.h
 * @brief Writes vtkPolyData to a MRIS file
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/07/30 23:38:35 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2007,
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

#ifndef vtkFSSurfaceWriter_h
#define vtkFSSurfaceWriter_h

#include "vtkWriter.h"
extern "C" {
#include "mrisurf.h"
}

class vtkPolyData;

class vtkFSSurfaceWriter : public vtkWriter {
  
 public:

  static vtkFSSurfaceWriter *New();
  vtkTypeRevisionMacro(vtkFSSurfaceWriter,vtkWriter);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description:
  // Get the input to this writer.
  vtkPolyData* GetInput();
  vtkPolyData* GetInput(int port);

  // Description:
  // Specify file name of the MRIS data file to write.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  MRIS* GetMRIS ();
  
 protected:
  vtkFSSurfaceWriter();
  ~vtkFSSurfaceWriter();
  
  void WriteData();
  
  virtual int FillInputPortInformation(int port, vtkInformation *info);
  
  MRIS* mMRIS;

  char* FileName;

 private:
  vtkFSSurfaceWriter(const vtkFSSurfaceWriter&);  // Not implemented.
  void operator=(const vtkFSSurfaceWriter&);  // Not implemented.
};

#endif


