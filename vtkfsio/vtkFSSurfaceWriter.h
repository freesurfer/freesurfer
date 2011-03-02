/**
 * @file  vtkFSSurfaceWriter.h
 * @brief Writes vtkPolyData to a MRIS file
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:56 $
 *    $Revision: 1.2 $
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


