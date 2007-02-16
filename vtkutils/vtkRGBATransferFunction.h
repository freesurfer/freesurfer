/**
 * @file  vtkRGBATransferFunction.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/02/16 21:14:01 $
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


/*=========================================================================

  Program:   Visualization Toolkit
  Original Module Name: vtkColorTransferFunction

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

  Modified by Kevin Teich

=========================================================================*/
// .NAME vtkRGBATransferFunction - Defines a transfer function for mapping a property to an RGBA color value.

// .SECTION Description
// This code is based on vtkColorTransferFunction. It was modified to add an alpha element to all color output. This required adding a 5th element to the Function units, and changing all Add***Point and Add***Segment functions to Add***APoint and Add***ASegment, adding an alpha parameter to each.

#ifndef __vtkRGBATransferFunction_h
#define __vtkRGBATransferFunction_h

#include "vtkScalarsToColors.h"

class vtkPiecewiseFunction;

#define VTK_CTF_RGB           0
#define VTK_CTF_HSV           1

class VTK_FILTERING_EXPORT vtkRGBATransferFunction : public vtkScalarsToColors {
public:
  static vtkRGBATransferFunction *New();
  vtkTypeRevisionMacro(vtkRGBATransferFunction,vtkScalarsToColors);
  void DeepCopy( vtkRGBATransferFunction *f );

  // Description:
  // Print method for vtkRGBATransferFunction
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // How many points are there defining this function?
  int GetSize() {
    return this->NumberOfPoints;
  };

  // Description:
  // Add/Remove a point to/from the function defined in RGB or HSV
  // Return the index of the point (0 based), or -1 on error.
  int AddRGBAPoint( double x, double r, double g, double b, double a );
  int AddHSVAPoint( double x, double h, double s, double v, double a );
  int RemovePoint( double x );

  // Description:
  // Add two points to the function and remove all the points
  // between them
  void AddRGBASegment( double x1, double r1, double g1, double b1, double a1,
                       double x2, double r2, double g2, double b2, double a12 );
  void AddHSVASegment( double x1, double h1, double s1, double v1, double a1,
                       double x2, double h2, double s2, double v2, double a2 );

  // Description:
  // Remove all points
  void RemoveAllPoints();

  // Description:
  // Returns an RGBA color for the specified scalar value
  double *GetColor(double x) {
    return vtkScalarsToColors::GetColor(x);
  }
  void GetColor(double x, double rgba[4]);

  // Description:
  // Get the color components individually.
  double GetRedValue( double x );
  double GetGreenValue( double x );
  double GetBlueValue( double x );
  double GetAlphaValue( double x );

  // Description:
  // Map one value through the lookup table.
  virtual unsigned char *MapValue(double v);

  // Description:
  // Returns min and max position of all function points.
  vtkGetVector2Macro( Range, double );

  // Description:
  // Remove all points out of the new range, and make sure there is a point
  // at each end of that range.
  // Return 1 on success, 0 otherwise.
  int AdjustRange(double range[2]);

  // Description:
  // Fills in a table of n function values between x1 and x2
  void GetTable( double x1, double x2, int n, double* table );
  void GetTable( double x1, double x2, int n, float* table );
  const unsigned char *GetTable( double x1, double x2, int n);

  // Description:
  // Construct a color transfer function from a table. Function range is
  // is set to [x1, x2], each function size is set to size, and function
  // points are regularly spaced between x1 and x2. Parameter "table" is
  // assumed to be a block of memory of size [3*size]
  void BuildFunctionFromTable( double x1, double x2, int size, double *table);

  // Description:
  // Sets and gets the clamping value for this transfer function.
  vtkSetClampMacro( Clamping, int, 0, 1 );
  vtkGetMacro( Clamping, int );
  vtkBooleanMacro( Clamping, int );

  // Description:
  // Set/Get the color space used for interpolation: RGB, or HSV.
  // In HSV mode, if HSVWrap is on, it  will take the shortest path in Hue
  // (going back through 0 if that is the shortest way around the hue circle)
  // whereas if HSVWrap is off it will not go through 0 (in order the match
  // the current functionality of vtkLookupTable)
  vtkSetClampMacro( ColorSpace, int, VTK_CTF_RGB, VTK_CTF_HSV );
  void SetColorSpaceToRGB() {
    this->SetColorSpace(VTK_CTF_RGB);
  };
  void SetColorSpaceToHSV() {
    this->SetColorSpace(VTK_CTF_HSV);
  };
  vtkGetMacro( ColorSpace, int );
  vtkSetMacro(HSVWrap, int);
  vtkGetMacro(HSVWrap, int);
  vtkBooleanMacro(HSVWrap, int);
  VTK_LEGACY(void SetColorSpaceToHSVNoWrap());

  // Description:
  // Returns a list of all nodes
  // Fills from a pointer to data stored in a similar list of nodes.
  double *GetDataPointer() {
    return this->Function;
  };
  void FillFromDataPointer(int, double*);

  // Description:
  // map a set of scalars through the lookup table
  virtual void MapScalarsThroughTable2(void *input, unsigned char *output,
                                       int inputDataType, int numberOfValues,
                                       int inputIncrement, int outputIncrement);

protected:
  vtkRGBATransferFunction();
  ~vtkRGBATransferFunction();

  // Determines the function value outside of defined points
  // Zero = always return 0.0 outside of defined points
  // One  = clamp to the lowest value below defined points and
  //        highest value above defined points
  int Clamping;

  // The color space in which interpolation is performed
  int ColorSpace;

  // Specify if HSW is warp or not
  int HSVWrap;

  // The color function
  double     *Function;
  int         FunctionSize;
  int         NumberOfPoints;

  // An evaluated color (0 to 255 RGBA)
  unsigned char UnsignedCharRGBAValue[4];

  // The min and max point locations for all three transfer functions
  double Range[2];

  vtkTimeStamp BuildTime;
  unsigned char *Table;
  int TableSize;

  // Description:
  // Set the range of scalars being mapped. The set has no functionality
  // in this subclass of vtkScalarsToColors.
  virtual void SetRange(double, double) {};
  void SetRange(double rng[2]) {
    this->SetRange(rng[0],rng[1]);
  };


private:
  vtkRGBATransferFunction(const vtkRGBATransferFunction&);  // Not implemented.
  void operator=(const vtkRGBATransferFunction&);  // Not implemented.
};

#endif

