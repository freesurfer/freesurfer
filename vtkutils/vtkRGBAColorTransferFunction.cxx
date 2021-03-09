/**
 * @brief Defines transfer function for mapping a property to an RGBA color value
 *
 * This code is based on vtkColorTransferFunction. It was modified to
 * add an alpha element to all color output. This required adding a
 * 5th element to the Function units, and changing all Add***Point and
 * Add***Segment functions to Add***APoint and Add***ASegment, adding
 * an alpha parameter to each.
 */
/*
 * Original Author: Kitware, Inc, modified by Ruopeng Wang
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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkRGBAColorTransferFunction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkRGBAColorTransferFunction.h"

#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <math.h>

vtkStandardNewMacro(vtkRGBAColorTransferFunction);

//=============================================================================
class vtkCTFNode
{
public:
  double X;
  double R;
  double G;
  double B;
  double A;
  double Sharpness;
  double Midpoint;
};

class vtkCTFCompareNodes
{
public:
  bool operator () ( const vtkCTFNode *node1,
                     const vtkCTFNode *node2 )
    {
      return node1->X < node2->X;
    }
};

class vtkCTFFindNodeEqual
{
public:
  double X;
  bool operator () ( const vtkCTFNode *node )
    {
      return node->X == this->X;
    }
};

class vtkCTFFindNodeInRange
{
public:
  double X1;
  double X2;
  bool operator () (const vtkCTFNode *node )
    {
      return ( node->X >= this->X1 &&
               node->X <= this->X2 );
    }
};

class vtkCTFFindNodeOutOfRange
{
public:
  double X1;
  double X2;
  bool operator () (const vtkCTFNode *node )
    {
      return ( node->X < this->X1 ||
               node->X > this->X2 );
    }
};
  
class vtkRGBAColorTransferFunctionInternals
{
public:
  std::vector<vtkCTFNode*> Nodes;
  vtkCTFCompareNodes          CompareNodes;
  vtkCTFFindNodeEqual         FindNodeEqual;
  vtkCTFFindNodeInRange       FindNodeInRange;
  vtkCTFFindNodeOutOfRange    FindNodeOutOfRange;
};

//=============================================================================
// Convert to and from a special polar version of CIELAB (useful for creating
// continuous diverging color maps).
inline void vtkRGBAColorTransferFunctionLabToMsh(const double lab[4], double msh[4])
{
  const double &L = lab[0];
  const double &a = lab[1];
  const double &b = lab[2];
  double &M = msh[0];
  double &s = msh[1];
  double &h = msh[2];

  M = sqrt(L*L + a*a + b*b);
  s = (M > 0.001) ? acos(L/M) : 0.0;
  h = (s > 0.001) ? atan2(b,a) : 0.0;
  msh[3] = lab[3];
}

inline void vtkRGBAColorTransferFunctionMshToLab(const double msh[4], double lab[4])
{
  const double &M = msh[0];
  const double &s = msh[1];
  const double &h = msh[2];
  double &L = lab[0];
  double &a = lab[1];
  double &b = lab[2];

  L = M*cos(s);
  a = M*sin(s)*cos(h);
  b = M*sin(s)*sin(h);
  lab[3] = msh[3];
}

// Given two angular orientations, returns the smallest angle between the two.
inline double vtkRGBAColorTransferFunctionAngleDiff(double a1, double a2)
{
  double adiff = a1 - a2;
  if (adiff < 0.0) adiff = -adiff;
  while (adiff >= 2*vtkMath::Pi()*2) adiff -= 2*vtkMath::Pi()*2;
  if (adiff > vtkMath::Pi()*2) adiff = 2*vtkMath::Pi()*2 - adiff;
  return adiff;
}

// For the case when interpolating from a saturated color to an unsaturated
// color, find a hue for the unsaturated color that makes sense.
inline double vtkRGBAColorTransferFunctionAdjustHue(const double msh[3],
                                                double unsatM)
{
  if (msh[0] >= unsatM - 0.1)
    {
    // The best we can do is hold hue constant.
    return msh[2];
    }
  else
    {
    // This equation is designed to make the perceptual change of the
    // interpolation to be close to constant.
    double hueSpin = (  msh[1]*sqrt(unsatM*unsatM - msh[0]*msh[0])
                      / (msh[0]*sin(msh[1])) );
    // Spin hue away from 0 except in purple hues.
    if (msh[2] > -0.3*vtkMath::Pi()*2)
      {
      return msh[2] + hueSpin;
      }
    else
      {
      return msh[2] - hueSpin;
      }
    }
}

// Interpolate a diverging color map.
inline void vtkRGBAColorTransferFunctionInterpolateDiverging(double s,
                                                         const double rgb1[4],
                                                         const double rgb2[4],
                                                         double result[4])
{
  double lab1[4], lab2[4];
  vtkMath::RGBToLab(rgb1, lab1);
  vtkMath::RGBToLab(rgb2, lab2);
  lab1[3] = rgb1[3];
  lab2[3] = rgb2[3];

  double msh1[4], msh2[4];
  vtkRGBAColorTransferFunctionLabToMsh(lab1, msh1);
  vtkRGBAColorTransferFunctionLabToMsh(lab2, msh2);

  // If the endpoints are distinct saturated colors, then place white in between
  // them.
  if (   (msh1[1] > 0.05) && (msh2[1] > 0.05)
      && (vtkRGBAColorTransferFunctionAngleDiff(msh1[2], msh2[2]) > 0.33*vtkMath::Pi()*2) )
    {
    // Insert the white midpoint by setting one end to white and adjusting the
    // scalar value.
    if (s < 0.5)
      {
      msh2[0] = 95.0;  msh2[1] = 0.0;  msh2[2] = 0.0;
      s = 2.0*s;
      }
    else
      {
      msh1[0] = 95.0;  msh1[1] = 0.0;  msh1[2] = 0.0;
      s = 2.0*s - 1.0;
      }
    }

  // If one color has no saturation, then its hue value is invalid.  In this
  // case, we want to set it to something logical so that the interpolation of
  // hue makes sense.
  if ((msh1[1] < 0.05) && (msh2[1] > 0.05))
    {
    msh1[2] = vtkRGBAColorTransferFunctionAdjustHue(msh2, msh1[0]);
    }
  else if ((msh2[1] < 0.05) && (msh1[1] > 0.05))
    {
    msh2[2] = vtkRGBAColorTransferFunctionAdjustHue(msh1, msh2[0]);
    }

  double mshTmp[4];
  mshTmp[0] = (1-s)*msh1[0] + s*msh2[0];
  mshTmp[1] = (1-s)*msh1[1] + s*msh2[1];
  mshTmp[2] = (1-s)*msh1[2] + s*msh2[2];  
  mshTmp[3] = (1-s)*msh1[3] + s*msh2[3];

  // Now convert back to RGB
  double labTmp[4];
  vtkRGBAColorTransferFunctionMshToLab(mshTmp, labTmp);
  vtkMath::LabToRGB(labTmp, result);
  result[3] = labTmp[3];
}

//----------------------------------------------------------------------------
// Construct a new vtkRGBAColorTransferFunction with default values
vtkRGBAColorTransferFunction::vtkRGBAColorTransferFunction()
{
  this->UnsignedCharRGBAValue[0] = 0;
  this->UnsignedCharRGBAValue[1] = 0;
  this->UnsignedCharRGBAValue[2] = 0;
  this->UnsignedCharRGBAValue[3] = 0;
  
  this->Range[0] = 0;
  this->Range[1] = 0;

  this->Clamping = 1;
  this->ColorSpace = VTK_CTF_RGB;
  this->HSVWrap = 1; //By default HSV will be wrap

  this->Scale = VTK_CTF_LINEAR;
  
  this->Function = NULL;

  this->Table = NULL;
  this->TableSize = 0;

  this->AllowDuplicateScalars = 1;

  this->Internal = new vtkRGBAColorTransferFunctionInternals;
}

//----------------------------------------------------------------------------
// Destruct a vtkRGBAColorTransferFunction
vtkRGBAColorTransferFunction::~vtkRGBAColorTransferFunction()
{
  delete [] this->Table;
  
  if ( this->Function )
    {
    delete [] this->Function;
    this->Function = NULL;
    }

  for(unsigned int i=0;i<this->Internal->Nodes.size();i++)
    {
    delete this->Internal->Nodes[i];
    }
  this->Internal->Nodes.clear(); 
  delete this->Internal;
}

// Return the number of points which specify this function
int vtkRGBAColorTransferFunction::GetSize()
{
  return static_cast<int>(this->Internal->Nodes.size());
}

// Since we no longer store the data in an array, we must
// copy out of the vector into an array. No modified check - 
// could be added if performance is a problem
double *vtkRGBAColorTransferFunction::GetDataPointer()
{
  int size = static_cast<int>(this->Internal->Nodes.size());
  
  if ( this->Function )
    {
    delete [] this->Function;
    this->Function = NULL;
    }
  
  if ( size > 0 )
    {
    this->Function = new double[size*5];
    for ( int i = 0; i < size; i++ )
      {
      this->Function[5*i  ] = this->Internal->Nodes[i]->X;
      this->Function[5*i+1] = this->Internal->Nodes[i]->R;
      this->Function[5*i+2] = this->Internal->Nodes[i]->G;
      this->Function[5*i+3] = this->Internal->Nodes[i]->B;
      this->Function[5*i+4] = this->Internal->Nodes[i]->A;
      }
    }
  
  return this->Function;
}

//----------------------------------------------------------------------------
// Add a point defined in RGB
int vtkRGBAColorTransferFunction::AddRGBAPoint( double x, double r,
                                           double g, double b, double a )
{
  return this->AddRGBAPoint( x, r, g, b, a, 0.5, 0.0 );
}

//----------------------------------------------------------------------------
// Add a point defined in RGB
int vtkRGBAColorTransferFunction::AddRGBAPoint( double x, double r,
                                           double g, double b, double a,
                                           double midpoint,
                                           double sharpness )
{
  // Error check
  if ( midpoint < 0.0 || midpoint > 1.0 )
    {
    vtkErrorMacro("Midpoint outside range [0.0, 1.0]");
    return -1;
    }

  if ( sharpness < 0.0 || sharpness > 1.0 )
    {
    vtkErrorMacro("Sharpness outside range [0.0, 1.0]");
    return -1;
    }
  
  // remove any node already at this X location
  if (!this->AllowDuplicateScalars)
    {
    this->RemovePoint( x );
    }
  
  // Create the new node
  vtkCTFNode *node = new vtkCTFNode;
  node->X         = x;
  node->R         = r;
  node->G         = g;
  node->B         = b;
  node->A         = a;
  node->Midpoint  = midpoint;
  node->Sharpness = sharpness;
  
  // Add it, then sort to get everyting in order
  this->Internal->Nodes.push_back(node);
  this->SortAndUpdateRange();
  
  // We need to find the index of the node we just added in order
  // to return this value
  unsigned int i;
  for ( i = 0; i < this->Internal->Nodes.size(); i++ )
    {
    if ( this->Internal->Nodes[i]->X == x )
      {
      break;
      }
    }
  
  int retVal;
  
  // If we didn't find it, something went horribly wrong so
  // return -1  
  if ( i < this->Internal->Nodes.size() )
    {
    retVal = i;
    }
  else
    {
    retVal = -1;
    }
  
  return retVal;  
}

//----------------------------------------------------------------------------
// Add a point defined in HSV
int vtkRGBAColorTransferFunction::AddHSVAPoint( double x, double h,
                                            double s, double v, double a )
{ 
  double r, b, g;
  
  vtkMath::HSVToRGB(h, s, v, &r, &g, &b);
  return this->AddRGBAPoint( x, r, g, b, a );
}

//----------------------------------------------------------------------------
// Add a point defined in HSV
int vtkRGBAColorTransferFunction::AddHSVAPoint( double x, double h,
                                           double s, double v, double a,
                                           double midpoint,
                                           double sharpness )
{ 
  double r, b, g;
  
  vtkMath::HSVToRGB(h, s, v, &r, &g, &b);
  return this->AddRGBAPoint( x, r, g, b, a, midpoint, sharpness );
}

// Sort the vector in increasing order, then fill in
// the Range
void vtkRGBAColorTransferFunction::SortAndUpdateRange()
{
  std::sort( this->Internal->Nodes.begin(),
                this->Internal->Nodes.end(),
                this->Internal->CompareNodes );
  
  int size = static_cast<int>(this->Internal->Nodes.size());
  if ( size )
    {
    this->Range[0] = this->Internal->Nodes[0]->X;
    this->Range[1] = this->Internal->Nodes[size-1]->X;
    }
  else
    {
    this->Range[0] = 0;
    this->Range[1] = 0;
    }
  
  this->Modified();  
}

//----------------------------------------------------------------------------
// Remove a point
int vtkRGBAColorTransferFunction::RemovePoint( double x )
{
  // First find the node since we need to know its
  // index as our return value
  unsigned int i;
  for ( i = 0; i < this->Internal->Nodes.size(); i++ )
    {
    if ( this->Internal->Nodes[i]->X == x )
      {
      break;
      }
    }
  
  int retVal;
  
  // If the node doesn't exist, we return -1
  if ( i < this->Internal->Nodes.size() )
    {
    retVal = i;
    }
  else
    {
    return -1;
    }
  
  // Now use STL to find it, so that we can remove it
  this->Internal->FindNodeEqual.X = x;
  
  std::vector<vtkCTFNode*>::iterator iter = 
    std::find_if(this->Internal->Nodes.begin(),
                    this->Internal->Nodes.end(),
                    this->Internal->FindNodeEqual );
  
  // Actually delete it
  if ( iter != this->Internal->Nodes.end() )
    {
    delete *iter;
    this->Internal->Nodes.erase(iter);
    this->Modified();
    }
  else
     {
     // This should never happen - we already returned if the node
     // didn't exist...
     return -1;
     }
   
  
  return retVal;
}


//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::MovePoint(double oldX, double newX)
{
  if (oldX == newX)
    {
    // Nothing to do.
    return;
    }

  this->RemovePoint(newX);
  for (unsigned int i = 0; i < this->Internal->Nodes.size(); i++ )
    {
    if ( this->Internal->Nodes[i]->X == oldX )
      {
      this->Internal->Nodes[i]->X = newX;
      this->SortAndUpdateRange();
      break;
      }
    }
}

//----------------------------------------------------------------------------
// Remove all points
void vtkRGBAColorTransferFunction::RemoveAllPoints()
{
  for(unsigned int i=0;i<this->Internal->Nodes.size();i++)
    {
    delete this->Internal->Nodes[i];
    }
  this->Internal->Nodes.clear(); 

  this->SortAndUpdateRange();
}

//----------------------------------------------------------------------------
// Add a line defined in RGB 
void vtkRGBAColorTransferFunction::AddRGBASegment( double x1, double r1, 
                                              double g1, double b1, double a1,
                                              double x2, double r2, 
                                              double g2, double b2, double a2 )
{
  int done;  
  
  // First, find all points in this range and remove them
  done = 0;
  while ( !done )
    {
    done = 1;
    
    this->Internal->FindNodeInRange.X1 = x1;
    this->Internal->FindNodeInRange.X2 = x2;
  
    std::vector<vtkCTFNode*>::iterator iter = 
      std::find_if(this->Internal->Nodes.begin(),
                      this->Internal->Nodes.end(),
                      this->Internal->FindNodeInRange );
  
    if ( iter != this->Internal->Nodes.end() )
      {
      delete *iter;
      this->Internal->Nodes.erase(iter);
      this->Modified();
      done = 0;
      }
    }

  // Now add the points
  this->AddRGBAPoint( x1, r1, g1, b1, a1, 0.5, 0.0 );
  this->AddRGBAPoint( x2, r2, g2, b2, a2, 0.5, 0.0 );
}

//----------------------------------------------------------------------------
// Add a line defined in HSV
void vtkRGBAColorTransferFunction::AddHSVASegment( double x1, double h1, 
                                              double s1, double v1, double a1,
                                              double x2, double h2, 
                                              double s2, double v2, double a2 )
{
  double r1, r2, b1, b2, g1, g2;
  
  vtkMath::HSVToRGB(h1, s1, v1, &r1, &g1, &b1);
  vtkMath::HSVToRGB(h2, s2, v2, &r2, &g2, &b2);
  this->AddRGBASegment( x1, r1, g1, b1, a1, x2, r2, g2, b2, a2 );
}

//----------------------------------------------------------------------------
// Returns the RGBA color evaluated at the specified location
#if VTK_MAJOR_VERSION >= 8 && VTK_MINOR_VERSION >= 2
const unsigned char *vtkRGBAColorTransferFunction::MapValue( double x )
#else
unsigned char *vtkRGBAColorTransferFunction::MapValue( double x )
#endif
{
  double rgba[4];
  this->GetColor( x, rgba );
  
  this->UnsignedCharRGBAValue[0] =
    static_cast<unsigned char>(255.0*rgba[0] + 0.5);
  this->UnsignedCharRGBAValue[1] =
    static_cast<unsigned char>(255.0*rgba[1] + 0.5);
  this->UnsignedCharRGBAValue[2] =
    static_cast<unsigned char>(255.0*rgba[2] + 0.5);
  this->UnsignedCharRGBAValue[3] = 
    static_cast<unsigned char>(255.0*rgba[3] + 0.5);
  return this->UnsignedCharRGBAValue;
}

//----------------------------------------------------------------------------
// Returns the RGBA color evaluated at the specified location
void vtkRGBAColorTransferFunction::GetColor(double x, double rgba[4])
{
  this->GetTable( x, x, 1, rgba );
}

//----------------------------------------------------------------------------
// Returns the red color evaluated at the specified location
double vtkRGBAColorTransferFunction::GetRedValue( double x )
{
  double rgba[4];
  this->GetColor( x, rgba );

  return rgba[0];
}

//----------------------------------------------------------------------------
// Returns the green color evaluated at the specified location
double vtkRGBAColorTransferFunction::GetGreenValue( double x )
{
  double rgba[4];
  this->GetColor( x, rgba );

  return rgba[1];
}

//----------------------------------------------------------------------------
// Returns the blue color evaluated at the specified location
double vtkRGBAColorTransferFunction::GetBlueValue( double x )
{
  double rgba[4];
  this->GetColor( x, rgba );

  return rgba[2];
}

//----------------------------------------------------------------------------
// Returns the alpha value evaluated at the specified location
double vtkRGBAColorTransferFunction::GetAlphaValue( double x )
{
  double rgba[4];
  this->GetColor( x, rgba );

  return rgba[3];
}
//----------------------------------------------------------------------------
// Returns a table of RGB colors at regular intervals along the function
void vtkRGBAColorTransferFunction::GetTable( double xStart, double xEnd, 
                                         int size, double* table )
{
  int i, j;
  int idx = 0;
  int numNodes = static_cast<int>(this->Internal->Nodes.size());
  
  // Need to keep track of the last value so that
  // we can fill in table locations past this with
  // this value if Clamping is On.
  double lastR = 0.0;
  double lastG = 0.0;
  double lastB = 0.0;
  double lastA = 0.0;
  if ( numNodes != 0 )
    {
    lastR = this->Internal->Nodes[numNodes-1]->R;
    lastG = this->Internal->Nodes[numNodes-1]->G;
    lastB = this->Internal->Nodes[numNodes-1]->B;
    lastA = this->Internal->Nodes[numNodes-1]->A;
    }
  
  double *tptr     = NULL;
  double x         = 0.0;  
  double x1        = 0.0;
  double x2        = 0.0;
  double rgba1[4]  = {0.0, 0.0, 0.0, 0.0};
  double rgba2[4]  = {0.0, 0.0, 0.0, 0.0};  
  double midpoint  = 0.0;
  double sharpness = 0.0;

  // If the scale is logarithmic, make sure the range is valid.
  bool usingLogScale = this->Scale == VTK_CTF_LOG10;
  if(usingLogScale)
    {
    // Note: This requires range[0] <= range[1].
    usingLogScale = this->Range[0] > 0.0;
    }

  double logStart = 0.0;
  double logEnd   = 0.0;
  double logX     = 0.0;
  if(usingLogScale)
    {
    logStart = log10(xStart);
    logEnd = log10(xEnd);
    }
    
  // For each table entry
  for ( i = 0; i < size; i++ )
    {
    
    // Find our location in the table
    tptr = table + 4*i;
    
    // Find our X location. If we are taking only 1 sample, make
    // it halfway between start and end (usually start and end will
    // be the same in this case)
    if ( size > 1 )
      {
      if(usingLogScale)
        {
        logX = logStart +
          (static_cast<double>(i)/static_cast<double>(size-1))
          *(logEnd-logStart);
        x = pow(static_cast<double>(10.0), logX);
        }
      else
        {
        x = xStart + (static_cast<double>(i)/static_cast<double>(size-1))
          *(xEnd-xStart);
        }
      }
    else
      {
      if(usingLogScale)
        {
        logX = 0.5*(logStart+logEnd);
        x = pow(static_cast<double>(10.0), logX);
        }
      else
        {
        x = 0.5*(xStart+xEnd);
        }
      }
    
    // Do we need to move to the next node?
    while ( idx < numNodes &&
            x > this->Internal->Nodes[idx]->X )
      {
      idx++;
      // If we are at a valid point index, fill in
      // the value at this node, and the one before (the
      // two that surround our current sample location)
      // idx cannot be 0 since we just incremented it.
      if ( idx < numNodes )
        {
        x1 = this->Internal->Nodes[idx-1]->X;
        x2 = this->Internal->Nodes[idx  ]->X;
        if(usingLogScale)
          {
          x1 = log10(x1);
          x2 = log10(x2);
          }
        
        rgba1[0] = this->Internal->Nodes[idx-1]->R;
        rgba2[0] = this->Internal->Nodes[idx  ]->R;

        rgba1[1] = this->Internal->Nodes[idx-1]->G;
        rgba2[1] = this->Internal->Nodes[idx  ]->G;
        
        rgba1[2] = this->Internal->Nodes[idx-1]->B;
        rgba2[2] = this->Internal->Nodes[idx  ]->B;
        
        rgba1[3] = this->Internal->Nodes[idx-1]->A;
        rgba2[3] = this->Internal->Nodes[idx  ]->A;
        
        // We only need the previous midpoint and sharpness
        // since these control this region
        midpoint  = this->Internal->Nodes[idx-1]->Midpoint;
        sharpness = this->Internal->Nodes[idx-1]->Sharpness;
        
        // Move midpoint away from extreme ends of range to avoid
        // degenerate math
        if ( midpoint < 0.00001 )
          {
          midpoint = 0.00001;
          }
        
        if ( midpoint > 0.99999 )
          {
          midpoint = 0.99999;
          }
        }
      }
    
    // Are we at the end? If so, just use the last value
    if ( idx >= numNodes )
      {
      tptr[0] = (this->Clamping)?(lastR):(0.0);
      tptr[1] = (this->Clamping)?(lastG):(0.0);
      tptr[2] = (this->Clamping)?(lastB):(0.0);
      tptr[3] = (this->Clamping)?(lastA):(0.0);
      }
    // Are we before the first node? If so, duplicate this nodes values
    else if ( idx == 0 )
      {
      tptr[0] = (this->Clamping)?(this->Internal->Nodes[0]->R):(0.0);
      tptr[1] = (this->Clamping)?(this->Internal->Nodes[0]->G):(0.0);
      tptr[2] = (this->Clamping)?(this->Internal->Nodes[0]->B):(0.0);
      tptr[3] = (this->Clamping)?(this->Internal->Nodes[0]->A):(0.0);
      }
    // Otherwise, we are between two nodes - interpolate
    else
      {
      // Our first attempt at a normalized location [0,1] - 
      // we will be modifying this based on midpoint and 
      // sharpness to get the curve shape we want and to have
      // it pass through (y1+y2)/2 at the midpoint.
      double s = 0.0;
      if(usingLogScale)
        {
        s = (logX - x1) / (x2 - x1);
        }
      else
        {
        s = (x - x1) / (x2 - x1);
        }
      
      // Readjust based on the midpoint - linear adjustment
      if ( s < midpoint )
        {
        s = 0.5 * s / midpoint;
        }
      else
        {
        s = 0.5 + 0.5*(s-midpoint)/(1.0-midpoint);
        }
      
      // override for sharpness > 0.99 
      // In this case we just want piecewise constant
      if ( sharpness > 0.99 )
        {
        // Use the first value since we are below the midpoint
        if ( s < 0.5 )
          {
          tptr[0] = rgba1[0];
          tptr[1] = rgba1[1];
          tptr[2] = rgba1[2];
          tptr[3] = rgba1[3];
          continue;
          }
        // Use the second value at or above the midpoint
        else
          {
          tptr[0] = rgba2[0];
          tptr[1] = rgba2[1];
          tptr[2] = rgba2[2];
          tptr[3] = rgba2[3];
          continue;
          }
        }
      
      // Override for sharpness < 0.01
      // In this case we want piecewise linear
      if ( sharpness < 0.01 )
        {
        // Simple linear interpolation 
        if ( this->ColorSpace == VTK_CTF_RGB )
          {
          tptr[0] = (1-s)*rgba1[0] + s*rgba2[0];
          tptr[1] = (1-s)*rgba1[1] + s*rgba2[1];
          tptr[2] = (1-s)*rgba1[2] + s*rgba2[2];
          tptr[3] = (1-s)*rgba1[3] + s*rgba2[3];
          }
        else if ( this->ColorSpace == VTK_CTF_HSV )
          {
          double hsv1[3], hsv2[3];
          vtkMath::RGBToHSV(rgba1, hsv1);
          vtkMath::RGBToHSV(rgba2, hsv2);
          
          if ( this->HSVWrap &&
               (hsv1[0] - hsv2[0] > 0.5 ||
                hsv2[0] - hsv1[0] > 0.5) )
            {
            if ( hsv1[0] > hsv2[0] )
              {
              hsv1[0] -= 1.0;
              }
            else
              {
              hsv2[0] -= 1.0;
              }
            }
          
          double hsvTmp[3];
          hsvTmp[0] = (1-s)*hsv1[0] + s*hsv2[0];
          if ( hsvTmp[0] < 0.0 )
            {
            hsvTmp[0] += 1.0;
            }
          hsvTmp[1] = (1-s)*hsv1[1] + s*hsv2[1];
          hsvTmp[2] = (1-s)*hsv1[2] + s*hsv2[2];   
          
          // Now convert this back to RGBA
          vtkMath::HSVToRGB( hsvTmp, tptr );
          tptr[3] = (1-s)*rgba1[3] + s*rgba2[3];
          }
        else if (this->ColorSpace == VTK_CTF_LAB)
          {
          double lab1[3], lab2[3];
          vtkMath::RGBToLab(rgba1, lab1);
          vtkMath::RGBToLab(rgba2, lab2);

          double labTmp[3];
          labTmp[0] = (1-s)*lab1[0] + s*lab2[0];
          labTmp[1] = (1-s)*lab1[1] + s*lab2[1];
          labTmp[2] = (1-s)*lab1[2] + s*lab2[2];

          // Now convert back to RGBA
          vtkMath::LabToRGB(labTmp, tptr);
          tptr[3] = (1-s)*rgba1[3] + s*rgba2[3];
          }
        else if (this->ColorSpace == VTK_CTF_DIVERGING)
          {
          vtkRGBAColorTransferFunctionInterpolateDiverging(s, rgba1, rgba2, tptr);
          }
        else
          {
          vtkErrorMacro("ColorSpace set to invalid value.");
          }
        continue;
        }    
      
      // We have a sharpness between [0.01, 0.99] - we will
      // used a modified hermite curve interpolation where we
      // derive the slope based on the sharpness, and we compress
      // the curve non-linearly based on the sharpness
      
      // First, we will adjust our position based on sharpness in 
      // order to make the curve sharper (closer to piecewise constant)
      if ( s < .5 )
        {
        s = 0.5 * pow(s*2,1.0 + 10*sharpness);
        }
      else if ( s > .5 )
        {
        s = 1.0 - 0.5 * pow((1.0-s)*2,1+10*sharpness);
        }
      
      // Compute some coefficients we will need for the hermite curve
      double ss = s*s;
      double sss = ss*s;
      
      double h1 =  2*sss - 3*ss + 1;
      double h2 = -2*sss + 3*ss;
      double h3 =    sss - 2*ss + s;
      double h4 =    sss -   ss;
      
      double slope;
      double t;

      if ( this->ColorSpace == VTK_CTF_RGB )
        {
        for ( j = 0; j < 4; j++ )
          {
          // Use one slope for both end points
          slope = rgba2[j] - rgba1[j];
          t = (1.0 - sharpness)*slope;
          
          // Compute the value
          tptr[j] = h1*rgba1[j] + h2*rgba2[j] + h3*t + h4*t;
          }
        }
      else if (this->ColorSpace == VTK_CTF_HSV)
        {
        double hsv1[4], hsv2[4];
        vtkMath::RGBToHSV(rgba1, hsv1);
        vtkMath::RGBToHSV(rgba2, hsv2);
        hsv1[3] = rgba1[3];
        hsv2[3] = rgba2[3];
        
        if ( this->HSVWrap &&
             (hsv1[0] - hsv2[0] > 0.5 ||
              hsv2[0] - hsv1[0] > 0.5) )
          {
          if ( hsv1[0] > hsv2[0] )
            {
            hsv1[0] -= 1.0;
            }
          else
            {
            hsv2[0] -= 1.0;
            }
          }
        
        double hsvTmp[4];
        
        for ( j = 0; j < 4; j++ )
          {
          // Use one slope for both end points
          slope = hsv2[j] - hsv1[j];
          t = (1.0 - sharpness)*slope;
          
          // Compute the value
          hsvTmp[j] = h1*hsv1[j] + h2*hsv2[j] + h3*t + h4*t;
          if ( j == 0 && hsvTmp[j] < 0.0 )
            {
            hsvTmp[j] += 1.0;
            }
          }
        // Now convert this back to RGB
        vtkMath::HSVToRGB( hsvTmp, tptr );
        tptr[3] = hsvTmp[3];
        }
      else if (this->ColorSpace == VTK_CTF_LAB)
        {
        double lab1[4], lab2[4];
        vtkMath::RGBToLab(rgba1, lab1);
        vtkMath::RGBToLab(rgba2, lab2);
        lab1[3] = rgba1[3];
        lab2[3] = rgba2[3];

        double labTmp[4];
        for (j = 0; j < 4; j++)
          {
          // Use one slope for both end points
          slope = lab2[j] - lab1[j];
          t = (1.0 - sharpness)*slope;
          
          // Compute the value
          labTmp[j] = h1*lab1[j] + h2*lab2[j] + h3*t + h4*t;
          }
        // Now convert this back to RGB
        vtkMath::LabToRGB(labTmp, tptr);
        tptr[3] = labTmp[3];
        }
      else if (this->ColorSpace == VTK_CTF_DIVERGING)
        {
        // I have not implemented proper interpolation by a hermite curve for
        // the diverging color map, but I cannot think of a good use case for
        // that anyway.
        vtkRGBAColorTransferFunctionInterpolateDiverging(s, rgba1, rgba2, tptr);
        }
      else
        {
        vtkErrorMacro("ColorSpace set to invalid value.");
        }
      
      // Final error check to make sure we don't go outside [0,1]
      for ( j = 0; j < 4; j++ )
        {
        tptr[j] = (tptr[j] < 0.0)?(0.0):(tptr[j]);
        tptr[j] = (tptr[j] > 1.0)?(1.0):(tptr[j]);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::GetTable( double xStart, double xEnd, 
                                         int size, float* table )
{
  double *tmpTable = new double [size*4];
  
  this->GetTable( xStart, xEnd, size, tmpTable );
  
  double *tmpPtr = tmpTable;
  float *tPtr = table;
  
  for ( int i = 0; i < size*4; i++ )
    {
    *tPtr = static_cast<float>(*tmpPtr);
    tPtr   ++;
    tmpPtr ++;
    }

  delete[] tmpTable;
}

//----------------------------------------------------------------------------
const unsigned char *vtkRGBAColorTransferFunction::GetTable( double xStart, double xEnd, 
                                                         int size)
{
  if (this->GetMTime() <= this->BuildTime &&
      this->TableSize == size)
    {
    return this->Table;
    }

  if ( this->Internal->Nodes.size() == 0 )
    {
    vtkErrorMacro( 
      "Attempting to lookup a value with no points in the function");
    return this->Table;
    }
  
  if (this->TableSize != size)
    {
    delete [] this->Table;
    this->Table = new unsigned char [size*4];
    this->TableSize = size;
    }
  
  double *tmpTable = new double [size*4];
  
  this->GetTable( xStart, xEnd, size, tmpTable );
  
  double *tmpPtr = tmpTable;
  unsigned char *tPtr = this->Table;
  
  for ( int i = 0; i < size*4; i++ )
    {
    *tPtr = static_cast<unsigned char>(*tmpPtr*255.0 + 0.5);
    tPtr   ++;
    tmpPtr ++;
    }

  delete[] tmpTable;

  this->BuildTime.Modified();

  return this->Table;
}

//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::BuildFunctionFromTable(double xStart,
                                                      double xEnd,
                                                      int size,
                                                      double *table)
{
  double inc = 0.0;
  double *tptr = table;

  this->RemoveAllPoints();
  
  if( size > 1 )
    {
    inc = (xEnd-xStart)/static_cast<double>(size-1);
    }
  
  int i;
  for (i=0; i < size; i++)
    {
    vtkCTFNode *node = new vtkCTFNode;
    node->X   = xStart + inc*i;
    node->R   = tptr[0];
    node->G   = tptr[1];
    node->B   = tptr[2];
    node->A   = tptr[3];
    node->Sharpness = 0.0;
    node->Midpoint  = 0.5;
  
    this->Internal->Nodes.push_back(node);
    tptr += 4;
    }
  
  this->SortAndUpdateRange();
}

//----------------------------------------------------------------------------
// For a specified index value, get the node parameters
int vtkRGBAColorTransferFunction::GetNodeValue( int index, double val[7] )
{
  int size = static_cast<int>(this->Internal->Nodes.size());
  
  if ( index < 0 || index >= size )
    {
    vtkErrorMacro("Index out of range!");
    val[0] = 0;	// initialize in case caller uses - detected as error by gcc4.5 on
    val[1] = 0;
    val[2] = 0;
    val[3] = 0;
    val[4] = 0;
    val[5] = 0;
    val[6] = 0;
    return -1;
    }
  
  val[0] = this->Internal->Nodes[index]->X;
  val[1] = this->Internal->Nodes[index]->R;
  val[2] = this->Internal->Nodes[index]->G;
  val[3] = this->Internal->Nodes[index]->B;
  val[4] = this->Internal->Nodes[index]->A;
  val[5] = this->Internal->Nodes[index]->Midpoint;
  val[6] = this->Internal->Nodes[index]->Sharpness;
  
  return 1;
}

//----------------------------------------------------------------------------
// For a specified index value, get the node parameters
int vtkRGBAColorTransferFunction::SetNodeValue( int index, double val[7] )
{
  int size = static_cast<int>(this->Internal->Nodes.size());

  if ( index < 0 || index >= size )
    {
    vtkErrorMacro("Index out of range!");
    return -1;
    }
  
  this->Internal->Nodes[index]->X = val[0];
  this->Internal->Nodes[index]->R = val[1];
  this->Internal->Nodes[index]->G = val[2];
  this->Internal->Nodes[index]->B = val[3];
  this->Internal->Nodes[index]->A = val[4];
  this->Internal->Nodes[index]->Midpoint = val[5];
  this->Internal->Nodes[index]->Sharpness = val[6];

  this->Modified();

  return 1;
}

//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::DeepCopy( vtkRGBAColorTransferFunction *f )
{
  if (f != NULL)
    {
    this->Clamping     = f->Clamping;
    this->ColorSpace   = f->ColorSpace;
    this->HSVWrap      = f->HSVWrap;
    this->Scale        = f->Scale;

    int i;
    this->RemoveAllPoints();
    for ( i = 0; i < f->GetSize(); i++ )
      {
      double val[7];
      f->GetNodeValue(i, val);
      this->AddRGBAPoint(val[0], val[1], val[2], val[3], val[4], val[5], val[6]);
      }
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::ShallowCopy( vtkRGBAColorTransferFunction *f )
{
  if (f != NULL)
    {
    this->Clamping     = f->Clamping;
    this->ColorSpace   = f->ColorSpace;
    this->HSVWrap      = f->HSVWrap;
    this->Scale        = f->Scale;

    int i;
    this->RemoveAllPoints();
    for ( i = 0; i < f->GetSize(); i++ )
      {
      double val[7];
      f->GetNodeValue(i, val);
      this->AddRGBAPoint(val[0], val[1], val[2], val[3], val[4], val[5], val[6]);
      }
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// Accelerate the mapping by copying the data in 32-bit chunks instead
// of 8-bit chunks.  The extra "long" argument is to help broken
// compilers select the non-templates below for unsigned char
// and unsigned short.
template <class T>
void vtkRGBAColorTransferFunctionMapData(vtkRGBAColorTransferFunction* self,
                                     T* input,
                                     unsigned char* output,
                                     int length, int inIncr,
                                     int outFormat, long)
{
  double          x;
  int            i = length;
  double          rgb[4];
  unsigned char  *optr = output;
  T              *iptr = input;
//  unsigned char   alpha = static_cast<unsigned char>(self->GetAlphaValue()*255.0);
  
  if(self->GetSize() == 0)
    {
    vtkGenericWarningMacro("Transfer Function Has No Points!");
    return;
    }

  while (--i >= 0) 
    {
    x = static_cast<double>(*iptr);
    self->GetColor(x, rgb);
    
    if (outFormat == VTK_RGB || outFormat == VTK_RGBA)
      {
      *(optr++) = static_cast<unsigned char>(rgb[0]*255.0 + 0.5);
      *(optr++) = static_cast<unsigned char>(rgb[1]*255.0 + 0.5);
      *(optr++) = static_cast<unsigned char>(rgb[2]*255.0 + 0.5);
      }
    else // LUMINANCE  use coeffs of (0.30  0.59  0.11)*255.0
      {
      *(optr++) = static_cast<unsigned char>(rgb[0]*76.5 + rgb[1]*150.45 +
                                             rgb[2]*28.05 + 0.5); 
      }
    
    if (outFormat == VTK_RGBA || outFormat == VTK_LUMINANCE_ALPHA)
      {
        *(optr++) = static_cast<unsigned char>(rgb[3]*255.0 + 0.5);
      }
    iptr += inIncr;
    }
}



//----------------------------------------------------------------------------
// Special implementation for unsigned char input.
void vtkRGBAColorTransferFunctionMapData(vtkRGBAColorTransferFunction* self,
                                     unsigned char* input,
                                     unsigned char* output,
                                     int length, int inIncr,
                                     int outFormat, int)
{
  int            x;
  int            i = length;
  unsigned char  *optr = output;
  unsigned char  *iptr = input;

  if(self->GetSize() == 0)
    {
    vtkGenericWarningMacro("Transfer Function Has No Points!");
    return;
    }
  
  const unsigned char *table = self->GetTable(0,255,256);
  switch (outFormat)
    {
    case VTK_RGB:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        *(optr++) = table[x+1];
        *(optr++) = table[x+2];
        iptr += inIncr;
        }
      break;
    case VTK_RGBA:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        *(optr++) = table[x+1];
        *(optr++) = table[x+2];
        *(optr++) = table[x+3];
        iptr += inIncr;
        }
      break;
    case VTK_LUMINANCE_ALPHA:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        *(optr++) = table[x+3];
        iptr += inIncr;
        }
      break;
    case VTK_LUMINANCE:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        iptr += inIncr;
        }
      break;
    }  
}

//----------------------------------------------------------------------------
// Special implementation for unsigned short input.
void vtkRGBAColorTransferFunctionMapData(vtkRGBAColorTransferFunction* self,
                                     unsigned short* input,
                                     unsigned char* output,
                                     int length, int inIncr,
                                     int outFormat, int)
{
  int            x;
  int            i = length;
  unsigned char  *optr = output;
  unsigned short *iptr = input;

  if(self->GetSize() == 0)
    {
    vtkGenericWarningMacro("Transfer Function Has No Points!");
    return;
    }
  

  const unsigned char *table = self->GetTable(0,65535,65536);
  switch (outFormat)
    {
    case VTK_RGB:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        *(optr++) = table[x+1];
        *(optr++) = table[x+2];
        iptr += inIncr;
        }
      break;
    case VTK_RGBA:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        *(optr++) = table[x+1];
        *(optr++) = table[x+2];
        *(optr++) = table[x+3];
        iptr += inIncr;
        }
      break;
    case VTK_LUMINANCE_ALPHA:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        *(optr++) = table[x+3];
        iptr += inIncr;
        }
      break;
    case VTK_LUMINANCE:
      while (--i >= 0) 
        {
        x = *iptr*4;
        *(optr++) = table[x];
        iptr += inIncr;
        }
      break;
    }  
}

//----------------------------------------------------------------------------
// An expensive magnitude calculation (similar to
// vtkLookupTable.cxx (vtkLookupTableMapMag).
template <class T>
void vtkRGBAColorTransferFunctionMagMapData(vtkRGBAColorTransferFunction* self,
                                     T* input,
                                     unsigned char* output,
                                     int length, int inIncr,
                                     int outFormat, int v)
{
  double tmp, sum;
  double *mag;
  int i, j;

  mag = new double[length];
  for (i = 0; i < length; ++i)
    {
    sum = 0;
    for (j = 0; j < inIncr; ++j)
      {
      tmp = static_cast<double>(*input);  
      sum += (tmp * tmp);
      ++input;
      }
    mag[i] = sqrt(sum);
    }

  vtkRGBAColorTransferFunctionMapData(self, mag, output, length, 1, outFormat, v);

  delete [] mag;
}

//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::MapScalarsThroughTable2(void *input, 
                                                       unsigned char *output,
                                                       int inputDataType, 
                                                       int numberOfValues,
                                                       int inputIncrement,
                                                       int outputFormat)
{
  if (this->UseMagnitude && inputIncrement > 1)
    {
    switch (inputDataType)
      {
      vtkTemplateMacro(
        vtkRGBAColorTransferFunctionMagMapData(this, static_cast<VTK_TT*>(input),
          output, numberOfValues, inputIncrement, outputFormat, 1);
        return
      );
    case VTK_BIT:
      vtkErrorMacro("Cannot compute magnitude of bit array.");
      break;
    default:
      vtkErrorMacro(<< "MapImageThroughTable: Unknown input ScalarType");
      }
    }

  switch (inputDataType)
    {
    vtkTemplateMacro(
      vtkRGBAColorTransferFunctionMapData(this, static_cast<VTK_TT*>(input),
                                      output, numberOfValues, inputIncrement,
                                      outputFormat, 1)
      );
    default:
      vtkErrorMacro(<< "MapImageThroughTable: Unknown input ScalarType");
      return;
    }
}

//----------------------------------------------------------------------------
void vtkRGBAColorTransferFunction::FillFromDataPointer(int nb, double *ptr)
{
  if (nb <= 0 || !ptr)
    {
    return;
    }

  this->RemoveAllPoints();

  while (nb)
    {
    this->AddRGBAPoint(ptr[0], ptr[1], ptr[2], ptr[3], ptr[4]);
    ptr += 5;
    nb--;
    }
}

//----------------------------------------------------------------------------
int vtkRGBAColorTransferFunction::AdjustRange(double range[2])
{
  if (!range)
    {
    return 0;
    }

  double *function_range = this->GetRange();
  
  // Make sure we have points at each end of the range

  double rgb[4];
  if (function_range[0] < range[0])
    {
    this->GetColor(range[0], rgb);
    this->AddRGBAPoint(range[0], rgb[0], rgb[1], rgb[2], rgb[3]);
    }
  else
    {
    this->GetColor(function_range[0], rgb);
    this->AddRGBAPoint(range[0], rgb[0], rgb[1], rgb[2], rgb[3]);
    }

  if (function_range[1] > range[1])
    {
    this->GetColor(range[1], rgb);
    this->AddRGBAPoint(range[1], rgb[0], rgb[1], rgb[2], rgb[3]);
    }
  else
    {
    this->GetColor(function_range[1], rgb);
    this->AddRGBAPoint(range[1], rgb[0], rgb[1], rgb[2], rgb[3]);
    }

  // Remove all points out-of-range
  int done;  
  
  done = 0;
  while ( !done )
    {
    done = 1;
    
    this->Internal->FindNodeOutOfRange.X1 = range[0];
    this->Internal->FindNodeOutOfRange.X2 = range[1];
  
    std::vector<vtkCTFNode*>::iterator iter = 
      std::find_if(this->Internal->Nodes.begin(),
                      this->Internal->Nodes.end(),
                      this->Internal->FindNodeOutOfRange );
  
    if ( iter != this->Internal->Nodes.end() )
      {
      delete *iter;
      this->Internal->Nodes.erase(iter);
      this->Modified();
      done = 0;
      }
    }

  this->SortAndUpdateRange();


  return 1;
}

//----------------------------------------------------------------------------
// Print method for vtkRGBAColorTransferFunction
void vtkRGBAColorTransferFunction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Size: " << this->Internal->Nodes.size() << endl;
  if ( this->Clamping )
    {
    os << indent << "Clamping: On\n";
    }
  else
    {
    os << indent << "Clamping: Off\n";
    }

  if ( this->ColorSpace == VTK_CTF_RGB )
    {
    os << indent << "Color Space: RGB\n";
    }
  else if ( this->ColorSpace == VTK_CTF_HSV && this->HSVWrap )
    {
    os << indent << "Color Space: HSV\n";
    }
  else if ( this->ColorSpace == VTK_CTF_HSV )
    {
    os << indent << "Color Space: HSV (No Wrap)\n";
    }
  else
    {
    os << indent << "Color Space: CIE-L*ab\n";
    }
  
  if ( this->Scale == VTK_CTF_LOG10 )
    {
    os << indent << "Scale: Log10\n";
    }
  else
    {
    os << indent << "Scale: Linear\n";
    }
  
  os << indent << "Range: " << this->Range[0] << " to " 
     << this->Range[1] << endl;

  os << indent << "AllowDuplicateScalars: " << this->AllowDuplicateScalars << endl;

  unsigned int i;
  for( i = 0; i < this->Internal->Nodes.size(); i++ )
    {
    os << indent << "  " << i 
       << " X: " << this->Internal->Nodes[i]->X 
       << " R: " << this->Internal->Nodes[i]->R 
       << " G: " << this->Internal->Nodes[i]->G 
       << " B: " << this->Internal->Nodes[i]->B 
       << " Sharpness: " << this->Internal->Nodes[i]->Sharpness 
       << " Midpoint: " << this->Internal->Nodes[i]->Midpoint << endl;
    }
}



