/**
 * @file  LayerPropertiesDTI.cxx
 * @brief Layer properties available to DTI layers
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
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


#include <assert.h>
#include "LayerPropertiesDTI.h"
#include "vtkLookupTable.h"
#include <wx/filename.h>

using namespace std;

LayerPropertiesDTI::LayerPropertiesDTI () : LayerPropertiesMRI(),
    m_nDirectionCode ( RGB )
{
  mDirectionCodedTable = vtkSmartPointer<vtkLookupTable>::New();
  SetColorMap( DirectionCoded );
  SetResliceInterpolation( 0 );
}

LayerPropertiesDTI::~LayerPropertiesDTI ()
{}

vtkLookupTable* LayerPropertiesDTI::GetDirectionCodedTable () const
{
  return mDirectionCodedTable;
}

void LayerPropertiesDTI::ColorMapChanged()
{
  if ( GetColorMap() != DirectionCoded )
    LayerPropertiesMRI::ColorMapChanged();
  else
  {
    int ns = 64;
    mDirectionCodedTable->SetNumberOfTableValues(ns*ns*ns);
    mDirectionCodedTable->ForceBuild();
    mDirectionCodedTable->SetTableRange(0, ns*ns*ns-1);
    switch ( m_nDirectionCode )
    {
    case RGB:
      for (int i = 0; i < ns; i++)
      {
        for (int j = 0; j < ns; j++)
        {
          for (int k = 0; k < ns; k++)
          {
            mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)i/ns, (float)j/ns, (float)k/ns);
          }
        }
      }
      break;
    case RBG:
      for (int i = 0; i < ns; i++)
      {
        for (int j = 0; j < ns; j++)
        {
          for (int k = 0; k < ns; k++)
          {
            mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)i/ns, (float)k/ns, (float)j/ns);
          }
        }
      }
      break;
    case GRB:
      for (int i = 0; i < ns; i++)
      {
        for (int j = 0; j < ns; j++)
        {
          for (int k = 0; k < ns; k++)
          {
            mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)j/ns, (float)i/ns, (float)k/ns);
          }
        }
      }
      break;
    case GBR:
      for (int i = 0; i < ns; i++)
      {
        for (int j = 0; j < ns; j++)
        {
          for (int k = 0; k < ns; k++)
          {
            mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)j/ns, (float)k/ns, (float)i/ns);
          }
        }
      }
      break;
    case BRG:
      for (int i = 0; i < ns; i++)
      {
        for (int j = 0; j < ns; j++)
        {
          for (int k = 0; k < ns; k++)
          {
            mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)k/ns, (float)i/ns, (float)j/ns);
          }
        }
      }
      break;
    case BGR:
      for (int i = 0; i < ns; i++)
      {
        for (int j = 0; j < ns; j++)
        {
          for (int k = 0; k < ns; k++)
          {
            mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)k/ns, (float)j/ns, (float)i/ns);
          }
        }
      }
      break;
    }
    this->SendBroadcast( "ColorMapChanged", NULL );
  }
}

void LayerPropertiesDTI::SetDirectionCode( DirectionCode nCode )
{
  if ( m_nDirectionCode != nCode )
  {
    m_nDirectionCode = nCode;
    this->ColorMapChanged();
  }
}
