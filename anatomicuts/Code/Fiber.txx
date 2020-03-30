/*=========================================================================

  Program:   Tensor ToolKit - TTK
  Module:    $URL: https://scm.gforge.inria.fr/svn/ttk/trunk/Common/itkFiber.txx $
  Language:  C++

  Copyright (c) INRIA 2010. All rights reserved.
  See LICENSE.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itk_Fiber_txx
#define _itk_Fiber_txx
#include "itkFiber.h"

namespace itk
{

  template <class T, unsigned int TDimension, class TTensorCoord>
  Fiber<T, TDimension, TTensorCoord>&
  Fiber<T, TDimension, TTensorCoord>::
  operator=(const Self& f)
  {
    this->SetPointList ( f.GetPointList() );
    return *this;
  }
  
  
  
  template <class T, unsigned int TDimension, class TTensorCoord>
  void
  Fiber<T, TDimension, TTensorCoord>  
  ::AddPoint (const FiberPoint &P)
  {
    m_FiberPointList.push_back(P);
  }



  template <class T, unsigned int TDimension, class TTensorCoord>
  typename Fiber<T, TDimension, TTensorCoord>::FiberPoint
  Fiber<T, TDimension, TTensorCoord>
  ::GetPoint (const int i) const
  {
    if(i<0 || i>(int)(m_FiberPointList.size()-1))
    {
      std::cerr << "Index i is out of range, returning 0" << std::endl;
      return FiberPoint();
    }

    return m_FiberPointList[i];
    
  }


  template <class T, unsigned int TDimension, class TTensorCoord>
  void
  Fiber<T, TDimension, TTensorCoord>
  ::MergeWith (const Self& F)
  {

    // first: check if both fibers are null
    if( this->GetNumberOfPoints()==0 && F.GetNumberOfPoints()==0 )
      return;

    if( F.GetNumberOfPoints()==0 )
      return;

    if( this->GetNumberOfPoints()==0 )
    {
      m_FiberPointList = F.GetPointList();
      return;
    }
    
    // normal cases:
    
    // second check whether the two first points are identical
    FiberPoint firstPoint1 = this->GetPoint (0);
    FiberPoint firstPoint2 = F.GetPoint (0);

    FiberPointListType  P;
    FiberPointListType  P1 = F.GetPointList();
    
    if(firstPoint1.Point == firstPoint2.Point)
    {
      for(int i=(int)(P1.size()-1);i>0;i--)
        P.push_back (P1[i]);
    }
    else
    {
      for(int i=(int)(P1.size()-1);i>=0;i--)
        P.push_back (P1[i]);
    }

    for(unsigned int i=0;i<m_FiberPointList.size();i++)
      P.push_back (m_FiberPointList[i]);

    m_FiberPointList  = P;
  }


  
  template <class T, unsigned int TDimension, class TTensorCoord>
  double
  Fiber<T, TDimension, TTensorCoord>
  ::GetLength (void) const
  {

    double length = 0.0;
    if( m_FiberPointList.size() < 2 )
      return length;
    
    for(unsigned int i=0;i<m_FiberPointList.size()-1;i++)
    {
      PointType current = m_FiberPointList[i].Point;
      PointType next = m_FiberPointList[i+1].Point;
      double dist = 0.0;
      
      for(unsigned int j=0;j<TDimension;j++)
        dist += (next[j] - current[j])*(next[j] - current[j]);

      length += sqrtf (dist);
      
    }

    return length;
  }

  

  template <class T, unsigned int TDimension, class TTensorCoord>
  double
  Fiber<T, TDimension, TTensorCoord>
  ::GetEuclideanLength (void) const
  {
    PointType Start = m_FiberPointList[0];
    PointType End = m_FiberPointList[ m_FiberPointList.size()-1 ];

    double dist = 0.0;

    for(unsigned int i=0;i<TDimension;i++)
      dist += ( End[i] - Start[i] )*( End[i] - Start[i] );

    return sqrt (dist);
    
  }

  template <class T, unsigned int TDimension, class TTensorCoord>
  void
  Fiber<T, TDimension, TTensorCoord>
  ::GetStatistics (StatisticsType type, double &mean, double &min, double &max, double &var) const
  {
      if (m_FiberPointList.size()==0)
      {
          mean = 0.0;
          min  = 0.0;
          max  = 0.0;
          var  = 0.0;
          return;
      }

      // first point
      PointType current        = m_FiberPointList[0].Point;
      TensorType currentTensor = m_FiberPointList[0].Tensor;

      if (m_FiberPointList.size()==1)
      {
          mean = this->GetTensorScalarValue (currentTensor, type);
          min  = mean;
          max  = mean;
          var  = 0.0;
          return;
      }

      unsigned int pointCount = m_FiberPointList.size();

      double minStat       = 99999999.9;
      double maxStat       = -1.0;
      double totalLength   = 0.0;
      double totalStat     = 0.0;
      double dist          = 0.0;

      PointType next        = m_FiberPointList[1].Point;
      TensorType nextTensor = m_FiberPointList[1].Tensor;

      dist = current.EuclideanDistanceTo (next) * 0.5;

      double stat = this->GetTensorScalarValue (currentTensor, type);

      if (stat<minStat)
          minStat = stat;
      if (stat>maxStat)
          maxStat = stat;

      totalStat += stat * dist;

      totalLength += dist;

      for(unsigned int i=1; i<m_FiberPointList.size()-1; i++)
      {
          current       = next;
          currentTensor = nextTensor;

          next       = m_FiberPointList[i+1].Point;
          nextTensor = m_FiberPointList[i+1].Tensor;

          double oldDist = dist;

          dist = current.EuclideanDistanceTo (next) * 0.5;

          stat = this->GetTensorScalarValue (currentTensor, type);;

          if (stat<minStat)
              minStat = stat;
          if (stat>maxStat)
              maxStat = stat;

          totalStat += stat * (dist + oldDist);

          totalLength += dist + oldDist;
      }

      // case of last point
      current = next;
      currentTensor = nextTensor;

      next       = m_FiberPointList[pointCount-1].Point;
      nextTensor = m_FiberPointList[pointCount-1].Tensor;

      dist = current.EuclideanDistanceTo (next) * 0.5;

      stat = this->GetTensorScalarValue (currentTensor, type);;

      if (stat<minStat)
          minStat = stat;
      if (stat>maxStat)
          maxStat = stat;

      totalStat += stat * dist;

      totalLength += dist;

      if (totalLength>0.0)
          totalStat /= totalLength;

      mean = totalStat;
      min  = minStat;
      max  = maxStat;
      var  = 0.0;
  }

  template <class T, unsigned int TDimension, class TTensorCoord>
  double
  Fiber<T, TDimension, TTensorCoord>
  ::GetTensorScalarValue(const TensorType &tensor, const StatisticsType &type) const
  {
      double scalar = 0.0;
      switch(type)
      {
      case ADC:
          scalar = tensor.GetTrace();
          break;

      case FA:
      default:
          scalar = tensor.GetFA();
          break;
      }

      return scalar;
  }

  template <class T, unsigned int TDimension, class TTensorCoord>
  void
  Fiber<T, TDimension, TTensorCoord>
  ::GetFAStatistics (double &mean, double &min, double &max, double &var) const
  {
      this->GetStatistics(FA, mean, min, max, var);
  }
  
  template <class T, unsigned int TDimension, class TTensorCoord>
  void
  Fiber<T, TDimension, TTensorCoord>
  ::GetADCStatistics(double &mean, double &min, double &max, double &var) const
  {
      this->GetStatistics(ADC, mean, min, max, var);
  }
  

} // end of namespace


#endif
