#ifndef _Fibers_h_
#define _Fibers_h_

#include <itkPoint.h>
#include <itkDiffusionTensor3D.h>
//#include <itkTensor.h>
#include <type_traits>

#include <typeinfo>
#include <ostream>


  template <  class T , unsigned int NDimension, class TTensorCoord = float >
    class Fiber
  {
  public:
    typedef T                                ScalarType;
    typedef Fiber                            Self;
    typedef Point<ScalarType, NDimension>    PointType;
    //typedef Tensor<TTensorCoord, NDimension> TensorType;
//typedef std::conditional<NDimension== 3, itk::DiffusionTensor3D<TTensorCoord>,itk::SymmetricSecondRankTensor<TTensorCoord>>::type TensorType;
//    typedef DiffusionTensor3D<TTensorCoord> TensorType;
typedef      itk::SymmetricSecondRankTensor<TTensorCoord,NDimension> TensorType;
    struct FiberPoint
    {
       PointType  Point;
       TensorType Tensor;
       FiberPoint(): Tensor(0.0) 
       {
         for(unsigned int i=0; i<NDimension; i++)
           Point[i] = 0.0;
       };
    };
     
    typedef std::vector<FiberPoint> FiberPointListType;
     

    itkStaticConstMacro (Dimension, unsigned int, NDimension);
    
    
    /** add a point to the tail of the fiber. */
    void AddPoint ( const FiberPoint & );

    /** set the list of points */
    void SetPointList ( const FiberPointListType &l)
    { m_FiberPointList = l; }
    
    /** get the list of points */
    FiberPointListType GetPointList (void) const
    { return m_FiberPointList; }
    
    /** return the ith point (if it exists)*/
    FiberPoint GetPoint (const int) const;

    /** merge two fibers */
    void MergeWith (const Self& );

    /** return the geodesic length of the fiber*/
    double GetLength (void) const;

    /** Get the Euclidean length*/
    double GetEuclideanLength (void) const;

    /** return the number of points */
    unsigned int GetNumberOfPoints (void) const
    { return m_FiberPointList.size(); }

    /** Empties the point list */
    void Clear (void)
    { m_FiberPointList.clear(); }
    
    /** Integrates the FA along the fiber */
    void GetFAStatistics (double &mean, double &min, double &max, double &var) const;

    /** Integrates the ADC along the fiber */
    void GetADCStatistics (double &mean, double &min, double &max, double &var) const;

    enum StatisticsType
    {
        FA,
        ADC
    };

    void GetStatistics (StatisticsType type, double &mean, double &min, double &max, double &var) const;


    Fiber(){};
    ~Fiber(){};
    Fiber (const Self& f)
    {
      m_FiberPointList  = f.GetPointList();
    }
    Self& operator=(const Self& f);

protected:
    inline double GetTensorScalarValue (const TensorType &tensor, const StatisticsType &type) const;
    
    
private:
    FiberPointListType  m_FiberPointList;

  };


  template <class T, unsigned int NDimension>
    std::ostream & operator<<(std::ostream &os, const Fiber<T,NDimension> &f)
  {
    for( unsigned int i=0; i<f.GetNumberOfPoints(); i++)
      os << f.GetPointList()[i].Point << " " << f.GetPointList()[i].Tensor << std::endl;

    return os;
  }
  

#endif
