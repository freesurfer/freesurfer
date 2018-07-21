#ifndef __HausdorffMembershipFunction_h
#define __HausdorffMembershipFunction_h

#include "LabelPerPointMembershipFunction.h"
//using namespace itk;
//#include "HausdorffMembershipFunction.txx"

template< class TVector >
class HausdorffMembershipFunction :  public LabelPerPointMembershipFunction< TVector >
{
	public:
		/** Standard class typedefs */
		typedef HausdorffMembershipFunction Self;
		typedef LabelPerPointMembershipFunction< TVector >    Superclass;
		typedef itk::SmartPointer<Self>                   Pointer;
		typedef itk::SmartPointer<const Self>             ConstPointer;

		/** Strandard macros */
		itkTypeMacro(HausdorffMembershipFunction,
				MembershipFunctionBase);
		itkNewMacro(Self);

		/** Typedef alias for the measurement vectors */
		typedef TVector MeasurementVectorType;
		typedef TVector CentroidType;

		/** Typedef to represent the length of measurement vectors */
		typedef typename Superclass::MeasurementVectorSizeType    MeasurementVectorSizeType;

		/**
		 * Method to get probability of an instance. The return value is the
		 * value of the density function, not probability. */
		double Evaluate(const MeasurementVectorType *m1,const MeasurementVectorType *m2)const ;
		//double Evaluate(const MeasurementVectorType &measurement) const{ std::cout << "not implemented " << std::endl;return -1;};

	protected:
		HausdorffMembershipFunction(void);
		virtual ~HausdorffMembershipFunction(void) {}

};
#include "HausdorffMembershipFunction.txx"
#endif
