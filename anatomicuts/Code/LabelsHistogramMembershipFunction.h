#ifndef __LabelsHistogramMembershipFunction_h
#define __LabelsHistogramMembershipFunction_h

#include "LabelPerPointMembershipFunction.h"

template< class TVector >
class  LabelsHistogramMembershipFunction :
	public LabelPerPointMembershipFunction< TVector >
{
	public:
		/** Standard class typedefs */
		typedef LabelsHistogramMembershipFunction Self;
		typedef LabelPerPointMembershipFunction< TVector >    Superclass;
		typedef itk::SmartPointer<Self>                   Pointer;
		typedef itk::SmartPointer<const Self>             ConstPointer;

		/** Strandard macros */
		itkTypeMacro(LabelsHistogramMembershipFunction,
				MembershipFunctionBase);
		itkNewMacro(Self);

		/** Typedef alias for the measurement vectors */
		typedef TVector MeasurementVectorType;
		typedef TVector CentroidType;

		/** Typedef to represent the length of measurement vectors */
		typedef typename Superclass::MeasurementVectorSizeType    MeasurementVectorSizeType;

		virtual double Evaluate(const MeasurementVectorType *m1,const MeasurementVectorType *m2)const ;
	protected:
		LabelsHistogramMembershipFunction(void):Superclass(){};
		virtual ~LabelsHistogramMembershipFunction(void) {}

	private:
		
};

#endif
