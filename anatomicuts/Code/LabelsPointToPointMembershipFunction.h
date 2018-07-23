#ifndef __LabelsPointToPointMembershipFunction_h
#define __LabelsPointToPointMembershipFunction_h

#include "LabelPerPointMembershipFunction.h"

template< class TVector >
class  LabelsPointToPointMembershipFunction :
	public LabelPerPointMembershipFunction< TVector >
{
	public:
		/** Standard class typedefs */
		typedef LabelsPointToPointMembershipFunction Self;
		typedef LabelPerPointMembershipFunction<TVector> Superclass; 
		typedef itk::SmartPointer<Self>                   Pointer;
		typedef itk::SmartPointer<const Self>             ConstPointer;

		/** Strandard macros */
		itkTypeMacro(LabelsPointToPointMembershipFunction,LabelPerPointMembershipFunction);
		itkNewMacro(Self);

		/** Typedef alias for the measurement vectors */
		typedef TVector MeasurementVectorType;
		typedef TVector CentroidType;

		/** Typedef to represent the length of measurement vectors */
		typedef typename Superclass::MeasurementVectorSizeType    MeasurementVectorSizeType;
		void SetLabelsCount(int count)
		{
			this->m_labelsCount = count;
		}
		virtual double Evaluate(const MeasurementVectorType *m1,const MeasurementVectorType *m2)const ;
	protected:
		LabelsPointToPointMembershipFunction():Superclass(){}
		virtual ~LabelsPointToPointMembershipFunction(void) {}

	private:
		int m_labelsCount;	
		
};
#include "LabelsPointToPointMembershipFunction.txx"
#endif
