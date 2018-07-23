#ifndef __LabelsEntropyAndIntersectionMembershipFunction_h
#define __LabelsEntropyAndIntersectionMembershipFunction_h

#include "LabelPerPointMembershipFunction.h"

template< class TVector >
class LabelsEntropyAndIntersectionMembershipFunction :
	public LabelPerPointMembershipFunction< TVector >
{
	public:
		/** Standard class typedefs */
		typedef LabelsEntropyAndIntersectionMembershipFunction Self;
		typedef LabelPerPointMembershipFunction< TVector >    Superclass;
		typedef itk::SmartPointer<Self>                   Pointer;
		typedef itk::SmartPointer<const Self>             ConstPointer;

		/** Strandard macros */
		itkTypeMacro(LabelsEntropyAndIntersectionMembershipFunction,
				MembershipFunctionBase);
		itkNewMacro(Self);

		/** Typedef alias for the measurement vectors */
		typedef TVector MeasurementVectorType;
		typedef TVector CentroidType;

		/** Typedef to represent the length of measurement vectors */
		typedef typename Superclass::MeasurementVectorSizeType    MeasurementVectorSizeType;

		virtual double Evaluate(const MeasurementVectorType *m1,const MeasurementVectorType *m2)const ;
		void SetIntersection(bool inter)
		{ this->m_intersection = inter; }
		void SetEntropy(bool entro)
		{ this->m_entropy= entro ; }
		void SetLabels(bool labels)
		{ this->m_labels = labels ; }
		void SetDice(bool dice)
		{ this->m_dice = dice ; }
		void SetEuclidean(bool euclidean)
		{ this->m_euclidean = euclidean ; }
		void SetKulczynskis(bool kulczynskis)
		{this->m_kulczynskis= kulczynskis;	}
		void SetRuzicka(bool js)
		{this->m_ruzicka= js;	}
		void SetJensenShannon(bool js)
		{this->m_jensenShannon= js;	}
		void SetLabelsAndEuclid(bool lae)
		{this->m_labelsAndEuclid= lae;	}
		void SetMeanEuclidean(bool me)
		{this->m_meanEuclidean= me;	}
		void SetGaussian(bool g)
		{this->m_gaussian= g;	}
		void SetMeanClosestPointGaussian(bool mcp)
		{this->m_meanClosestPointGaussian= mcp;	}
		void SetMeanClosestPointInvert(bool mcp)
		{this->m_meanClosestPointInvert= mcp;	}
		void SetMeanAndCovInvert(bool m)
		{this->m_meanAndCovInvert= m;	}
		void SetMeanAndCovGaussian(bool m)
		{this->m_meanAndCovGaussian= m;	}

	protected:
		LabelsEntropyAndIntersectionMembershipFunction(void):Superclass(){
			this->m_intersection = false;
			this->m_entropy = false;
			this->m_labels=false; 
			this->m_labelsAndEuclid=false; 
			this->m_euclidean=false; 
			this->m_dice=false; 
			this->m_kulczynskis=false;
			this->m_jensenShannon=false;
			this->m_ruzicka=false;
			this->m_meanEuclidean=false;
			this->m_gaussian=false;
			this->m_meanClosestPointInvert=false;
			this->m_meanAndCovInvert=false;
			this->m_meanClosestPointGaussian=false;
			this->m_meanAndCovGaussian=false;
		}
		virtual ~LabelsEntropyAndIntersectionMembershipFunction(void) {}

	private:
		bool m_intersection;
		bool m_entropy;		
		bool m_labels;		
		bool m_euclidean;		
		bool m_dice;		
		bool m_kulczynskis;
		bool m_jensenShannon;
		bool m_ruzicka;
		bool m_labelsAndEuclid;
		bool m_meanEuclidean;
		bool m_gaussian;
		bool m_meanClosestPointInvert;
		bool m_meanAndCovInvert;
		bool m_meanClosestPointGaussian;
		bool m_meanAndCovGaussian;
};
#include "LabelsEntropyAndIntersectionMembershipFunction.txx"
#endif
