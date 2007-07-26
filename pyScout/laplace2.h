
#ifndef _h_laplace_2h
#define _h_laplace_2h

#include <vector>

#include "WrapperPolyData.h"
#include "tracer.h"


// typedefs

typedef std::vector<double> DoubleVector;


typedef std::vector<DoubleVector> WrapperLineType;
typedef std::vector<WrapperLineType> WrapperLineSetType;

/*

compute profiles at a given spacing, after an initial offset

*/
WrapperLineSetType
ComputeProfiles(int offset,
		double dspacing,
		WrapperLineType& referenceLine,
		const Tracer& tracer);

/*

compute isolines for each of the values inserted

in order to return profiles ordered consistently
with the first point, the coordinates of the first
point are provided

*/
WrapperLineSetType
ComputeIsolines(const DoubleVector& vec,
		const Tracer& tracer,
		double x0,
		double y0);


Tracer
ComputeSolution(WrapperPolyData* wrapper,
		int pad1, int pad2,
		double dresolution,
		int convergenceCriterion);

int InitializePackage();
int Finalize();

#endif
