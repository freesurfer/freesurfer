/*
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


#ifndef GaussianDistribution_h
#define GaussianDistribution_h

#include "mri_ca_util.h"
#include "math.h"
using namespace std;


const float fNormalizationConstant=1/sqrt(2*M_PI);
const float fMinAllowedVariance=0.001;
const float fReplacementVariance=25.0;
// This class should have a fixed size so that it can be binary written/read out easily by CSparse2DGausDistMatrix
class CGaussianDistribution
{

public:
  float density(const TypeVectorFloat& measure)
  {
    float fDensity;
    // [[]] Need to enhance this routine to compute
    // the density of a given measurement vector occurring
    //
    //   _             _    T   -1 _
    // p(x) = exp[-1/2(x-Mx)  Vx  (x-Mx)]
    //        ---------------------------
    //                 N/2     1/2
    //            (2Pi)    |Vx|
    //
    //
    // where:
    // _
    // x  = measurement vector, measure
    // Mx = mean vector
    // Vx = CovMatrix
    if (!bMeanVectorIsUpToDate)
    {
      meanVector();
    }

    if (!bCovMatrixIsUpToDate)
    {
      covMatrix();
    }

    //
    // PreCompute fixed denominator term
    // PreCompute inverse of Vx

    //                                 _
    // compute vector substraction for x-Mx


    float fVariance=covariance[0][0];
    if (fVariance<fMinAllowedVariance)
    {
      fVariance=fReplacementVariance;
    }
    float fStdDev=sqrt(fVariance);
    float fTerm1=fNormalizationConstant/fStdDev;
    float fTerm3=exp(-0.5*(measure[0]-mean[0])*(measure[0]-mean[0])/fVariance);
    fDensity=fTerm1*fTerm3;
    return(fDensity);
  }


  CGaussianDistribution()
  {
    bInitialized=false;
  }


  CGaussianDistribution(int nMeasureDim)
  {
    bInitialized=true;
    floatNumMeasuresEntered=0.0;
    nDim=nMeasureDim;
    sumX.clear();
    sumXOuterProdX.clear();

    for (int i=0; i<nDim; i++)
    {
      sumX.push_back(0.0);
      mean.push_back(0);

      TypeVectorFloat vectFloatRow;
      for (int j=0; j<nDim; j++)
      {
        vectFloatRow.push_back(0.0);
      }
      covariance.push_back(vectFloatRow);
      sumXOuterProdX.push_back(vectFloatRow);
    }

    bMeanVectorIsUpToDate=true;
    bCovMatrixIsUpToDate=true;

  }



  void init(TypeVectorFloat& vectFloatMeasureMean,TypeMatrixFloat& matrixFloatMeasureVariance, float  floatNewNumMeasuresEntered)
  {
    bInitialized=true;
    floatNumMeasuresEntered=floatNewNumMeasuresEntered;
    nDim=vectFloatMeasureMean.size();

    bMeanVectorIsUpToDate=true;
    bCovMatrixIsUpToDate=true;
    mean=vectFloatMeasureMean;
    covariance=matrixFloatMeasureVariance;

    sumX.clear();
    sumXOuterProdX.clear();
    for (int i=0; i<nDim; i++)
    {
      sumX.push_back(floatNumMeasuresEntered*vectFloatMeasureMean[i]);

      TypeVectorFloat vectFloatRow;
      for (int j=0; j<nDim; j++)
      {
        vectFloatRow.push_back(floatNumMeasuresEntered*matrixFloatMeasureVariance[i][j]);
        // reconstruct the sum[xy] so we can add new measures to the distribution and recompute the parameters
        // Cov(x,y)=E[xy]-E[x]E[y]=(1/numMeasEntered)*wtdsum[xy]  - E[x] * E[y]
        // thus sum[xy]=numMeasEntered*(Cov(x,y) + E[x] * E[y]);
        vectFloatRow.push_back(floatNumMeasuresEntered*(mean[i]*mean[j] + covariance[i][j]));

      }
      sumXOuterProdX.push_back(vectFloatRow);
    }
  }

  TypeVectorFloat meanVector()
  {
    if (!bMeanVectorIsUpToDate)
    {
      bMeanVectorIsUpToDate=true;
      mean.clear();
      for (int i=0; i<nDim; i++)
      {
        mean.push_back(sumX[i]/floatNumMeasuresEntered);
      }

    }
    return(mean);
  }

  TypeMatrixFloat covMatrix()
  {
    meanVector();

    if (!bCovMatrixIsUpToDate)
    {
      bCovMatrixIsUpToDate=true;
      covariance.clear();
      for (int i=0; i<nDim; i++)
      {
        TypeVectorFloat vectFloatRow;
        for (int j=0; j<nDim; j++)
        { // Cov[x,y] = E[xy] - E[x]E[y]
          vectFloatRow.push_back(  (sumXOuterProdX[i][j]/floatNumMeasuresEntered  -  mean[i]*mean[j])  );
        }
        covariance.push_back(vectFloatRow);
      }
    }
    return(covariance);
  }


  void insertMeasure(const TypeVectorFloat& measure, float fFractionOfAMeasure=1.0)
  {
    bMeanVectorIsUpToDate=false;
    bCovMatrixIsUpToDate=false;

    if (bInitialized==false)
    {
      bInitialized=true;
      floatNumMeasuresEntered=fFractionOfAMeasure;
      nDim=measure.size();

      sumX.clear();
      sumXOuterProdX.clear();
      for (int i=0; i<nDim; i++)
      {
        float fTerm1=fFractionOfAMeasure*measure[i];
        sumX.push_back(fTerm1);

        TypeVectorFloat vectFloatRow;
        for (int j=0; j<nDim; j++)
        {
          vectFloatRow.push_back(fTerm1*fFractionOfAMeasure*measure[j]);
        }
        sumXOuterProdX.push_back(vectFloatRow);
      }
    }
    else
    {
      if (measure.size()==(unsigned int) nDim)
      {
        floatNumMeasuresEntered+=fFractionOfAMeasure;
        for (int i=0; i<nDim; i++)
        { // add mean to sum of means
          float fTerm1=fFractionOfAMeasure*measure[i];
          sumX[i]=sumX[i]+fTerm1;

          // add outer product to the sum of outer products
          for (int j=0; j<nDim; j++)
          {     // E[X^2] = wtd sum of Xi^2  / number of measurements
            sumXOuterProdX[i][j]=sumXOuterProdX[i][j] + (fTerm1*measure[j]);
          }

        }
      }
    }

  }

  float numMeasuresEntered()
  {
    return(floatNumMeasuresEntered);
  }


  int getDim()
  {
    return(nDim);
  }

  bool write(CConfigFile& configfile)
  {
    string strSectionName="GaussianDistribution";
    // bring the stats up to date
    meanVector();
    covMatrix();

    configfile.write(floatNumMeasuresEntered,  strSectionName, "NumberOfMeasuresEntered");
    configfile.write(mean,strSectionName,"MeasureMean");
    configfile.write(covariance, strSectionName, "MeasureVariance");

    return(true);
  }
  // [[]] This function is wrong and needs to be discarded or corrected
  bool get(CConfigFile& configfile)
  {
    string strSectionName="GaussianDistribution";
    float  dNumMeasuresEntered;
    configfile.get(dNumMeasuresEntered,  strSectionName, "NumberOfMeasuresEntered");
    floatNumMeasuresEntered=dNumMeasuresEntered;
    configfile.get(mean, strSectionName,"MeasureMean");
    configfile.get(covariance, strSectionName, "MeasureVariance");

    sumX.clear();
    sumXOuterProdX.clear();
    for (int i=0; i<nDim; i++)
    {
      sumX.push_back(floatNumMeasuresEntered*mean[i]);

      TypeVectorFloat vectFloatRow;
      for (int j=0; j<nDim; j++)
      {
        vectFloatRow.push_back(floatNumMeasuresEntered*covariance[i][j]);
      }
      sumXOuterProdX.push_back(vectFloatRow);
    }

    bInitialized=true;
    bMeanVectorIsUpToDate=true;
    bCovMatrixIsUpToDate=true;

    return(true);
  }



private:  // size is critical. Must conserve space!
  bool bInitialized;

  int nDim;
  float  floatNumMeasuresEntered;
  TypeVectorFloat sumX;  // sum of the meas
  TypeMatrixFloat sumXOuterProdX; // sum of the outer product of (meas w itself)

  bool bMeanVectorIsUpToDate;
  bool bCovMatrixIsUpToDate;
  TypeVectorFloat mean;
  TypeMatrixFloat covariance;

private:
  friend istream& operator>>(istream&,CGaussianDistribution&);
  friend ostream& operator<<(ostream&,CGaussianDistribution&);
  friend ostream& print(ostream&,CGaussianDistribution&);
};

ostream& print(ostream& os, CGaussianDistribution& gausDist)
{
  os << "mean[0]= " << gausDist.mean[0];
  os << " , covariance[0,0]= " << gausDist.covariance[0][0];
  return(os);
}


istream& operator>>(istream& is, CGaussianDistribution& gausDist)
{

  // Read in from the stream as binary data
  gausDist.bInitialized=true;
  gausDist.mean.clear();
  gausDist.covariance.clear();
  gausDist.sumX.clear();
  gausDist.sumXOuterProdX.clear();
  is.read(&gausDist.floatNumMeasuresEntered, sizeof(gausDist.floatNumMeasuresEntered));

  // Dont read or write this because it wastes space, instead assume that the input gausDist has the appropriate size
  // is.read(&gausDist.nDim, sizeof(gausDist.nDim));


  for (int i=0; i<gausDist.nDim; i++)
  {
    float d;
    is.read(&d, sizeof(d));
    gausDist.mean.push_back(d);
  }

  for (int i=0; i<gausDist.nDim; i++)
  {
    TypeVectorFloat vectFloat;
    for (int j=0; j<gausDist.nDim; j++)
    {

      float d;
      is.read(&d, sizeof(d));
      vectFloat.push_back(d);
    }
    gausDist.covariance.push_back(vectFloat);
  }


  for (int i=0; i<gausDist.nDim; i++)
  {
    gausDist.sumX.push_back(gausDist.floatNumMeasuresEntered*gausDist.mean[i]);

    TypeVectorFloat vectFloatRow;
    for (int j=0; j<gausDist.nDim; j++)
    {
      // reconstruct the sum[xy] so we can add new measures to the distribution and recompute the parameters
      // Cov(x,y)=E[xy]-E[x]E[y]=(1/numMeasEntered)*sum[xy]  - E[x] * E[y]
      // thus sum[xy]=numMeasEntered*(Cov(x,y) + E[x] * E[y]);
      vectFloatRow.push_back(gausDist.floatNumMeasuresEntered*(gausDist.mean[i]*gausDist.mean[j] + gausDist.covariance[i][j]));
    }
    gausDist.sumXOuterProdX.push_back(vectFloatRow);
  }

  gausDist.bMeanVectorIsUpToDate=true;
  gausDist.bCovMatrixIsUpToDate=true;


  return is;
}





ostream& operator<<(ostream& os, CGaussianDistribution& gausDist)
{
  // write out the dist as binary data
  os.write(&gausDist.floatNumMeasuresEntered, sizeof(gausDist.floatNumMeasuresEntered));
  // Dont write this out bc it wastes space since it is the same for every gaus dist in the sparse 2D matrix in which this class is contained
  // os.write(&gausDist.nDim, sizeof(gausDist.nDim));

  // bring mean and cov up to date
  gausDist.meanVector();
  gausDist.covMatrix();


  for (int i=0; i<gausDist.nDim; i++)
  {
    float d=gausDist.mean[i];
    os.write(&d, sizeof(d));
  }

  for (int i=0; i<gausDist.nDim; i++)
  {
    TypeVectorFloat vectFloat;
    for (int j=0; j<gausDist.nDim; j++)
    {


      float d=gausDist.covariance[i][j];


      if (d<0)
      {
        //                cout << "Error negative variance\n ";
      }


      os.write(&d, sizeof(d));
    }
  }

  return os;
}

#endif















