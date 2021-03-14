/**
 * @brief Stores all data associated with a subject.
 *
 * This is one row from the input data table file (qdec.table.dat).
 */
/*
 * Original Author: Nick Schmansky
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

#include <cstring>
#include <sstream>

#include "QdecSubject.h"

// Constructors/Destructors
//

QdecSubject::QdecSubject(std::string isId, std::vector<QdecFactor *> iFactors) {
  msId     = isId;
  mFactors = iFactors;
}

QdecSubject::~QdecSubject() {
  while (mFactors.size() != 0) {
    delete mFactors.back();
    mFactors.pop_back();
  }
}

/**
 * @return string
 * @param  isFactorName
 */
std::string QdecSubject::GetDiscreteFactorValue(const char *isFactorName) {
  for (unsigned int i = 0; i < mFactors.size(); i++) {
    if (mFactors[i]->IsDiscrete()) {
      if (0 == strcmp(mFactors[i]->GetFactorName().c_str(), isFactorName)) {
        return mFactors[i]->GetDiscreteValue();
      }
    }
  }

  std::stringstream ssErr;
  ssErr << "ERROR: QdecSubject::GetDiscreteFactor failure: could not find "
           "factor name: "
        << isFactorName << " for subject " << this->GetId();
  // cerr << ssErr.str() << std::endl;
  throw std::runtime_error(ssErr.str().c_str());
  return nullptr;
}

/**
 * @return double
 * @param  isFactorName
 */
double QdecSubject::GetContinuousFactorValue(const char *isFactorName) {
  for (unsigned int i = 0; i < mFactors.size(); i++) {
    if (0 == strcmp(mFactors[i]->GetFactorName().c_str(), isFactorName)) {
      return mFactors[i]->GetContinuousValue();
    }
  }

  std::stringstream ssErr;
  ssErr << "ERROR: QdecSubject::GetContinuousFactor failure: could not find "
           "factor name: "
        << isFactorName << " for subject " << this->GetId();
  // cerr << ssErr.str() << std::endl;
  throw std::runtime_error(ssErr.str().c_str());
  return 0.0;
}

/**
 * @return vector < QdecFactor* >
 */
std::vector<QdecFactor *> QdecSubject::GetContinuousFactors() {
  std::vector<QdecFactor *> factors;
  for (unsigned int i = 0; i < mFactors.size(); i++) {
    if (mFactors[i]->IsContinuous()) {
      factors.push_back(mFactors[i]);
    }
  }

  return factors;
}

/**
 * @return QdecFactor
 */
QdecFactor *QdecSubject::GetFactor(const char *isFactorName) {
  for (unsigned int i = 0; i < mFactors.size(); i++) {
    if (0 == strcmp(mFactors[i]->GetFactorName().c_str(), isFactorName)) {
      return mFactors[i];
    }
  }

  std::stringstream ssErr;
  ssErr << "ERROR: QdecSubject::GetFactor failure: could not find "
           "factor name: "
        << isFactorName << " for subject " << this->GetId();
  // cerr << ssErr.str() << std::endl;
  throw std::runtime_error(ssErr.str().c_str());
  return nullptr;
}

/**
 * @param  isFactorName
 */
void QdecSubject::DeleteFactor(const char *isFactorName) {
  std::vector<QdecFactor *>::iterator iter = mFactors.begin();
  while (iter != mFactors.end()) {
    QdecFactor *factor     = *iter;
    std::string factorName = factor->GetFactorName();
    if (0 == strcmp(factorName.c_str(), isFactorName)) {
      mFactors.erase(iter);
      return;
    }
    ++iter;
  }
}
