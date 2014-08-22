/**
 * @file  testMultiplication.cpp  
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Benjamin Lewin 
 * CVS Revision Info:
 *    $Author: blewin $
 *    $Date: 2014/08/22 21:22:50 $
 *    $Revision: 1.2 $
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
//testMultiplication.cpp
//
//Created 7/28/2014
//By: Benjamin Lewin
//
//A test file to establish the best method for multiplication
//with or without reshaping
//


#include <iostream>
#include <vnl/vnl_matrix.h>
#include "matrix3d.h"
#include <iomanip>
#include <time.h>

const int ROWS = 400;
const int COLS = 400;
const int SLICES = 400;
const int NUM_BASIS = 15;

const int TEST_ITERATIONS = 1;

int main()
{
  

  Matrix3d data(ROWS, COLS, SLICES);
  vnl_matrix<float> Bxt(NUM_BASIS, ROWS);
  vnl_matrix<float> result(NUM_BASIS, COLS * SLICES);
  vnl_matrix<float> result2(NUM_BASIS, COLS * SLICES);
  
  int count = 0;
  for(int i = 0; i < SLICES; i++)
    for(int j = 0; j < ROWS; j++)
      for(int k = 0; k < COLS; k++)
        data.setVal(j, k, i, count++ % 11);

  count = 0;
  for(int i = 0; i < NUM_BASIS; i++)
    for(int j = 0; j < ROWS; j++)
      Bxt(i, j) = count++ % 7 ;


// create 2d matrix to contain all elements then multiply
std::cout << "Begin test with full reshape:\n";
clock_t begin1 = clock();
vnl_matrix<float> reshaped(ROWS, COLS * SLICES);
for(int j = 0; j < TEST_ITERATIONS; j++) {
  for(int i = 0; i < SLICES; i++)
    reshaped.set_columns(i * COLS, data.getSlice(i));
  result2 = Bxt * reshaped;
}
clock_t end1 = clock();
double time1 = difftime(end1, begin1);
std::cout << "Time elapsed: " << time1 << " ms\n\n";

// multiply each slice into the result matrix seperately
std::cout << "Begin test with multiplication by slice:\n";
clock_t begin2 = clock();
for(int j = 0; j < TEST_ITERATIONS; j++)
  for(int i = 0; i < SLICES; i++)
    result.set_columns(i * COLS, Bxt * data.getSlice(i));
clock_t end2 = clock();
double time2 = difftime(end2, begin2);
std::cout << "Time elapsed: " << time2 << " ms\n\n";

std::cout << "Full reshape took " << time1 / time2 << " times longer\n";



  if(result == result2)
    std::cout << "METHODS EQUIVALENT\n";
  else
    std::cout << "SOMETHING HAS GONE HORRIBLY WRONG :(\n";

/*
//print the resulting matrix
  int rWidth = COLS * SLICES;
  for(int i = 0; i < NUM_BASIS; i++) {
    for(int j = 0; j < rWidth; j++) {
      std::cout << std::setw(4) << result(i, j) << ' ';
    }
    std::cout << std::endl;
  }
*/

}

