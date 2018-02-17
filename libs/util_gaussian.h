
#ifndef __UTIL_GAUSSIAN_H__
#define __UTIL_GAUSSIAN_H__

#include <cv.h>
#include "util.h"

// マハラノビス距離を返す関数
inline double _mahala(CvMat* _x, CvMat* _mu, CvMat* _sigma, int _cov_type)
{
  int dimension = cvGetSize(_x).width;
  int numframes = cvGetSize(_x).height;
  double mahala = 0;

  if (dimension != cvGetSize(_mu).width)
    Error(1111, "Invalid Dimension: %d; %d", 
	  dimension, cvGetSize(_mu).width);
  if (dimension != cvGetSize(_sigma).width)
    Error(1111, "Invalid Dimension: %d; %d", 
	  dimension, cvGetSize(_sigma).width);

  // マハラノビス距離の算出
  if (_cov_type == DIAGC)
  {
    for (int i = 0; i < numframes; ++i)
      for (int j = 0; j < dimension; ++j)
	mahala += square(cvmGet(_x, i, j) - cvmGet(_mu, 0, j)) 
	  / cvmGet(_sigma, 0, j);  
  }
  return mahala;
}




#endif
