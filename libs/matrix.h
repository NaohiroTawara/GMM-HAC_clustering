#ifndef matrix_h
#define matrix_h

#include <Eigen/Dense>


// 表示用マクロ:
#define print(var)  \
  std::cout<<#var"= "<<std::endl<<(var)<<std::endl

#ifdef DOUBLE // Execute each calcration with double precision
typedef Eigen::MatrixXd MATRIX;
typedef Eigen::VectorXd VECTOR;
typedef double SCALAR;
#else
typedef Eigen::MatrixXf MATRIX;
typedef Eigen::VectorXf VECTOR;
typedef float SCALAR;
#endif

#endif

