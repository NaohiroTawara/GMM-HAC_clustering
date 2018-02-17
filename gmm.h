//
//  gmm.h
//  spkr_cl_GMM-HAC
//
//  GMM のパラメタ管理と EM アルゴリズムによるパラメタの最尤推定を行うクラス
//  Created by 直弘 俵 on 14/07/29.
//  Copyright 2012年-2014年 早稲田大学. All rights reserved.
//

#ifndef spkr_cl_GMM_HAC_gmm_h
#define spkr_cl_GMM_HAC_gmm_h

#include <iostream>
#include <vector>
#include <float.h>

#include "util.h"
#include "matrix.h"
#include "omp.h"

using namespace std;
using namespace Eigen;


#ifndef log2pi
#define log2pi 1.837877066409
#endif


inline int isinfX(SCALAR a)
{ return isinf(a);}


inline double LAddS2(double &logL1, double &logL2)
{
  if (logL1==-DBL_MAX || isinf(logL1))
    return logL2;
  if (logL2==-DBL_MAX || isinf(logL2))
    return logL1;
  if (logL1 > logL2)
    return logL1 + log(1 + exp(logL2 - logL1));
  else
    return logL2 + log(1 + exp(logL1 - logL2)); 
}



inline VECTOR LAddS_mat(VECTOR& logL1, VECTOR& logL2)
{
  int w1 = logL1.rows();
  //  int w2 = logL2.rows();
  //  int h2 = logL2.cols();

  VECTOR out = VECTOR(w1);

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
  for (int i=0; i<w1; ++i){
    double l1 = logL1.coeff(i);
    double l2 = logL2.coeff(i);
    out(i) = LAddS2(l1,l2);
  }
  return out;
}


inline VECTOR LAddS_mat2(VECTOR& logL1, VECTOR& logL2)
{
  VECTOR tmp = (logL2.array()>logL1.array())
    .select((logL2.array() + (1 + (logL1.array() - logL2.array()).exp()).log()),
	    (logL1.array() + (1 + (logL2.array() - logL1.array()).exp()).log()));
  VECTOR tmp2 = 
    ((logL2.array() == -DBL_MAX) || (logL2.unaryExpr(ptr_fun(isinfX)).array() == 1)).select(logL1,tmp);
  VECTOR out = 
    (((logL1.array() == -DBL_MAX) || (logL1.unaryExpr(ptr_fun(isinfX)).array() == 1)).select(logL2, tmp2));

    return out;
}

/*
inline double mahala(MATRIX _x, MATRIX _mu, MATRIX _sigma, int _cov_type)
{
  int num = _x.rows();
  
  return   (_cov_type==FULLC) ?
    // full covariance
    (((_x.array() - _mu.colwise().replicate(num).array()).MATRIX() * 
      _sigma.inverse()).array() *
     (_x.array() - _mu.colwise().replicate(num).array())).sum() : 
    // diag
    ((_x.array() - _mu.colwise().replicate(num).array()) *
      _sigma.array().inverse().colwise().replicate(num).array()*
       (_x.array() - _mu.colwise().replicate(num).array())).sum();
}
*/

inline MATRIX mahala(const MATRIX& _x, const MATRIX& _mu, 
		     const MATRIX& _sigma, int _cov_type)
{
  int num = _x.rows();
  
  if (_cov_type==FULLC)
    return   
      // full covariance    
      ((-0.5 * (_x.array() - _mu.colwise().replicate(num).array()).matrix() * 
	_sigma.inverse()).array() *
       (_x.array() - _mu.colwise().replicate(num).array())).rowwise().sum();
  else
    return
      // diag
      ( -(_x.array() - _mu.colwise().replicate(num).array()).pow(2) *
	(_sigma.array().inverse().colwise().replicate(num).array() * 0.5) ).rowwise().sum();
}




inline MATRIX mahala2(const MATRIX& _x, const MATRIX& _mu, 
		      const MATRIX& _sigma, int _cov_type)
{
  int num = _x.rows();
  int dim = _x.cols();
  MATRIX lik = MATRIX(num,1);

  if (_cov_type==FULLC)
    lik = 
      ( (_x.array() - _mu.colwise().replicate(num).array()).matrix() * 
       _sigma.inverse()).array() *
      (_x.array() - _mu.colwise().replicate(num).array());
  else
  {
    for (int i=0; i<num; ++i)
    {
      SCALAR l = 0;
      for (int j=0; j<dim; ++j)
      {
	l -= (_x(i,j)-_mu(j))*(_x(i,j)-_mu(j)) *
	  (1.0/_sigma(j)) * 0.5;
      }
      lik(i)= l;
    }
  }
  return lik;
}


inline VECTOR mahala3(const MATRIX& _x, const MATRIX& _mu, 
		      const MATRIX& _sigma, int _cov_type)
{
  int num = _x.rows();
  int dim = _x.cols();
  VECTOR lik = VECTOR(num);
  int i;

  if (_cov_type==FULLC)
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < num; ++i)
    {
      lik(i) = 
	- 0.5 *( 
	( (_x.block(i,0,1,dim).array() - _mu.array()).matrix() * _sigma
	  ).array() *
	(_x.block(i,0,1,dim).array() - _mu.array())).sum();
    }
  }else{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i<num; ++i)
      lik(i) = 
	-0.5 * ((_x.block(i,0,1,dim).array() - _mu.array()).pow(2) *
		(_sigma.array())).sum();
  }
  return lik;
}



inline SCALAR floorX(SCALAR a)
{
  return floor(a);
}

inline MATRIX randperm(int _n, int _s)
{
  MATRIX m = MATRIX::Random(1,_n);
  return (((m.array()+1)*0.5*_s).unaryExpr(ptr_fun(floorX)));
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

class CGMM
{
public:
  int m_nummixtures;  // # of mixtures
  int m_dimension;    // # of dimensions
  int m_covtype;      // type of Covariance MATRIX
  int m_numData;      // # of data
  int m_trace;        // 標準出力の有無
    
  //VectorXd m_const; // Normalize terms  [1 x # of mixtures]
  //VectorXd m_pi;    // Mixture weights  [1 x # of mixtures]

  vector <SCALAR> m_pi;
  vector <SCALAR> m_const;

  // Mean vector  [1 x dimension] x # of mixtures
  vector <MATRIX> m_mu;    
  // Covariance MATRIX  [dimension x dimension] x # of mixtures
  //        or
  // Covariance MATRIX  [1 x dimension] x # of mixtures
  vector <MATRIX> m_Sigma; 

public:
    
 CGMM(const int _nummixtures, const int _dimension, const int _covtype)
   : m_nummixtures(_nummixtures), m_dimension(_dimension), 
     m_covtype(_covtype), m_trace(0)
  { 
    malloc(_nummixtures, _dimension, m_covtype); 
  };
    
  ~CGMM()
    {};
    
private:    
  void malloc(const int _nummixtures, const int _dimension, const int _covtype);
  void free();

public:
  int getNumMixtures() { return m_nummixtures;}
  int getDimension()   { return m_dimension;}
  int getCovType()     { return m_covtype;}

  SCALAR getPi(const int _i)
  {
    if (_i < 0 || _i >m_nummixtures)
      Error(1111, "invalid index of mixtures¥n");
    return m_pi.at(_i);
  }
    
  SCALAR getConst(const int _i)
  {
    if (_i < 0 || _i >m_nummixtures)
      Error(1111, "invalid index of mixtures¥n");
    return m_const.at(_i);
  }
    
  const MATRIX getMu(const int _i)
  {
    if (_i < 0 || _i >m_nummixtures)
      Error(1111, "invalid index of mixtures¥n");
    return m_mu.at(_i);
  }
  
  const MATRIX getSigma(const int _i)
  {
    if (_i < 0 || _i >m_nummixtures)
      Error(1111, "invalid index of mixtures¥n");
    return m_Sigma.at(_i);
  }

  void setTrace(const int _i) { m_trace = _i;};
    
    
public:

  /// Initialize GMM 
  void init(vector<const MATRIX*>& _data);

  /// ML estimation with EM algorithm
  void runEM(vector<const MATRIX*> & _data, const int _numIter, const int _minDiff);

  /// Calcurate logarithmic likelihood
  SCALAR calcLogLikelihood(const MATRIX& _data);
  SCALAR calcLogLikelihood(vector<const MATRIX*>& _data);

  /// Show information of GMM
  string printInfo();


 private:

  VECTOR calcLogLikelihood(const MATRIX& _data, const int _j);
  
};


#endif
