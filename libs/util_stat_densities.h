/*
** --------------------------------------------------
**  util_stat_densities.h
**
**  KL Divergence 系ユーティリティ関数群
**  copy right N.TAWARA
**             2011/6/23
** --------------------------------------------------
*/

#ifndef __UTIL_STAT_DENSITIES_H__
#define __UTIL_STAT_DENSITIES_H__

#include <cv.h>
#include <math.h>
#include "util_mat.h"
#include "util.h"




// -------------------------------------------------------
// １変量ガウス分布
inline double _KLGaussian(const double mu_q, const double Sigma_q,
                         const double mu_p, const double Sigma_p)
{
  return
    0.5 * (log(Sigma_p) -log(Sigma_q))
      + (mu_q * mu_q + mu_p * mu_p + Sigma_q - 2 * mu_q * mu_p)
          / ( 2 * Sigma_p)
      - 0.5;
}

// -------------------------------------------------------
// 多変量ガウス分布
static double _KLGaussian(CvMat* mu_q, CvMat* Sigma_q,
                         CvMat* mu_p, CvMat* Sigma_p,
                         const int covtype)
{
  int dimension = cvGetSize(mu_q).width;
  if (dimension != cvGetSize(mu_p).width)
  {
    cerr << "ERROR[util_stat_densities:KLGaussian()] dimension is not match "
	 << endl;
    exit(-1);
  }
  CvMat* rowvec  = cvCreateMat(1, dimension, CV_32F);
  CvMat* colvec  = cvCreateMat(dimension, 1, CV_32F);
  CvMat* matrix  = (covtype == FULLC) ?
    cvCreateMat(dimension, dimension, CV_32F) :   /* full */
    cvCreateMat(1, dimension, CV_32F);            /* diagonal */
  CvMat* scalar = (covtype == FULLC) ? 
    cvCreateMat(1, 1, CV_32F) : /* full */
    NULL;                       /* diagonal */

  double logdet_p = 0.0;
  double logdet_q = 0.0;
  double trace    = 0.0;
  double mahala   = 0.0;
  
  if (covtype == FULLC)
  {
    CvMat* Lambda_q = cvCreateMat(dimension, dimension, CV_32F); /* full */
    cvSVD(Sigma_q, rowvec);
    cvLog(rowvec, rowvec);

    logdet_q = cvSum(rowvec).val[0];
    cvSVD(Sigma_p, rowvec);
    cvLog(rowvec, rowvec);
    logdet_p = cvSum(rowvec).val[0];

    cvInvert(Sigma_q, Lambda_q);
    cvMatMul(Lambda_q, Sigma_p, matrix);
    trace = cvTrace(matrix).val[0];

    cvSub(mu_q, mu_p, rowvec);
    cvmTranspose(rowvec, colvec);

    cvMatMul(rowvec, Sigma_p, rowvec);
    cvMatMul(rowvec, colvec, scalar);

    mahala = cvmGet(scalar, 0, 0);
    cvReleaseMat(&Lambda_q);
  }
  else
  {
    cvLog(Sigma_q, rowvec);
    logdet_q = cvSum(rowvec).val[0];
    cvLog(Sigma_p, rowvec);
    logdet_p = cvSum(rowvec).val[0];

    cvSub(mu_q, mu_p, rowvec);
    cvPow(rowvec, rowvec, 2);
    for (int d = 0; d < dimension; ++d)
    {
      trace  += cvmGet(Sigma_p, 0, d) / cvmGet(Sigma_q, 0, d);
      mahala += cvmGet(rowvec, 0, d)  * cvmGet(Sigma_p, 0, d);
    }
  }

  cvReleaseMat(&rowvec);
  cvReleaseMat(&colvec);
  cvReleaseMat(&matrix);

  if (scalar) cvReleaseMat(&scalar);
  
  return
    0.5 * (logdet_q - logdet_p
           + trace  + mahala - dimension);

}

// -------------------------------------------------------
// Wishart分布
static double _KLWishart(const double _q, CvMat* __Q,
                         const double _p, CvMat* __P,
			 const int _covtype)
{
  int dimension = cvGetSize(__Q).width;
  if (dimension != cvGetSize(__P).width)
  {
    cerr << "ERROR[util_stat_densities:KLWishart()] dimension is not match "
	 << endl;
    exit(-1);
  }
  CvMat* rowvec  = cvCreateMat(1, dimension, CV_32F);
  double glog2   = dimension * log(2);
  //CvMat* inv_P = (_covtype == FULLC) ?
  //  cvCreateMat(dimension, dimension, CV_32F) :   /* full */
  //  cvCreateMat(1, dimension, CV_32F);            /* diagonal */ 
  //CvMat* inv_Q = (_covtype == FULLC) ?
  //  cvCreateMat(dimension, dimension, CV_32F) :   /* full */
  //  cvCreateMat(1, dimension, CV_32F);            /* diagonal */ 
  
  //  cvInvert(__P, inv_P);
  //  cvInvert(__Q, inv_Q);

  double trace = 0;
  if (_covtype == FULLC)
  {
    CvMat* matrix  = cvCreateMat(dimension, dimension, CV_32F);
    cvInvert(__Q, matrix);
    cvMatMul(matrix, __P, matrix);
    trace = cvTrace(matrix).val[0];
    cvReleaseMat(&matrix);
  }
  else
  {
    for (int i = 0; i < dimension; ++i)
      trace += cvmGet(__P, 0, i) / cvmGet(__Q, 0, i);
  }
  if (_covtype == FULLC)
  {
    cvSVD(__Q, rowvec);
    cvLog(rowvec, rowvec);
  }
  else
  {
    cvLog(__Q, rowvec);
  }
  double logdet_Q = cvSum(rowvec).val[0];
  if (_covtype == FULLC)
  {
    cvSVD(__P, rowvec);
    cvLog(rowvec, rowvec);
  }
  else
  {
    cvLog(__P, rowvec);
  }
  double logdet_P = cvSum(rowvec).val[0];
  
  double v = 0.25 * dimension * (dimension - 1 ) * logPI;
  double logGammad_q = v;
  double logGammad_p = v;
  for (int d = 0; d < dimension; ++d)
  {
    logGammad_q += lgamma(0.5 * (_q - d + 1));
    logGammad_p += lgamma(0.5 * (_p - d + 1));
  }
  /// 正規化項の算出
  double logZ_q = _q * 0.5 * (-logdet_Q) + logGammad_q;
  double logZ_p = _p * 0.5 * (-logdet_P) + logGammad_p;
  if (_covtype == FULLC)
  { // glog2 は対角の場合必要ない
    logZ_q += _q * 0.5 * glog2;
    logZ_p += _p * 0.5 * glog2;
  }
  /// エントロピー項の算出
  double L_q = -logdet_Q;
  double L_p = -logdet_P;
  for (int d = 0; d < dimension; ++d)
  {
    L_q += digamma((_q + 1 - d) * 0.5);
    L_p += digamma((_p + 1 - d) * 0.5);
  }
  if (_covtype == FULLC)
  { // glog2 は対角の場合必要ない
    L_q += glog2;
    L_p += glog2;
  }  
  //  cvReleaseMat(&inv_P);
  //  cvReleaseMat(&inv_Q);

  return
        0.5 * (_q - dimension - 1) * L_q
      - 0.5 * (_p - dimension - 1) * L_p
      - 0.5 * _q * dimension
      + 0.5 * _q * trace
      + logZ_p - logZ_q;
}

// -------------------------------------------------------
// ガンマ分布
static double _KLGamma(const double b_q, const double c_q,
		       const double b_p, const double c_p)
{
  double inv_c_q = 1.0 / c_q;
  double inv_c_p = 1.0 / c_p;
  return
    (b_q - 1) * digamma(b_q)
    - log(inv_c_q) - b_q - lgamma(b_q)
    + lgamma(b_p) + b_p * log(inv_c_p)
    - (b_p - 1) * (digamma(b_q) + log(inv_c_q))
    + inv_c_q * b_q / inv_c_p;
}

// -------------------------------------------------------
// ディリクレ分布
static double _KLDir(DoubleVect& l_q, DoubleVect& l_p)
{
  double sumloggamma_lq = 0.0;
  double sumloggamma_lp = 0.0;
  double sumdif         = 0.0;
  double sum_lq         = _sum(l_q);
  double sum_lp         = _sum(l_p);
  DoubleVect::iterator iter_lq = l_q.begin();
  DoubleVect::iterator iter_lp = l_p.begin();

  for (; iter_lq != l_q.end(); ++iter_lq, ++iter_lp)
  {
    sumloggamma_lq += lgamma(*iter_lq);
    sumloggamma_lp += lgamma(*iter_lp);

    sumdif += (*iter_lq - *iter_lp) *
      (digamma(*iter_lq) - digamma(sum_lq));
  }
  
  return
        lgamma(sum_lq) - lgamma(sum_lp)
      + sumloggamma_lp - sumloggamma_lq
      + sumdif;
}

#endif
