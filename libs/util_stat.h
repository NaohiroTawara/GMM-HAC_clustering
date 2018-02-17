/*
** --------------------------------------------------
**  util_stat.h
**
**  確率関係関数群
**  copy right N.TAWARA
**             2010/8/4
** --------------------------------------------------
*/

#ifndef __UTIL_STAT_H__
#define __UTIL_STAT_H__

#include "util_mat.h"
#include "util.h"
#include <cv.h>

#include <boost/random.hpp>

#define  RAND_TYPE  minstd_rand



class Rand
{
 private:
  
  boost::RAND_TYPE gen_;
  boost::uniform_real<> dst_;
  boost::variate_generator<boost::RAND_TYPE, boost::uniform_real<> > rand_;
  
 public:
  
 Rand(int _seed) : gen_(_seed), dst_(0, 1), rand_(gen_, dst_) {}
  
  boost::uniform_real<>::result_type operator()() { return rand_(); }
};


/*
 * @brief マハラノビス距離( (x-mu)^2/sigma )を算出
 * @param _x     入力データ [1 x dim] 
 * @param _mu    平均       [dim x 1]
 * @param _sigma 共分散行列 [dim x 1]
 * @memo  対角共分散行列専用
 */
inline double _mahala(const CvMat* _x, CvMat* _mu, CvMat* _sigma, int _cov_type)
{
  double f = 0.0;
  int dim  = cvGetSize(_x).width;
  if (_cov_type == DIAGC)
  {
    for (int d = 0; d < dim; ++d)
    {
      float x = cvmGet(_x, 0, d);
      float u = cvmGet(_mu, 0, d);
      f += -square(x - u) / (2 * cvmGet(_sigma, 0, d));
    }
  }
  else
  {
  }    
  return f;
}


// -------------------------------------------------------
// Student-t 分布の対数確率密度関数
inline double _distStudent(const CvMat* _x, CvMat*  _mu,
                           CvMat* _Sigma, const double _r)
{
  const int covtype      = cvGetSize(_Sigma).height > 1? FULLC: DIAGC;
  const int dimension    = cvGetSize(_Sigma).width;
  
  CvMat* colvec = cvCreateMat(1, dimension, CV_32F);
  double logdet; // _x に依存しない項
  double mahala = 0; // _x に依存する項

  if (covtype == FULLC)
  {
    logdet = lgamma((_r + dimension) * 0.5) - lgamma(_r * 0.5) - log(_r * PI) * dimension * 0.5;
    cvSVD(_Sigma, colvec);
    cvLog(colvec, colvec);
  }
  else
  {
    logdet =   dimension * (lgamma(_r + 0.5) - lgamma(_r))
             - dimension * 0.5 * log(PI * _r);
    cvLog(_Sigma, colvec);
  }
  logdet -= 0.5 * cvSum(colvec).val[0];

  if (covtype == DIAGC)
  {
    for (int d = 0; d < dimension; ++d)
    {
      float x = cvmGet(_x,  0, d);
      float u = cvmGet(_mu, 0, d);
//      f += - square(x - u) / (cvmGet(_Sigma, 0, d));
      mahala -=  (_r + 0.5) * log(1 + square(x - u)  / (cvmGet(_Sigma, 0, d) * _r));
    }
  }

  cvReleaseMat(&colvec);
  
  return logdet + mahala;
}



/**
 * @brief  Dirichlet分布からのサンプリング
 * @param  _vec 
 */
/*
static vector <float> _randdir(int _num_class, vector <float>& _alpha)
{
  vector <float> pi(_num_class);
  float sum = 0.0;
  for (int i = 0; i < _num_class; ++i) {
    pi[i] = gengam(1, _alpha[i]);
    sum += pi[i];
  }
  for (int i = 0; i < _num_class; ++i) {
    pi[i] = pi[i] / sum;
  }
  return pi;
}

*/
 /**
  * @brief  多項分布からのサンプリング
  * @param  _vec 合計が１になるように正規化されたパラメータ
 * @return サンプリング値. エラー時は-1
 */
 /*
inline int _randMult(vector <double> vec, Rand<boost::uniform_int<> > _rand)
{
  double u = _rand();

  cout <<u<<endl;
  double ss = 0;
  vector<double>::iterator iter = vec.begin();
  for (int i = 0; iter != vec.end(); ++iter, ++i)
  {
    ss += (*iter);
    if (u <= ss)
      return i;
  }
  return -1;
}
*/
 
 inline int _randMult(DoubleVect vec, Rand* _rand)
{
  double u = (*_rand)(); 
  double ss = 0;
  int num = vec.size();

  for (int i = 0; i< num; ++i) 
  {
    ss += vec[i];
    if (u <= ss)
      return i;
  }
  return num-1;
}

/**
 * @brief 順列をランダムに並び替える
 * @param _num  順列の要素数
 * @param _perm 並び替える順列．変更される
 */
inline void _randperm(IntVect& _perm, Rand* _rand)
{
  int num = _perm.size();
  for (int i = 0; i < num; ++i)
  {
    int n = static_cast<int>((*_rand)() * static_cast<double>(num));
    int tmp = _perm[i];
    _perm[i] = _perm[n];
    _perm[n] = tmp;
  }
}


/**
 * @brief 順列をランダムに並び替える
 * @param _num  順列の要素数
 * @param _perm 並び替える順列．変更される
 */
/*
inline void _randperm(vector <int>& _perm, boost::minstd_rand _gen)
{
  int num = _perm.size();
  boost::uniform_smallint<> dst(0, num - 1);
  boost::variate_generator< boost::minstd_rand&, boost::uniform_smallint<> > 
    rand( _gen, dst );
 

  for (int i=0; i<num; ++i)
  {
    int n = rand();
    int tmp = _perm[i];
    _perm[i] = _perm[n];
    _perm[n] = tmp;
  }
}
*/

/*
template <class T> inline int _randMult(vector <T> vec, boost::RAND_TYPE _gen)
{

#if 0
  #if 0
    boost::uniform_smallint<> dst(0, RAND_MAX);
    boost::variate_generator< boost::RAND_TYPE&, boost::uniform_smallint<> > 
      rand( _gen, dst );
    T u = (double)rand() / (double) RAND_MAX;
  #else
    //for (int i=0;i<2;++i){
      boost::uniform_real<> dst(0, 1);
      boost::variate_generator< boost::RAND_TYPE&, boost::uniform_real<> > 
	rand( _gen, dst );
      //cout <<"!"<<rand()<<endl;
      // }
    T u = rand();

  #endif
#else
  T u = (double)rand() / (double) RAND_MAX;
#endif
  //    cout << "u: " << u <<endl;
  T ss = 0;
  int num = vec.size();
  for (int i = 0; i< num; ++i) 
  {
    ss += vec[i];
    if (u <= ss)
      return i;
  }
  return -1;
}
*/
/*
template <class T> inline int _randMult(vector <T> vec, boost::RAND_TYPE* _gen)
{
  exit(3);
  boost::uniform_real<> dst(0, 1);
  boost::variate_generator< boost::RAND_TYPE&, boost::uniform_real<> > 
    rand( (*_gen), dst );
  
  T u = rand();
  cout << "u: "<<u<<endl;  
  T ss = 0;
  int num = vec.size();
  for (int i = 0; i< num; ++i) 
  {
    ss += vec[i];
    if (u <= ss)
      return i;
  }
  return -1;
}
  */ 



#endif
