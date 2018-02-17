/*
* ---------------------------------------------------------------------   
*   util_mat.h
*
*   < Author >��N.TAWARA  2010/8/4
*         ����������ʹ���ؿ���
*
*   < Note > ����������ѥ������-D"_TR"���ץ�����̵ͭ�ˤ�����ϥǡ��������ž�֤��ڤ��ؤ�����
*                  �դ���  �� X=[x_1 | x_2 | ...]^t  [num * dim]
*                  �Ĥ��ʤ���X=[x_1 | x_2 | ...]     [dim * num]
*                  ���������Ū�ˤϽĥ٥��ȥ�򲣤��¤٤��ۤ����������ɤ��������ꥢ��������Ψ�ϲ��٥��ȥ��Ĥ��¤٤������ɤ��Ϥ�
*            �ѹ�����
*
*            ToDo:
*  
*            ���͡�
*---------------------------------------------------------------------
*
*/

#ifndef __UTIL_MAT_H__
#define __UTIL_MAT_H__

#include <cv.h>
#include <highgui.h>
#include <iostream>
#include <vector>

#include "util.h"
using namespace std;

typedef vector<CvMat*> CvMatVect;



/**
 * @brief �������Ȥ���Ϳ����줿�٥��ȥ�ˤĤ��ơ���������Ȥ���Ϳ����줿�����顼�ͤ�
 *        ���פ������ǤΥ���ǥå�����٥��ȥ�������֤�
 *
 *
 */
template <class T> inline IntVect _findFactors(CvMat* _vect, T _sch)
{
  IntVect index(0);
  int vectsize = cvGetSize(_vect).width;
  for (int i = 0; i < vectsize; ++i)
    if ( cvmGet(_vect, 0, i) == _sch)
      index.push_back(i);

  return index;
}

/**
 * @brief Ϳ����줿��Ĥι�������Ǥ�Ʊ�����ɤ���Ĵ�٤�
 *
 */
inline int _isIdenticalMatrix(const CvMat* _mat1, const CvMat* _mat2)
{
  int cols1 = cvGetSize(_mat1).width;
  int rows1 = cvGetSize(_mat1).height;
  int cols2 = cvGetSize(_mat2).width;
  int rows2 = cvGetSize(_mat2).height;
  
  if ((cols1 != cols2) || (rows1 != rows2))
    return 0;

  for (int j = 0; j < rows1; ++j) 
    for (int i = 0; i < cols1; ++i)
      if (cvmGet(_mat1, j, i) != cvmGet(_mat2, j, i))
	return 0;

  return 1;
}

/**
 * @brief �������Ȥ���Ϥ���
 * @param _mat  ��Ȥ�ɽ������������
 */
inline void _print_matrix(const CvMat* _mat)
{
  int cols = cvGetSize(_mat).width;
  int rows = cvGetSize(_mat).height;
  cout <<" ["<<rows<<" x "<<cols<<"] = "<<endl;
  for (int j = 0; j < rows; ++j){
    cout << "  ";
    for (int i = 0; i < cols; ++i)
      cout << cvmGet(_mat, j,i) << ",";
    cout << endl;
  }
}
/**
 * @brief �������Ȥ���Ϥ���
 * @param _os  ���ϥ��ȥ꡼��
 * @param _mat ��Ȥ�ɽ������������
 */
inline void _print_matrix(ostream & _os, CvMat* _mat)
{
  int cols = cvGetSize(_mat).width;
  int rows = cvGetSize(_mat).height;
  for (int j = 0; j<rows; ++j)
  {
    for (int i = 0; i<cols; ++i)
      _os << cvmGet(_mat, j,i) << ",";
    _os << endl;
  }
}

template <class T> inline void _print_matrix(vector <T> vec)
{
  int cols = vec.size();
  for (int i=0; i<cols;++i)
    cout << vec[i] << ",";
  cout << endl;
}


inline void _getMeanMat(CvMat* _src_mat, CvMat* _dst_mat)
{
  int num_samples = cvGetSize(_src_mat).height;
  
  cvReduce(_src_mat, _dst_mat, CV_REDUCE_SUM);
  cvConvertScale(_dst_mat, _dst_mat, 1.0 / static_cast<double>(num_samples));
}

inline void _getCentrizedMat(CvMat* _src_mat, CvMat* _mean_vec, CvMat* _dst_mat)
{
  int dim = cvGetSize(_src_mat).width;
  int num_samples = cvGetSize(_src_mat).height;
  for (int f = 0; f < num_samples; ++f)
    for (int d = 0; d < dim; ++d)
      cvmSet(_dst_mat, f, d, cvmGet(_src_mat, f, d) - cvmGet(_mean_vec, 0, d));
}

/**
 * @brief ��ʬ��������֤��ؿ�
   @param _src_mat ���Ϲ���
   @param _src_vec ���Ϲ����ʿ�ѥ٥��ȥ�
   @param _dst ���Ϲ���
   @memo  ����ʬ��
  */
inline void _getCovMat(CvMat* _src_mat, CvMat* _mean_vec, CvMat* _dst_mat, const int _covtype)
{
  int dim = cvGetSize(_src_mat).width;
  int num_samples = cvGetSize(_src_mat).height;

  CvMat* rowvec = cvCreateMat(1, dim, CV_32F);
  CvMat* colvec = (_covtype == FULLC) ? 
    cvCreateMat(dim, 1, CV_32F) : /* full */
    NULL;                         /* diagonal */
  CvMat* matrix = (_covtype == FULLC ) ? 
    cvCreateMat(dim, dim, CV_32F) : /* full */
    NULL;                           /* diagonal */
  
  CvMat* cent_mat = cvCreateMat(num_samples, dim, CV_32F);
  _getCentrizedMat(_src_mat, _mean_vec, cent_mat);
  
  // ��ʬ������λ���
  if (_covtype == FULLC) 
  { // full covariance matix
    CvMat* trans_cent_mat =  cvCreateMat(dim, num_samples, CV_32F);
    cvmTranspose(cent_mat, trans_cent_mat);
    cvMatMul(trans_cent_mat, cent_mat, _dst_mat);
    cvReleaseMat(&trans_cent_mat);
  } 
  else 
  {                  // diag covariance matix
    cvPow(cent_mat, cent_mat, 2);
    cvReduce(cent_mat, _dst_mat, CV_REDUCE_SUM);
  }
  cvConvertScale(_dst_mat, _dst_mat, 1.0 / static_cast<double>(num_samples - 1));

  if (rowvec) cvReleaseMat(&rowvec);
  if (colvec) cvReleaseMat(&colvec);
  if (matrix) cvReleaseMat(&matrix);
  cvReleaseMat(&cent_mat);
}


/**
 * @brief  �٥��ȥ�����Ǥι���ͤ��֤�
 * @param  _vec �٥��ȥ�
 */
template <class T>  inline T _sum(vector <T> _vec)
{
  T sum = 0;
  int num = _vec.size();
  for (int i=0;i < num; ++i)
    sum += _vec[i];
  return sum;
}
/**
 * @brief  �٥��ȥ�����Ǥκ����ͤ��֤�
 * @param  _vec �٥��ȥ�
 */
template <class T>  inline T _max(vector <T> _vec)
{
  T max = (T)(-RAND_MAX);
  int num = _vec.size();
  for (int i=0;i < num; ++i)
    if (max < _vec[i]) max = _vec[i];
  return max;
}

/**
 * @brief  �����������
   @param _src_mat ���Ϲ���
   @param _dst ���Ϲ���
  */
inline void _AbsMat(CvMat* _src_mat)
{
  int dim = cvGetSize(_src_mat).width;
  int num_samples = cvGetSize(_src_mat).height;
  for (int i=0; i<dim;++i)
    for (int j=0; j<num_samples;++j)
      cvmSet(_src_mat, i, j, _abs(cvmGet(_src_mat,i,j)));
  
}

inline string _mat2str(CvMat* _mat)
{
  ostringstream os;
  int cols = cvGetSize(_mat).width;
  int rows = cvGetSize(_mat).height;
  os <<" ["<<rows<<" x "<<cols<<"] = "<<endl;
  for (int j = 0; j < rows; ++j){
    os << "  ";
    for (int i = 0; i < cols; ++i)
      os << cvmGet(_mat, j,i) << ",";
    os << endl;
  }
  return os.str();
}


#endif
