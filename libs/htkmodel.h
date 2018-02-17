/*
* ---------------------------------------------------------------------   
*   htktools::htkmodel.h
*
*   < Author >　N.TAWARA  2011/8/3
*         htk 形式の GMM モデルファイルのパーサー
*
*
*---------------------------------------------------------------------
*
*/

#ifndef __CHTKMODELREADER_H__
#define __CHTKMODELREADER_H__

#include <fstream>
#include <cv.h>
#include <highgui.h>
#include "util.h"

using namespace std;

namespace htktools {

class CHTKModel
{
private:
  ifstream m_ifs;     // ファイルストリーム

  string m_name;
  int m_dimension;    // 次元
  int m_num_mixtures; // 混合数
  int m_cov_type;     // 共分散行列のタイプ { FULLC, DIAGC}

  int m_cnt;          // 最後に読み込んだ混合要素のインデクス

private:

  /* @brief htk モデルファイルのヘッダ（~o から <STATE> まで）を読み込む
   */
  void readHeader();
  
public:
  /*
   * @brief  コンストラクタ
   * @param _filename: htk モデルファイル名
   */
  CHTKModel(const char* _filename)
       : m_name(""),        m_dimension(0),
         m_num_mixtures(0), m_cov_type(0),
         m_cnt(0)
    { m_ifs.open(_filename); readHeader(); };
  
  ~CHTKModel()
    { m_ifs.close(); }
  
  /*
   * @brief: 次のガウス分布を読み込む
   * @param _weight: 混合重み
   * @param _mean:   平均ベクトル
   * @param _Cov:    共分散行列
   * @return: 読み込んだ混合要素のインデックス
   */
  int readGaussian(double* weight, CvMat* mean, CvMat* Cov);
  int readGaussian(double* weight, CvMat* mean, CvMat* Cov, 
		   double* _gconst);

  int getNumDimension() { return m_dimension; };
  int getNumMixtures()  { return m_num_mixtures; };
  int getCovType()      { return m_cov_type; };
  string getName()      { return m_name; };
  
};

}

#endif

