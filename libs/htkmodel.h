/*
* ---------------------------------------------------------------------   
*   htktools::htkmodel.h
*
*   < Author >��N.TAWARA  2011/8/3
*         htk ������ GMM ��ǥ�ե�����Υѡ�����
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
  ifstream m_ifs;     // �ե����륹�ȥ꡼��

  string m_name;
  int m_dimension;    // ����
  int m_num_mixtures; // �����
  int m_cov_type;     // ��ʬ������Υ����� { FULLC, DIAGC}

  int m_cnt;          // �Ǹ���ɤ߹�����������ǤΥ���ǥ���

private:

  /* @brief htk ��ǥ�ե�����Υإå���~o ���� <STATE> �ޤǡˤ��ɤ߹���
   */
  void readHeader();
  
public:
  /*
   * @brief  ���󥹥ȥ饯��
   * @param _filename: htk ��ǥ�ե�����̾
   */
  CHTKModel(const char* _filename)
       : m_name(""),        m_dimension(0),
         m_num_mixtures(0), m_cov_type(0),
         m_cnt(0)
    { m_ifs.open(_filename); readHeader(); };
  
  ~CHTKModel()
    { m_ifs.close(); }
  
  /*
   * @brief: ���Υ�����ʬ�ۤ��ɤ߹���
   * @param _weight: ����Ť�
   * @param _mean:   ʿ�ѥ٥��ȥ�
   * @param _Cov:    ��ʬ������
   * @return: �ɤ߹�����������ǤΥ���ǥå���
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

