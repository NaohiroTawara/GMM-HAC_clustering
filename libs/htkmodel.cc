#include "htkmodel.h"

using namespace std;

namespace htktools {
  
void CHTKModel::readHeader()
{
  string token;
  int i_value;
  
  /* 1 mixture ��GMM���ɤ߹�����  */
  m_ifs >> token;
  checkFormat(token, "~o");
  m_ifs >> token;
  checkFormat(token, "<STREAMINFO>");/* <STREAMINFO> ���ɤ����Ф� */
  m_ifs >> token;
  m_ifs >> token;
    
  m_ifs >> token;/* �����ɤ߹��� */
  checkFormat(token, "<VECSIZE>");
  m_ifs >> m_dimension;
    
  m_ifs >> token; /* ʬ���μ�����ɤ߹��� */
  if (token == "<FULLC>")
    Error(1111, "readHeader: FULL covariance is not supported");
    //m_cov_type = FULLC;
  else
    m_cov_type = DIAGC;
  m_ifs >> token;
  
  checkFormat(token, "~h");
  m_ifs >> token;
  token.erase(0, 1);
  token.erase(token.size() - 1, 1); // �������ä�
  m_name = token; /* ���饹̾�Τ����� */
  m_ifs >> token;
  checkFormat(token, "<BEGINHMM>");
  m_ifs >> token;
  checkFormat(token, "<NUMSTATES>");
  m_ifs >> i_value;
  m_ifs >> token;
  checkFormat(token, "<STATE>");
  m_ifs >> i_value;

  ifstream::pos_type beg = m_ifs.tellg();
  m_ifs >> token;
  if (token == "<NUMMIXES>")
  {
    checkFormat(token, "<NUMMIXES>");
    m_ifs >> m_num_mixtures; /* ������ɤ߹���*/
    //    m_ifs >> token;
  }
  else
  {
    m_num_mixtures = 1;
    m_ifs.seekg(beg);
  }
}

int CHTKModel::readGaussian(double* _weight, CvMat* _mean, CvMat* _Cov)
{
  double tmp;
  return readGaussian(_weight, _mean, _Cov, &tmp);
}
  
int CHTKModel::readGaussian(double* _weight, CvMat* _mean, 
			    CvMat* _Cov,     double* _gconst)
{
  int mean_dim = cvGetSize(_mean).width;
  int cov_width  = cvGetSize(_Cov).width;
//  int cov_height = cvGetSize(_Cov).height;

  if (mean_dim != m_dimension)
    Error(1111, 
	  "readGaussian: Dimension is invalid: [1 x %d], Required [1 x %d]",
          mean_dim, m_dimension);
  if (cov_width!= m_dimension)
    Error(1111, 
	  "readGaussian: Dimension is invalid: [1 x %d],  Required [1 x %d]",
          cov_width, m_dimension);
  
  string token;  // string ���ΰ���ѿ�
  float f_value; // float ���ΰ���ѿ�
  int i_value;   // int ���ΰ���ѿ�

  m_ifs >> token;

  if (token == "<TRANSP>")
    Error(1111, "readGaussian: The end of file");
  
  if (token == "<MIXTURE>")
  {
    /* 1 mixture ��GMM���ɤ߹����硤������̵�뤵��� */
    checkFormat(token, "<MIXTURE>");
    m_ifs >> i_value;
    if (i_value != (m_cnt + 1))
      Error(1111, "readGaussian: Error in reading <MIXTURE>");
    if (_weight != null)
      m_ifs >> (*_weight); /* ����Ť��ɤ߹��� */
    else
      m_ifs >> f_value;
    m_ifs >> token;
  }
  else
  {
    if (_weight != null)
      (*_weight) = 1.0;
  }
  /* ʿ�ѥ٥��ȥ��ɤ߹��� */
  checkFormat(token, "<MEAN>");

  m_ifs >> i_value;
  if (i_value != m_dimension)
    Error(1111, "readGaussian: Dimension is invalid in <MEAN>: %d Required %d",
          i_value, m_dimension);
  for (int d = 0; d < m_dimension; ++d)
  {
    m_ifs >> f_value;
    cvmSet(_mean, 0, d, f_value);
  }

  /* ��ʬ�������ɤ߹��� */
  m_ifs >> token;
  if (m_cov_type == DIAGC)
  {
    checkFormat(token, "<VARIANCE>");
    m_ifs >> i_value;
    if (i_value != m_dimension)
      Error(1111, "readGaussian: Dimension is invalid in <MEAN>: %d Required %d",
            i_value, m_dimension);

    for (int d = 0; d < m_dimension; ++d)
    {
      m_ifs >> f_value;
      cvmSet(_Cov, 0, d, f_value);
    }
  }
  else
    Error(1111, "readGaussian: Full covariant matrix is not supported");

  m_ifs >> token;

  checkFormat(token, "<GCONST>");
  m_ifs >> (*_gconst);
  
  ++m_cnt;
  
  return m_cnt;  
}

}

