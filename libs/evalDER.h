/* ------------------------------------------------------
* 
* evalDER.h
* 
* < Author > N. TAWARA 2012/2/17
*
* < Note> Diarization Error Rate (DER) �򻻽Ф��뤿��Υ��饹
*         A* õ���ˤ�� DER ���Ǿ��Ȥʤ��üԤȥ��饹�����Ȥߤ�׻������λ���DER�ͤ��֤�   
*
*
* ------------------------------------------------------
*/
#ifndef _SPKR_CL_EVALUATION_DER__H_
#define _SPKR_CL_EVALUATION_DER__H_

#include "eval.h"
#include <list>
#include <algorithm>



using namespace std;

class CSpkrClEvaluation_DER: public ISpkrClEvaluation
{
public:

  // ���󥹥ȥ饯��


  CSpkrClEvaluation_DER(const int _num_data, 
			const int _num_class, 
			const int _true_num_spk)
    : m_min_DER(DBL_MAX),
    ISpkrClEvaluation(_num_data, _num_class, _true_num_spk)
  {} 
  
  // �ǥ��󥹥ȥ饯��
 ~CSpkrClEvaluation_DER()
    {} 
  
  // ɾ��
  Result evaluate(IntVect& _datacc);

  /*  
  // MLF�ե������ɤ߹���
  int readMlf(string _filename); 

  // ���饹��������
  void setNumClass(const int _n) { m_num_class = _n; }

  int getNumData()         { return m_num_data;  }
  int getNumClass()        { return m_num_class; }
  int getNumTrueSpeakers() { return m_true_num_spk; }
  */
 private:

  double m_min_DER;

  DoubleMat m_DER_matrix;

  DoubleVect m_spkr_utt_length;

  // �ü� i �� ���饹�� j �˳�����Ƥ��Ȥ��Υ��顼Ĺ���֤�
  double getErrorLength(const int _i, const int _j)
    { return 
	getNumClass() > getNumTrueSpeakers()
	? m_DER_matrix[_i][_j]
	: m_DER_matrix[_j][_i]; 
    }

  // �ü� i ��ȯ�ä����ƥ��顼�Ȥߤʤ����Ȥ��Υ��顼Ĺ���֤�
  double getErrorLength(const int _i)
  { return m_spkr_utt_length[_i];}

  // DER ���ФΤ���κƵ��ؿ�
  void calcDER(const int _i, const double _DER, list<int> pibot);

  void setDER(IntVect& _datacc);

  void sorting(const int _i, list<int>* pibot);

  struct pair_less {
    bool operator()(const pair<double, int>& x, 
		    const pair<double, int>& y) const 
    { return x.first < y.first; }
  };
};

#endif

