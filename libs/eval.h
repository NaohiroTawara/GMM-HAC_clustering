/* ------------------------------------------------------
* 
* eval.h
* 
* < Author > N. TAWARA 2012/2/17
*
* < Note> ���饹��������ٻ��Х��饹�Υ��󥿡��ե��������饹
*
*
*
* ------------------------------------------------------
*/

#ifndef _SPKR_CL_EVALUATION__H_
#define _SPKR_CL_EVALUATION__H_

// Spekaer Diarization�η�̤�ɾ���򤹤륯�饹
  
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "util.h"
using namespace std;

//#define square(X) (X)*(X)
/* ɾ����̤��ݻ����빽¤�� */
struct Result
{
  map<string, double>result;
};

/* �������Ȥ˴ؤ��������ݻ����빽¤�� */
struct MlfInfo
{
  string speaker; // �ü�̾
  int seg_id;     // ���������ֹ�
  double length;  // ȯ��Ĺ
  MlfInfo(string _speaker, int _id, double _length)
  : speaker(_speaker), seg_id(_id) , length(_length) {};
  ~MlfInfo(void){};
};

typedef vector<MlfInfo*> MlfVect;

class ISpkrClEvaluation
{
private:
  int m_num_data;       // �������ȿ�
  int m_num_class;      // ���饹����
  int m_true_num_spk;   // �����üԿ�

  MlfVect m_mlf_vec; // �������Ⱦ�������
  map<string, int> m_spk_map; //�ü�̾���ü��ֹ���Ϣ�դ��륳��ƥ�
public:

  // ���󥹥ȥ饯��
  ISpkrClEvaluation(const int _num_data, 
		    const int _num_class, 
		    const int _true_num_spk)
    : m_num_data(_num_data), 
      m_num_class(_num_class),
      m_true_num_spk(_true_num_spk) 
      {} ;
  
  // �ǥ��󥹥ȥ饯��
  ~ISpkrClEvaluation()
    {} 
  
  // ɾ��
  virtual Result evaluate(IntVect& _datacc) = 0;

  // MLF�ե������ɤ߹���
  int setMlf(vector<string> _spkr_label, DoubleVect _length_vec);

  // ���饹��������
  void setNumClass(const int _n) { m_num_class = _n; }


  double getFullLength()   
  {
    double full_length = 0;
    MlfVect::iterator iter_m = m_mlf_vec.begin();
    for (;iter_m != m_mlf_vec.end(); ++iter_m)
      full_length += (*iter_m)->length;
    return full_length;
  }

  int getNumData()         { return m_num_data;  }
  int getNumClass()        { return m_num_class; }
  int getNumTrueSpeakers() { return m_true_num_spk; }
  MlfVect* getMlfVect()    { return &m_mlf_vec;}
  map<string, int>* getSpkrMap() { return &m_spk_map;}

 protected:

  int checkDataCc(IntVect& _datacc);
};

#endif

