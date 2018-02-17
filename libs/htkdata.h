/*
*
* ---------------------------------------------------------------------
*
*  htktools::htkdata.h
*  
* ---------------------------------------------------------------------  
*   < Author >��N.TAWARA  2011/08/04
*
*   < Note > �������ǡ����ɤ߹��߷ϴؿ�
*            �ѹ���: 2011/08/04: ����
*            ToDo: 
*            ���͡� 64bit�Ķ��Ǥ�ư���ʤ����� ->�����Ѥ�
*/

#ifndef __HTKDATA_H__
#define __HTKDATA_H__

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <sstream>
#include <math.h>

using namespace std;

namespace htktools{
 
// HTK �Υإå���������
const int HEADER_SIZE = 3;

// HTK��ħ�̥ǡ����Υإå���������Ǽ���빽¤��
struct _htkHeader {
  unsigned int nSamples;        // �ե졼���
  unsigned int sampPeriod;      // �ե졼�����
  unsigned short int sampSize;  // ���ե졼������ΥХ��ȿ�
  unsigned short int paramKind; // �ǡ����μ������ꤹ�륳�����ֹ�
};

class CHTKData
{
public:

  int m_dimension;
  int m_num_samples;
  vector<vector <float> > m_data;

public:

  CHTKData(const char* _filename)
    { if (readData(_filename)) exit(-1);};

  // �ɤ߹�����ե������ _m �� _n ���ܤ��ͤ��֤�
  float getData(const int _m, const int _n)
    { return m_data[_m][_n];}

  // �ɤ߹�����ե�����μ���������ˤ��֤�
  int getDimension()
    { return m_dimension;}

  // �ǡ��������֤�
  int getNumSamples()
    { return m_num_samples;}

  // HTK��ħ�̥ǡ����Υإå������ɤ߹���
  struct _htkHeader* readHtkHeader(float *tmp);



private:

  unsigned swapint(unsigned int);
  float swapfloat(float*);
  int getLength(ifstream&, const int, const int);

  // htk �ե�������ɤ߹���
  int readData(const char* _filename);
  
};

}


#endif
