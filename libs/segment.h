/*
* ---------------------------------------------------------------------   
*   segment.h
*
*   < Author >��N.TAWARA  2014/7/31
*               ȯ�å��åȤȤ���ȯ�üԥ�٥��������륯�饹
*
*   < Note > ����: ȯ�åǡ����μ��Τ� CvMat �������ݻ�����
*
*            �ѹ�����2011/11/24 �ɤ߹��ߵ�ǽ������
*                   2014/05/09 diag�����Υ�����ץȥե�����˲ä���
*                              scp�����Υ�����ץȥե�������б�
*                   2014/07/31 �ǡ����η�����CvMat����MatrixXd (eigen)���ѹ�
*
*            ToDo:
*  
*            ���͡� 
*---------------------------------------------------------------------
*
*/



#ifndef __SEGMENT_H__
#define __SEGMENT_H__


#include <string.h>
#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include "fileList.h"
#include "htkdata.h"
#include "matrix.h"
#include "util.h"

class CSegment
{
private:
  int m_dimension;    // ��ħ�̤μ���
  int m_numsegments;  // ȯ�ÿ�
  int m_numallframes; // ��ե졼���

  vector<int> m_numframes;    // [1 x # segments]

  vector <MATRIX*> m_data;  // [# segments x [# frames x dimension] ]
  vector <string> m_speakerlabel; // ���򥯥饹������

public:
  // ========================================
  // ���󥹥ȥ饯�����ǥ��󥹥ȥ饯��
  // ========================================
  CSegment(const char* _filename)
    : m_dimension(0), 
      m_numsegments(0),
      m_numallframes(0),
      m_numframes(0),
      m_speakerlabel(0)
    {
      //       readFeatureFromScript(_filename);
       readFeatureFromSCP(_filename);
    }
  
  CSegment(int _dimension)
    : m_dimension(_dimension), 
      m_numsegments(0),
      m_numallframes(0),
      m_numframes(0),
      m_speakerlabel(0)
  {}

  ~CSegment()
    { 
      free();
    };

 private:
      
  void free();

  /*
   * @brief           : diag �����ե����뤫��ȯ�åǡ������ɤ߹���
   * @param _filename : ������ץȥե�����̾
   */
  void readFeatureFromScript(const char* _filename);

  /*
   * @brief            : scp �����ե����뤫��ȯ�åǡ������ɤ߹���
   * @param  _filename : ������ץȥե�����̾
   */
  void readFeatureFromSCP(const char* _filename);

  /*
   * @brief            : ���ꤵ�줿�ե������¸�߳�ǧ�ȥե����ޥåȥ����å�����ɤ߹���
   * @param  _filename : �ե�����̾
   * @exception        : �ե����뤬¸�ߤ��ʤ������㳰
   */
  MATRIX* readFeature(const char* _filename) throw (string);

  /*
   * @brief     : m_data �κǸ����� ȯ�� _ss ���ɲä���
   * @param _ss    : �ɲä���ȯ�ù���
   * @param _stime : ����frame
   * @param _etime : ��λframe
   * @param _s     : �üԥ�٥�
   */
  void addSegment(const MATRIX& _ss, const int _stime, const int _etime, 
		  const string _s);

  void addSegment(const MATRIX& _ss, string _s);

  //// �����ե������csv �����ˤ��� 1 ȯ���ɤ߹���
  MATRIX* readCSVFeature(const char* _filename);
  //// �����ե������htk �����ˤ��� 1 ȯ���ɤ߹���
  MATRIX* readHTKFeature(const char* _filename);

 public:

  //// ȯ�åǡ����ΰ�Τߤ��������ؿ�
  void freeSegment(const int _i);

 public:
  // ========================================
  // private ���Х�������
  // ========================================

  string getSpeakerLabel(const int _i)
  { return m_speakerlabel[_i];}
  
  int getDimension()
  { return m_dimension;}

  int getNumSegments()
  { return m_numsegments; }

  int getNumAllFrames()
  { return m_numallframes;}

  int getNumFrames(const int _i);

  /*
   * @brief   : ȯ���ֹ����ꤷ�ơ�ȯ�ý�������뤿��δؿ�
   * @param _i: ȯ���ֹ�
   * @return  : [# frames x m_dimension]
   */
  const MATRIX* getSegment(const int _i);
  
  /*
   * @brief    : ȯ���ֹ�ȥե졼���ֹ����ꤷ��ľ����ħ�̤����뤿��δؿ�
   * @param _i : ȯ���ֹ�
   * @param _j : �ե졼���ֹ�
   * @return   : [1 x m_dimension] �Υ������������ʬ����Υݥ���
   */
  MATRIX getFrame(const int _i, const int _j);

  /*
   * @brief     : m_data �� ��Ƭ���� _i ���ܤ�ȯ�� _ss ���ɲä���
   * @param _ss : �ɲä���ȯ��
   * @param _i  : �ɲä�����
   */
  //  void setSegment(CvMat* _ss, const int _i);
  
  vector<int> getNumFramesVector()
  { return m_numframes; }
  
  //  void compression(const int _dim);
  
  void dimensionSelection(const int _sid, const int _new_dim);
};

#endif
