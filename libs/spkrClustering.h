
/*
* ---------------------------------------------------------------------   
*   spkrClustering.h
*
*   < Author >��N.TAWARA  2011/11/24
*               �üԥ��饹����󥰤�Ԥ����֥������ȤΥ��󥿥ե�����
*
*   < Note > ����: ���дؿ����ѿ���ޤि�ḷ̩�ˤϥ��󥿥ե������ǤϤʤ�
*
*            �ѹ�����
*
*            ToDo:
*  
*            ���͡� 
*---------------------------------------------------------------------
*
*/



#ifndef __SPKRCLUSTERING_H__
#define __SPKRCLUSTERING_H__

#include "segment.h"

class ISpkrClustering
{
  
protected:
  
  ostream* const m_ros;
  
 private:
  int m_num_clusters;  //// # of clusters
  int m_dimension;     //// dimension
  int m_covtype;       //// type of covariant matrix

  CSegment* m_datass;  //// features for all segments
  

 public:

  // -----------------------------------------------------------------------
  // ���󥹥ȥ饯�����ǥ��󥹥ȥ饯��
  // -----------------------------------------------------------------------

 ISpkrClustering(const int _num_clusters, const int _covtype, ostream* const _ros)
   : m_dimension(0),
     m_num_clusters(_num_clusters),
     m_covtype(_covtype),
     m_datass(NULL),
     m_ros(_ros)
    {};

  ~ISpkrClustering()
    {};
  
 public:
  
  // -----------------------------------------------------------------------
  // �ȥåץ�٥�����
  // -----------------------------------------------------------------------

  /// (1) ��ħ�̥��֥������Ȥ򥻥åȤ���ؿ�
  void setFeature(CSegment* _segment)
  { 
    m_datass = _segment; 
    m_dimension = _segment->getDimension();
  }

  /// (2) �������Ԥ��ؿ�
  virtual void init() = 0;
  
  /// (3) ���饹����󥰼¹Դؿ�
  virtual int run() = 0;
 
  /// (4-1) ���饹����󥰷�̤��֤��ؿ��ʳƥ������Ȥ��ɤΥ��饹���˳����Ƥ�줿����
  virtual int getClusteringResultSeg(IntVect* _segment) = 0;

  /// (4-2) ���饹����󥰷�̤��֤��ؿ��ʳƥ��饹���ˤɤΥ������Ȥ������Ƥ�줿����
  virtual int getClusteringResultCl(vector<IntVect>* _cluster) = 0;

//----------------------------------------------------------------------

 public:

  // -----------------------------------------------------------------------
  // ��������
  // -----------------------------------------------------------------------

    /// �üԥ��饹�������֤��ؿ�
  int getNumClusters() const
  { return m_num_clusters; }
    
  /// ��ħ�̤μ��������֤��ؿ�
  int getDimension(void) const
  { return m_dimension; }
  
  /// ȯ��������֤��ؿ�
  int getNumSegments() const
  { return m_datass->getNumSegments(); }

  /// ���ꤵ�줿ȯ�äΥե졼������֤��ؿ�
  int getNumFrames(const int _i)
  { return m_datass->getNumFrames(_i); }

  /// ��ե졼������֤��ؿ�
  int getNumAllFrames() const
  { return m_datass->getNumAllFrames();}

  /// ��ʬ������Υ����פ��֤��ؿ�
  int getCovType() const
  { return m_covtype; }

  /// ���ꤵ�줿ȯ�ä�ȯ�ý�����֤��ؿ�
  const MATRIX* getSegment(const int _i)
  { return m_datass->getSegment(_i); }

  /// ���ꤵ�줿ȯ�äλ��ꤵ�줿�ե졼����֤��ؿ�
  const MATRIX getFrame(const int _i, const int _j)
  { return m_datass->getFrame(_i, _j);}

 protected:

  void setNumClusters(const int _numclusters)
  { m_num_clusters = _numclusters;}
   
  IntVect getNumFramesVector()
  { return m_datass->getNumFramesVector(); }
  
};



#endif

