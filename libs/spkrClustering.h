
/*
* ---------------------------------------------------------------------   
*   spkrClustering.h
*
*   < Author >　N.TAWARA  2011/11/24
*               話者クラスタリングを行うオブジェクトのインタフェース
*
*   < Note > 現状: メンバ関数・変数を含むため厳密にはインタフェースではない
*
*            変更点：
*
*            ToDo:
*  
*            備考： 
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
  // コンストラクタ・デコンストラクタ
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
  // トップレベル制御
  // -----------------------------------------------------------------------

  /// (1) 特徴量オブジェクトをセットする関数
  void setFeature(CSegment* _segment)
  { 
    m_datass = _segment; 
    m_dimension = _segment->getDimension();
  }

  /// (2) 初期化を行う関数
  virtual void init() = 0;
  
  /// (3) クラスタリング実行関数
  virtual int run() = 0;
 
  /// (4-1) クラスタリング結果を返す関数（各セグメントがどのクラスタに割当てられたか）
  virtual int getClusteringResultSeg(IntVect* _segment) = 0;

  /// (4-2) クラスタリング結果を返す関数（各クラスタにどのセグメントが割当てられたか）
  virtual int getClusteringResultCl(vector<IntVect>* _cluster) = 0;

//----------------------------------------------------------------------

 public:

  // -----------------------------------------------------------------------
  // アクセサ
  // -----------------------------------------------------------------------

    /// 話者クラスタ数を返す関数
  int getNumClusters() const
  { return m_num_clusters; }
    
  /// 特徴量の次元数を返す関数
  int getDimension(void) const
  { return m_dimension; }
  
  /// 発話総数を返す関数
  int getNumSegments() const
  { return m_datass->getNumSegments(); }

  /// 指定された発話のフレーム数を返す関数
  int getNumFrames(const int _i)
  { return m_datass->getNumFrames(_i); }

  /// 総フレーム数を返す関数
  int getNumAllFrames() const
  { return m_datass->getNumAllFrames();}

  /// 共分散行列のタイプを返す関数
  int getCovType() const
  { return m_covtype; }

  /// 指定された発話の発話集合を返す関数
  const MATRIX* getSegment(const int _i)
  { return m_datass->getSegment(_i); }

  /// 指定された発話の指定されたフレームを返す関数
  const MATRIX getFrame(const int _i, const int _j)
  { return m_datass->getFrame(_i, _j);}

 protected:

  void setNumClusters(const int _numclusters)
  { m_num_clusters = _numclusters;}
   
  IntVect getNumFramesVector()
  { return m_datass->getNumFramesVector(); }
  
};



#endif

