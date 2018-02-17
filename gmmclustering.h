//
//  gmm-hac_clustering.h
//  spkr_cl_GMM-HAC
//
//  クラスタリング処理を管理するクラス
//  Created by 直弘 俵 on 12/07/30.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#ifndef spkr_cl_GMM_HAC_gmm_hac_clustering_h
#define spkr_cl_GMM_HAC_gmm_hac_clustering_h

#include "spkrClustering.h"
#include "gmmmlcluster.h"     // GMM で各発話をモデル化する
#include "segment.h"
#include "matrix.h"

#define CL_TYPE CGMMMLCluster // 各クラスタのモデル

using namespace std;

class CGMMClustering: public ISpkrClustering
{
    
private:
    
  float m_alpha;         // マージを行う最小距離    
  list<CL_TYPE*> m_clusters; // クラスタ集合のリスト
  
  vector<DoubleVect> m_distmat;
  
  int last_updated_i;
  int last_updated_i_new;
  int updatedFlag;
    
public:
    
  /// コンストラクタ
 CGMMClustering(const int _init_numclusters,
		const int _covtype,
		ostream* const _ros)
   : ISpkrClustering(_init_numclusters, _covtype, _ros)
    { malloc(_init_numclusters, _covtype); };
  
  /// デストラクタ
  ~CGMMClustering()
    { free(); }

private:
    
  // メモリ管理関数
  void malloc(const int _numclusters, const int _covtype);
  
  void free();
  
  /// cl1 と cl2 をマージした時の距離を計算する関数
  double calcDistance(CL_TYPE* _cl1, CL_TYPE* _cl2);
    
  /// cl1 と cl2 をマージする関数
  void mergeClusters(CL_TYPE* _cl1, CL_TYPE* _cl2);
    
public:

  string getDistMat();

  /// パラメタを設定する関数
  void setAlpha(const float _alpha)
  { m_alpha = _alpha; }
    
  /// パラメタを返す関数
  float getAlpha()
  { return m_alpha; }
    
  /// クラスタの prior を設定する
  void setBasisFromFile(char* _filename);
  
  void setTrace(const int _t)
  {
    list<CGMMMLCluster*>::iterator iter_c = m_clusters.begin();
    for (; iter_c != m_clusters.end(); ++iter_c)
      (*iter_c)->setTrace(_t);
  }
    
  string showClusteringResult();
    
  /// クラスタの初期化を行う
  void initClusters(const int _nummixtures, 
		    const int _maxiter, 
		    const int _mindiff);

  void initDistMat();

  // -------------------------------------------------------
  // トップレベル制御関数（インタフェースで宣言済み）
  // -------------------------------------------------------
  
  /// (1) 初期化を行う関数
  void init();
  
  /// (2) クラスタリング実行関数
  int run();
  
  /// (3) クラスタリング結果を返す関数
  int getClusteringResultSeg(IntVect* _segment_cc);
  int getClusteringResultCl(vector<IntVect>* _cluster);
  
  string printInfo();
};

#endif
