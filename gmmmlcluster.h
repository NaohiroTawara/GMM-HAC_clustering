//
//  gmmMLcluster.h
//  spkr_cl_GMM-HAC
//
//  各クラスタに相当するクラス
//  Created by 直弘 俵 on 12/07/30.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#ifndef spkr_cl_GMM_HAC_gmmmlcluster_h
#define spkr_cl_GMM_HAC_gmmmlcluster_h

#include <iostream>

#include "cluster.h"
#include "matrix.h"
#include "gmm.h"
#include <list>

using namespace std;

class CGMMMLCluster: public CCluster
{
private:
    
  CGMM* m_gmm;       // GMM をモデルとして用いる
    
  int m_nummixtures; // 混合数
  int m_maxiter;     // EM の最大イタレーション   
  int m_mindiff;     // EM の最小更新量
public:
    
    /// コンストラクタ
    CGMMMLCluster(const int _covtype)
    : CCluster(_covtype), m_nummixtures(0), m_maxiter(0), m_mindiff(0)
    {};
    
    ~CGMMMLCluster() 
    { free(); };
     
public:
    
    void initCluster(vector <const MATRIX*>& _feature, 
		     vector<int>& _seg_id,
                     const int _nummixtures,
                     const int _maxiter,
                     const int _mindiff);
    
    /// クラスタ間距離 (CLR) を定義
    double calcDistance(CGMMMLCluster* _cl);
    
    /// クラスタ間のマージ処理を定義
    void mergeCluster(CGMMMLCluster* _cl);
    
    void setTrace(const int _t) { m_gmm->setTrace(_t);}
    
    int getNumMixtures() { return m_nummixtures;}
    
    string printInfo();
    
    double calcLogLikelihood(const MATRIX& _data);

    void updateCluster();

private:
    
    void malloc(const int _nummixtures, const int _dimension);
    void free();

};

#endif
