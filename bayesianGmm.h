//
//  bayesianGmm.h
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 12/08/01.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#ifndef spkr_cl_GMM_HAC_bayesianGmm_h
#define spkr_cl_GMM_HAC_bayesianGmm_h

#include "gmm.h"
using namespace std;

class CBayesianGMM: public CGMM
{
private:
    DoubleVect m_pi0;         // Dirichlet のハイパーパラメタ
    DoubleVect m_xi0;         // Gaussian のハイパーパラメタ
    DoubleVect m_eta0;        // Gamma(Wishart) のハイパーパラメタ
    vector <CvMat*> m_mu0;    // 平均
    vector <CvMat*> m_Sigma0; // 共分散行列
    
public:
    CBayesianGMM(const int _numMixtures, const int _dimension, const int _covtype)
    : CGMM(_numMixtures, _dimension, _covtype), m_pi0(0), m_xi0(0), m_eta0(0), m_mu0(0), m_Sigma0(0)
    { malloc(_numMixtures, _dimension, _covtype); };
    
    ~CBayesianGMM()
    { free();};
    
    void readPriorFromHtkFormat(const char* _filename);
    
    /// =======================================
    /// オーバーライドメソッド
    
    /// GMM を初期化する
    void init(vector<CvMat*> _data);
    
    /// EM アルゴリズムにより分布を最尤推定する
    void runEM(vector<CvMat*> _data, const int _numIter, const int _minDiff);
    
private:    
    void malloc(const int _nummixtures, const int _dimension, const int _covtype);
    void free();
    
};

#endif
