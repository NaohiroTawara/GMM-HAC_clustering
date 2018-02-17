//
//  gmm.h
//  spkr_cl_GMM-HAC
//
//  GMM のパラメタ管理と EM アルゴリズムによるパラメタの最尤推定を行うクラス
//  Created by 直弘 俵 on 12/07/31.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#ifndef spkr_cl_GMM_HAC_gmm_h
#define spkr_cl_GMM_HAC_gmm_h

#include "util.h"
#include "util_mat.h"
#include "htkmodel.h"
#include <cv.h>
#include <list>

using namespace std;

inline double mahala(CvMat* _x, const CvMat* _mu,
                     CvMat* _sigma, int _cov_type)
{
    double f = 0.0;
    int dim  = cvGetSize(_x).width;
    if (_cov_type == DIAGC)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int d = 0; d < dim; ++d)
        {
            float x = cvmGet(_x, 0, d);
            float u = cvmGet(_mu, 0, d);
            float s = cvmGet(_sigma, 0, d);
            //float x = _x->data.fl[d];
            //                float s = _sigma->data.fl[d];
            //float u = _mu->data.fl[d];
            f += -square(x - u) /(2*s);
        }
    }
    else
    {
    }    
    return f;
}
inline double mahala(float* _x, float* _mu, float* _sigma, int _cov_type, int _dim)
{
    double f = 0.0;
    if (_cov_type == DIAGC)
    {
        for (int d = 0; d < _dim; ++d)
        {
            float x = *(_x+d);
            float u = *(_mu + d);
            float s = *(_sigma + d);
            f += -square(x - u) /(2*s);
        }
    }
    else
    {
    }
    return f;
}

class CGMM
{
public:
    int m_nummixtures;  // 混合数
    int m_dimension;    // 特徴量の次元
    int m_covtype;      // 共分散行列の種類
    int m_numData;      // データ数
    int m_trace;        // 標準出力の有無
    
    DoubleVect m_const;      // 正規化項   [1 x # of mixtures]
    DoubleVect m_pi;         // 混合重み   [1 x # of mixtures]
    vector <CvMat*> m_mu;    // 平均      [1 x dimension] x # of mixtures
    vector <CvMat*> m_Sigma; // 共分散行列 [dimension x dimension] x # of mixtures

public:
    
    CGMM(const int _nummixtures, const int _dimension, const int _covtype)
    : m_nummixtures(_nummixtures), m_dimension(_dimension), m_covtype(_covtype), m_trace(0)
    { 
        malloc(_nummixtures, _dimension, m_covtype); 
    };
    
    ~CGMM()
    {};
    
private:    
    void malloc(const int _nummixtures, const int _dimension, const int _covtype);
    void free();

public:
    int getNumMixtures() { return m_nummixtures;}
    int getDimension()   { return m_dimension;}
    int getCovType()     { return m_covtype;}
    
    double getPi(const int _i)
     {
         if (_i < 0 || _i >m_nummixtures)
             Error(1111, "invalid index of mixtures¥n");
         return m_pi.at(_i);
     }
    
    double getConst(const int _i)
    {
        if (_i < 0 || _i >m_nummixtures)
            Error(1111, "invalid index of mixtures¥n");
        return m_const.at(_i);
    }
    
    CvMat* getMu(const int _i)
    {
        if (_i < 0 || _i >m_nummixtures)
            Error(1111, "invalid index of mixtures¥n");
        return m_mu.at(_i);
    }
    
    CvMat* getSigma(const int _i)
    {
        if (_i < 0 || _i >m_nummixtures)
            Error(1111, "invalid index of mixtures¥n");
        return m_Sigma.at(_i);
    }
    
    void setTrace(const int _i) { m_trace = _i;};
    
    
public:
    
    /// GMM を初期化する
    void init(vector<CvMat*> _data);

    /// EM アルゴリズムにより分布を最尤推定する
    void runEM(vector<CvMat*> _data, const int _numIter, const int _minDiff);
    
    string printInfo();

    double calcLogLikelihood(CvMat* _data);
    double calcLogLikelihood(float* _data);

    /// 他の GMM との KL Divergence (Sigma point approximation) を計算する
    double calcKL(CGMM* _gmm);
private:
    
    double calcLogLikelihood(CvMat* _data, const int _j);
    double calcLogLikelihood(float* _data, const int _j);
    double calcLogLikelihood(vector<CvMat*> _data);


};


#endif
