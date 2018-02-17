//
//  bayesianGmm.cc
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 12/08/01.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#include <iostream>
#include "bayesianGmm.h"

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

void CBayesianGMM::malloc(const int _nummixtures, const int _dimension, const int _covtype)
{
    vector<CvMat*>::iterator iter_mu    = m_mu0.begin();
    vector<CvMat*>::iterator iter_Sigma = m_Sigma0.begin();
    for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
    {
        if (*iter_mu) cvReleaseMat(&(*iter_mu));
        if (*iter_Sigma) cvReleaseMat(&(*iter_Sigma));
        *iter_mu = cvCreateMat(1, _dimension, CV_32F);
        *iter_Sigma = (_covtype == FULLC) ?
        cvCreateMat(_dimension, _dimension, CV_32F):
        cvCreateMat(1, _dimension, CV_32F);
    }
    m_pi0.resize(_nummixtures);
    m_xi0.resize(_nummixtures);
    m_eta0.resize(_nummixtures);
    m_mu0.resize(_nummixtures);
    m_Sigma0.resize(_nummixtures);
}

//----------------------------------------------------------------------

void CBayesianGMM::free()
{
    vector<CvMat*>::iterator iter_mu    = m_mu0.begin();
    vector<CvMat*>::iterator iter_Sigma = m_Sigma0.begin();
    for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
    {
        if (*iter_mu) cvReleaseMat(&(*iter_mu));
        if (*iter_Sigma) cvReleaseMat(&(*iter_Sigma));
    }  
}

//----------------------------------------------------------------------

void CBayesianGMM::readPriorFromHtkFormat(const char* _filename)
{
    htktools::CHTKModel* htkmodel
    = new htktools::CHTKModel(_filename);
    
    int dimension    = htkmodel->getNumDimension();
    int num_mixtures = htkmodel->getNumMixtures();
    int covtype      = htkmodel->getCovType();
    string name      = htkmodel->getName();
    
    
    if (dimension != m_dimension)
        Error(1111, "Dimension does not match, expect %d, have %d", m_dimension, dimension);
    if (num_mixtures != m_nummixtures)
        Error(1111, "Num of mixtures does not match, expect %d, have %d", m_nummixtures, num_mixtures);
    if (covtype != m_covtype)
        Error(1111, "Covtype does not match, expect %d, have %d", m_covtype, covtype);
    
    for (int i = 0; i < num_mixtures; ++i)
        htkmodel->readGaussian((&m_pi0[i]), m_mu0[i], m_Sigma0[i]);

    delete htkmodel;
}

//----------------------------------------------------------------------

