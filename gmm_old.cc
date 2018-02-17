//
//  gmm.cc
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 13/07/31.
//  revised by 直弘 俵 on 14/07/14: KL divergence を追加．ただし動かず.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#include <iostream>
#include "gmm_old.h"

#include <cv.h>
#include <highgui.h>

void CGMM::malloc(const int _nummixtures, const int _dimension, const int _covtype)
{
    vector<CvMat*>::iterator iter_mu    = m_mu.begin();
    vector<CvMat*>::iterator iter_Sigma = m_Sigma.begin();
    for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
    {
        if (*iter_mu) cvReleaseMat(&(*iter_mu));
        if (*iter_Sigma) cvReleaseMat(&(*iter_Sigma));
    }
    
    m_pi.resize(_nummixtures);
    m_mu.resize(_nummixtures);
    m_Sigma.resize(_nummixtures);
    m_const.resize(_nummixtures);
    
    iter_mu    = m_mu.begin();
    iter_Sigma = m_Sigma.begin();
    for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
    {
        *iter_mu    = cvCreateMat(1, _dimension, CV_32F);
        *iter_Sigma = (_covtype == FULLC) ?
        cvCreateMat(_dimension, _dimension, CV_32F):
        cvCreateMat(1, _dimension, CV_32F);        
    }
}

//----------------------------------------------------------------------

void CGMM::free()
{
    vector<CvMat*>::iterator iter_mu    = m_mu.begin();
    vector<CvMat*>::iterator iter_Sigma = m_Sigma.begin();
    for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
    {
        if (*iter_mu) cvReleaseMat(&(*iter_mu));
        if (*iter_Sigma) cvReleaseMat(&(*iter_Sigma));
    }  
}

//----------------------------------------------------------------------
/* データをランダムに分割してパラメタを初期化する */
void CGMM::init(vector<CvMat*> _data)
{
    int numSegments  = static_cast<int>(_data.size());
    int numallframes = 0;

    IntVect stat(m_nummixtures, 0);
    CvMat* tmp_matrix = (m_covtype == FULLC) ?
        cvCreateMat(m_dimension, m_dimension, CV_32F):
        cvCreateMat(1, m_dimension, CV_32F);
    CvMat* tmp_colvec     = cvCreateMat(1, m_dimension, CV_32F);

    cvSetZero(tmp_matrix);
    cvSetZero(tmp_colvec);

    vector<CvMat*>::iterator iter_mu = m_mu.begin();
    vector<CvMat*>::iterator iter_Sg = m_Sigma.begin();
    DoubleVect::iterator     iter_pi = m_pi.begin();
    DoubleVect::iterator     iter_ct = m_const.begin();
    IntVect::iterator        iter_st = stat.begin();
    for (; iter_mu != m_mu.end()
         ; ++iter_mu, ++iter_Sg, ++iter_pi, ++iter_ct, ++iter_st)
    {
        cvSetZero(*iter_mu);
        cvSetZero(*iter_Sg);
        (*iter_pi) = 0;
        (*iter_st) = 0;
        (*iter_ct) = 0;
    }

    /* 初期割当てを決定 */
    vector<IntVect> assignment(numSegments);
    vector<CvMat*>::iterator iter_d = _data.begin();
    vector<IntVect>::iterator iter_a = assignment.begin();
    for (; iter_d != _data.end(); ++iter_d, ++iter_a)
    {
        iter_a->resize(cvGetSize(*iter_d).height);
        IntVect::iterator iter_ai = iter_a->begin();
        for (;  iter_ai != iter_a->end(); ++iter_ai)
            (*iter_ai) = rand() % m_nummixtures;
    }

    // パラメタの更新
    iter_d = _data.begin();
    iter_a = assignment.begin();
    for (; iter_d != _data.end(); ++iter_d, ++iter_a)
    {
        numallframes += static_cast<int>(cvGetSize(*iter_d).height);
        IntVect::iterator iter_ai = iter_a->begin();
        for (int i = 0; iter_ai != iter_a->end(); ++iter_ai, ++i)
        {
            CvMat stub;
            CvMat* submat = cvGetRow(*iter_d, &stub, i);
            int k = (*iter_ai);
            /* 零次統計量の更新 */
            stat[k] += 1;
            /* 一次統計量の更新 */
            cvAdd(submat, m_mu[k], m_mu[k]);
            /* 二次統計量の更新 */
            if (m_covtype == FULLC)
                cvMatMul(submat, submat, tmp_matrix);
            else
                cvPow(submat, tmp_matrix, 2);
            cvAdd(tmp_matrix, m_Sigma[k], m_Sigma[k]);
        }
    }

    iter_mu = m_mu.begin();
    iter_Sg = m_Sigma.begin();
    iter_pi = m_pi.begin();
    iter_ct = m_const.begin();
    iter_st = stat.begin();
    for (; iter_mu != m_mu.end()
         ; ++iter_mu, ++iter_Sg, ++iter_pi, ++iter_ct, ++iter_st)
    {

        (*iter_pi) = (*iter_st) / static_cast<double>(numallframes);
        cvConvertScale(*iter_mu, *iter_mu, 1.0 / static_cast<double> (*iter_st));
        cvConvertScale(*iter_Sg, *iter_Sg, 1.0 / static_cast<double> (*iter_st));
        if (m_covtype == FULLC)
            cvMatMul(*iter_mu, *iter_mu, tmp_matrix);
        else
            cvPow(*iter_mu, tmp_matrix, 2);
        cvSub(*iter_Sg, tmp_matrix, *iter_Sg);

        if (m_covtype == FULLC)
            (*iter_ct) = - 0.5 * (m_dimension * log2PI + log(cvDet(*iter_Sg)));
        else
        {
            cvLog(*iter_Sg, tmp_matrix);
            (*iter_ct) = - 0.5 * (m_dimension * log2PI + cvSum(tmp_matrix).val[0]);
        }
    }
    cvReleaseMat(&tmp_matrix);
    cvReleaseMat(&tmp_colvec);
}

//----------------------------------------------------------------------

void CGMM::runEM(vector<CvMat*> _data,
                 const int _numiter,
                 const int _minDiff)
{
    int endFlag = 0;
    int iter    = 0;
 
    int dimension    = m_dimension;
    int numMixtures  = m_nummixtures;
    int numAllframes = 0;
    vector<CvMat*>::iterator iter_d = _data.begin();
    for (; iter_d != _data.end(); ++iter_d)
        numAllframes += static_cast<int>(cvGetSize(*iter_d).height);
    
    DoubleVect loglik(numMixtures);
    DoubleVect::iterator iter_l;
    
    DoubleVect     s0_stat(numMixtures, 0); // 零次統計量
    vector<CvMat*> s1_stat(numMixtures);    // 一次統計量
    vector<CvMat*> s2_stat(numMixtures);    // 二次統計量 
    DoubleVect::iterator     iter_s0;
    vector<CvMat*>::iterator iter_s1;
    vector<CvMat*>::iterator iter_s2;
    
    DoubleVect::iterator     iter_pi;
    vector<CvMat*>::iterator iter_mu;
    vector<CvMat*>::iterator iter_Sg;
    DoubleVect::iterator     iter_ct;

    CvMat* tmp_matrix = (m_covtype == FULLC) ?
        cvCreateMat(dimension, dimension, CV_32F):
        cvCreateMat(1, dimension, CV_32F);

    iter_s1 = s1_stat.begin();
    iter_s2 = s2_stat.begin();
    for (; iter_s1 != s1_stat.end(); ++iter_s1, ++iter_s2)
    {
        (*iter_s1) = cvCreateMat(1, dimension, CV_32F);
        (*iter_s2) = (m_covtype == FULLC) ?
            cvCreateMat(dimension, dimension, CV_32F):
            cvCreateMat(1, dimension, CV_32F);
    }
    
    double p_loglik = -DBL_MAX;

    clock_t start,end;
    while (!endFlag && iter < _numiter)
    {

      iter_s0 = s0_stat.begin();
      iter_s1 = s1_stat.begin();
      iter_s2 = s2_stat.begin();
        for (; iter_s0 != s0_stat.end(); ++iter_s0, ++iter_s1, ++iter_s2)
        {
            (*iter_s0) = 0;
            cvSetZero(*iter_s1);
            cvSetZero(*iter_s2);
        }
        /* E-step */
	double t1=0;
	double t2=0;
        DoubleVect::iterator iter_l;
        iter_d = _data.begin();
        for (; iter_d != _data.end(); ++iter_d)
        {

            int numFrames = static_cast<int>(cvGetSize(*iter_d).height);
            /* 1 frame づつ潜在変数の期待値を算出 */
            for (int i = 0; i < numFrames; ++i)
            {
#if 0
                CvMat stub;
                CvMat* submat = cvGetRow(*iter_d, &stub, i);
#else
                float* submat = ((*iter_d)->data.fl + i * m_dimension);
#endif
                double sumLogLik = -DBL_MAX;
                iter_l = loglik.begin();
		start = clock();
                for (int j = 0; iter_l != loglik.end(); ++iter_l, ++j)
                {
                    (*iter_l) = calcLogLikelihood(submat, j);
                    sumLogLik = LAddS(sumLogLik, (*iter_l));
                }
		end = clock();
		t1 += ((double)(end-start)/CLOCKS_PER_SEC);
                iter_l  = loglik.begin();
                iter_s0 = s0_stat.begin();
                iter_s1 = s1_stat.begin();
                iter_s2 = s2_stat.begin();


		start = clock();
                for (; iter_l != loglik.end(); 
		     ++iter_l, ++iter_s0, ++iter_s1, ++iter_s2)
                {
                    /* 潜在変数の期待値を算出 (r(z_{nk}))*/
                    (*iter_l) = exp((*iter_l) - sumLogLik);
                    /* 零次統計量の更新 */
                    (*iter_s0) += (*iter_l);
#if 0
                    /* 一次統計量の更新 */
                    cvAddWeighted(*iter_s1, 1, submat, *iter_l, 0, *iter_s1);
                    /* 二次統計量の更新 */
                    if (m_covtype == FULLC)
                        cvMatMul(submat, submat, tmp_matrix);
                    else
                        cvPow(submat, tmp_matrix, 2);
#else
                    for (int d = 0; d < m_dimension; ++d)
                    {
                        /* 一次統計量の更新 */
                        (*iter_s1)->data.fl[d] += submat[d] * (*iter_l);
                        /* 二次統計量の更新 */
                        tmp_matrix->data.fl[d] = square(submat[d]);
                    }
#endif
		    cvAddWeighted(*iter_s2, 1, tmp_matrix, *iter_l, 0, *iter_s2);
                }
		end = clock();
		t2 += ((double)(end-start)/CLOCKS_PER_SEC);
            } // # frames
        } // # segments
	printf("  E-step: Elapsed time is %.4f seconds\n",t1+t2);
	printf("\t(EXP)Elapsed time is %.4f seconds\n",t1);
	printf("\t(STA)Elapsed time is %.4f seconds\n",t2);



	start = clock();
        /* M-Step */
        iter_s0 = s0_stat.begin();
        iter_s1 = s1_stat.begin();
        iter_s2 = s2_stat.begin();
        iter_pi = m_pi.begin();
        iter_mu = m_mu.begin();
        iter_Sg = m_Sigma.begin();
        iter_ct = m_const.begin();
        for (; iter_s0 != s0_stat.end(); 
             ++iter_s0, ++iter_s1, ++iter_s2,
             ++iter_pi, ++iter_mu, ++iter_Sg, ++iter_ct)
        {
            (*iter_pi) = (*iter_s0) / static_cast<double>(numAllframes);
            cvConvertScale(*iter_s1, *iter_mu, 1.0 / (*iter_s0));
            if (m_covtype == FULLC)
                cvMatMul(*iter_mu, *iter_mu, tmp_matrix);
            else
                cvPow(*iter_mu, tmp_matrix, 2);
            cvAddWeighted(*iter_s2, 1.0 / (*iter_s0), tmp_matrix, - 1.0, 0, *iter_Sg);
            if (m_covtype == FULLC)
                (*iter_ct) = - 0.5 * (m_dimension * log2PI + log(cvDet(*iter_Sg)));
            else
            {
                cvLog(*iter_Sg, tmp_matrix);
                (*iter_ct) = - 0.5 * (m_dimension * log2PI + cvSum(tmp_matrix).val[0]);
            }
        }
    end = clock();
    printf("  M-step:\n");
    printf("\tElapsed time is %.4f seconds\n",(double)(end-start)/CLOCKS_PER_SEC);


        /* 尤度計算 */
        double new_loglik = calcLogLikelihood(_data);
        double diff       = new_loglik - p_loglik;

	//        if (diff <= _minDiff)
	//   endFlag = 1;
        
        if (m_trace)
	  printf("iter %d: LogLikelihood %.4f, delta %.4f\n", 
		 iter, 
		 new_loglik, (iter==0)?new_loglik:diff);

        p_loglik = new_loglik;
        iter++;
    }
    if (m_trace)
      cout << "stop." <<endl;
    iter_s1 = s1_stat.begin();
    iter_s2 = s2_stat.begin();
    for (; iter_s1 != s1_stat.end(); ++iter_s1, ++iter_s2)
    {
        cvReleaseMat(&(*iter_s1));
        cvReleaseMat(&(*iter_s2));
    }
    cvReleaseMat(&tmp_matrix);
}

//----------------------------------------------------------------------

double CGMM::calcLogLikelihood(CvMat* _frame, const int _j)
{
    return m_pi[_j] + m_const[_j] + mahala(_frame, m_mu[_j], m_Sigma[_j], m_covtype);
}

//----------------------------------------------------------------------

double CGMM::calcLogLikelihood(float* _frame, const int _j)
{
    return m_pi[_j] + m_const[_j] + 
           mahala(_frame, m_mu[_j]->data.fl, 
		  m_Sigma[_j]->data.fl, m_covtype, m_dimension);
}

//----------------------------------------------------------------------

double CGMM::calcLogLikelihood(vector<CvMat*> _data)
{
    double loglik  = 0;
    vector<CvMat*>::iterator iter_d = _data.begin();
    for (; iter_d != _data.end(); ++iter_d)
        loglik += calcLogLikelihood(*iter_d);
    return loglik;
}

//----------------------------------------------------------------------

double CGMM::calcLogLikelihood(CvMat* _data)
{
    double loglik  = 0;
    int numdata = static_cast<int>(cvGetSize(_data).height);
    for (int i = 0; i < numdata; ++i)
    {
        double lg = -DBL_MAX;
#if 0
        CvMat stub;
        CvMat* submat = cvGetRow(_data, &stub, i);
        for (int j = 0; j < m_nummixtures; ++j)
        {
            double l = calcLogLikelihood(submat, j);
            lg = LAddS(lg, l);
        }
#else

        float* submat = (_data->data.fl + i * m_dimension);
        for (int j = 0; j < m_nummixtures; ++j)
        {
	  //	  double l = log(m_pi[j]) + calcLogLikelihood(submat, j);
	  	  double l =  calcLogLikelihood(submat, j);
	  lg = LAddS(lg, l);
        }
#endif
        loglik += lg;
    }
    return loglik;
}

//----------------------------------------------------------------------

double CGMM::calcLogLikelihood(float* _data)
{
    double loglik  = -DBL_MAX;
    for (int j = 0; j < m_nummixtures; ++j)
    {
      double l = log(m_pi[j]) + calcLogLikelihood(_data, j);
      loglik = LAddS(loglik, l);
    }
    return loglik;
}

//----------------------------------------------------------------------
double CGMM::calcKL(CGMM* _gmm)
{

  const int dimension = m_dimension;
  double KL = 0;
  for (int j = 0; j < m_nummixtures; j++)
  {
    const CvMat* mu    = getMu(j);
    const CvMat* sigma = getSigma(j);

    float x[m_dimension];
    for (int d = 0; d < m_dimension; ++d)
      x[d] = cvmGet(mu, 0, d);
    
    double loglik = 0;
    for (int d = 0; d < m_dimension; ++d)
    {
      x[d]   +=      sqrt(m_dimension * cvmGet(sigma, 0, d));
      loglik += _gmm->calcLogLikelihood(x);
      x[d]   -=  2 * sqrt(m_dimension * cvmGet(sigma, 0, d));
      loglik += _gmm->calcLogLikelihood(x);
      x[d]   +=      sqrt(m_dimension * cvmGet(sigma, 0, d));
      //      for (int d = 0; d < m_dimension; ++d) cout<<x[d]<<",";cout<<endl;
    }
    KL = KL + m_pi[j] * loglik;
  }

  return KL / (2 * m_dimension);
}

//----------------------------------------------------------------------

string CGMM::printInfo()
{
    ostringstream os;
    int numMixtures = getNumMixtures();

    os << "dimension    = " << getDimension()   << endl;
    os << "num mixtures = " << numMixtures << endl;
    os << "covtype      = " << ((getCovType()==FULLC)?"FULLC":"DIAGC") << endl;
    os << "pi           = " << endl;
    for(int i = 0; i< numMixtures; ++i )
        os << m_pi[i] << ",";
    os << endl;
    os << "const        = " << endl;
    for(int i = 0; i< numMixtures; ++i )
        os << m_const[i] << ",";
    os << endl;
    os << "mu           = " << endl;
    for(int i = 0; i< numMixtures; ++i )
        os << _mat2str(m_mu[i]);
    os << endl;
    os << "Sigma        = " << endl;
    for(int i = 0; i< numMixtures; ++i )
        os << _mat2str(m_Sigma[i]);

    return os.str();
}



//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------



int main(int _argc, char* _argv[])
{
  int nummixtures = atoi(_argv[1]);
  int dimension   = atoi(_argv[2]);
  int covtype     = atoi(_argv[3]);
  int numdata     = 1000;
  int numsegments = 1000;
  CGMM* gmm = new CGMM(nummixtures, dimension, covtype);
  clock_t start,end;  
  srand(0);

  cout << "Creating data"<<endl;
  start = clock();
  vector<CvMat*> segments(numsegments);
  for (int i=0; i<numsegments; ++i)
  {
    segments[i] = cvCreateMat(numdata,dimension,CV_32F);
    for (int ii=0; ii<numdata;++ii)
      for (int jj=0; jj<dimension;++jj)
	cvmSet(segments[i], ii,jj, rand()%10);
    
  }
  end = clock();
  printf("\tElapsed time is %.4f seconds\n",(double)(end-start)/CLOCKS_PER_SEC);



  cout << "Initializing"<<endl;
  start = clock();
  gmm->init(segments);
  end = clock();
  printf("\tElapsed time is %.4f seconds\n",(double)(end-start)/CLOCKS_PER_SEC);



 gmm->setTrace(1);
 // cout<<gmm->printInfo()<<endl;

 cout << "CalcLogLikelihood"<<endl;
  start = clock();
 cout << gmm->calcLogLikelihood(segments[0])<<endl;
  end = clock();
  printf("\tElapsed time is %.4f seconds\n",(double)(end-start)/CLOCKS_PER_SEC);

  


  /*
  CvMat* x = cvCreateMat(numdata,dimension, CV_32F);
  for (int ii=0; ii<numdata;++ii)
    for (int jj=0; jj<dimension;++jj)
      cvmSet(x, ii,jj, rand()%10);
  CvMat* mu = cvCreateMat(1,dimension, CV_32F);
  for (int jj=0; jj<dimension;++jj)
    cvmSet(mu, 0,jj, rand()%10);
  CvMat* sigma = cvCreateMat(1,dimension, CV_32F);
  for (int jj=0; jj<dimension;++jj)
    cvmSet(sigma, 0,jj, rand()%10);
  start = clock();
  double l=0;
  for (int i = 0; i < numdata; ++i)
  {
    float* submat = (x->data.fl + i * dimension);
    l += mahala(submat, mu->data.fl, 
		sigma->data.fl, covtype, dimension);
  }
  cout << l<<endl;
  end = clock();
  printf("%.4f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);

  start = clock();
  l=0;
  for (int i = 0; i < numdata; ++i)
  {
    CvMat stub;
    CvMat* submat = cvGetRow(x, &stub, i);
    l += mahala(submat, mu, 
		sigma, covtype);
  }
  cout << l<<endl;
  end = clock();
  printf("%.4f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
;exit(3);
  */
  cout << "Execute EM algorithm"<<endl;
  gmm->runEM(segments, 1000, 0);

  delete gmm;
}

