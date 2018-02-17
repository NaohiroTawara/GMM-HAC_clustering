//
//  gmm.cc
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 31/07/12.
//  Last updated by 直弘 俵 on 25/07/14.
//  Copyright 2012年-2014年 早稲田大学. All rights reserved.
//
//  変更点：
//    OpenMP を用いた並列化に対応
//    Eigen  を用いた行列計算
//    共分散を逆行列の形で保持することにより呼び出し回数の多い calcLikelihood の計算負荷軽減
//    単精度・倍精度をコンパイル時オプション(-DDOUBLE)で切り替え可

#include "gmm.h"

//----------------------------------------------------------------------



void CGMM::malloc(const int _nummixtures, 
		  const int _dimension, 
		  const int _covtype)
{
  m_pi.resize(_nummixtures);
  m_mu.resize(_nummixtures);
  m_Sigma.resize(_nummixtures);
  m_const.resize(_nummixtures);

  vector<MATRIX>::iterator iter_mu    = m_mu.begin();
  vector<MATRIX>::iterator iter_Sigma = m_Sigma.begin();
  for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
  {
    iter_mu->resize(1,_dimension);
  
    if (_covtype == FULLC)
      iter_Sigma->resize(_dimension, _dimension);
    else
      iter_Sigma->resize(1, _dimension);
  }
}

//----------------------------------------------------------------------

void CGMM::free()
{
  /*
  vector<_MAT>::iterator iter_mu    = m_mu.begin();
  vector<_MAT>::iterator iter_Sigma = m_Sigma.begin();
  for (;iter_mu != m_mu.end(); ++iter_mu, ++iter_Sigma)
  {
    if (*iter_mu) cvReleaseMat(&(*iter_mu));
    if (*iter_Sigma) cvReleaseMat(&(*iter_Sigma));
  }  
  */
}


//----------------------------------------------------------------------
/* データをランダムに分割してパラメタを初期化する */
void CGMM::init(vector<const MATRIX*>& _data)
{
  int numallframes = 0;


  vector<SCALAR> stat(m_nummixtures,0);

  vector<SCALAR>::iterator iter_pi = m_pi.begin();
  vector<MATRIX>::iterator iter_mu = m_mu.begin();
  vector<MATRIX>::iterator iter_Sg = m_Sigma.begin();
  vector<SCALAR>::iterator iter_ct = m_const.begin();
  vector<SCALAR>::iterator iter_st = stat.begin();

  for (; iter_mu != m_mu.end(); ++iter_mu, ++iter_Sg)
  {
    (*iter_mu) = MATRIX::Zero(1,m_dimension);
    if(m_covtype == FULLC)
      (*iter_Sg) = MATRIX::Zero(m_dimension, m_dimension);
    else
      (*iter_Sg) = MATRIX::Zero(1, m_dimension);
  }

  /* GMM の要素に対する各フレームの初期割当てを決定 */
  vector<const MATRIX*>::iterator iter_d = _data.begin();  
  for (; iter_d != _data.end(); ++iter_d)
  {
    int numframes = (*iter_d)->rows();
    numallframes += numframes;
    MATRIX assignment = randperm(numframes, m_nummixtures);

    for(int i=0;i < numframes; ++i)
    {
      int k = assignment(0, i);
      /* 零次統計量の更新 */
      stat[k] += 1;
      /* 一次統計量の更新 */
      m_mu[k] += (*iter_d)->row(i);
      /* 二次統計量の更新 */
      if (m_covtype == FULLC)
	m_Sigma[k] += (*iter_d)->row(i).transpose() * (*iter_d)->row(i);
      else
	m_Sigma[k].array() += (*iter_d)->row(i).array().pow(2);
    }
  }

  iter_pi = m_pi.begin();
  iter_st = stat.begin();
  iter_mu = m_mu.begin();
  iter_Sg = m_Sigma.begin();
  iter_ct = m_const.begin();
  for (; iter_mu != m_mu.end(); 
       ++iter_st, ++iter_pi, ++iter_mu, ++iter_Sg, ++iter_ct)
  {
    (*iter_pi) = (*iter_st)  / numallframes;
    (*iter_mu).array() *= (1.0 / (*iter_st));
    (*iter_Sg).array() *= (1.0 / (*iter_st));


    if (m_covtype == FULLC)
      (*iter_Sg)         -= (*iter_mu).transpose() * (*iter_mu);
    else
      (*iter_Sg).array() -= (*iter_mu).array().pow(2);

    if (m_covtype == FULLC)
      (*iter_ct) = - 0.5 * (m_dimension * log2pi + log((*iter_Sg).determinant()));
    else
      (*iter_ct) = - 0.5 * (m_dimension * log2pi + (*iter_Sg).array().log().sum());

    // 繰り返し呼び出される calcLogLikelihood で計算をしなくて良いようにここで逆行列にする
    if (m_covtype == FULLC)
      (*iter_Sg) = (*iter_Sg).inverse();
    else
      (*iter_Sg) = (*iter_Sg).array().inverse();
  }
}


//----------------------------------------------------------------------

void CGMM::runEM(vector<const MATRIX*>& _data,
                 const int _numiter,
                 const int _minDiff)
{
  int endFlag = 0;
  int iter    = 0;
 
  int dimension    = m_dimension;
  int numMixtures  = m_nummixtures;
  int numallframes = 0;
  SCALAR p_loglik  = -DBL_MAX; // Loglikelihood in previous iteration

  vector<const MATRIX*>::iterator iter_d = _data.begin();
  for (; iter_d != _data.end(); ++iter_d)
    numallframes += (*iter_d)->rows();

  vector<SCALAR> s0_stat(numMixtures); // 0th statistics
  vector<MATRIX> s1_stat(numMixtures); // 1st statistics
  vector<MATRIX> s2_stat(numMixtures); // 2nd statistics
  vector<VECTOR> loglik(numMixtures);  // loglikelihood MATRIX

  vector<SCALAR>::iterator iter_s0;
  vector<MATRIX>::iterator iter_s1;
  vector<MATRIX>::iterator iter_s2;

  vector<SCALAR>::iterator iter_pi;
  vector<MATRIX>::iterator iter_mu;
  vector<MATRIX>::iterator iter_Sg;
  vector<VECTOR>::iterator iter_ll;
  vector<SCALAR>::iterator iter_ct;

#ifdef _OPENMP 
  //    double start, end;
#else
  //    clock_t start,end;
#endif


  while (!endFlag && iter < _numiter)
  {

    /* Clear all statistics */
    iter_s0 = s0_stat.begin();
    iter_s1 = s1_stat.begin();
    iter_s2 = s2_stat.begin();
    for (; iter_s1 != s1_stat.end(); ++iter_s0, ++iter_s1, ++iter_s2)
    {
      (*iter_s0) = 0;
      (*iter_s1) = MATRIX::Zero(1, dimension);
      if(m_covtype == FULLC)
	(*iter_s2) = MATRIX::Zero(dimension, dimension);
      else
	(*iter_s2) = MATRIX::Zero(1, dimension);
    }

    /* E-step */
    double t1=0, t2=0;
    iter_d = _data.begin();
    for (int s=0; iter_d != _data.end(); ++iter_d,++s)
    {

#ifdef _OPENMP
      //      start = omp_get_wtime();
#else
      //      start = clock();
#endif
      /* Calcurate log expectation for each frame */
      VECTOR logsummat = VECTOR::Constant((*iter_d)->rows(), -DBL_MAX);
      iter_ll = loglik.begin();
      for (int j = 0; iter_ll != loglik.end(); ++iter_ll, ++j)
      {
	(*iter_ll) = calcLogLikelihood(**iter_d, j);
	logsummat  = LAddS_mat(logsummat, *iter_ll);  
      }
#ifdef _OPENMP
      //            end = omp_get_wtime();
      //            t1 += end-start;
#else
      //      end = clock();
      //      t1 += ((double)(end-start)/CLOCKS_PER_SEC);
#endif

      iter_s0 = s0_stat.begin();
      iter_s1 = s1_stat.begin();
      iter_s2 = s2_stat.begin();
      iter_ll = loglik.begin();

#ifdef _OPENMP
      //      start = omp_get_wtime();
#else
      //      start = clock();
#endif
      for (int s=0; iter_s1 != s1_stat.end(); 
	   ++iter_ll, ++iter_s0, ++iter_s1, ++iter_s2, ++s)
      {
	/* Calcurate the expectation value of latent variables (r(z_{nk}))*/
	(*iter_ll).array() = ((*iter_ll).array() - logsummat.array()).exp();

	/* Update 0th order of statistics */
	(*iter_s0) += (*iter_ll).array().sum();

#if 0
	(*iter_s1).array() += 
	  ((*iter_d)->array() * 
	   (*iter_ll).replicate(1,dimension).array()).colwise().sum();
	if (m_covtype == FULLC)
	  (*iter_s2) += 
	    ( (*iter_d)->array() * 
	      (*iter_ll).replicate(1, dimension).array()
	    ).matrix().transpose() * (**iter_d);
	else
	  (*iter_s2).array() += 
	    ( (*iter_d)->array().pow(2) * 
	      (*iter_ll).replicate(1,dimension).array()
	    ).colwise().sum();	
#else
	/* Update 1st order of statistics */
	/* Update 2nd order of statistics */
	int f, d;
	int numframes = (*iter_d)->rows();
	if (m_covtype == FULLC)
	{
	  (*iter_s1).array() += 
	    ((*iter_d)->array() * 
	     (*iter_ll).replicate(1,dimension).array()).colwise().sum();
	  (*iter_s2) += 
	    ( (*iter_d)->array() * 
	      (*iter_ll).replicate(1, dimension).array()
	      ).matrix().transpose() * (**iter_d);
	}else{  
	  for (d=0; d < dimension;++d)
	  {
	    SCALAR l1 = (*iter_s1)(0,d), l2 = (*iter_s2)(0,d);
	    //#ifdef _OPENMP
	    //#pragma omp parallel for reduction(+:l1,l2)
	    //#endif
	    for (f=0; f<numframes;++f)
	    {
	      l1 += 
		((**iter_d)(f,d) * (*iter_ll)(f));
	      l2 += 
		( (**iter_d)(f,d) * (**iter_d)(f,d) *(*iter_ll)(f));
	    }
	    (*iter_s1)(0,d) = l1;
	    (*iter_s2)(0,d) = l2;
	  }
	  /*
	    for (f=0; f<numframes;++f)
	    for (d=0; d < dimension;++d)
	    {
	      (*iter_s1)(0,d) += 
		((**iter_d)(f,d) * (*iter_ll)(f,0));
	      (*iter_s2)(0,d) += 
		( (**iter_d)(f,d) * (**iter_d)(f,d) *(*iter_ll)(f,0));
	    }
	  */
	}
#endif
      }
#ifdef _OPENMP
      //      end = omp_get_wtime();
      //      t2 += end-start;
#else
      //      end = clock();
      //      t2 += ((double)(end-start)/CLOCKS_PER_SEC);
#endif      
    } // # segments
    //        printf("  E-step: Elapsed time is %.4f seconds\n",t1+t2);
    //    printf("\t(EXP)Elapsed time is %.4f seconds\n",t1);
    //    printf("\t(STA)Elapsed time is %.4f seconds\n",t2);
    //    start = clock();

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
      /* Update weights */

      (*iter_pi) =  (*iter_s0) / numallframes;

      /* Update mean vector */
      (*iter_mu).array() = (*iter_s1).array() * (1.0 / (*iter_s0));

      /* Update covariant MATRIX */
      if (m_covtype == FULLC)
	(*iter_Sg) = 
	  (*iter_mu).transpose() * (*iter_mu);
      else
	(*iter_Sg).array() = 
	  (*iter_mu).array() * (*iter_mu).array();
      (*iter_Sg) = (*iter_s2).array() * (1.0 / (*iter_s0)) - (*iter_Sg).array();


      /* Update normalization term */
      if (m_covtype == FULLC)
	(*iter_ct) = - 0.5 * (dimension * log2pi + log((*iter_Sg).determinant()));
      else
	(*iter_ct) = - 0.5 * (dimension * log2pi + (*iter_Sg).array().log().sum());

      // 繰り返し呼び出される calcLogLikelihood で計算をしなくて良いようにここで逆行列にする
      if (m_covtype == FULLC)
	(*iter_Sg) = (*iter_Sg).inverse();
      else
	(*iter_Sg) = (*iter_Sg).array().inverse();
      
    }
    //    end = clock();
    //    printf("  M-step:\n");
    //    printf("\tElapsed time is %.4f seconds\n",(SCALAR)(end-start)/CLOCKS_PER_SEC);
    

    /* 尤度計算 */
    SCALAR new_loglik = calcLogLikelihood(_data);
    SCALAR diff       = new_loglik - p_loglik;
    
    
    if (diff <= _minDiff)
      endFlag = 0;
    
    if (m_trace)
      printf("iter %d: LogLikelihood %.4f, delta %.4f\n", 
	     iter, 
	     new_loglik, (iter==0)?new_loglik:diff);

    p_loglik = new_loglik;
    
    iter++;
   }
  if (m_trace)
    cout << "stop." <<endl;
}

//----------------------------------------------------------------------

VECTOR CGMM::calcLogLikelihood(const MATRIX& _data, const int _j)
{
  return 
    m_pi[_j] +
    m_const[_j] + 
    //  mahala(_data, m_mu[_j], m_Sigma[_j], m_covtype).array();
    mahala3(_data, m_mu[_j], m_Sigma[_j], m_covtype).array();
}

//----------------------------------------------------------------------

SCALAR CGMM::calcLogLikelihood(const MATRIX& _data)
{
  VECTOR logmat  = VECTOR::Constant(_data.rows(), -DBL_MAX);
  for (int j = 0; j < m_nummixtures; ++j)
  {
    VECTOR tmp = calcLogLikelihood(_data, j);
    logmat = LAddS_mat(logmat, tmp);
  }
  return logmat.sum();
}

//----------------------------------------------------------------------

SCALAR CGMM::calcLogLikelihood(vector<const MATRIX*>& _data)
{
    SCALAR loglik  = 0;
    vector<const MATRIX*>::iterator iter_d = _data.begin();
    for (; iter_d != _data.end(); ++iter_d)
      loglik += calcLogLikelihood(**iter_d);
    return loglik;
}

//----------------------------------------------------------------------

string CGMM::printInfo()
{
  ostringstream os;
  int numMixtures = getNumMixtures();

  os << "dimension    = " << getDimension() << endl;
  os << "num mixtures = " << getNumMixtures() << endl;
  os << "covtype      = " << ((getCovType()==FULLC)?"FULLC":"DIAGC") << endl;
  os << "pi           = " << endl;
  for (int j=0; j<numMixtures;++j)
    os << m_pi[j] << ",";
  os << endl;
  os << "const        = " << endl;
  for (int j=0; j<numMixtures;++j)
    os << m_const[j] << ",";
  os << endl;
  os << "mu           = " << endl;
  for(int i = 0; i< numMixtures; ++i )
    os << m_mu[i]<<endl;
  os << endl;
  os << "Sigma        = " << endl;
  for(int i = 0; i< numMixtures; ++i )
    os << m_Sigma[i]<<endl;
  
  return os.str();
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------


#ifdef DEBUG
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
  //  Eigen::MATRIX randM = randperm(100,5);

  cout << "Creating test data"<<end;
  start = clock();  
  vector<const MATRIX*> segments(numsegments);
  for (int i=0; i<numsegments; ++i)
  {
    //    MATRIX* seg = new MATRIX(4,3);
    MATRIX* seg = new MATRIX(numdata,dimension);
    for (int ii=0; ii<numdata;++ii)
      for (int jj=0; jj<dimension;++jj)
        (*seg)(ii,jj) = rand()%10;
    segments[i] = seg;
  }
  end = clock();
  printf("\tElapsed time is %.4f seconds\n",(SCALAR)(end-start)/CLOCKS_PER_SEC);

  cout << "Initializing"<<endl;
  start = clock();
  gmm->init(segments);
  end = clock();
  printf("\tElapsed time is %.4f seconds\n",(SCALAR)(end-start)/CLOCKS_PER_SEC);


  gmm->setTrace(1);
  // cout<<gmm->printInfo()<<endl;
  
  cout << "CalcLogLikelihood"<<endl;
  start = clock();
  cout << gmm->calcLogLikelihood(*segments[0])<<endl;
  end = clock();
  printf("\tElapsed time is %.4f seconds\n",(SCALAR)(end-start)/CLOCKS_PER_SEC);



  /*
  MATRIX x = MATRIX(numdata,dimension);
  for (int ii=0; ii<numdata;++ii)
    for (int jj=0; jj<dimension;++jj)
      x(ii,jj) = rand()%10;
  MATRIX mu = MATRIX(1,dimension);
  for (int jj=0; jj<dimension;++jj)
    mu(0,jj) = rand()%10;
  MATRIX sigma = MATRIX(1,dimension);
  for (int jj=0; jj<dimension;++jj)
    sigma(0,jj) = rand()%10;
  clock_t start,end;
  cout << mahala(x,mu,sigma,DIAGC).sum()<<endl;
  end = clock();
  printf("%.4f秒かかりました\n",(SCALAR)(end-start)/CLOCKS_PER_SEC);
 ;exit(3);
  */
  cout << "Execute EM algorithm"<<endl;
   gmm->runEM(segments, 1000, 0);


 /*
 print(m.block(0,0,3,3));


 print(m.row(1).array()+1);
 print(m.colwise().mean());

 print(randM);
 print(randM.rowwise().sum());

 print(randM.col(1).array() *  randM.rowwise().sum().array().inverse());
 */
  delete gmm;
}
#endif

