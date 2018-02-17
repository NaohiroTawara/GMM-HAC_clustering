//
//  gmmcluster.cc
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 12/07/30.
//  Copyright 2012年 早稲田大学. All rights reserved.
//


#include "gmmmlcluster.h"

//----------------------------------------------------------------------

void CGMMMLCluster::malloc(const int _nummixtures, 
			   const int _dimension)
{
  m_gmm = new CGMM(_nummixtures, _dimension, getCovType());
}

//----------------------------------------------------------------------

void CGMMMLCluster::free()
{
  delete m_gmm;
}

//----------------------------------------------------------------------

void CGMMMLCluster::mergeCluster(CGMMMLCluster* _cl)
{
  CCluster::mergeCluster(_cl);
  updateCluster();
}

//----------------------------------------------------------------------

void CGMMMLCluster::updateCluster()
{
  vector<const MATRIX*> data(0);
  list<const MATRIX*>::iterator iter_d  = getSegListIterator();
  list<const MATRIX*>::iterator iter_dt = getSegListIteratorTail();
  for (;iter_d != iter_dt; ++iter_d)
    data.push_back(*iter_d);

  m_gmm->runEM(data, m_maxiter, m_mindiff);
}

//----------------------------------------------------------------------

double CGMMMLCluster::calcDistance(CGMMMLCluster* cl2)
{
  
  double CLR = 0;
    
  CGMMMLCluster* cl1 = this;
  list<const MATRIX*>::iterator iter_d1  = cl1->getSegListIterator();
  list<const MATRIX*>::iterator iter_d2  = cl2->getSegListIterator();
  list<const MATRIX*>::iterator iter_dt1 = cl1->getSegListIteratorTail();
  list<const MATRIX*>::iterator iter_dt2 = cl2->getSegListIteratorTail();
  
  for (;iter_d1 != iter_dt1; ++iter_d1)
  {
    CLR += cl1->calcLogLikelihood(**iter_d1) - 
           cl2->calcLogLikelihood(**iter_d1);
  }
  for (;iter_d2 != iter_dt2; ++iter_d2)
  {
    CLR += cl2->calcLogLikelihood(**iter_d2) - 
           cl1->calcLogLikelihood(**iter_d2);
  }
  return CLR;
  
}

//----------------------------------------------------------------------

void CGMMMLCluster::initCluster(vector <const MATRIX*>& _data, 
				vector <int>& _seg_id,
                                const int _nummixtures,
                                const int _maxiter,
                                const int _mindiff)
{
  
  m_maxiter = _maxiter;
  m_mindiff = _mindiff;
  int dimension = CCluster::initCluster(_data, _seg_id);
  malloc(_nummixtures, dimension);
  m_gmm->init(_data);
  updateCluster();
  
}

//----------------------------------------------------------------------

double CGMMMLCluster::calcLogLikelihood(const MATRIX& _data)
{
  return m_gmm->calcLogLikelihood(_data);
}

//----------------------------------------------------------------------

string CGMMMLCluster::printInfo()
{
  ostringstream os;
  os << m_gmm->printInfo();
  return os.str();
 
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifdef DEBUG
int main()
{
  int nummixtures = 10;
  int covtype     = DIAGC;
  int dimension   = 10;
  int numdata     = 1000;
  int numsegments = 1000;

  int maxiter = 100;
  int mindiff = 0;

  list<CGMMMLCluster*> gmmclusters(2);
  list<CGMMMLCluster*>::iterator iter_cl = gmmclusters.begin();

  cout << "Creating test data"<<endl;
  vector<const MATRIX*> segments(numsegments);
  for (int i=0; i<numsegments; ++i)
  {
    MATRIX* seg = new MATRIX(numdata,dimension);
    for (int ii=0; ii<numdata;++ii)
      for (int jj=0; jj<dimension;++jj)
        (*seg)(ii,jj) = rand()%10;
    segments[i] = seg;
  }
  cout << "done"<<endl;

  iter_cl = gmmclusters.begin();
  for(; iter_cl != gmmclusters.end(); ++iter_cl)
    (*iter_cl) = new CGMMMLCluster(covtype);
  

  iter_cl = gmmclusters.begin();
  for(int i=0; iter_cl != gmmclusters.end(); ++iter_cl, ++i)
  {
    vector<int> id(1, i);
    (*iter_cl)->initCluster(segments, id, nummixtures, maxiter, mindiff);
  }

  iter_cl = gmmclusters.begin();
  for(; iter_cl != gmmclusters.end(); ++iter_cl)
    delete (*iter_cl);
   
}
#endif

