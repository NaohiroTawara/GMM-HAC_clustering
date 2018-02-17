//
//  gmm-hac_clustering.cc
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 12/07/31.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#include <iostream>
#include "gmmclustering.h"

//----------------------------------------------------------------------

void CGMMClustering::malloc(const int _numclusters, const int _covtype)
{
  last_updated_i     = -1;
  last_updated_i_new = -1;
  updatedFlag        = 1;
  m_clusters.resize(_numclusters, NULL);
  // 空のクラスタを生成    
  list<CGMMMLCluster*>::iterator iter_cl = m_clusters.begin();    
  for(; iter_cl != m_clusters.end(); ++iter_cl)
    (*iter_cl) = new CGMMMLCluster(_covtype);
  
  m_distmat.resize(_numclusters - 1);
  vector<DoubleVect>::iterator iter_B_i = m_distmat.begin();
  for (int i = 0; iter_B_i != m_distmat.end(); ++iter_B_i, ++i)
    iter_B_i->resize(_numclusters - i - 1, 0);
}

//----------------------------------------------------------------------

void CGMMClustering::free()
{
  list<CGMMMLCluster*>::iterator iter_cl = m_clusters.begin();
  
  for(; iter_cl != m_clusters.end(); ++iter_cl)
    if (*iter_cl) delete (*iter_cl);
}

//----------------------------------------------------------------------

void CGMMClustering::initClusters(const int _nummixtures, 
				  const int _maxiter, 
				  const int _mindiff)
{
  list<CGMMMLCluster*>::iterator iter_cl = m_clusters.begin();    
  for(int i = 0; iter_cl != m_clusters.end(); ++iter_cl, ++i)
  {
    const MATRIX *feature =  getSegment(i);
    vector <const MATRIX*> features(1, feature);
    IntVect id(1, i);
    (*iter_cl)->initCluster(features, id, _nummixtures, _maxiter, _mindiff);
  }
  initDistMat();
}

//----------------------------------------------------------------------

void CGMMClustering::initDistMat(void)
{
  vector<DoubleVect>::iterator iter_B_i    = m_distmat.begin();
  list<CGMMMLCluster*>::iterator iter_cl_i = m_clusters.begin();
  for (; iter_cl_i != m_clusters.end(); ++iter_cl_i, ++iter_B_i)
  {
    DoubleVect::iterator iter_B_j = iter_B_i->begin();
    list<CGMMMLCluster*>::iterator iter_cl_j = iter_cl_i;
    iter_cl_j++;
    for (; iter_cl_j != m_clusters.end(); ++iter_cl_j, ++iter_B_j)
      (*iter_B_j) = calcDistance((*iter_cl_i), (*iter_cl_j));
  }
}

//----------------------------------------------------------------------

string CGMMClustering::getDistMat()
{
  ostringstream os;

  vector<DoubleVect>::iterator iter_B_i = m_distmat.begin();
  for (int i=0;iter_B_i != m_distmat.end(); ++iter_B_i, ++i)
  {
    if (iter_B_i->at(0) == -DBL_MAX) continue;
    DoubleVect::iterator iter_B_j = iter_B_i->begin();
    for (;iter_B_j != iter_B_i->end(); ++iter_B_j)
      if (*iter_B_j != DBL_MAX) os << (*iter_B_j)<<",";
    os <<endl;
  }
  return os.str();
}

//----------------------------------------------------------------------
/*
Comment:
m_clusters は類似度（距離）行列の上（下）三角行列の値を保持する２次元配列
1行目がデータ数-1の要素を持ち，2行目はデータ数-2, 3行目はデータ数-3, ..., 
N-1 行目は1 の要素を持つ．
CGMMClusteringインスタンス生成後初めてrun()を実行すると
まず全てのデータの組み合わせについて類似度を計算し m_clusters の全要素を埋める．
次にその中で類似度が最大となるクラスタをマージ（mergeClusters(・, ・)）する．
このとき以下の処理で m_clusters に pivot を埋め込む
 (1) m_clusters マージされる側のクラスタに対応する 行の1列目に DBL_MAX をpivotとして代入する．
 (2) DB_MAX を埋め込んだ行以外の全ての行についてマージされる側のクラスタに対応する列の要素に -DBL_MAX を代入する．
従って以上の処理完了後の m_clusters を見た時に，
ある行の 1 列目が DBL_MAX が格納されていれば
その行に対応するクラスタは既にマージされて消滅されていることがわかり，
ある列の要素に -DBL_MAX が格納されていれば，
その列に対応するクラスタも既にマージされて消滅されていることがわかる．
*/
int CGMMClustering::run()
{
  double min_dist = -DBL_MAX;
  int numsegments = getNumSegments();
  
  vector<DoubleVect>::iterator iter_B_i;
  DoubleVect::iterator iter_B_j;
  
  while(min_dist <= m_alpha)
  {
    min_dist = DBL_MAX;
    
    list<CGMMMLCluster*>::iterator min_i;
    list<CGMMMLCluster*>::iterator min_j;
    int min_jj;
    
    iter_B_i = m_distmat.begin();
    list<CGMMMLCluster*>::iterator iter_cl_i = m_clusters.begin();
    for (int i=0, ii = 0; iter_cl_i != m_clusters.end() && ii < numsegments - 1;
	 ++iter_cl_i, ++iter_B_i, ++i, ++ii)
    {
      list<CGMMMLCluster*>::iterator iter_cl_j = iter_cl_i;
      iter_cl_j++;
      while (iter_B_i->at(0) == -DBL_MAX)
      {
	++ii; ++iter_B_i;
	if (ii == getNumSegments() - 1) break;
      }
      if (ii == getNumSegments() - 1) break;
      
      iter_B_j = iter_B_i->begin();
      for (int j = i + 1, jj = ii + 1; iter_cl_j != m_clusters.end();
	   ++iter_cl_j, ++iter_B_j, ++j, ++jj)
      {
	while ((*iter_B_j) == DBL_MAX) {++iter_B_j; ++jj;}
                
	double dist;
	if (updatedFlag &&
	    (ii == last_updated_i || jj == last_updated_i))
	  { // 前のiterationで更新されたクラスタが係るクラスタのみ類似度を更新
	  dist = calcDistance((*iter_cl_i), (*iter_cl_j));
	  (*iter_B_j) = dist;
	}
	else
	  dist = (*iter_B_j);

	if (dist < min_dist)
        {
	  min_dist = dist;
	  min_i    = iter_cl_i;
	  min_j    = iter_cl_j;
	  min_jj = jj;
	  last_updated_i_new = ii;
	}
      }
    }

    if (min_dist <= m_alpha)
    { // Merge clustesers which has a minimam distance
      last_updated_i = last_updated_i_new;
      mergeClusters((*min_i), (*min_j));
      iter_B_i = m_distmat.begin();
      for (int i = 0; iter_B_i != m_distmat.end() 
	     && i < min_jj && i < numsegments - 2; ++iter_B_i,++i)
      { // 前述のルールに従って pivot を追加
	if (iter_B_i->at(min_jj - i -1) != -DBL_MAX)
	  iter_B_i->at(min_jj - i -1) = DBL_MAX;
      }
      if (min_jj<numsegments-1)
	m_distmat.at(min_jj).at(0) = -DBL_MAX;
      updatedFlag = 1;
    }
    else
      updatedFlag = 0;    
  }
  return 0;
}

//----------------------------------------------------------------------

void CGMMClustering::init()
{
    
}

//----------------------------------------------------------------------

int CGMMClustering::getClusteringResultSeg(IntVect* _data_cc)
{
  int numclusters = getNumClusters();
  _data_cc->resize(getNumSegments());
  
  list<CL_TYPE*>::iterator iter_cl = m_clusters.begin();
  for (int speaker_id = 0; iter_cl != m_clusters.end(); 
       ++iter_cl, ++speaker_id)
  {
    list<int>::iterator iter_id  
      = (*iter_cl)->getSegIdListIterator();
    list<int>::iterator iter_idt 
      = (*iter_cl)->getSegIdListIteratorTail();
    for (; iter_id != iter_idt; ++iter_id)
      _data_cc->at(*iter_id) = speaker_id;
  }
  return numclusters;
}

//----------------------------------------------------------------------

int CGMMClustering::getClusteringResultCl(vector<IntVect>* _cluster)
{  
  int numclusters = getNumClusters();

  _cluster->resize(numclusters,IntVect());

  vector<IntVect>::iterator iter_cc = _cluster->begin();
  list<CL_TYPE*>::iterator iter_cl  = m_clusters.begin();

  for(; iter_cl != m_clusters.end(); ++iter_cl, ++iter_cc)
  {
    iter_cc->resize((*iter_cl)->getNumSegments());
    list<int>::iterator iter_id  
      = (*iter_cl)->getSegIdListIterator();
    list<int>::iterator iter_idt 
      = (*iter_cl)->getSegIdListIteratorTail();
    copy(iter_id, iter_idt, iter_cc->begin());
  }

  return numclusters;
}

//----------------------------------------------------------------------

double CGMMClustering::calcDistance(CL_TYPE* _cl1, CL_TYPE* _cl2)
{
#ifdef _OPENMP
  //  double start, end;    
  //  start = omp_get_wtime();
#else
  //  clock_t start,end;
  //  start = clock();
#endif
  double distance = _cl1->calcDistance(_cl2);
#ifdef _OPENMP
  //  end = omp_get_wtime();
  //  printf("Dist: Elapsed time is %.4f seconds\n",(double)(end-start));
#else
  //  end = clock();
  //  printf("Dist: Elapsed time is %.4f seconds\n",(double)(end-start)/CLOCKS_PER_SEC);
#endif
  return distance;
}

//----------------------------------------------------------------------

void CGMMClustering::mergeClusters(CL_TYPE* _cl1, CL_TYPE* _cl2)
{
  // クラスタの統合
  _cl1->mergeCluster(_cl2);
  // リストの更新
  m_clusters.remove(_cl2);
  // 片方のクラスタのインスタンスを削除
  delete _cl2;    
  // クラスタ数の更新
  setNumClusters(getNumClusters() - 1);
}

//----------------------------------------------------------------------

void CGMMClustering::setBasisFromFile(char* _filename)
{
    
}

//----------------------------------------------------------------------

string CGMMClustering::showClusteringResult()
{
  ostringstream os;
    
  IntVect data_cc(getNumSegments());
  getClusteringResultSeg(&data_cc);
  vector<IntVect> cluster_cc;
  getClusteringResultCl(&cluster_cc);

  vector<IntVect>::iterator iter_c = cluster_cc.begin();
  for (int i = 0; iter_c != cluster_cc.end(); ++iter_c, ++i)
  {
    os << "cluster" << i << ": ";
    IntVect::iterator iter_c_c = iter_c->begin();
    for (; iter_c_c != iter_c->end(); ++iter_c_c)
      os << (*iter_c_c) << ", ";
    os << endl;
  }

  os << endl;
  IntVect::iterator iter_cc = data_cc.begin();
  os << "Allignment: ";
  for (; iter_cc != data_cc.end(); ++iter_cc)
      os << (*iter_cc) << ", ";
  os << endl;

  return os.str();
}

//----------------------------------------------------------------------

string CGMMClustering::printInfo()
{
  ostringstream os;
  list<CGMMMLCluster*>::iterator iter_c = m_clusters.begin();
  for (int i = 0; iter_c != m_clusters.end(); ++iter_c, ++i)
    os << "Cluster: "<<i<<endl<<(*iter_c)->printInfo()<<endl;
  return os.str();
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

