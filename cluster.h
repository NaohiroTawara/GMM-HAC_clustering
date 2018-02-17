//
//  cluster.h
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 on 12/08/01.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#ifndef spkr_cl_GMM_HAC_cluster_h
#define spkr_cl_GMM_HAC_cluster_h

#include "util.h"
#include <vector>
#include "matrix.h"
//#include "cv.h"

#include <list>


using namespace std;

class CCluster
{
private:
  int m_numsegments;      // このクラスタに割当てられた発話総数
  int m_numallframes;     // このクラスタに割り当てられた発話の総フレーム数
  int m_dimension;        // 特徴量の次元
  int m_covtype;          // 共分散行列の種類
  
  list<int>           m_segidlist; // このクラスタに割当てられた発話番号のリスト
  list<const MATRIX*> m_seglist;   // このクラスタに割当てられた発話のポインタリスト
    
public:
    
  // コンストラクタ
 CCluster(const int _covtype)
   : m_numsegments(0), m_numallframes(0), m_dimension(0), m_covtype(_covtype)
  {};
 
  // デコンストラクタ
  ~CCluster() 
    {};
  
  int initCluster(vector <const MATRIX*> _feature, vector<int> _seg_id)
  {
    m_segidlist.resize(0);
    m_seglist.resize(0);
        
    vector<const MATRIX*>::iterator iter_f = _feature.begin();
    vector<int>::iterator iter_id          = _seg_id.begin();
    m_dimension    = (*iter_f)->cols();
    for (;iter_id != _seg_id.end(); ++iter_id, ++iter_f)
    {
      m_segidlist.push_back(*iter_id);
      m_seglist.push_back(*iter_f);
      m_numsegments  += 1;
      m_numallframes += (*iter_f)->rows();
      int new_dim = (*iter_f)->cols();
      if (m_dimension != new_dim)
      	Error(1111, "Dimension does not match, expected %d, have %d", 
      	      m_dimension, new_dim);
    }
   
    return m_dimension;
  }
    
    int getNumSegments() 
    { return m_numsegments; }
    
    int getNumAllFrames()
    { return m_numallframes; }
        
    int getDimension()
    { return m_dimension; }
    
    int getCovType()
    { return m_covtype;}
        
    void mergeCluster(CCluster* _cl)
    {
      
      m_segidlist.merge(*(_cl->getSegIdList()));
      m_seglist.merge(*(_cl->getSegList()));
      m_numsegments += _cl->getNumSegments();
      
    }
    
    int calcDistance(CCluster* _cl){return 0;};

    list<int>::iterator getSegIdListIterator()
    { return m_segidlist.begin();}

    list<int>::iterator getSegIdListIteratorTail()
    { return m_segidlist.end();}
    
    list<const MATRIX*>::iterator getSegListIterator()
    { return m_seglist.begin();}
    
    list<const MATRIX*>::iterator getSegListIteratorTail()
    { return m_seglist.end();}
    
    list<int>* getSegIdList()
    { return &m_segidlist; }

    list<const MATRIX*>* getSegList()
    { return &m_seglist; }


};

#endif
