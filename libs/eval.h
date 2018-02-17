/* ------------------------------------------------------
* 
* eval.h
* 
* < Author > N. TAWARA 2012/2/17
*
* < Note> クラスタリング精度算出クラスのインターフェースクラス
*
*
*
* ------------------------------------------------------
*/

#ifndef _SPKR_CL_EVALUATION__H_
#define _SPKR_CL_EVALUATION__H_

// Spekaer Diarizationの結果を評価をするクラス
  
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "util.h"
using namespace std;

//#define square(X) (X)*(X)
/* 評価結果を保持する構造体 */
struct Result
{
  map<string, double>result;
};

/* セグメントに関する情報を保持する構造体 */
struct MlfInfo
{
  string speaker; // 話者名
  int seg_id;     // セグメント番号
  double length;  // 発話長
  MlfInfo(string _speaker, int _id, double _length)
  : speaker(_speaker), seg_id(_id) , length(_length) {};
  ~MlfInfo(void){};
};

typedef vector<MlfInfo*> MlfVect;

class ISpkrClEvaluation
{
private:
  int m_num_data;       // セグメント数
  int m_num_class;      // クラスタ数
  int m_true_num_spk;   // 真の話者数

  MlfVect m_mlf_vec; // セグメント情報配列
  map<string, int> m_spk_map; //話者名と話者番号を関連付けるコンテナ
public:

  // コンストラクタ
  ISpkrClEvaluation(const int _num_data, 
		    const int _num_class, 
		    const int _true_num_spk)
    : m_num_data(_num_data), 
      m_num_class(_num_class),
      m_true_num_spk(_true_num_spk) 
      {} ;
  
  // デコンストラクタ
  ~ISpkrClEvaluation()
    {} 
  
  // 評価
  virtual Result evaluate(IntVect& _datacc) = 0;

  // MLFファイル読み込み
  int setMlf(vector<string> _spkr_label, DoubleVect _length_vec);

  // クラスタ数設定
  void setNumClass(const int _n) { m_num_class = _n; }


  double getFullLength()   
  {
    double full_length = 0;
    MlfVect::iterator iter_m = m_mlf_vec.begin();
    for (;iter_m != m_mlf_vec.end(); ++iter_m)
      full_length += (*iter_m)->length;
    return full_length;
  }

  int getNumData()         { return m_num_data;  }
  int getNumClass()        { return m_num_class; }
  int getNumTrueSpeakers() { return m_true_num_spk; }
  MlfVect* getMlfVect()    { return &m_mlf_vec;}
  map<string, int>* getSpkrMap() { return &m_spk_map;}

 protected:

  int checkDataCc(IntVect& _datacc);
};

#endif

