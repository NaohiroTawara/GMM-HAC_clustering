/* ------------------------------------------------------
* 
* evalK.h
* 
* < Author > N. TAWARA 2012/2/17
*
* < Note> K 値 を算出するためのクラス
*
*
*
* ------------------------------------------------------
*/

#ifndef _SPKR_CL_EVALUATION_K__H_
#define _SPKR_CL_EVALUATION_K__H_

#include "eval.h"

// K 値に基づく評価


using namespace std;


class CSpkrClEvaluation_K: public ISpkrClEvaluation
{
public:

  // コンストラクタ
  CSpkrClEvaluation_K(const int _num_data, 
		      const int _num_class, 
		      const int _true_num_spk)
    : ISpkrClEvaluation(_num_data, _num_class, _true_num_spk)
  {} 
  
  // デコンストラクタ
  ~CSpkrClEvaluation_K()
    {} 
  
  // 評価
  Result evaluate(IntVect& _datacc);

};

#endif

