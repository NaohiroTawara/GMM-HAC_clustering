/* ------------------------------------------------------
* 
* evalK.h
* 
* < Author > N. TAWARA 2012/2/17
*
* < Note> K �� �򻻽Ф��뤿��Υ��饹
*
*
*
* ------------------------------------------------------
*/

#ifndef _SPKR_CL_EVALUATION_K__H_
#define _SPKR_CL_EVALUATION_K__H_

#include "eval.h"

// K �ͤ˴�Ť�ɾ��


using namespace std;


class CSpkrClEvaluation_K: public ISpkrClEvaluation
{
public:

  // ���󥹥ȥ饯��
  CSpkrClEvaluation_K(const int _num_data, 
		      const int _num_class, 
		      const int _true_num_spk)
    : ISpkrClEvaluation(_num_data, _num_class, _true_num_spk)
  {} 
  
  // �ǥ��󥹥ȥ饯��
  ~CSpkrClEvaluation_K()
    {} 
  
  // ɾ��
  Result evaluate(IntVect& _datacc);

};

#endif

