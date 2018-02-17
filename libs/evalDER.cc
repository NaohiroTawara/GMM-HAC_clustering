#include "evalDER.h"

void CSpkrClEvaluation_DER::setDER(IntVect& _datacc)
{
  int true_num_spk = getNumTrueSpeakers();
  int num_class = getNumClass();

  m_DER_matrix.resize(true_num_spk);
  m_spkr_utt_length.resize(true_num_spk, 0);

  DoubleMat::iterator iter_d = m_DER_matrix.begin();
  for (;iter_d != m_DER_matrix.end(); ++iter_d)
    iter_d->resize(num_class, 0);

  // [true_num_spkr x m_num_class]
  DoubleVect::iterator iter_l   = m_spkr_utt_length.begin();
  iter_d = m_DER_matrix.begin();
  for (int i = 0; iter_d != m_DER_matrix.end(); ++iter_d, ++iter_l, ++i)
  { // i: 話者番号
    MlfVect::iterator iter_m = getMlfVect()->begin();
    for (; iter_m != getMlfVect()->end(); ++iter_m)
    {
      int spkr_id = getSpkrMap()->find((*iter_m)->speaker)->second;
      if (i == spkr_id)
	(*iter_l) += (*iter_m)->length;
    }

    DoubleVect::iterator iter_d_i = iter_d->begin();
    for (int j = 0; iter_d_i != iter_d->end(); ++iter_d_i, ++j)
    { // j: クラスタ番号
      IntVect::iterator iter_c = _datacc.begin();
      iter_m = getMlfVect()->begin();
      for (; iter_m != getMlfVect()->end(); ++iter_m, ++iter_c)
      {
	int spkr_id = getSpkrMap()->find((*iter_m)->speaker)->second;
	if (i == spkr_id && j != (*iter_c))
	  (*iter_d_i) += (*iter_m)->length;
      }
    }
  }
}

void CSpkrClEvaluation_DER::sorting(const int _i,
				 list<int>* pibot)
{
  vector<pair<double, int> > m;
  list<int>::iterator iter_p    = pibot->begin();
  for (; iter_p != pibot->end(); ++iter_p)
    m.push_back(pair<int,int>(getErrorLength(_i, (*iter_p)), (*iter_p)));

  sort(m.begin(), m.end(), pair_less());

  vector<pair<double, int> >::iterator iter_m = m.begin();
  iter_p    = pibot->begin();
  for (; iter_m != m.end(); ++iter_m, ++iter_p)
    (*iter_p) = iter_m->second;
}

void CSpkrClEvaluation_DER::calcDER(const int _i,
				    const double _DER,
				    list<int> pibot)
{
  int num_true_spk = getNumTrueSpeakers();
  int num_class    = getNumClass();
  int DER_matrix_height = min(num_true_spk, num_class);
 
  /* エラー長が昇順になるように pibot のソート */
  if ((_i + 1) < DER_matrix_height)
  sorting(_i, &pibot);

  for (int cnt = 0; cnt < pibot.size(); ++cnt)
  {
    if ((_i + 1) == DER_matrix_height)
    {
      if (num_true_spk <= num_class)
      { // 推定クラスタ数のほうが多い場合
	double err = DBL_MAX;
	list<int>::iterator itr = pibot.begin();
	for(;itr != pibot.end(); ++itr)
	{
	  double val = getErrorLength(_i, (*itr));
	  if (val < err) err = val;
	}
	if (_DER + err < m_min_DER)
	  m_min_DER = _DER + err;
      }
      else
      { //真のクラスタ数のほうが多い場合
	list<int>::iterator itr = pibot.begin();
	for(;itr != pibot.end(); ++itr)
	{
	  double err = getErrorLength(_i, (*itr));
	  list<int>::iterator itr2 = pibot.begin();
	  for(;itr2 != pibot.end(); ++itr2)
	  {
	    if ((*itr) == (*itr2)) continue;
	    err += getErrorLength(*itr2);
	  }
	  if (_DER + err < m_min_DER)
	    m_min_DER = _DER + err;
	}
      }
      //      cout <<"  "<<m_min_DER<<endl;
      return;
    }

    int j = pibot.front();
    /*
    cout <<"("<<_i+1<<", "<<j+1<<")"<<endl<<" ";
    list<int>::iterator iter_p = pibot.begin();
    for(;iter_p != pibot.end(); ++iter_p)
      cout <<(*iter_p)<<", ";
    cout<<endl;
    */    
    double err = getErrorLength(_i, j);
    //    cout << "  "<<err<<" "<<endl;
    list<int>::iterator itr = pibot.begin();

    pibot.pop_front();

    if (_DER + err < m_min_DER)
      calcDER(_i + 1, _DER+ err, pibot);

    pibot.push_back(j);
  }
}

/*
  ** @brief 評価
  ** @param  _datacc クラスタリング結果配列
  ** @return 結果を格納したResultクラスのインスタンス
  */
Result CSpkrClEvaluation_DER::evaluate(IntVect& _datacc)
{
  checkDataCc(_datacc);
  Result result;
  
  list <int> pibot;
  if (getNumClass() > getNumTrueSpeakers())
    pibot.resize(getNumClass());
  else
    pibot.resize(getNumTrueSpeakers());
  list <int>::iterator iter_p   = pibot.begin();

  for (int j = 0; iter_p != pibot.end(); ++iter_p, ++j)
    (*iter_p) = j;

  setDER(_datacc);

  /*
  DoubleMat::iterator iter_d = m_DER_matrix.begin();
  DoubleVect::iterator iter_l   = m_spkr_utt_length.begin();
  
  for (; iter_d != m_DER_matrix.end(); ++iter_d, ++iter_l)
  {
    cout << (*iter_l) << ": ";
    DoubleVect::iterator iter_d_i = iter_d->begin();
    for (; iter_d_i != iter_d->end(); ++iter_d_i)
      cout << (*iter_d_i) << ",";
    cout << endl;
  }
  */  
  calcDER(0, 0, pibot);

  m_min_DER /= static_cast<double>(getFullLength());
  result.result.insert(pair<string, double>("DER", m_min_DER));

  return result;
}

/*
  ==================================================
   クラステスト用main関数
   ==================================================
*/
/*
int main(int _argc, char** _argv)
{ 
  int num_data           = atoi(_argv[1]); // 発話総数
  int num_class          = atoi(_argv[2]); // クラスタ数
  int true_num_spk       = atoi(_argv[3]); // 真のクラスタ数
  string mlf_filename    = _argv[4];       // 正解ファイル名
  string datacc_filename = _argv[5];     // クラスタリング結果ファイル

  IntVect data_cc;

  DoubleVect length_vec(num_data);

  DoubleVect::iterator iter_l = length_vec.begin();

  cout << "num_Data: "<<num_data<<endl;
  cout << "num_class: "<<num_class <<endl;
  cout << "true_num_spk: "<< true_num_spk <<endl;

  for(; iter_l != length_vec.end(); ++iter_l)
    (*iter_l) = 1;

  ifstream ifs(datacc_filename.c_str());
  if (ifs == NULL)
    Error(1111, "cannot read file");
  
  int cc;
  
  while (ifs >> cc)
    data_cc.push_back(cc);

  CSpkrClEvaluation_DER* spkr_eval = 
    new CSpkrClEvaluation_DER(num_data, num_class, true_num_spk);

  spkr_eval->readMlf(mlf_filename, length_vec);

  Result result = spkr_eval->evaluate(data_cc);
  cout << "min_DER: " << result.result["DER"] << endl;
  delete spkr_eval;
  return 0;
}
*/

