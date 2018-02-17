#include "evalK.h"

/*
  ** @brief 評価
  ** @param  _datacc クラスタリング結果配列
  ** @return 結果を格納したResultクラスのインスタンス
  */
Result CSpkrClEvaluation_K::evaluate(IntVect& _datacc)
{
  checkDataCc(_datacc);
  
  Result result;
  int num_data     = getNumData();  //_dpm->getNumData();
  int num_class    = getNumClass(); // _dpm->getNumClass();
  int true_num_spk = getNumTrueSpeakers();
  
  /* 話者jが話したクラスタi内のセグメント数 */
  vector <IntVect> n(num_class); // [class x speaker]
  vector <IntVect>::iterator iter_n;
  IntVect::iterator          iter_n_ij;
  iter_n = n.begin();
  for (; iter_n != n.end(); ++iter_n)
    iter_n->resize(true_num_spk, 0);

  //  cout << "cluster:" << num_class << ", speaker:"
  //     << m_true_num_spk << endl;
  
  // for (int j = 0; j < num_data; ++j)
  //  n[_datacc[j]][m_spk_map[m_mlf_vec[j]->speaker]]++;

  MlfVect::iterator iter_m = getMlfVect()->begin();    // [numSegment]
  IntVect::iterator iter_c = _datacc.begin();
  for (; iter_m != getMlfVect()->end(); ++iter_m, ++iter_c)
    n[(*iter_c)][(*getSpkrMap())[(*iter_m)->speaker]]++;

  /*
  for (int i=0; i<num_class; ++i){
    for (int j=0; j<true_num_spk; ++j)
      cout <<n[i][j]<<",";
    cout <<endl;
  }
  */    

  /* n_{j}: 話者 j が発話したセグメント総数 */
  IntVect n_j(true_num_spk, 0);
  /* n_{i}: クラスタ i 内のセグメント総数 */
  IntVect n_i(num_class, 0);
  IntVect::iterator iter_n_i;
  IntVect::iterator iter_n_j;

  iter_n = n.begin();
  iter_n_i = n_i.begin();
  for (; iter_n_i != n_i.end(); ++iter_n_i, ++iter_n)
  {
    iter_n_ij = iter_n->begin();
    iter_n_j = n_j.begin();
    for (; iter_n_j != n_j.end(); ++iter_n_j, ++iter_n_ij)
    {
      (*iter_n_i) += (*iter_n_ij);
      (*iter_n_j) += (*iter_n_ij);
    }
  }
  
  iter_n_i = n_i.begin();
  for (int i = 0; iter_n_i != n_i.end(); ++iter_n_i, ++i)
    if ((*iter_n_i) == 0)
      Error(1111, "%There is no utterance in %d-th speaker", i);
  iter_n_j = n_j.begin();
  for (int j = 0; iter_n_j != n_j.end(); ++iter_n_j, ++j)
    if ((*iter_n_j) == 0)
      Error(1111, "%There is no utterance in %d-th speaker", j);
  
  /* p_{i}: cluster purityの算出 */
  DoubleVect p(num_class, 0);
  DoubleVect::iterator iter_p = p.begin();
  iter_n = n.begin();
  iter_n_i = n_i.begin();
  for (; iter_n != n.end(); ++iter_n, ++iter_n_i, ++iter_p)
  {
    iter_n_ij = iter_n->begin();
    for (; iter_n_ij != iter_n->end(); ++iter_n_ij)
      (*iter_p) += (*iter_n_ij)*(*iter_n_ij);
    (*iter_p) /= static_cast<double>((*iter_n_i)*(*iter_n_i));
  }

  /* I_{BBN}: BBN尺度の算出 */
  double I_BBN = 0;
  float Q = 0.5;
  for (int i=0; i < num_class; ++i)
    I_BBN += n_i[i] * p[i];
  I_BBN -= Q * num_class;

  /* I_{rand}: Rand Indexの算出 */
  double I_rand = 0;
  for (int i=0; i < num_class; ++i)
    I_rand += (n_i[i]) * (n_i[i]) * (0.5 - p[i]);
  for (int j=0; j < true_num_spk; ++j)
    I_rand += 0.5 * n_j[j];
  
  /* acp: average cluster purityの算出 */
  double acp = 0;
  iter_p = p.begin();
  iter_n_i = n_i.begin();
  for (; iter_p != p.end(); ++iter_p, ++iter_n_i)
    acp += (*iter_p) * (*iter_n_i);
  acp /= static_cast<double>(num_data);

  /* asp: average speaker purityの算出 */
  DoubleVect q(true_num_spk, 0);
  DoubleVect::iterator iter_q = q.begin();
  iter_n_j = n_j.begin();
  for (int j = 0; iter_q != q.end(); ++iter_q, ++iter_n, ++iter_n_j, ++j)
  {
    iter_n = n.begin();
    for (; iter_n != n.end(); ++iter_n)
    {
      (*iter_q) += iter_n->at(j)*iter_n->at(j);
    }
    (*iter_q) /= static_cast<double>((*iter_n_j)*(*iter_n_j));
  }
  double asp=0;
  iter_q = q.begin();
  iter_n_j = n_j.begin();
  for (; iter_q != q.end(); ++iter_q, ++iter_n_j)
    asp += (*iter_q) * (*iter_n_j);
  asp /= static_cast<double>(num_data);
  
  result.result.insert(pair<string, double>("acp", acp));
  result.result.insert(pair<string, double>("asp", asp));
  result.result.insert(pair<string, double>("K", sqrt(acp*asp)));
  result.result.insert(pair<string, double>("I_BBN", I_BBN));
  result.result.insert(pair<string, double>("I_rand", I_rand));
  return result;
}

/*
 ** クラステスト用main関数
 */
/*
int main(int _argc, char** _argv)
{
  int num_data        = atoi(_argv[1]); // 発話総数
  int num_class       = atoi(_argv[2]); // クラスタ数
  int true_num_spk    = atoi(_argv[3]); // 真のクラスタ数
  string mlf_filename = _argv[4];       // 正解ファイル名
  string datacc_filename= _argv[5];     // クラスタリング結果ファイル

  vector <int> data_cc;
  
  DoubleVect length_vec(num_data);
  DoubleVect::iterator iter_l = length_vec.begin();

  cout << "num_Data: "<<num_data<<endl;
  cout << "num_class: "<<num_class <<endl;
  cout << "true_num_spk: "<< true_num_spk <<endl;

  for(; iter_l != length_vec.end(); ++iter_l)
    (*iter_l) = 1;

  ifstream ifs(datacc_filename.c_str());
  if (ifs==NULL)
  {
    cerr << "error: cannot read file: " << datacc_filename << endl;
    exit(-1);
  }
  int cc;
  while (ifs >> cc)
    data_cc.push_back(cc);

  CSpkrClEvaluation_K* spkr_eval = new CSpkrClEvaluation_K(num_data,
							   num_class,
							   true_num_spk);

    spkr_eval->readMlf(mlf_filename, length_vec);
  
  Result result = spkr_eval->evaluate(data_cc);
  cout << "acp: "   << result.result["acp"] 
       << ", asp: " << result.result["asp"] 
       << ", K: "   << result.result["K"]
       << endl;
  delete spkr_eval;
  return 0;

}
*/


