#include "eval.h"

int ISpkrClEvaluation::checkDataCc(IntVect& _datacc)
{
  IntVect::iterator iter_c = _datacc.begin();
  IntVect cnt_vec(0);
  int max_id = -1;
  for (;iter_c != _datacc.end(); ++iter_c)
  {
    if ((*iter_c) >= max_id)
    {
      max_id = (*iter_c) + 1;
      cnt_vec.resize(max_id, 0);
    }
    cnt_vec.at(*iter_c)++;
  }
  /*
  IntVect::iterator iter_v = cnt_vec.begin();
  for (int i = 0; iter_v != cnt_vec.end(); ++iter_v, ++i)
    if ((*iter_v) == 0)
      Error(1111, "There is no utterances of speaker [%d] in label vector", i);
  */
  if (cnt_vec.size() != getNumClass())
    Error(1111, "The number of clustersin label vector [%d] doesn't match with the given number of clusters [%d]", 
	  cnt_vec.size(), getNumClass());
  
  return 0;
}

/*
 ** @brief  MLFファイル（正解ラベルが列挙されたファイル）を読み込む
 ** @param  _filename:MLFファイル名
 ** @return 真の話者数
 */
int ISpkrClEvaluation::setMlf(vector<string> _spkr_label, DoubleVect _length_vec)
{
  int id      = 0; // 発話 id
  int num_spk = 0; // 話者 id

  DoubleVect::iterator iter_l     = _length_vec.begin();
  vector<string>::iterator iter_s = _spkr_label.begin();
  for (; iter_s != _spkr_label.end(); ++iter_s, ++iter_l)
  {
    string speaker = (*iter_s);
    m_mlf_vec.push_back(new MlfInfo(speaker, id, (*iter_l)));
    if (m_spk_map.find(speaker) == m_spk_map.end())
    {
      m_spk_map.insert(pair<string, int>(speaker, num_spk));
      ++num_spk;
    }
    ++id;
  }

  if (id != m_num_data)
  {
    Error(1111, "Number of data in MLF file is invalid: %d!!\n", id);
  }
  if (num_spk != m_true_num_spk)
  {
    Error(1111, "Number of speakers in MLF file is invalid: %d!!\n", num_spk);
  }
  return m_true_num_spk;
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
  
  ifstream ifs(datacc_filename.c_str());
  if (ifs==NULL)
  {
    cerr << "error: cannot read file: " << datacc_filename << endl;
    exit(-1);
  }
  int cc;
  while (ifs >> cc)
    data_cc.push_back(cc);

  Spkr_cl_evaluation* spkr_eval = new Spkr_cl_evaluation(num_data,
                                                         num_class,
                                                         true_num_spk);
  spkr_eval->readMlf(mlf_filename);
  Result result = spkr_eval->evaluate(data_cc);
  cout << "acp: "<<result.acp<<", asp: "<<result.asp << ", K: "<<result.K<<endl;
  delete spkr_eval;
  return 0;

}
  */


