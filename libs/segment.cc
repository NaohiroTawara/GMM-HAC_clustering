#include "segment.h"

//----------------------------------------------------------------------

void CSegment::free()
{
  vector<MATRIX*>::iterator iter_dat = m_data.begin();
  m_speakerlabel.resize(0);
  
  for (; iter_dat != m_data.end(); ++iter_dat)
    if (*iter_dat) delete (*iter_dat);
}

//----------------------------------------------------------------------

void CSegment::freeSegment(const int _i)
{
  if (_i >= m_numsegments)
    ERROR(1111, "freeSegment: Invalid index of segment");
  m_speakerlabel.erase(m_speakerlabel.begin() + _i);

  if (m_data[_i]) delete m_data[_i];
}

//----------------------------------------------------------------------

MATRIX* CSegment::readHTKFeature(const char* _filename)
{
  htktools::CHTKData* htk = new htktools::CHTKData(_filename);

  int dimension    = htk->getDimension();
  int num_dframes  = htk->getNumSamples();

  int num_frames = num_dframes;

  MATRIX* feature = new MATRIX(num_frames, dimension);
  for (int t = 0; t < num_frames; ++t)
    for (int d = 0; d < dimension; ++d)
      (*feature)(t, d) =  htk->getData(t, d);
  
  delete htk;
  
  return feature;
}
//----------------------------------------------------------------------

MATRIX* CSegment::readCSVFeature(const char* _filename)
{
  ifstream ifs(_filename);
  if (ifs)
    ERROR(1111, "Cannot read CSV file: %s", _filename);

  vector<vector<double> > f_mat;

  string str;
  int numframes = 0;
  int dimension = 0;
  while(getline(ifs, str))
  {
    /* 一列 (frame) 読み込み */
    int dim = 0;
    string sub = str;
    DoubleVect f_vec;
    string::size_type index = 0;
    while (sub.size() > 0 && index != string::npos)
    {
      /* 各次元の特徴量を読み込む */
      index = sub.find_first_of(',');
      double fvalue = atof(sub.substr(0, index).c_str());
      f_vec.push_back(fvalue);
      sub = sub.substr(sub.find_first_of(',') + 1);
      dim++;
    }
    if (numframes > 0)
      if (dim != (f_mat.back()).size())
        ERROR(1111, "Invalid dimension: %d", dimension);

    f_mat.push_back(f_vec);
    numframes++;
    dimension = dim;
  }

  MATRIX* feature = new MATRIX(numframes, dimension);

  vector<vector <double> >::iterator iter_f = f_mat.begin();
  for (int i = 0; iter_f != f_mat.end(); ++iter_f, ++i)
  {
    DoubleVect::iterator iter_ft = iter_f->begin();
    for (int j = 0; iter_ft != iter_f->end(); ++iter_ft, ++j)
      (*feature)(i, j) =  (*iter_ft);
  }
  
  return feature;
}

//----------------------------------------------------------------------

void CSegment::readFeatureFromSCP(const char* _filename)
{

  // スクリプトファイル(.scp)からファイル名と話者idを読み込む
  ifstream ifs(_filename);
  if (!ifs)
    ERROR(1111, "Cannot read scp file: %s", _filename);

  int isfirst          = true;
  int pre_dimension    = 0;
  int seg_index        = 0; // 読み込み済みのデータ数
  string pre_filename  = "";
  MATRIX* feature = 0;
  string str;

  while (ifs >> str)
  {
    int s_time = 0, e_time = 0;
    string speaker, filename;

    vector<string> v;
    boost::algorithm::split( v, str, boost::is_any_of("="));
    speaker  = v[0];
    filename = v[1];
    boost::algorithm::split( v, v[1], boost::is_any_of("["));
    if (v.size()>1)
    {
      filename = v[0];
      boost::algorithm::split( v, v[1], boost::is_any_of(","));
      s_time = atoi(v[0].c_str());
      boost::algorithm::split( v, v[1], boost::is_any_of("]"));
      e_time = atoi(v[0].c_str());
    }

   /*
    boost::regex re("=|,|\\]|\\[");
    vector<string> v;
    boost::regex_split( back_inserter(v), str, re);  
    string speaker  = v[0];
    string filename = v[1];
    */

    try{
      // 前に読み込んだファイルと違うファイルの場合新たに読み込み直す
      if (filename != pre_filename){
	if (feature) delete feature;
	feature  = readFeature(filename.c_str());
      }
    } catch (string errmsg) {
      ERROR(1111, "%s (at line %d in %s)", errmsg.c_str(), seg_index+1, _filename);
    }

    int numframes  = (*feature).rows();
    int dimension  = (*feature).cols();

    if (e_time==0) 
      e_time = numframes - 1;
	
    if (!isfirst && pre_dimension != dimension)
      ERROR(1111, "Invalid dimension: %d (at line %d in %s)", 
	    dimension, seg_index+1, _filename);

    if (s_time >= numframes  || e_time >= numframes)
      ERROR(1111, "Out of index: [%d, %d] (at line %d in %s), frame-length: %d", 
	    s_time,e_time, seg_index+1, _filename, numframes);
    // 読み込んだ行列の指定されたframeを作業領域にコピー
    addSegment(*feature, s_time, e_time, speaker); 

    pre_dimension = dimension;
    pre_filename  = filename;    
    isfirst       = false;
    seg_index++;
  } // end of while (ifs >> str)
  if (feature) delete feature;
  m_dimension  = pre_dimension;
}

//----------------------------------------------------------------------

void CSegment::readFeatureFromScript(const char* _filename)
{
  // スクリプトファイル(.diag)からファイル名と話者idを読み込む
  FileList script(_filename);

  int first            = 1;
  int pre_dimension    = 0;
  int seg_index        = 0; // 読み込み済みのデータ数
  while (!script.isEndList())
  {
    /* ファイル名とラベル名を分離 */
    char del[] = "=";
    char file[(int)script.get().size()];
    strcpy(file, script.get().c_str());
    char *sub      = strstr(file, del);
    if (!sub)
      ERROR(1111, "Invalid file format at line %d: lack of '=label'", seg_index+1);
    char *speaker  = sub + 1;
    (*sub) = '\0';
    
    MATRIX* feature = readFeature(file);
    int dimension   = (*feature).cols();
    if (!first && pre_dimension != dimension)
      ERROR(1111,"Invalid dimension: %d", dimension);
    
    addSegment(*feature, speaker);    // 読み込んだデータを作業領域にコピー
    delete feature;
    
    pre_dimension = dimension;
    (script)++;
    
    first = 0;
    seg_index++;
  }
  
  m_dimension  = pre_dimension;
  cout << "  dimension: "         << getDimension()    << "\n";
  cout << "  # of segments: "     << getNumSegments()  << "\n";
  cout << "  total # of frames: " << getNumAllFrames() << "\n";
  
}

//----------------------------------------------------------------------

MATRIX* CSegment::readFeature(const char* _filename) throw (string)
{
  ifstream ifs_temp(_filename, ios::binary);
  if (!ifs_temp)
    throw ("Cannot read feature file: " + string(_filename));

  // 先頭の 1byte * 16 = 16byte（HTK 形式の header の長さ）を先行読み込み
  char *tmp = (char*)calloc(16,sizeof(char)); 
  for (int i = 0; i < 16; ++i) tmp[i] = 0xFF;
  ifs_temp.read( reinterpret_cast<char *>(tmp),  sizeof(char)*16);
  ifs_temp.close();
  int isbin = 0;
  // 読み込んだ先頭の 16bit でフォーマットチェック
  for (int i = 0; i < 16; ++i)
  {
    if (tmp[i] == '\377') break;
    if (tmp[i] != 'e' && tmp[i] != '-' && tmp[i] != ',' &&
	tmp[i] != ' ' && tmp[i] != '.' && tmp[i] != '\n' && 
	tmp[i] != '\r' && !isdigit(tmp[i]))
      isbin = 1;
  }
  delete tmp;

  /* 発話データ読み込み */
  MATRIX* feature;
  if (isbin)
    feature = readHTKFeature(_filename);
  else
    feature = readCSVFeature(_filename);

  return feature;
}

//----------------------------------------------------------------------
// 指定されたフレーム区間を追加
void CSegment::addSegment(const MATRIX& _ss, const int _sframe, const int _eframe,
			  const string _spk)
{
  int dimension = _ss.cols();
  int numframes = _ss.rows();
  int duration  = _eframe - _sframe + 1;
  if (_eframe > numframes)
    ERROR(1111, "Out of index: [%d, %d]", _sframe, _eframe);

  MATRIX* tmp = new MATRIX(duration, dimension);
  (*tmp)      = _ss.block(_sframe, 0, duration, dimension);
  m_data.push_back(tmp);
  m_numframes.push_back(duration);
  m_numallframes += duration;
  m_numsegments++;

  m_speakerlabel.push_back(_spk);
}

//----------------------------------------------------------------------
// 全フレームを追加
void CSegment::addSegment(const MATRIX& _ss, string _s)
{
  addSegment(_ss, 0, _ss.rows(), _s);
}

//----------------------------------------------------------------------

MATRIX CSegment::getFrame(const int _i, const int _j)
{
  if (_i >= m_numsegments)
    ERROR(1111, "Out of index: %d / %d", _i, m_numsegments);
  if (_j >= m_numframes[_i])
    ERROR(1111, "Out of frame index: %d / %d in utterance %d", 
	  _j, m_numframes[_i], _i);
  return m_data[_i]->row(_j);
}

//----------------------------------------------------------------------

const MATRIX* CSegment::getSegment(const int _i)
{
  if (_i >= m_numsegments)
    ERROR(1111, "Out of index: %d", _i);
  return m_data[_i];
}

//----------------------------------------------------------------------

int CSegment::getNumFrames(const int _i)
{
  if (_i >= m_numsegments)
    ERROR(1111, "Out of index: %d", _i);
  return m_numframes[_i];
}

//----------------------------------------------------------------------

void CSegment::dimensionSelection(const int _sid, const int _new_dim)
{
  // 0 ~ new_dim の次元のみを取り出す
  int numsegments   = getNumSegments();
  int numallframes  = getNumAllFrames();
  int org_dimension = getDimension();

  if (_sid < 0 || _new_dim <= 0 || _sid + _new_dim - 1 >= org_dimension)
    ERROR(1111, "Invalid index of segment");

  vector<MATRIX*>::iterator iter_dat = m_data.begin();
  for (int u = 0; iter_dat != m_data.end(); ++u, ++iter_dat)
  {
    MATRIX* sub_data = new MATRIX(getNumFrames(u), _new_dim);
    for (int t = 0; t < getNumFrames(u); ++t)
      for (int d = 0; d < _new_dim; ++d)
	(*sub_data)(t, d) = (**iter_dat)(t, d);	
    delete (*iter_dat);
    (*iter_dat) = sub_data;
  }
  m_dimension  = _new_dim;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifdef DEBUG
int main(int _argc, char* _argv[])
{
  char* filename =_argv[1];
  cout << filename << endl;
  
  CSegment* segment = new CSegment(filename);
  int numsegments = segment->getNumSegments();
  int dimension   = segment->getDimension();

  cout << "Num segments: " << numsegments << endl;
  cout << "dimension: "  <<  dimension  <<endl;
  for (int i = 0; i < numsegments; ++i)
  {
    cout << i << ": numframes = " << segment->getNumFrames(i) <<endl;
  }
  cout << (segment->getFrame(0,0)) << endl;
  delete segment;
  return 0;
}

#endif


