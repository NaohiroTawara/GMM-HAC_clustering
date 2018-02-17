#include "htkdata.h"

using namespace std;

namespace htktools {

struct _htkHeader* CHTKData::readHtkHeader(float *_tmp)
{
  // short int : 16bit (2byte)
  // int       : 32bit (4byte)
  static struct _htkHeader header;
  
  header.nSamples   = *((int *)_tmp);              // フレーム数(4byte)
  header.sampPeriod = (int) *((int *)_tmp+1);      // フレーム周期(4byte)
  header.sampSize   = (short int) *((short int*)_tmp+4); // １フレーム当りのバイト数(2byte)
  header.paramKind  = (short int) *((short int*)_tmp+5); // データの種類を指定するコード番号(2byte)
#if 0
  cout << "header.nSamples   = " << header.nSamples<<endl;
  cout << "header.sampPeriod = " << header.sampPeriod<<endl;
  cout << "header.sampSize   = " << header.sampSize<<endl;
  cout << "header.paramKind  = " << header.paramKind<<endl;
#endif
  return &header;
}

int CHTKData::getLength(ifstream &_ifs, const int channels, const int size)
{
  //seek．
  _ifs.seekg(0,ios::end);
  int len = _ifs.tellg() / ( size*channels);
  
  _ifs.seekg(0,ios::beg);
  return len;
}

unsigned int CHTKData::swapint(unsigned int _x)
{
  // x & 0000FF00 は上位16ビットと下位8ビットをクリアする行為
  return ((_x << 24) | ((_x & 0x0000FF00) << 8) |
            ((_x & 0x00FF0000) >> 8) | (_x >> 24));
}
/*
   0A 0B 0C 0D
        v
   0D 0C 0B 0A
  
  */
float CHTKData::swapfloat(float* _x)
{
  // float型をうまくキャストするために整数型のポインタベースにする
  unsigned int *temp_i = reinterpret_cast<unsigned int*>(_x);
  volatile unsigned int value = swapint(*temp_i);
  volatile float temp_f = *( reinterpret_cast<volatile float*>(&value));
  *_x = temp_f;
  return 0;
}

int CHTKData::readData(const char* _filename)
{
  ifstream ifs_temp( _filename, ios::binary);
  if (!ifs_temp) {
    cerr << "ERROR[CHTKData::readData()]: cannot read HTK file: " << _filename << endl;
    return -1;
  }
//  cout <<_filename<<endl;
  int size = getLength(ifs_temp, 1, sizeof(float));
  float *tmp = new float[size];
//  cout <<"size: "<<size * sizeof(float)<<endl;

  ifs_temp.read( reinterpret_cast<char *>(tmp), sizeof(float) * size);
  ifs_temp.close();

  // エンディアン変換
  // 4byte (32bits) の順番を反転
  for(int t=0; t<size; ++t)
    swapfloat(&tmp[t]);
  // 先頭の12byte はヘッダー

  struct _htkHeader *header =readHtkHeader(tmp);
  
  int frame_num = header->nSamples;
  int dim = (size - HEADER_SIZE) / frame_num;
  /*
  cout <<"frame_num = "<<frame_num<<endl;
  cout <<"dim       = "<<dim<<endl;
    */
  m_data.resize(frame_num, vector <float>(dim));

  for(int i = 0; i < frame_num; ++i)
    for(int j = 0; j < dim; ++j)
      m_data[i][j] = tmp[i*dim+3+j];

  m_dimension   = dim;
  m_num_samples = frame_num;

  delete tmp;
  return 0;
}

}

/*
int main(int _argc, char** _argv)
{
  CHTKData* htk = new CHTKData(_argv[1]);
  delete htk;
}
*/



/* end of file */

