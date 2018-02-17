/*
*
* ---------------------------------------------------------------------
*
*  htktools::htkdata.h
*  
* ---------------------------------------------------------------------  
*   < Author >　N.TAWARA  2011/08/04
*
*   < Note > 現状：データ読み込み系関数
*            変更点: 2011/08/04: 完成
*            ToDo: 
*            備考： 64bit環境では動かないかも ->修正済み
*/

#ifndef __HTKDATA_H__
#define __HTKDATA_H__

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <sstream>
#include <math.h>

using namespace std;

namespace htktools{
 
// HTK のヘッダーサイズ
const int HEADER_SIZE = 3;

// HTK特徴量データのヘッダー情報を格納する構造体
struct _htkHeader {
  unsigned int nSamples;        // フレーム数
  unsigned int sampPeriod;      // フレーム周期
  unsigned short int sampSize;  // １フレーム当りのバイト数
  unsigned short int paramKind; // データの種類を指定するコード番号
};

class CHTKData
{
public:

  int m_dimension;
  int m_num_samples;
  vector<vector <float> > m_data;

public:

  CHTKData(const char* _filename)
    { if (readData(_filename)) exit(-1);};

  // 読み込んだファイルの _m 行 _n 列目の値を返す
  float getData(const int _m, const int _n)
    { return m_data[_m][_n];}

  // 読み込んだファイルの次元（列数）を返す
  int getDimension()
    { return m_dimension;}

  // データ数を返す
  int getNumSamples()
    { return m_num_samples;}

  // HTK特徴量データのヘッダーを読み込む
  struct _htkHeader* readHtkHeader(float *tmp);



private:

  unsigned swapint(unsigned int);
  float swapfloat(float*);
  int getLength(ifstream&, const int, const int);

  // htk ファイルを読み込む
  int readData(const char* _filename);
  
};

}


#endif
