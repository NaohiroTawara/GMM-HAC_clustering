//
//  lvcontainer.h
//  VB
//  latentvariables をバイナリ形式で保存・読み込むためのクラス
//
//  utterance レベル潜在変数：各発話の，各クラスタに対する，割当て尤度が最大となるクラスタ番号を保存
//    データ容量：発話数 x integer(クラスタ番号を保持)
//  frame レベル潜在変数：各発話に含まれる各フレームの，各クラスタの各コンポーネントに対する，
//                     割当て尤度が最大となるコンポーネント番号を保存
//    データ容量：フレーム数 x クラスタ数 x integer(コンポーネント番号を保持)
//
//  Created by 俵 直弘 on 12/12/29.
//  Copyright (c) 2012年 俵 直弘. All rights reserved.
//

#ifndef VB_lvcontainer_h
#define VB_lvcontainer_h

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

/*
 format:

 Offset Size Name Description
 
 0  1   Endian      0:Big endian 1: Little endian
 1  4   version     Version of .lv file
 5  4   NumSeg      Number of segments
 9  4   FileSize    Size of file
 13 4   NumCluster  Number of clusters
 17 1   DesType     Description type of each variables
                    (0: index 1: likelihood)
 18 4   NumFrame    Number of frames in each segment
 
 -   -   NumMixtures number of mixtures in each clustersœ
 for j=0:SegNum-1
 0   4     SegSize     numbef of samples in this segment (4 byte)
 4   4     u-LV        segment level latent variables (4 byte)
 ================= for i=0:SegSize-1 ============================
 8+i   4   f-LV        sample  level latent variables (4 byte)}
 ================================================================
 number of utterances(4byte)
 number of frames (
 
 */

typedef vector<double> DoubleVect;


const int HEADER_SIZE = 3;

struct _lvHeader {
  unsigned char endian;        // 1byte
  unsigned short int  version; // 2byte
  unsigned short int NumSeg;   // 2byte
  unsigned short int DesType;  // 2byte
};



class CLVContainer
{
public:

  CLVContainer() {};
  
  ~CLVContainer() {};
  

  /* 潜在変数を保存する */
  void saveULV(const char* _filename, DoubleVect& _ulv);
  
  void saveFLV(const char* filename, vector<DoubleVect>& _flv);

  void saveULV(const char* filename, vector<DoubleVect>& _ulv);
  void saveFLV(const char* filename, vector<vector<DoubleVect> >& _flv);

  /* 潜在変数を読み込む */
  void loadULV(const char* filename, DoubleVect& _ulv);
  void loadFLV(const char* filename, vector<DoubleVect>& _flv);
  
  void loadULV(const char* filename, vector<DoubleVect>& _ulv);
  void loadFLV(const char* filename, vector<vector<DoubleVect> >& _flv);

  float* writeLVHeader(struct _lvHeader _header);
  
private:

  struct _lvHeader* readLVHeader(float *tmp);
  int getLength(ifstream &_ifs, const int channels, const int size);

  
  unsigned swapint(unsigned int);
  float swapfloat(float*);
  
  
};


#endif
