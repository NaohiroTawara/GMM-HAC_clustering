//
//  lvcontainer.cc
//  VB
//
//  Created by 俵 直弘 on 12/12/29.
//  Copyright (c) 2012年 俵 直弘. All rights reserved.
//

#include "lvcontainer.h"

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

unsigned int CLVContainer::swapint(unsigned int _x)
{
  // x & 0000FF00 §œæÂ∞Ã16•”•√•»§»≤º∞Ã8•”•√•»§Ú•Ø•Í•¢§π§Îπ‘∞Ÿ
  return ((_x << 24) | ((_x & 0x0000FF00) << 8) |
          ((_x & 0x00FF0000) >> 8) | (_x >> 24));
}

//----------------------------------------------------------------------

/*
 0A 0B 0C 0D
 v
 0D 0C 0B 0A
 
 */
float CLVContainer::swapfloat(float* _x)
{
  // float∑ø§Ú§¶§ﬁ§Ø•≠•„•π•»§π§Î§ø§·§À¿∞øÙ∑ø§Œ•›•§•Û•ø•Ÿ°º•π§À§π§Î
  unsigned int *temp_i = reinterpret_cast<unsigned int*>(_x);
  volatile unsigned int value = swapint(*temp_i);
  volatile float temp_f = *( reinterpret_cast<volatile float*>(&value));
  *_x = temp_f;
  return 0;
}

//----------------------------------------------------------------------

int CLVContainer::getLength(ifstream &_ifs, const int channels, const int size)
{
  _ifs.seekg(0,ios::end);
  int len = _ifs.tellg() / ( size*channels);
  
  _ifs.seekg(0,ios::beg);
  return len;
}

//----------------------------------------------------------------------

struct _lvHeader* CLVContainer::readLVHeader(float *_tmp)
{
  // char      :  8bit (1byte)
  // short int : 16bit (2byte)
  // int       : 32bit (4byte)
  static struct _lvHeader header;
  
  header.endian   = *((char *)_tmp);
  header.version  = (int) *((int *)_tmp+1);
  header.NumSeg   = (short int) *((short int*)_tmp+4);
  header.DesType  = (short int) *((short int*)_tmp+5);
  return &header;
}

//----------------------------------------------------------------------

float* CLVContainer::writeLVHeader(struct _lvHeader _header)
{
  float* tmp = (float*)malloc(HEADER_SIZE);
  
//  tmp = (void*)_header.endian;
  tmp[0] = 10;

  //_header.version;
  //_header.NumSeg;
  //_header.DesType;
  return tmp;
}

//----------------------------------------------------------------------

void CLVContainer::saveULV(const char* _filename, DoubleVect& _ulv)
{
  ifstream ifs_temp( _filename, ios::binary);
  if (ifs_temp == NULL)
  {
    cerr << "ERROR[CLVContainer::saveULV()]: cannot read lv file: "
         << _filename << endl;
    exit(-1);
  }
  
  int size   = getLength(ifs_temp, 1, sizeof(float));
  float *tmp = new float[size];

  // •®•Û•«•£•¢•Û —¥π
  // 4byte (32bits) §ŒΩÁ»÷§Ú»ø≈æ
  for(int t = 0; t < size; ++t)
    swapfloat(&tmp[t]);

  
  struct _lvHeader *header =readLVHeader(tmp);
  
}

//----------------------------------------------------------------------

int main(int _argc, char* _argv[])
{
  char* filename = _argv[1];
  int numsegments = 10; // number of segments
  int nummixtures = 8;   // number of mixtures;
  int numframes   = 2;  // number of frames;
  int numclusters = 8;  // number of clusters;
  
  DoubleVect testulv_1(numsegments);          // ULV (index)
  vector<DoubleVect> testulv_2(numsegments);  // ULV (likelihood)
  
  DoubleVect::iterator         iter_lv1 = testulv_1.begin();
  vector<DoubleVect>::iterator iter_lv2 = testulv_2.begin();
  for(; iter_lv1 != testulv_1.end(); ++iter_lv1, ++iter_lv2)
  {
    (*iter_lv1) = rand() % numclusters;
    iter_lv2->resize(numclusters);
    DoubleVect::iterator       iter_lv21 = iter_lv2->begin();
    double sum = 0;
    for(; iter_lv21 != iter_lv2->end(); ++iter_lv21)
    {
      (*iter_lv21) = (double)rand()/INT_MAX;
      sum += (*iter_lv21);
    }

    iter_lv21 = iter_lv2->begin();
    for(; iter_lv21 != iter_lv2->end(); ++iter_lv21)
      (*iter_lv21) /= sum;
  }
  
  iter_lv1 = testulv_1.begin();
  for(; iter_lv1 != testulv_1.end(); ++iter_lv1, ++iter_lv2)
    cout << (*iter_lv1) << ",";
  cout <<endl;
  
  iter_lv2 = testulv_2.begin();
  for(; iter_lv2 != testulv_2.end(); ++iter_lv2)
  {
    DoubleVect::iterator       iter_lv21 = iter_lv2->begin();
    for(; iter_lv21 != iter_lv2->end(); ++iter_lv21)
      printf("%.2f,",*iter_lv21);
    cout << endl;
  }
  
  cout << filename <<endl;
  
  CLVContainer* clv = new CLVContainer();

  static struct _lvHeader header;
  void* vo = clv->writeLVHeader(header);
  cout << vo <<endl;
  delete(clv);
}

