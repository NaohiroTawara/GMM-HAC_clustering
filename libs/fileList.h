/* ###############################################
** FileList
** スクリプトファイルを管理するクラス
**                copy right N.Tawara
**                           2010/8/4
** ###############################################
*/


#ifndef __FILELIST_H__
#define __FILELIST_H__
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>

using namespace std;

class FileList
{
  
 public:
  
  FileList(string _filename)
       : m_list_filename(_filename), m_ite(0)
    {
      ifstream ifs(_filename.c_str());
      if (!ifs)
      {
        cerr << "error[fileList]: cannot read file: " << _filename << endl;
        exit(-1);
      }
      string name;
      while (ifs >> name)
        m_data_filename.push_back(name);
      m_length = m_data_filename.size();
    }
  
  ~FileList() {};

  int isEndList() {return (m_ite==m_length);}
  
  string get()
  {
    if (m_length<=0 || m_ite>=m_length)
    {
      cerr << "ERROR [fileList.h]: invalid vector size:"<< m_ite<<endl;
      exit(-1);
    }
    else
      return m_data_filename[m_ite];
  }

  void operator ++(int n)
  { m_ite++; }

  FileList operator ++()
  { ++m_ite; return *this; }
  
  int m_ite;

private:
  string m_list_filename;
  vector <string> m_data_filename;
  int m_length;

};

#endif
