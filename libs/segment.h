/*
* ---------------------------------------------------------------------   
*   segment.h
*
*   < Author >　N.TAWARA  2014/7/31
*               発話セットとその発話者ラベルを管理するクラス
*
*   < Note > 現状: 発話データの実体は CvMat 形式で保持する
*
*            変更点：2011/11/24 読み込み機能も内包
*                   2014/05/09 diag形式のスクリプトファイルに加え，
*                              scp形式のスクリプトファイルに対応
*                   2014/07/31 データの形式をCvMatからMatrixXd (eigen)に変更
*
*            ToDo:
*  
*            備考： 
*---------------------------------------------------------------------
*
*/



#ifndef __SEGMENT_H__
#define __SEGMENT_H__


#include <string.h>
#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include "fileList.h"
#include "htkdata.h"
#include "matrix.h"
#include "util.h"

class CSegment
{
private:
  int m_dimension;    // 特徴量の次元
  int m_numsegments;  // 発話数
  int m_numallframes; // 総フレーム数

  vector<int> m_numframes;    // [1 x # segments]

  vector <MATRIX*> m_data;  // [# segments x [# frames x dimension] ]
  vector <string> m_speakerlabel; // 正解クラスタ系列

public:
  // ========================================
  // コンストラクタ・デコンストラクタ
  // ========================================
  CSegment(const char* _filename)
    : m_dimension(0), 
      m_numsegments(0),
      m_numallframes(0),
      m_numframes(0),
      m_speakerlabel(0)
    {
      //       readFeatureFromScript(_filename);
       readFeatureFromSCP(_filename);
    }
  
  CSegment(int _dimension)
    : m_dimension(_dimension), 
      m_numsegments(0),
      m_numallframes(0),
      m_numframes(0),
      m_speakerlabel(0)
  {}

  ~CSegment()
    { 
      free();
    };

 private:
      
  void free();

  /*
   * @brief           : diag 形式ファイルから発話データを読み込む
   * @param _filename : スクリプトファイル名
   */
  void readFeatureFromScript(const char* _filename);

  /*
   * @brief            : scp 形式ファイルから発話データを読み込む
   * @param  _filename : スクリプトファイル名
   */
  void readFeatureFromSCP(const char* _filename);

  /*
   * @brief            : 指定されたファイルの存在確認とフォーマットチェック後に読み込む
   * @param  _filename : ファイル名
   * @exception        : ファイルが存在しない場合に例外
   */
  MATRIX* readFeature(const char* _filename) throw (string);

  /*
   * @brief     : m_data の最後尾に 発話 _ss を追加する
   * @param _ss    : 追加する発話行列
   * @param _stime : 開始frame
   * @param _etime : 終了frame
   * @param _s     : 話者ラベル
   */
  void addSegment(const MATRIX& _ss, const int _stime, const int _etime, 
		  const string _s);

  void addSegment(const MATRIX& _ss, string _s);

  //// 外部ファイル（csv 形式）から 1 発話読み込む
  MATRIX* readCSVFeature(const char* _filename);
  //// 外部ファイル（htk 形式）から 1 発話読み込む
  MATRIX* readHTKFeature(const char* _filename);

 public:

  //// 発話データ領域のみを解放する関数
  void freeSegment(const int _i);

 public:
  // ========================================
  // private メンバアクセサ
  // ========================================

  string getSpeakerLabel(const int _i)
  { return m_speakerlabel[_i];}
  
  int getDimension()
  { return m_dimension;}

  int getNumSegments()
  { return m_numsegments; }

  int getNumAllFrames()
  { return m_numallframes;}

  int getNumFrames(const int _i);

  /*
   * @brief   : 発話番号を指定して，発話集合を得るための関数
   * @param _i: 発話番号
   * @return  : [# frames x m_dimension]
   */
  const MATRIX* getSegment(const int _i);
  
  /*
   * @brief    : 発話番号とフレーム番号を指定して直接特徴量を得るための関数
   * @param _i : 発話番号
   * @param _j : フレーム番号
   * @return   : [1 x m_dimension] のサイズを持つ部分行列のポインタ
   */
  MATRIX getFrame(const int _i, const int _j);

  /*
   * @brief     : m_data の 先頭から _i 番目に発話 _ss を追加する
   * @param _ss : 追加する発話
   * @param _i  : 追加する場所
   */
  //  void setSegment(CvMat* _ss, const int _i);
  
  vector<int> getNumFramesVector()
  { return m_numframes; }
  
  //  void compression(const int _dim);
  
  void dimensionSelection(const int _sid, const int _new_dim);
};

#endif
