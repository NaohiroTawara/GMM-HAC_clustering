//
//  main.cpp
//  spkr_cl_GMM-HAC
//
//  Created by 直弘 俵 onzcf 12/07/30.
//  Copyright 2012年 早稲田大学. All rights reserved.
//

#include "gmmclustering.h"
#include "segment.h"
#include "evalK.h"

using namespace std;

float const DEF_INIT_MAX_ALPHA = 1.0;
float const DEF_INIT_MIN_ALPHA = 100;
float const DEF_INIT_D_ALPHA   = 1.0;

string g_script_filename = "";     // 発話データのパスが記載されたスクリプトファイル
string g_gmm_filename    = "";     // GMM モデルファイル
string g_out_dirname     = "";     // 結果の出力先のディレクトリ名
int    g_num_clusters    = 1;      // 話者クラスタ数
int    g_num_mixtures    = 1;      // 混合数
int    g_covtype         = DIAGC;  // covariance type (full:1, diagonal:0)
int    g_max_iter        = 10;     // EM の最大イタレーション数
int    g_min_diff        = 0;      // EM の最小更新尤度差
int    g_trace           = 0;      // debug 出力の有無
int    g_true_num_spkr   = 1;      // 真の話者数
float  g_alpha_max       = DEF_INIT_MAX_ALPHA; // 最大閾値
float  g_alpha_min       = DEF_INIT_MIN_ALPHA; // 最小閾値
float  g_alpha_step      = DEF_INIT_D_ALPHA; // 閾値のステップ幅

vector<string> g_spkr_label;

static const char* USAGE_MESSAGE =
"\nUSAGE: ./spkr_cl_GMM-HAC [option...]"
"\n"
" Option                                         Default\n"
"\n"
" -a f      Set tuning parameter to f (minimum)     1.0\n"
" -b f      Set tuning parameter to f (maximum)     100\n"
" -d f      Set step of tuning parameter to f       1.0\n"
" -f N      Set covariance type flag to N           0\n"
" -h        show this message\n"
" -i N      Set EM iterations to N                  10\n"
" -m N      Set number of mixtures to N             1\n"
" -s N      *Set true number of speakers            1\n"
" -M mmf    Set GMM model file to f                 none\n"
" -O mmf    *Set output dirname to mmf              none\n"
" -S f      *Set script filename to f               none\n"
" -T N      Set trace flags to N                    0\n"
"\n"
"*: required\n"
;

void Usage(void);
void GetOption(int _argc, char* _argv[]);

//----------------------------------------------------------------------

void Usage(void)
{
    cerr << USAGE_MESSAGE;
    exit(1);
}


//----------------------------------------------------------------------

void GetOption(int _argc, char* _argv[])
{
    extern char* optarg;
    
    int opt;
    while ((opt = getopt(_argc, _argv, "a:b:d:f:h:i:m:s:M:O:S:T:")) != -1)
    {
        switch (opt)
        {
            case 'a':
                g_alpha_min  = atof(optarg);
                break;
            case 'b':
                g_alpha_max  = atof(optarg);
                break;
            case 'd':
                g_alpha_step = atof(optarg);
                break;
            case 'f':
                g_covtype = atoi(optarg);
                break;
            case 'i':
                g_max_iter = atoi(optarg);
                break;
            case 'm':
                g_num_mixtures = atoi(optarg);
                break;
            case 's':
                g_true_num_spkr = atoi(optarg);
                break;
            case 'M':
                g_gmm_filename = optarg;
                break;
            case 'O':
                g_out_dirname = optarg;
                break;
            case 'S':
                g_script_filename = optarg;
                break;
            case 'T':
                g_trace = atoi(optarg);
                break;
            case 'h':
            default:
                Usage();
        }
    }
  
  if (g_true_num_spkr <= 0)
    Error(1111, "main.cc: please set the valid value for option '-s'");
  if (g_out_dirname == "")
    Error(1111, "main.cc: please set the valid value for option '-O'");
  if (g_script_filename == "")
    Error(1111, "main.cc: please set the valid value for option '-S'");
}

//----------------------------------------------------------------------

int main (int _argc, char * _argv[])
{
    double t1, t2;
    GetOption(_argc, _argv);
  
  
    /* ================================ */
    /* 特徴量の読み込み */
    /* ================================ */
  
    cout << "Reading feature files " << " << " <<endl
         << "   " << g_script_filename << endl;
    CSegment* feature = new CSegment(g_script_filename.c_str());
    cout << "done. " << endl;

    int numsegments = feature->getNumSegments();
    for (int i = 0; i < numsegments; ++i)
      g_spkr_label.push_back(feature->getSpeakerLabel(i));

    /* ================================ */
    /* インスタンス生成 */
    /* ================================ */
    CGMMClustering* gmm_hac
      = new CGMMClustering(numsegments, g_covtype, &cout);

    /* ================================ */
    /* 読み込んだ特徴量のポインタをセットする */
    /* ================================ */
    gmm_hac->setFeature(feature);
    
    /* ================================ */
    /* クラスタを初期化する */
    /* ================================ */

    cout << "Num of mixtures: " << g_num_mixtures <<endl;
    gmm_hac->initClusters(g_num_mixtures, g_max_iter, g_min_diff);

    /* ================================ */
    /* 初期類似度行列を出力する */
    /* ================================ */
    {
      stringstream out_filename;
      out_filename << g_out_dirname << "/initCM.csv";
      ofstream ofs_out((out_filename.str()).c_str());
      ofs_out << gmm_hac->getDistMat();
      ofs_out.close();
    }
    /* ================================ */
    /* GMM 初期モデルファイルの読み込み */
    /* ================================ */
    /*
    cout << "Reading basis model file " << " << " << endl
         << "   " << g_gmm_filename << endl;
    gmm_hac->setBasisFromFile(g_basis_filename.c_str());
    cout << "done. " << endl;
*/    
    
   gmm_hac->setTrace(g_trace);
    /* ================================ */
    /* しきい値を変えながらクラスタリングの実行 */
    /* ================================ */

    int pre_num_clusters = -1;
    t1 = gettimeofday_sec();
    for (float alpha = g_alpha_min; alpha < g_alpha_max; alpha += g_alpha_step)
    {
        cout << "alpha: " << alpha << endl;
        gmm_hac->setAlpha(alpha);
        
        gmm_hac->run();
        
        /* ================================ */
        /* クラスタリング結果の取得 */
        /* ================================ */
        
        IntVect data_cc(gmm_hac->getNumSegments());
        gmm_hac->getClusteringResultSeg(&data_cc);
        
        /* ================================ */
        /* 評価 */
        /* ================================ */
        int num_clusters = gmm_hac->getNumClusters();
        int num_segments = gmm_hac->getNumSegments();
        
        if (num_clusters != pre_num_clusters)
        {
            t2 = gettimeofday_sec();
            
            DoubleVect length_vec(num_segments);
            DoubleVect::iterator iter_l = length_vec.begin();
            for (int i = 0; iter_l != length_vec.end(); ++iter_l, ++i)
                (*iter_l) = gmm_hac->getNumFrames(i);

            cout << "class: "  << num_clusters;
            CSpkrClEvaluation_K* evalK =
                new CSpkrClEvaluation_K(num_segments, 
                                        num_clusters, 
                                        g_true_num_spkr);
          
            evalK->setMlf(g_spkr_label, length_vec);
          
            Result resultK = evalK->evaluate(data_cc);
            cout << ", acp: "  << resultK.result["acp"] 
                 << ", asp: "  << resultK.result["asp"] 
                 << ", K: "    << resultK.result["K"] 
                 << ", time: " << t2 - t1
                 << endl;
            delete evalK;
             
            /* ================================ */
            /* 結果を保存 */
            /* ================================ */
          
            stringstream out_filename;
            out_filename << g_out_dirname << "/" << alpha << ".txt";
//            out_filename << "./out/" << alpha << ".txt";
            ofstream ofs_out((out_filename.str()).c_str());
            ofs_out << "class: " << num_clusters << endl;
            ofs_out << "acp:   " << resultK.result["acp"]   << endl;
            ofs_out << "asp:   " << resultK.result["asp"]   << endl;
            ofs_out << "K:     " << resultK.result["K"]     << endl;
            ofs_out << "time:  " << t2 - t1                 << endl;
            ofs_out << endl;

            ofs_out << gmm_hac->showClusteringResult();
	    ofs_out << endl;
	    ofs_out << "Similality matrix" << endl;
	    ofs_out << gmm_hac->getDistMat();
            ofs_out.close();
            
            pre_num_clusters = num_clusters;
	    

        }
        
    }
    

    delete gmm_hac;
    delete feature;
    
    return 0;
}

