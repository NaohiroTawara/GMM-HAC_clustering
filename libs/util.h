/*
** --------------------------------------------------
**  util.h
**
**  あると便利な関数群
**  copy right N.TAWARA
**             2011/8/4
** --------------------------------------------------
*/

#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

static const double PI = 6 * asin( 0.5 );
static const double logPI  = log(6 * asin( 0.5 ));
static const double log2PI = log(2) + log(6 * asin( 0.5 ));

static const bool abortOnError = false;

static const int null = 0;
//static const char *defargs[2]  = { "<Uninitialised>", "" };
//static const char **arglist    = defargs;

#define ERROR(id, fmt, ...) Error(id, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

#define _abs(X) (X>=0)?(X):(-X)

enum {DIAGC, FULLC}; /* 0: diagonal, 1: full covariance */


typedef vector<double> DoubleVect;
typedef vector<int> IntVect;
//typedef vector<CvMat*> MatVect;
typedef vector<DoubleVect> DoubleMat;

inline double gettimeofday_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

static void Error(const int errcode, const char* file, 
		  const char* function, const int line,
		  const char *message, ...)
{
   va_list ap;             /* Pointer to unnamed args */
   FILE *f;

   fflush(stdout);        /* Flush any pending output */
   va_start(ap, message);
   if (errcode <= 0)
   {
     fprintf(stdout," WARNING [%+d] in %s::  ", errcode, function);
      f = stdout;
      vfprintf(f, message, ap);
      va_end(ap);
//      fprintf(f," in %s\n", arglist[0]);
   }
   else
   {
     fprintf(stderr, "  ERROR [%+d] in %s::%s() at %d  ", 
	     errcode, file, function, line);
      f = stderr;
      vfprintf(f, message, ap);
      va_end(ap);
      fprintf(f, "\n");
//      fprintf(f,"\n FATAL ERROR - Terminating program %s\n", arglist[0]);
   }
   fflush(f);
   if (errcode > 0)
   {
      if (abortOnError) abort();
      else exit(errcode);
   }
}


static void Error(const int errcode, const char *message, ...)
{
   va_list ap;             /* Pointer to unnamed args */
   FILE *f;

   fflush(stdout);        /* Flush any pending output */
   va_start(ap, message);
   if (errcode <= 0)
   {
      fprintf(stdout," WARNING [%+d]  ",errcode);
      f = stdout;
      vfprintf(f, message, ap);
      va_end(ap);
//      fprintf(f," in %s\n", arglist[0]);
   }
   else
   {
      fprintf(stderr, "  ERROR [%+d]  ", errcode);
      f = stderr;
      vfprintf(f, message, ap);
      va_end(ap);
      fprintf(f, "\n");
//      fprintf(f,"\n FATAL ERROR - Terminating program %s\n", arglist[0]);
   }
   fflush(f);
   if (errcode > 0)
   {
      if (abortOnError) abort();
      else exit(errcode);
   }
}

// token と target が一致しているか調べる
inline void checkFormat(string _token, string _target)
{
  if (_token != _target)
    Error(1, "checkFormat: Invalid format: \"%s\" Requested: \"%s\"",
          _token.c_str(), _target.c_str());
}


/*
** @brief log同士の足し算 sumLogL = log(A+B)
** @param logL1  = log(A)
** @param logL2  = log(B)
*/
inline double LAddS(double &logL1, double &logL2)
{
  double sumLogL; 
  if (logL1==-DBL_MAX || isinf(logL1))
    return logL2;
  if (logL2==-DBL_MAX || isinf(logL2))
    return logL1;
  if (logL1 > logL2)
    sumLogL = logL1 + log(1 + exp(logL2 - logL1));
  else
    sumLogL = logL2 + log(1 + exp(logL1 - logL2));
  return sumLogL;
}
/*
** @brief log同士の引き算 sumLogL = log(A-B)
** @param logL1  = log(A)
** @param logL2  = log(B)
*/
inline double LSubS(double &logL1, double &logL2)
{
  
  if (logL1 < logL2)
    Error(1111, "LSubS: Result is negative: %f - %f = %f",
          logL1, logL2, logL1 - logL2);
  
  if (logL1 == logL2)
    return -DBL_MAX;
  
//  double subLogL;
  if (logL1 == -DBL_MAX || isinf(logL1))
  {
    if (logL2 == -DBL_MAX || isinf(logL2))
      return logL2;
    else
      Error(1111, "LSubS: Result is negative: %f - %f = %f",
            logL1, logL2, logL1 - logL2);
  }  
  if (logL2 == -DBL_MAX || isinf(logL2))
    return logL1;
  
  return logL1 + log(1 - exp(logL2 - logL1));
}

/**
 * @brief  ベクトルの log 要素の合計値を返す
 * @param  _vec ベクトル
 */
inline double  _logsum(vector <double> _vec)
{
  double sum = -DBL_MAX;
  int num = _vec.size();
  for (int i = 0;i < num; ++i)
    sum = LAddS(sum, _vec[i]);
  return sum;
}

inline double digamma(double _value)
{
  const double dx = 0.01;
  return (lgamma(_value + dx) - lgamma(_value)) / dx;
}

#endif
