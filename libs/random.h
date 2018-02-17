//
//  random.h
//  MCMC
//
//  Created by 俵 直弘 on 13/03/04.
//  Copyright (c) 2013年 俵 直弘. All rights reserved.
// References:
// [1]  L. Devroye, "Non-Uniform Random Variate Generation",
// Springer-Verlag, 1986
// http://cgm.cs.mcgill.ca/~luc/rnbookindex.html


#ifndef __MCMC__random__
#define __MCMC__random__

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define RANDOM ((double)rand()/(double)RAND_MAX)

#if 0
inline double normrand( double mu, double sigma )
{
  double z = sqrt( -2.0*log(RANDOM)) * sin( 2.0 * M_PI * RANDOM );
  return mu + sigma * z;
}
#endif


inline double exprand (void)
{
  return (- log(RANDOM));
}


inline double gamrand (double a)
{
  double x, y, z;
	double u, v, w, b, c, e;
	int accept = 0;
	if (a < 1)
	{
		/* Johnk's generator. Devroye (1986) p.418 */
		e = exprand();
		do {
			x = pow(RANDOM, 1 / a);
			y = pow(RANDOM, 1 / (1 - a));
		} while (x + y > 1);
		return (e * x / (x + y));
	} else {
		/* Best's rejection algorithm. Devroye (1986) p.410 */
		b = a - 1;
		c = 3 * a - 0.75;
		do {
			/* generate */
			u = RANDOM;
			v = RANDOM;
			w = u * (1 - u);
			y = sqrt(c / w) * (u - 0.5);
			x = b + y;
			if (x >= 0)
			{
				z = 64 * w * w * w * v * v;
				if (z <= 1 - (2 * y * y) / x)
				{
					accept = 1;
				} else {
					if (log(z) < 2 * (b * log(x / b) - y))
						accept = 1;
				}
			}
		} while (accept != 1);
		return x;
	}
}

inline void dirrand (vector<double>* theta, vector<double>* alpha, double prec = 0)
{
	double z = 0;
  vector<double>::iterator iter_t = theta->begin();
  vector<double>::iterator iter_a = alpha->begin();
  for (; iter_t != theta->end(); ++iter_t, ++iter_a)
  { /* theta must have been allocated */
		if (prec != 0)
			(*iter_t)= gamrand((*iter_a) * prec);
		else
			(*iter_t) = gamrand(*iter_a);
		z += (*iter_t);
  }
  for (iter_t = theta->begin(); iter_t != theta->end(); ++iter_t)
		(*iter_t) /= z;
}





#endif /* defined(__MCMC__random__) */
