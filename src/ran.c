// ***********************************************************
//
// file:	ran.cpp
//
// author:
//
// created:
//
// purpose:
//
// modified:
//
// ***********************************************************

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define INITDUM -179

extern int RandSeed;

/* Table of constant values */

double ran()
{
    /* Initialized data */

    static long int iff = 0;
    static long int idum = INITDUM;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    static long int j;
    static double r__[97];
    static long int ix1, ix2, ix3;

    if (idum < 0 || iff == 0) {
/*********************************************************
	FILE *fp;
	fp=fopen("INITRANF","r");
	if(fp) {
	   if(fscanf(fp,"%d",&idum)!=1) {
		perror("fscanf: ");
		fprintf(stderr,"Random seed will be set to %d\n",INITDUM);
		idum=INITDUM;
	   } else {
		if(idum > -1) {
		   fprintf(stderr,"Invalid random seed: %d\n",idum);
		   fprintf(stderr,"Random seed will be set to: %d\n",INITDUM);
		   idum=INITDUM;
		}
	   }
	}
***********************************************************/
	idum = iff = 1;

if ((RandSeed%2)==0) RandSeed+=1;
if (RandSeed>0) RandSeed=-RandSeed;

//	ix1 = (54773 - idum) % 259200;
	ix1 = (54773 - RandSeed) % 259200;
	ix1 = (ix1 * 7141 + 54773) % 259200;
	ix2 = ix1 % 134456;
	ix1 = (ix1 * 7141 + 54773) % 259200;
	ix3 = ix1 % 243000;
	for (j = 1; j <= 97; ++j) {
	    ix1 = (ix1 * 7141 + 54773) % 259200;
	    ix2 = (ix2 * 8121 + 28411) % 134456;
	    r__[j - 1] = ((double) ix1 + (double) ix2 * (double)
		    7.4373772832748258e-6) * (double)3.8580246913580248e-6;
	}
    }


    ix1 = (ix1 * 7141 + 54773) % 259200;
    ix2 = (ix2 * 8121 + 28411) % 134456;
    ix3 = (ix3 * 4561 + 51349) % 243000;
    j = ix3 * 97 / 243000 + 1;
    if (j > 97 || j < 1) {
	fprintf(stderr,"Error in ran: invalid index: %d\n",j);
	exit(-1);
    }
    ret_val = r__[j - 1];
    r__[j - 1] = ((double) ix1 + (double) ix2 * 7.4373772832748258e-6) *
		   3.8580246913580248e-6;
    return ret_val;
} /* ran_ */

int intran(int *iseed)
//int *iseed;
{
    /* System generated locals */
    int ret_val;

    *iseed = (*iseed * 9301 + 49297) % 233280;
    ret_val = *iseed * 2147483647 / 233280 + 1;
    return ret_val;
} /* intran */



double gaussian (const double sigma)
{
  double x, y, w;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = 2.0 * ran() - 1.0;
      y = 2.0 * ran() - 1.0;

      /* see if it is in the unit circle */
      w = x * x + y * y;
    }
  while (w >= 1.0 || w == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt ( (-2.0 * log(w)) / w );
}
