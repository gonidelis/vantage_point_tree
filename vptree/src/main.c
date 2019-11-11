/* Project no.1 on Parallel and Destributed Systems @AUTh - ECE Dept.
*  GONIDELIS GIANNIS
*  DATE: 03/11/2019
*  ----------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../inc/vptree.h"

#include <sys/time.h>
#include <sys/times.h>

struct timeval startwtime, endwtime;
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

double seq_time;
double p_time;


#define N 1000000
#define D 20         //works for any d


int main()
{

    srand(time(NULL));

    vptree *T ;

    //allocate the proper space for X
    double *X = malloc(D*N*(sizeof(double)));

    //fill X with numbers in space (0,1)
    for(int i = 0 ; i<N ; i++)
    {
        for(int j=0; j<D; j++)
        {
            X[i+j*N] = (double) (rand()  / (RAND_MAX + 2.0));
        }

    }


    /*
    double arr[] = {0.34 , 0.01 , 0.567, 0.21, 0.67, 0.09, 0.99,
        0.1, 0.23, 0.54, 0.002, 0.78, 0.21, 0.35};
    double *X = arr;

    */

    /*
    //-------print X as mapped in memory------
    printf("N = %d\n", N);
    //check X array format
    for(int i = 0 ; i<D*N ; i++)
    {
        if(i%(N)==0 && i!=0 )
        {
            printf("-----\n");
        }
        printf("%f\n", X[i]);
    }
    printf("\n");
    //---------------
    */

    /*
    //------print X as a n-by-d array-------
    for(int i = 0 ; i<N ; i++)
    {
        for(int j=0; j<D; j++)
        {

            printf("%f ",X[i+j*N] );
        }
        printf("\n");
    }
    printf("\n");
    //------------------
    */

    gettimeofday (&startwtime, NULL);
    st_time = times(&st_cpu);

    T =buildvp(X, N, D);


    en_time = times(&en_cpu);
    gettimeofday (&endwtime, NULL);
    p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
  		      + endwtime.tv_sec - startwtime.tv_sec);

    printf(" time = %fsec\n",p_time);


    destroy(T);
    free(X);


    return 0;
}
