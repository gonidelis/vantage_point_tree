/* Project no.1 on Parallel and Destributed Systems @AUTh - ECE Dept.
*  GONIDELIS GIANNIS
*  DATE: 03/11/2019
*  ----------------------------------
* Version no.3 : parallel implementation utilizing openmp interface
*/

#include <stdio.h>
#include <stdlib.h>
//#include "quickselect.h"
#include "../inc/vptree.h"
#include <math.h>
#include <omp.h>


#define SWAP(x, y) { double temp = x; x = y; y = temp; }
#define N (sizeof(A)/sizeof(A[0]))

// Partition using Lomuto partition scheme
double partition(double a[], int left, int right, int pivotIndex)
{
	// Pick pivotIndex as pivot from the array
	double pivot = a[pivotIndex];

	// Move pivot to end
	SWAP(a[pivotIndex], a[right]);

	// elements less than pivot will be pushed to the left of pIndex
	// elements more than pivot will be pushed to the right of pIndex
	// equal elements can go either way
	int pIndex = left;
	int i;

	// each time we finds an element less than or equal to pivot, pIndex
	// is incremented and that element would be placed before the pivot.
	for (i = left; i < right; i++)
	{
		if (a[i] <= pivot)
		{
			SWAP(a[i], a[pIndex]);
			pIndex++;
		}
	}

	// Move pivot to its final place
	SWAP(a[pIndex], a[right]);

	// return pIndex (index of pivot element)
	return pIndex;
}

// Returns the k-th smallest element of list within left..right
// (i.e. left <= k <= right). The search space within the array is
// changing for each round - but the list is still the same size.
// Thus, k does not need to be updated with each round.
double quickselect(double A[], int left, int right, int k)
{
	// If the array contains only one element, return that element
	if (left == right)
		return A[left];

	// select a pivotIndex between left and right
	int pivotIndex = (left + right)/2;

	pivotIndex = partition(A, left, right, pivotIndex);

	// The pivot is in its final sorted position
	if (k == pivotIndex)
		return A[k];

	// if k is less than the pivot index
	else if (k < pivotIndex)
		return quickselect(A, left, pivotIndex - 1, k);

	// if k is more than the pivot index
	else
		return quickselect(A, pivotIndex + 1, right, k);
}


double * calcDistance(double * X ,int n, int d)
{
    //calculate and return an array[n-1] of all the distances
    //from the last point
    double *distances = calloc(n,sizeof(double));
    omp_set_nested(1);
	//in parallel for n>10.000 else -> sequentially
    #pragma omp parallel for if(n>25000)
    for(int i=0 ; i<n-1; i++)
    {
        for (int j=0; j< d; j++)
        {
            double dist = X[(j+1)*n-1]-X[j*n+i];
            distances[i] += dist*dist;

        }
        distances[i] = sqrt(distances[i]);
    }
    return distances;
}

vptree * buildvp_recur(double *X, int n, int d, int *xID)
{

    //initializations
    double * dist ;
    double *distances=malloc(n*sizeof(double));
    vptree *T = NULL;
    int innerSize= 0;
    int outerSize = 0;
    //allocate for zero elements (no allocation) and realloc
    //for the additional ones, as innerSize is not predefined
    double *innerX = malloc(innerSize*d * sizeof(double));
    int *innerID = malloc(innerSize * sizeof(int));
    double *outerX = malloc(outerSize*d * sizeof(double));
    int *outerID = malloc(outerSize * sizeof(int));

    if(n==0)
    {
        free(innerX);
        free(outerX);
        free(innerID);
        free(outerID);
        return T;
    }

    T=calloc(1, sizeof (vptree));

    //define vantage-point index (idx) and coefficients (vp)
    T->vp = malloc(d * sizeof(double));

    //id = n-1 as C is zero based
    T->idx = xID[n-1];
    for (int i=0 ; i < d ; i++)
    {
        T->vp[i] = X[n-1+i*n];
    }

    if(n==1)
    {
        free(innerX);
        free(outerX);
        free(innerID);
        free(outerID);
        return T;
    }



    //define distance from each other point
    dist= calcDistance(X , n , d);
	//keep distances[]-array twice as quickselect alters its contents
    for(int i = 0; i<n;i++)
    {
        distances[i]=dist[i];
    }

    //find median in case n -> odd
    if ((n-1)%2!=0)
    {
        T->md = quickselect(dist, 0, n-1, n/2  );
    }
    //find median in case n -> even
    else
    {
        T->md = quickselect(dist,0, n-1 , n/2 + 1) + quickselect(dist,0,  n-1, n/2 );
        T->md /= 2.0;
    }


    /*Define inner - outer
    *Split X into innerX and outerX:
    *innerX is used as an array of pointers at first
    *as it only keeps the X positions (i) of its future
    *elements just for the 1st dimension. As soon as
    *innerSize has been determined the pointers resolve
    *to the real values of innerX. Same goes for outerX.
    */
    for(int i=0 ; i<n-1 ; i++)
    {
        if(distances[i]<= T->md)
        {
            innerSize ++;
            //for each element allocate d blocks
            innerX = realloc(innerX , innerSize*d * sizeof(double));
            innerID = realloc(innerID , innerSize * sizeof(int));

            innerX[innerSize-1] = i;
            //xID maintains the original elements' id's
            innerID[innerSize-1] = xID[i];

        }
        else
        {
            outerSize++;
            outerX = realloc(outerX, outerSize*d * sizeof(double));
            outerID = realloc(outerID , outerSize * sizeof(int));

            outerX[outerSize-1] = i;
            outerID[outerSize-1] = xID[i];
        }
    }

    //resolve innerX[] to real coefficients
    for(int i=0 ; i<innerSize; i++)
    {
        for(int j=d-1; j>=0; j--)
        {
            /*
            fill the coefficients backwards in order for the
            *pointer information not to be lost,
            *as the first dimenstion holds the position information
            */
            innerX[i+j*innerSize]=X[(int)innerX[i]+j*n];
        }
    }
    //resolve outerX[] to real coefficients
    for(int i=0 ; i<outerSize; i++)
    {
        for(int j=d-1; j>=0; j--)
        {
            //resolve outerX[i] to real coefficients
            outerX[i+j*outerSize]=X[(int)outerX[i]+j*n];
        }
    }

    //call buildvp() recursively - Pre-order tree traversal
    int activeness = omp_get_active_level();
	//define (low) level of parallelism due to hardware limitations
    #pragma omp parallel sections if(activeness<2)
    {
        #pragma omp section
        T->inner = buildvp_recur(innerX, innerSize, d,innerID);
        #pragma omp section
        T->outer = buildvp_recur(outerX, outerSize, d, outerID);
    }


    free(innerX);
    free(outerX);
    free(innerID);
    free(outerID);
    free(distances);
    free(dist);

    return T;

}

vptree * buildvp(double *X, int n, int d)
{
    //keep pointer-array of indexes
    //update array on each row swap/change
    int *Xid = malloc(n*sizeof(int));
    vptree *T;
    for(int i=0;i<n;i++)
    {
        Xid[i] = i;
    }

    T=buildvp_recur(X, n, d, Xid);

    free(Xid);
    return T;
}


vptree * getInner(vptree * T)
{
    return T->inner;
}

vptree * getOuter(vptree * T)
{
    return T->outer;
}

double getMD(vptree * T)
{
    return T->md;
}

double *getVP(vptree *T)
{
    return T->vp;
}

int getIDX(vptree *T)
{
    return T->idx;
}

void destroy(vptree *T)
{
    if(T==NULL)
    {
        return;
    }

    destroy(T->inner);
    destroy(T->outer);

    free(T->vp);
    free(T);
}
