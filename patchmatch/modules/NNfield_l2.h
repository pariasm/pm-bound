#include "NNfield.h"
#include "dImage.h"

/*
 *
 * This is an example of an NNF module.
 *
 */

/* data needed for the computation of the patch distance */
typedef struct {
	dImage im1;         /* images */
	dImage im2;
	double *window;     /* weight function for patch metric */
	int psz;            /* size of the patch */
	int *idx1toPos;     /* compute image coordinates of a point given */
	int *idx2toPos;     /* its index in the order of D_1 (see README file) */
} dataStruct;

int isinside(int ncol, int nrow, int x, int y) {
	return (x>=0 && y>=0 && x<ncol && y<nrow);
}

/* l2 distance function, only works for grayscale images */
DOUBLE l2_distance (int idx1,int idx2, void *data)
{
	dataStruct* d = (dataStruct*) data;

	/* coordinates of the patch centers c1 and c2 */
	int c1x = d->idx1toPos[2*idx1 + 0];
	int c1y = d->idx1toPos[2*idx1 + 1];
	int c2x = d->idx2toPos[2*idx2 + 0];
	int c2y = d->idx2toPos[2*idx2 + 1];

	int psz  = d->psz;
	int hpsz = psz /2;

	dImage im1 = d->im1;
	dImage im2 = d->im2;
	double *p_window = d->window;

	int ncol1 = im1->ncol;
	int nrow1 = im1->nrow;
	int ncol2 = im2->ncol;
	int nrow2 = im2->nrow;

	int kx,ky;
	DOUBLE dista = 0;
	for(ky = -hpsz; ky <= hpsz; ky++)
	for(kx = -hpsz; kx <= hpsz; kx++)
	{
		int p1x = c1x + kx;
		int p1y = c1y + ky;
		int p2x = c2x + kx;
		int p2y = c2y + ky;

		double m1 = isinside(ncol1, nrow1, p1x, p1y);
		double m2 = isinside(ncol2, nrow2, p2x, p2y);

		int idxp1 = p1x + ncol1*p1y;
		int idxp2 = p2x + ncol2*p2y;

		double tmp;
		dista += (m1 && m2) ? (tmp = im1->val[idxp1] - im2->val[idxp2] )*(tmp)* (*p_window++) : 1e10;
	}

	return dista;
}

extern NNfield( int *mask1, dImage im1, int *mask2, dImage im2,
                double *window, int psz, int maxit, list nnfinout)
{
	/*
	 * mask1 and mask2 encode two things at the same time.
	 *   - the source and target domains for the NNF
	 *   - the ordering in which pixels in those domains are traversed
	 *
	 * If mask1 at (i,j) is negative, (i,j) does NOT belong to D_1.
	 * If mask1 at (i,j) is positive, then it does belong, and
	 * mask1(i,j) indicates the index of (i,j) in the ordering of D_1.
	 *
	 */
	int i,j;

	dataStruct NNF_data;
	NNF_data.im1 = im1;
	NNF_data.im2 = im2;
	NNF_data.window = window;
	NNF_data.psz = psz;

	int ncol1 = im1->ncol;
	int nrow1 = im1->nrow;
	int ncol2 = im2->ncol;
	int nrow2 = im2->nrow;

	NNF_data.idx1toPos = mxCalloc( 2*ncol1*nrow1, sizeof(int));
	for (i=0;i<ncol1;i++) 
	for (j=0;j<nrow1;j++) {
		int idx = mask1[i+j*ncol1]; 
		if( idx>=0 ) {
			NNF_data.idx1toPos[2*idx+0] = i;
			NNF_data.idx1toPos[2*idx+1] = j;
		}
	}

	NNF_data.idx2toPos = mxCalloc( 2*ncol2*nrow2, sizeof(int));
	for (i=0;i<ncol2;i++)
	for (j=0;j<nrow2;j++) {
		int idx = mask2[i+j*ncol2]; 
		if( idx>=0 ) {
			NNF_data.idx2toPos[2*idx+0] = i;
			NNF_data.idx2toPos[2*idx+1] = j;
		}
	}

#ifndef EXHNNF
	/* NNfield_decoupled(mask1,ncol1, nrow1, mask2,ncol2,nrow2,
		maxit ,nnfinout, &NNF_data, l2_distance);*/
	NNfield_basic(mask1, ncol1, nrow1, mask2, ncol2, nrow2,
		maxit ,nnfinout, &NNF_data, l2_distance);
#else
	/* this is included just for comparison purposes */
	printf(" --> exhaustive computation of nearest neighbor <-- ");
	NNfield_exhaustive(mask1, ncol1, nrow1, mask2, ncol2, nrow2, 
		maxit ,nnfinout, &NNF_data, l2_distance);
#endif

	mxFree(NNF_data.idx2toPos);
	mxFree(NNF_data.idx1toPos);
}
