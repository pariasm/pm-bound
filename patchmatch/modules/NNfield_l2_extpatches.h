#include "NNfield.h"

/*
 *
 * This is an example of an NNF module.
 *
 * In this case, patches are previously extracted and stored
 * in a double array of psz*psz*|D_1|*channels, where psz is 
 * the side of the patch and |D_1| is the number of pixels 
 * in D_1 and channels is the number of image channels.
 *
 * Note that this setting can be used to compute matchings
 * based on any set of pixel-wise features (color, Gabor 
 * coefficients, etc).
 *
 */
typedef struct {
   double *features1;  /* patches in image1 */
   double *features2;  /* patches in image2 */
   int     feat_dim;   /* dimension of a patch (psz*psz)*channels */
   int     channels;   /* number of channels of image */
} dataStruct;


/*
 * User defined patch distance function.
 */
DOUBLE l2_distance (int idx1,int idx2, void *data) 
{

   int jj;
   dataStruct* d = (dataStruct*) data;

   double *p_patch1 = d->features1 + d->feat_dim*idx1;
   double *p_patch2 = d->features2 + d->feat_dim*idx2;

   DOUBLE dist = 0, tmp;
   int Ccount = d->channels;
   for (jj = 0; jj < d->feat_dim; jj ++) 
      dist += (tmp = *p_patch1++ - *p_patch2++)*tmp;

   return dist ;

}


/*
 * Your NNF interface.
 */
extern NNfield ( int *mask1, int ncol1, int nrow1, double* features1_in, 
	              int *mask2, int ncol2, int nrow2, double* features2_in, 
 					  int patch_dim_in, int channels, int maxit, list nnfinout)
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

   dataStruct NNF_data ;
   NNF_data.features1 = features1_in;
   NNF_data.features2 = features2_in;
   NNF_data.feat_dim  = patch_dim_in*channels;
   NNF_data.channels = channels;

#ifndef EXHNNF
   NNfield_basic(mask1,ncol1,nrow1,mask2,ncol2,nrow2,
                 maxit,nnfinout,&NNF_data,l2_distance);
#else 
   NNfield_exhaustive(mask1,ncol1,nrow1,mask2,ncol2,nrow2,
        		          maxit,nnfinout,&NNF_data,l2_distance);
#endif

}
