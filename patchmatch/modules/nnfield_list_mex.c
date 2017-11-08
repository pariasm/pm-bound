#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "NNfield_l2.h"
#include "dImage.h"

/* funciones de verificacion de parametros */
int mxIsRealMatrix(const mxArray *data,int dim) {
	return !(mxIsChar               (data) ||
	         mxIsSparse             (data) ||
	         mxIsComplex            (data) ||
	         mxGetNumberOfDimensions(data) > dim);
}

double mxExtractStructScalarField(const mxArray *STR,
		int str_idx, char *field, double def_val) {
	if (mxGetField(STR,str_idx,field) == NULL) return def_val;
	else if (!mxIsDouble(mxGetField(STR,str_idx,field)) ||
		(mxGetM(mxGetField(STR,str_idx,field)) != 1) ||
		(mxGetN(mxGetField(STR,str_idx,field)) != 1))
	{
		char errmsg[100];
		sprintf(errmsg,"Field %s is not a double scalar.",field);
		mexErrMsgTxt(errmsg);
	}
	else
		return (*mxGetPr(mxGetField(STR,str_idx,field)));
}

#define IM1_IN    prhs[0]  /* image 1                             */
#define MASK1_IN  prhs[1]  /* mask 1 (1 = ROI)                    */
#define IM2_IN    prhs[2]  /* image 2                             */
#define MASK2_IN  prhs[3]  /* mask 2 (1 = ROI)                    */
#define WINDOW_IN prhs[4]  /* weighting function for patch metric */
#define PARAMS_IN prhs[5]  /* parameters: mxit, psz               */
#define NNFx_IN   prhs[6]
#define NNFy_IN   prhs[7]
#define NNFx_OUT  plhs[0]  /* OUTPUT : nnf                        */
#define NNFy_OUT  plhs[1]  /* OUTPUT : nnf                        */
#define SCRS_OUT  plhs[2]  /* OUTPUT : nnf scores                 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* get input data ------------------------------------------------------------------  */

	/* sizes */
	int sz1[2] = {mxGetM(MASK1_IN), mxGetN(MASK1_IN) };
	int sz2[2] = {mxGetM(MASK2_IN), mxGetN(MASK2_IN) };

	int ndims = mxGetNumberOfDimensions(IM1_IN);
	int *dims = (int *)mxGetDimensions (IM1_IN);
	int nch = (ndims == 3) ? dims[2] : 1;
	double **im_channels = (double **)mxCalloc(nch,sizeof(double *));

	int ii;
	for (ii = 0; ii < nch; ii++) im_channels[ii] = mxGetPr(IM1_IN) + sz1[0]*sz1[1]*ii;
	dImage im1 = copy_image_components_to_dimage_with_frame(sz1[0],sz1[1],nch,0,0,0,im_channels);

	for (ii = 0; ii < nch; ii++) im_channels[ii] = mxGetPr(IM2_IN) + sz2[0]*sz2[1]*ii;
	dImage im2 = copy_image_components_to_dimage_with_frame(sz2[0],sz2[1],nch,0,0,0,im_channels);

	double *msk1 = mxGetPr(MASK1_IN);
	double *msk2 = mxGetPr(MASK2_IN);

	/* length of candidate lists */
   int list_length = (int)mxExtractStructScalarField(PARAMS_IN,0,"list",5); 

	/* number of iterations for the NNF algorithm */
   int mxit = (int)mxExtractStructScalarField(PARAMS_IN,0,"mxit",5);

	/* patch size */
	int tmp[2];
	tmp[0] = (int)(mxGetPr(mxGetField(PARAMS_IN,0,"psz"))[0]);
	tmp[1] = (int)(mxGetPr(mxGetField(PARAMS_IN,0,"psz"))[1]);
	int psz = tmp[0];

	/* weighting function for patch distance */
	double *window = mxGetPr(WINDOW_IN);





	/* parse binary input masks and create indexed masks for NNF ------------------------ */
	int *idx1 = (int *)mxCalloc(sz1[0]*sz1[1],sizeof(int));
	int l1 = 0;
	for (ii = 0; ii < sz1[0]*sz1[1]; ii++) if (msk1[ii]) idx1[l1++] = ii;
	/*
	 * idx1 and idx2 impose an order in the pixels in regions D_1 and
	 * D_2 ~ see README file.
	 *
	 */

	int *idx2 = (int *)mxCalloc(sz2[0]*sz2[1],sizeof(int));
	int l2 = 0;
	for (ii = 0; ii < sz2[0]*sz2[1]; ii++) if (msk2[ii]) idx2[l2++] = ii;


	/*
	 * pos2idx1 and pos2idx2 are the 'indexed masks'
	 *   - a negative value indicates that the position does not belong to
	 *     the region
	 *   - a positive value specifies the index in the region ordering
	 *
	 * See README file.
	 *
	 */
	int *pos2idx1 = (int *)mxCalloc(sz1[0]*sz1[1],sizeof(int));
	int *pos2idx2 = (int *)mxCalloc(sz2[0]*sz2[1],sizeof(int));

	for (ii = 0; ii < sz1[0]*sz1[1] ; ii++) pos2idx1 [ii] = -1;
	for (ii = 0; ii < sz2[0]*sz2[1] ; ii++) pos2idx2 [ii] = -1;

	for (ii = 0; ii < l1 ; ii++) pos2idx1 [idx1 [ii]] = ii;
	for (ii = 0; ii < l2 ; ii++) pos2idx2 [idx2 [ii]] = ii;

	/* compute nnf ---------------------------------------------------------------------- */

	/* create nnf */
	list_max_global = list_length;  
	list nnf = NNfield_create_nnflists(l1);
	for (ii = 0; ii < l1; ii++) {
		nnf[ii].quant = 1;
		nnf[ii].l[0].ox  = INT_MAX;
		nnf[ii].l[0].oy  = 0;
		nnf[ii].l[0].pos = -1; 
	}

	/* if an nnf is provided, copy it in nnf */
	if (nrhs >= 8) {
		int ndims = mxGetNumberOfDimensions(NNFx_IN);
		int *dims = (int *)mxGetDimensions (NNFx_IN);
		int quant = (ndims == 3) ? dims[2] : 1;

		double *nnfx_in = mxGetPr(NNFx_IN);
		double *nnfy_in = mxGetPr(NNFy_IN);
		for (ii = 0; ii < sz1[0]*sz1[1]; ++ii) {
			int idx = pos2idx1[ii];
			if (idx >= 0) {
				int jj;
				nnf[idx].quant = quant;
				for (jj = 0; jj < nnf[idx].quant; jj++) {
					nnf[idx].l[jj].ox    = nnfx_in[ ii + jj*sz1[0]*sz1[1] ];
					nnf[idx].l[jj].oy    = nnfy_in[ ii + jj*sz1[0]*sz1[1] ];
/*					nnf[idx].l[jj].score = scrs_in[ ii + jj*sz1[0]*sz1[1] ]; */
				}
			}
		}
	}

	NNfield(pos2idx1,im1,pos2idx2,im2,window,psz,mxit,nnf);




	/* store output as a MATLAB array --------------------------------------------------- */
	int nnfo_ndims = 3;
	int nnfo_dims[3];
	nnfo_dims[0] = sz1[0];
	nnfo_dims[1] = sz1[1];
	nnfo_dims[2] = list_length;
	NNFx_OUT = mxCreateNumericArray(nnfo_ndims,nnfo_dims,mxDOUBLE_CLASS,mxREAL);
	NNFy_OUT = mxCreateNumericArray(nnfo_ndims,nnfo_dims,mxDOUBLE_CLASS,mxREAL);
	SCRS_OUT = mxCreateNumericArray(nnfo_ndims,nnfo_dims,mxDOUBLE_CLASS,mxREAL);
	double *nnfx_out = mxGetPr(NNFx_OUT);
	double *nnfy_out = mxGetPr(NNFy_OUT);
	double *scrs_out = mxGetPr(SCRS_OUT);

	for (ii = 0; ii < l1; ii++) {
		int jj;
		for (jj = 0; jj < nnf[ii].quant; jj++) {
			nnfx_out[ idx1[ii] + jj*sz1[0]*sz1[1] ] = nnf[ii].l[jj].ox;
			nnfy_out[ idx1[ii] + jj*sz1[0]*sz1[1] ] = nnf[ii].l[jj].oy;
			scrs_out[ idx1[ii] + jj*sz1[0]*sz1[1] ] = nnf[ii].l[jj].score;
		}
	}



	/* free! */
	mxFree(pos2idx1); 
	mxFree(pos2idx2); 
	mxFree(idx1); 
	mxFree(idx2); 
	del_dimage(im1);
	del_dimage(im2);

	NNfield_destroy_nnflists(nnf, l1);

}
