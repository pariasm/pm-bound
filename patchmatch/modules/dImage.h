#ifndef DIMAGE_FACCIOLO_H
#define DIMAGE_FACCIOLO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* declarations */
typedef struct dimage {
	int nrow;
	int ncol;
	int channels;
	double *val;
} *dImage;

static dImage new_dimage(void);


static dImage new_dimage2(int nx, int ny, int channels);

static void del_dimage(dImage i) ;

static dImage copy_dimage(dImage dest, dImage src);


/* definitions */
dImage new_dimage(void)
{
	dImage im;

	if ( !(im = (dImage) (mxMalloc(sizeof(struct dimage))) ) ) {
		mexPrintf("[new_dimage] Not enough memory\n");
		exit(1);
		return(NULL);
	}

	im->nrow = im->ncol = 0;
	im->channels = 0;
	im->val = NULL;
	return(im);
}


dImage new_dimage2(int nx, int ny, int channels)
{
	dImage im = new_dimage();
	im->ncol = nx;
	im->nrow = ny;
	im->channels = channels;
	im->val = (double*) mxCalloc(im->ncol*im->nrow*im->channels,sizeof(double));
	return (im);
}


void del_dimage(dImage im) 
{
	if (im->val) mxFree(im->val);
	mxFree(im);
}

dImage copy_dimage(dImage dest, dImage src)
{
	if ((dest->ncol == src->ncol) && (dest->nrow == src->nrow) && (dest->channels == src->channels)){
		memcpy(dest->val,src->val,src->channels*src->ncol*src->nrow*sizeof(double)); 
	} else {
		mexPrintf("[copy_dimage] dImage sizes does not match\n");
	}
	return dest;
}

dImage copy_data_to_dimage_with_frame ( int ncol, int nrow, int channels, int framesize_col, int framesize_row, int fill_value, double* data )  
{
	int i,j,k;

	dImage im = new_dimage2(ncol + 2*framesize_col , nrow + 2*framesize_row, channels);
	for ( i =0; i<im->ncol * im->nrow * channels; i++ ) { im->val[i] = fill_value; }

	for ( i =0; i<ncol; i++ ) {
		for ( j =0; j<nrow; j++ ) {
			for (k=0;k<channels;k++) 
				im->val[(i+framesize_col + im->ncol*(j+framesize_row))*channels + k ] =  data [(i + ncol*j)*channels + k];
		}
	}
	return im;
}

void copy_dimage_with_frame_to_data ( dImage im , int framesize_col, int framesize_row, double* data )  
{
	int i,j,k;
	int ncol = im->ncol - framesize_col*2;
	int nrow = im->nrow - framesize_row*2;
	int channels = im->channels;

	for ( i =0; i<ncol; i++ ) {
		for ( j =0; j<nrow; j++ ) {
			for (k=0;k<channels;k++) 
				data [(i + ncol*j)*channels + k] = im->val[(i+framesize_col + im->ncol*(j+framesize_row))*channels + k ];
		}
	}
}

dImage copy_image_components_to_dimage_with_frame (int ncol, int nrow, int channels, 
		int framesize_col, int framesize_row, int fill_value, 
		double** pdata)  
{
	int i,j,k;

	dImage im = new_dimage2(ncol + 2*framesize_col, nrow + 2*framesize_row, channels);

	for ( i =0; i < im->ncol * im->nrow * channels; i++ ) { im->val[i] = fill_value; }

	for ( i = 0; i < ncol    ; i++ )
	for ( j = 0; j < nrow    ; j++ )
	for ( k = 0; k < channels; k++ ) 
	{
				double tmp = pdata[k][i + ncol*j];
				im->val[ (i+framesize_col + im->ncol*(j+framesize_row))*channels + k ] = tmp;
	}

	return im;
}

dImage copy_image_components_to_dimage_with_frame_inplace (int ncol, int nrow, int channels, 
		int framesize_col, int framesize_row, int fill_value, 
		double** pdata, dImage im)  
{
	int i,j,k;

	if(im->ncol != ncol || im->nrow != nrow || im->channels != channels ) {
		mexPrintf("ERROR!!!! : copy_image_components_to_dimage_with_frame_inplace of different sizes\n");
		return 0;
	}
	for( i =0; i<im->ncol * im->nrow * channels; i++ ) { im->val[i] = fill_value; }

	for ( i =0; i<ncol; i++ ) {
		for ( j =0; j<nrow; j++ ) {
			for (k=0;k<channels;k++) 
				im->val[(i+framesize_col + im->ncol*(j+framesize_row))*channels + k ] = pdata[k][i + ncol*j];
		}
	}
	return im;
}

void copy_dimage_with_frame_into_components(dImage im, int framesize_col, int framesize_row, double** pdata)  
{
	int i,j,k;
	int ncol = im->ncol - framesize_col*2;
	int nrow = im->nrow - framesize_row*2;
	int channels = im->channels;

	for ( i =0; i<ncol; i++ ) {
		for ( j =0; j<nrow; j++ ) {
			for (k=0;k<channels;k++) 
				pdata[k][(i + ncol*j)] = im->val[(i+framesize_col + im->ncol*(j+framesize_row))*channels + k ];
		}
	}
}

void im_gradient_dImage(dImage im,char *num_scheme,dImage gim)
{
	int p,ncol = im->ncol,nrow = im->nrow;
	int nchl = im->channels,k;

	for(p = 0; p < ncol*nrow; p++) {

		int pc = p % ncol;
		int pr = p / ncol;

		/* derivative in the direction of the rows (column is fixed) */
		switch (num_scheme[0]) {

			case 'f':
				if (pr < nrow - 1) 
					for (k = 0; k < nchl; k++)
						gim->val[2*nchl*p + k] = im->val[nchl*(p + ncol) + k] - im->val[nchl*p + k];
				else for (k = 0; k < nchl; k++) gim->val[2*nchl*p + k] = 0;
				break;

			case 'b':
				if (pr > 0)
					for (k = 0; k < nchl; k++)
						gim->val[2*nchl*p + k] = im->val[nchl*p + k] - im->val[nchl*(p - ncol) + k];
				else for (k = 0; k < nchl; k++) gim->val[2*nchl*p + k] = 0;
				break;

			default:
				mexErrMsgTxt("Numeric scheme must be 'f' of 'b'");
		}

		/* derivative in the direction of the columns (row is fixed) */
		switch (num_scheme[1]) {

			case 'f':
				if (pc < ncol - 1) 
					for (k = 0; k < nchl; k++)
						gim->val[2*nchl*p + nchl + k] = im->val[nchl*(p + 1) + k] - im->val[nchl*p + k];
				else for (k = 0; k < nchl; k++) gim->val[2*nchl*p + nchl + k] = 0;
				break;

			case 'b':
				if (pc > 0)
					for (k = 0; k < nchl; k++)
						gim->val[2*nchl*p + nchl + k] = im->val[nchl*p + k] - im->val[nchl*(p - 1) + k];
				else for (k = 0; k < nchl; k++) gim->val[2*nchl*p + nchl + k] = 0;
				break;

			default:
				mexErrMsgTxt("Numeric scheme must be 'f' of 'b'");
		}
	}
}

/* |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| */
/* This function computes the divergence of a vector field in a region of the image.
	The divergence can be computed according to four numerical schemes (ff,bb,fb,bf) where f
	stands for forward differences, b for backwards. Therefore fb means that the differences
	in the derivative aint  x is computed with forward differences, whereas the derivative 
	aint  direction y with backward differences. The x direction increases with the column, 
	and the y with the rows.

	The boundary of the regions has two parts: the one due to the limits of the image domain
	(image boundary) and the one due to the region (region boundary). For both, the boudary 
	conditions are zero derivatives.                                                       

NOTE: The first component increases with the columns.                                 

NOTE: The difference between im_region_divergence_matrix and im_region_divergence is 
that im_region_divergence returns a list div of length rl, such that div[i] is the 
divergence in pixels p = ridx[i]. im_region_divergence_matrix returns a matrix 
div such that div[p] is the divergence at p if p is in the region or zero if not.*/
void im_region_divergence_dImage(
		dImage  vim,                                 /* field (first component)                */
		int     *ridx,                               /* indices of positions in the region     */
		int     rl,                                  /* length of indices vector               */
		double  *msk,                                /* mask with rlabel in the region         */
		double  rlabel,                              /* label of the region                    */
		char    *num_scheme,                         /* numeric scheme                         */
		dImage  div)                                 /* divergence [size sz]                   */
{
	int i,p;                                              
	int isbound,j;
	int ncol = vim->ncol,nrow = vim->nrow;
	int nchl = vim->channels,k;

	/* the 4-neighborhood: p + neigh[i] is the position of neighbor i */
	int neigh[4] = {-1,ncol,1,-ncol};

	for(i = 0;i < rl;i++ ) {

		p = ridx[i];                               /* i-th pixel of the region */

		/* check if p is in the image boundary     */
		if((p % ncol > 0) && (p % ncol < ncol - 1) &&
				(p / ncol > 0) && (p / ncol < nrow - 1)) {

			/* check if p is in the region boundary */
			isbound = 0; 
			for (j = 0;j < 4;j++)
				isbound = isbound + (msk[p + neigh[j]] != rlabel);

			if (isbound == 0) {
				switch (num_scheme[0])  {     /* derivative across columns (d/dx) */ 
					case 'f':
						for (k = 0; k < nchl; k++)
							div->val[nchl*p + k] = vim->val[2*nchl*(p + ncol) + k] - vim->val[2*nchl*p + k];
						break;
					case 'b':
						for (k = 0; k < nchl; k++)
							div->val[nchl*p + k] = vim->val[2*nchl*p + k] - vim->val[2*nchl*(p - ncol) + k];
						break;
					default:
						mexErrMsgTxt("Numeric scheme must be 'f' of 'b'");
				}

				switch (num_scheme[1]) {        /* derivative across rows (d/dy) */
					case 'f':
						for (k = 0; k < nchl; k++)
							div->val[nchl*p + k] += vim->val[2*nchl*(p + 1) + nchl + k] - vim->val[2*nchl*p + nchl + k];
						break;
					case 'b':
						for (k = 0; k < nchl; k++)
							div->val[nchl*p + k] += vim->val[2*nchl*p + nchl + k] - vim->val[2*nchl*(p - 1) + nchl + k];
						break;
					default:
						mexErrMsgTxt("Numeric scheme must be 'f' of 'b'");
				}
			}
			else
				for (k = 0; k < nchl; k++)
					div->val[nchl*p + k] = 0;
		}
		else
			for (k = 0; k < nchl; k++)
				div->val[nchl*p + k] = 0;
		/* if-else image boundary */
	}                           /* for through the image  */
}                              /* function               */

/* This version of extract_patches_mirror_inplace reveices a pointer to previously 
 * allocated memory location where the patches are copied.                                */
void extract_patches_mirror_dImage_inplace(
		dImage im,
		int    *idx,
		int    numel,
		int    *psz,
		double *patches)
{
	int  sz[2],hpsz[2],pdim;
	int  col,row,i,k,i1,i2,nchannels;
	int  pos[2],idx_tmp[2];

	sz[0] = im->ncol;
	sz[1] = im->nrow;
	nchannels = im->channels;
	pdim = psz[0]*psz[1]*nchannels;
	hpsz[0] = psz[0]/2;
	hpsz[1] = psz[1]/2;

	for (i=0; i < numel; i++)
	{
		idx_tmp[0] = idx[i]%sz[0];                    /* indices to images coordinates */
		idx_tmp[1] = idx[i]/sz[0];

		/* If the center is NOT in the domain of the image, the image has to be mirrorred */
		/* This case is sepparated from the normal case because the normal case can be
		 * implemented faster.                                                            */
		if (idx_tmp[0] < hpsz[0] || idx_tmp[1] < hpsz[1] || idx_tmp[0] >= sz[0] - hpsz[0] ||
				idx_tmp[1] >= sz[1] - hpsz[1])
		{
			pos[0] = idx_tmp[0] - hpsz[0];
			pos[1] = idx_tmp[1] - hpsz[1];

			/* Go through the patch */
			for (col=0; col < psz[1]; col++)
				for (row=0; row < psz[0]; row++)
				{
					((pos[0] + row) < 0)?(i1 = - pos[0]-row):(i1 = pos[0]+row);
					((pos[1] + col) < 0)?(i2 = - pos[1]-col):(i2 = pos[1]+col);

					(i1 > (sz[0]-1))?(i1 = 2*(sz[0]-1)-pos[0]-row):(i1 = i1);
					(i2 > (sz[1]-1))?(i2 = 2*(sz[1]-1)-pos[1]-col):(i2 = i2);

					for (k = 0; k < nchannels; k++)
						patches[i*pdim + (psz[0]*col + row)*nchannels + k] = im->val[nchannels*(i2*sz[0] + i1) + k];
				}
		}
		else                                  /* plain and simple extraction of the patch */
		{
			pos[0] = idx_tmp[0] - hpsz[0];
			pos[1] = idx_tmp[1] - hpsz[1];

			/* Go through the patch */
			for (col=0; col < psz[1]; col++)
				memcpy(patches + psz[0]*nchannels*col + i*pdim,
						&(im->val[nchannels*(pos[0] + sz[0]*(pos[1] + col))]),
						psz[0]*nchannels*sizeof(double));
		}
	}
}

#endif

