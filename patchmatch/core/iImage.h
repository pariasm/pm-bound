#ifndef IIMAGE_FACCIOLO_H
#define IIMAGE_FACCIOLO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* declarations */
typedef struct iimage {
	int nrow;
	int ncol;
	int channels;
	int *val;
} *iImage;

static iImage new_iimage(void);


static iImage new_iimage2(int nx, int ny, int channels);

static void del_iimage(iImage i) ;

static iImage copy_iimage(iImage dest, iImage src);


/* definitions */
iImage new_iimage(void)
{
	iImage im;

	if ( !(im = (iImage) (malloc(sizeof(struct iimage))) ) )
	{
		printf("[new_iimage] Not enough memory\n");
		exit(1);
		return(NULL);
	}

	im->nrow = im->ncol = 0;
	im->channels = 0;
	im->val = NULL;
	return(im);
}


iImage new_iimage2(int nx, int ny, int channels){
	iImage im = new_iimage();
	im->ncol = nx;
	im->nrow = ny;
	im->channels = channels;
	im->val = (int*) calloc(im->ncol*im->nrow*im->channels,sizeof(int));
	return (im);
}


void del_iimage(iImage im) {
	if (im->val) free(im->val);
	free(im);
}

iImage copy_iimage(iImage dest, iImage src){
	if ((dest->ncol == src->ncol) && (dest->nrow == src->nrow) && (dest->channels == src->channels)){
		memcpy(dest->val,src->val,src->channels*src->ncol*src->nrow*sizeof(int)); 
	} else {
		printf("[copy_iimage] Image sizes does not match\n");
	}
	return dest;
}

iImage copy_data_to_iimage_with_frame ( int ncol, int nrow, int channels, int framesize_col, int framesize_row, int fill_value, int* data )  {
	int i,j,k;

	iImage im = new_iimage2(ncol + 2*framesize_col , nrow + 2*framesize_row, channels);
	for ( i =0; i<im->ncol * im->nrow * channels; i++ ) { im->val[i] = fill_value; }

	for ( i =0; i<ncol; i++ ) {
		for ( j =0; j<nrow; j++ ) {
			for (k=0;k<channels;k++) 
				im->val[(i+framesize_col + im->ncol*(j+framesize_row))*channels + k ] =  data [(i + ncol*j)*channels + k];
		}
	}
	return im;
}
void copy_iimage_with_frame_to_data ( iImage im , int framesize_col, int framesize_row, int* data )  {
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
#endif

