#include "iImage.h" 
#include "NNfield.h" 
//#include <omp.h>

/*#define DEBUG_STATS*/

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

/* LOCALLY GLOBAL VARIABLES */
int NUM_NEIGH = 2;                    /* neighbors for the forward/backward propagation */
int NEIGH_FWD[2*2] = {-1, 0, 0, -1};
int NEIGH_BWD[2*2] = { 1, 0, 0,  1};
int SUCCESFUL_RANDOM_SEARCH; /* these are just for statistics */
int FAILED_RANDOM_SEARCH;

int VERBOSE = 0;

/* a position in the image */
typedef struct {
	int x;
	int y;
	int pos;
} idxelement ;


/* Auxiliary functions ------------------------------------------------------------------ */

void write_double_image (double *im, int ncol, int nrow, char *filename)
{
	FILE *file = fopen(filename,"w");

	fprintf(file,"%d\n%d\n",ncol,nrow);

	int ii;
	for (ii = 0; ii < ncol*nrow; ii++)
		fprintf(file,"%f\n",im[ii]);

	fclose(file);
}

idxelement sample_inside_circle (idxelement center, int radius, iImage msk)
{
	int mxtrials = 5;       /* parameter */

	idxelement rnd;

	/* repeat: pick a random point in box of side (2*radius[kk] + 1) 
	 * until : the picked point is in msk                           */
	int stop = 0;
	int ntrials = 0;

	/* search box */
	int box_x[2] = { max(0, center.x - radius), min(msk->ncol, center.x + radius) };
	int box_y[2] = { max(0, center.y - radius), min(msk->nrow, center.y + radius) };
	int box_ncol = box_x[1] - box_x[0];
	int box_nrow = box_y[1] - box_y[0];
	while (!stop) {

		int tmpos = (int) (box_ncol*box_nrow* ((double)rand() / ((double)RAND_MAX + 1.0)));
		rnd.x = box_x[0] + tmpos % box_ncol;
		rnd.y = box_y[0] + tmpos / box_ncol;
		rnd.pos = rnd.x + msk->ncol*rnd.y;

		/* by construction, this shouldn't happen TODO: remove this check //
		if(rnd.x < 0 || rnd.y < 0 || rnd.x >= msk->ncol || rnd.y >= msk->nrow) 
		{
			printf("WFT??\n");
			rnd.pos = -1;
			rnd.x = -1;
			rnd.y = -1;
			FAILED_RANDOM_SEARCH++;
			return rnd;
		} //*/

		/* check if rnd is in image and in msk2 */
		if (msk->val[rnd.pos] >= 0) {
			stop = 1;
			SUCCESFUL_RANDOM_SEARCH++;
		} else FAILED_RANDOM_SEARCH++; 

		/* K K K A A A H H H ! ! ! ! */
		if (ntrials > mxtrials) { 
			rnd.pos = -1;
			rnd.x   = -1;
			rnd.y   = -1;
			return rnd;
			if (VERBOSE) printf("KKKAAAHHH ! ! ! > no candidate found \n");
		}

		ntrials++;
	}

	return rnd;
}

int update_score_n_merge (list dst, idxelement idst, list src, iImage msk2)
{
	/* list to store offsets with updated scores */
	list candidate_list = NewList();
	candidate_list->quant 		= 0;
	candidate_list->bestScore  = 0;
	candidate_list->worstScore = 0;

	int j;
	for (j = 0; j < src->quant; j++) 
	{ 
		/* check if proposed offset is valid for dst */
		int cand_x   = idst.x + src->l[j].ox;
		int cand_y   = idst.y + src->l[j].oy;

		if ((cand_x >= 0) && (cand_x < msk2->ncol) &&
		    (cand_y >= 0) && (cand_y < msk2->nrow)) 
		{
			int cand_pos = msk2->val[cand_x + cand_y * msk2->ncol];
			if (cand_pos >= 0)
			{
				int qq = candidate_list->quant;

				/* update score */
				candidate_list->l[qq].ox = src->l[j].ox;
				candidate_list->l[qq].oy = src->l[j].oy;
				candidate_list->l[qq].score = PATCH_DISTANCE (idst.pos, cand_pos, globaldata) ;
				candidate_list->quant++;
			}
		}
	}

	/* merge candidate list with current list */
	SortList(candidate_list);
	int changes = MergeListsIntoA(dst,candidate_list);

	DelList(candidate_list);

	return changes;
}

int update_nnf_score ( idxelement *idx1, int lidx1, idxelement *idx2, int lidx2,
		iImage msk2, list nnf) 
{
	int i,j;

	/* fill the scores */
	for (i = 0; i < lidx1; i++) {
		idxelement *p1 = idx1+i;

		for (j = 0; j < nnf[ p1->pos ].quant; j++) {
			idxelement  p2;
			p2.x = p1->x + nnf[ p1->pos ].l[j].ox;
			p2.y = p1->y + nnf[ p1->pos ].l[j].oy;

			/* check that everything is ok */
			if( p2.x >= 0 && p2.x < msk2->ncol &&
			    p2.y >= 0 && p2.y < msk2->nrow ) {

				p2.pos = msk2->val[p2.x + p2.y * msk2->ncol];
				if (p2.pos >= 0) {
					nnf[p1->pos].l[j].score = PATCH_DISTANCE (p1->pos, p2.pos, globaldata);

				} else {
					/* the field is not correctly defined */
					if (VERBOSE)
						printf("the field is not correctly defined: p2.pos = %d\n", p2.pos);
					return 0;
					/* FIXME : lo correcto seria eliminar de la lista este elemento */
					nnf[p1->pos].l[j].score = (double)INT_MAX;
				}
			} else {
				/* the field is not correctly defined */
				if (VERBOSE) {
					printf("the field is not correctly defined: \n"
						"p1 = [%d,%d] ~ idx = %d ~ j = %d ~ ox = [%d,%d]\n",
						p1->x, p1->y, p1->pos, j,
						nnf[ p1->pos ].l[j].ox, nnf[ p1->pos ].l[j].oy);
					printf("p2 = [%d,%d] ~ mask-size = [%d,%d]\n", p2.x, p2.y,
						msk2->ncol, msk2->nrow );
				}
				return 0;
			}
		}

		SortList(&nnf[ p1->pos ]);
	}
	return 1;
}

void random_initialize_nnf (idxelement *idx1, int lidx1, idxelement *idx2, int lidx2, iImage msk2, list nnf) 
{
	int j,i;

	/* fill the random offset and scores */
	for (i = 0; i < lidx1; i++) {
		idxelement p1 = idx1[i];
		nnf[ p1.pos ].quant = LIST_MAX;
		for (j = 0; j < LIST_MAX; j++) {
			int repeated;
			do { 
				int pos2 = (int) ((((double)rand()) / ((double)RAND_MAX + 1.0)) * lidx2 );
				idxelement p2 = idx2[pos2];
				nnf[ p1.pos ].l[j].ox = p2.x-p1.x;
				nnf[ p1.pos ].l[j].oy = p2.y-p1.y;

				repeated = 0;
				int k = 0;
				while (k < j && !repeated)
				{
					repeated = (nnf[p1.pos].l[j].ox == nnf[p1.pos].l[k].ox &&
					            nnf[p1.pos].l[j].oy == nnf[p1.pos].l[k].oy );
					++k;
				}
			} while (repeated);
		}
	}

	/* update the cores for this field */
	update_nnf_score (idx1, lidx1, idx2, lidx2, msk2, nnf);

	/* for debugging
	for (i = 0; i < lidx1; i++) {
		idxelement p1 = idx1[i];
		nnf[ p1.pos ].quant = LIST_MAX;
		for (j = 0; j < LIST_MAX; j++) {
			if (nnf[ p1.pos ].l[j].score == 0) {
				printf("[%d,%d,%d] ~ [%d,%d]\n", p1.x, p1.y, j, 
						nnf[p1.pos].l[j].ox, nnf[p1.pos].l[j].oy);
			}
		}
	}//*/
}

idxelement* compute_idxelement(iImage msk, int *lidx) 
{
	int i,j;
	idxelement * idx = calloc(msk->ncol*msk->nrow,sizeof(idxelement)) ;  

	*lidx = 0;
	for (j=0;j<msk->nrow;j++) { 
		for (i=0;i<msk->ncol;i++) { 
			int pos = i+j*msk->ncol;
			if (msk->val[pos] >= 0 ) {
				idx[*lidx].x = i;
				idx[*lidx].y = j;
				idx[*lidx].pos = msk->val[pos];
				(*lidx)++;
			} 
		}
	}
	return idx;
}

extern list NNfield_create_nnflists ( int size ) 
{
	int i;
	list nnf = (list)calloc(size,sizeof(struct list_s));
	for (i = 0; i < size ; i++) {
		nnf[i].l = (listElement*)calloc(LIST_MAX,sizeof(listElement));
		nnf[i].quant = 1;
	}
	return (nnf);
}

extern void NNfield_destroy_nnflists ( list nnf, int size ) 
{
	int i=0;
	for (i = 0; i < size ; i++) { free(nnf[i].l); }
	free(nnf);
}


/* Joint propagation and seach ---------------------------------------------------------- */

int forward_pass(iImage msk1, iImage msk2,
		idxelement * idx1, int lidx1, idxelement * idx2, int lidx2,
		double alpha, int w, list nnf)
{
	int changes = 0;

	/* radii of concentric regions for random search */
	int nradii = (int)( log( 1.0 / ((double)w) ) / log(alpha) ); 
	int *radii = calloc(nradii,sizeof(int));
	double radii_tmp = w;
	radii[0] = (int)w;
	for (int i = 1; i < nradii; ++i) {
		radii_tmp = alpha*radii_tmp;
		radii[i] = (int)(radii_tmp);
	}

	/* alloc list for search candidates */
	list candidate_list = NewList();
	idxelement i2rnd;
//	int samples_per_radii = LIST_MAX/nradii + 1;
	int samples_per_radii = 1;

	/* loop through maks 1 */
	for (int i = 0; i < lidx1; i++) 
	{
		int idx = idx1[i].pos;
		int px  = idx1[i].x;
		int py  = idx1[i].y;

		list current_list = &nnf[idx];

		/* propagation of neighbor's offsets */
		for(int t = 0; t < NUM_NEIGH; t++) 
		{
			int neighx = px + NEIGH_FWD[ 2*t + 0 ];
			int neighy = py + NEIGH_FWD[ 2*t + 1 ];
			int idxneigh = msk1->val[neighx + neighy * msk1->ncol];

			/* check if neighbor is in msk1 before retriving its offset */
			if( idxneigh >= 0) 
				changes += update_score_n_merge(current_list, idx1[i], &nnf[idxneigh], msk2);
		}

		/* random search */
		for (int ll = 0; ll < current_list->quant; ll++)
		{
			idxelement i2;
			i2.x = px + current_list->l[ll].ox;
			i2.y = py + current_list->l[ll].oy;

			/* for each radii */
			candidate_list->quant = 0;
			candidate_list->bestScore = 0;
			candidate_list->worstScore = 0;
			for (int kk = 0; kk < nradii ; kk++) 
			for (int jj = 0; jj < samples_per_radii; jj++)
			{
				i2rnd = sample_inside_circle (i2,radii[kk],msk2);
				if (i2rnd.pos < 0) continue;

				listElement tmp_el;
				tmp_el.ox = i2rnd.x - px;
				tmp_el.oy = i2rnd.y - py;
				tmp_el.score = PATCH_DISTANCE(idx, msk2->val[i2rnd.pos], globaldata);
				AddToList(candidate_list,tmp_el);
			}

			changes += MergeListsIntoA(current_list,candidate_list);
		}
	}

	free(radii);
	DelList(candidate_list);

	return changes;
}

int backward_pass(iImage msk1, iImage msk2,
		idxelement * idx1, int lidx1, idxelement * idx2, int lidx2,
		double alpha, int w, list nnf)
{ 
	int changes = 0;

	/* radii of concentric regions for random search */
	int nradii = (int)( log( 1.0 / ((double)w) ) / log(alpha) ); 
	int *radii = calloc(nradii,sizeof(int));
	double radii_tmp = w;
	radii[0] = (int)w;
	for (int i = 1; i < nradii; ++i) {
		radii_tmp = alpha*radii_tmp;
		radii[i] = (int)(radii_tmp);
	}

	/* alloc list for search candidates */
	list candidate_list = NewList();
	idxelement i2rnd;
	int samples_per_radii = LIST_MAX/nradii + 1;

	/* loop through maks 1 */
	for (int i = lidx1 - 1; i >= 0; --i)
	{
		int idx = idx1[i].pos;
		int px  = idx1[i].x;
		int py  = idx1[i].y;

		list current_list = &nnf[idx];

/*		if (px == 90 && py == 41)
		{
				printf("BEF PROP: p1 = [%d,%d] ~ idx = %d:\n", px, py, idx);
				for (int j = 0; j < current_list->quant; ++j)
					printf("\tj = %d ~ o = [%d,%d] ~ s = %f\n", j, 
					current_list->l[j].ox, current_list->l[j].oy,
					current_list->l[j].score);
		} //*/

		/* propagation of neighbor's offsets */
		for(int t = 0; t < NUM_NEIGH; t++)
		{
			int neighx = px + NEIGH_BWD[ 2*t + 0 ];
			int neighy = py + NEIGH_BWD[ 2*t + 1 ];
			int idxneigh = msk1->val[neighx + neighy * msk1->ncol];

			/* check if neighbor is in msk1 before retriving its offset */
			if( idxneigh >= 0) 
				changes += update_score_n_merge(current_list, idx1[i], &nnf[idxneigh], msk2);
		}

		/* random search */
		for (int ll = 0; ll < current_list->quant; ll++)
		{
			idxelement i2;
			i2.x = px + current_list->l[ll].ox;
			i2.y = py + current_list->l[ll].oy;

			/* for each radii */
			candidate_list->quant = 0;
			candidate_list->bestScore = 0;
			candidate_list->worstScore = 0;
			for (int kk = 0; kk < nradii ; kk++) 
			for (int jj = 0; jj < samples_per_radii; jj++)
			{
				i2rnd = sample_inside_circle (i2,radii[kk],msk2);
				if (i2rnd.pos < 0) continue;

				listElement tmp_el;
				tmp_el.ox = i2rnd.x - px;
				tmp_el.oy = i2rnd.y - py;
				tmp_el.score = PATCH_DISTANCE(idx, msk2->val[i2rnd.pos], globaldata);

				AddToList(candidate_list,tmp_el);
			}

			changes += MergeListsIntoA(current_list,candidate_list);
		}
	}

	free(radii);
	DelList(candidate_list);

	return changes;
}


/* Propagation of candidates in domain 1 (decoupled from search) ------------------------ */

int forward_propagation (iImage msk1, iImage msk2, idxelement * idx1, int lidx1, idxelement * idx2, int lidx2, list nnf) 
{ 

	int i,t;
	int changes = 0;

	/* loop through maks 1 */
	for (i = 0; i < lidx1; i++) 
	{ 

		int idx = idx1[i].pos;
		int px  = idx1[i].x;
		int py  = idx1[i].y;

		list current_list = &nnf[idx];

		/* test neighborhood offsets */
		for(t=0;t<NUM_NEIGH;t++) 
		{
			
			int neighx = px + NEIGH_FWD[ 2*t + 0 ];
			int neighy = py + NEIGH_FWD[ 2*t + 1 ];
			int idxneigh = msk1->val[neighx + neighy * msk1->ncol];

			/* check if neighbor is in msk1 before retriving its offset */
			if( idxneigh >= 0) 
				changes += update_score_n_merge (current_list,idx1[i],&nnf[idxneigh],msk2);

		}
	}

	return changes;
}

int backward_propagation (iImage msk1, iImage msk2, idxelement * idx1, int lidx1, idxelement * idx2, int lidx2, list nnf) 
{ 

	int i,t;
	int changes = 0;

	/* loop through maks 1 */
	for (i = lidx1 - 1; i >= 0; i--) 
	{ 

		int idx = idx1[i].pos;
		int px  = idx1[i].x;
		int py  = idx1[i].y;

		list current_list = &nnf[idx];

		/* test neighborhood offsets */
		for(t=0;t<NUM_NEIGH;t++) 
		{

			int neighx = px + NEIGH_BWD[ 2*t + 0 ];
			int neighy = py + NEIGH_BWD[ 2*t + 1 ];
			int idxneigh = msk1->val[neighx + neighy * msk1->ncol];

			/* check if neighbor is in msk1 before retriving its offset */
			if( idxneigh >= 0) 
				changes += update_score_n_merge (current_list,idx1[i],&nnf[idxneigh],msk2);

		}
	}

	return changes;
}


/* Search for candidates in domain 2 (decoupled from propagation) ----------------------- */

int random_search_all ( idxelement *idx1, int lidx1, idxelement *idx2, int lidx2, iImage msk2, double alpha, int w, list nnf) 
{

	int    ii,jj,kk,ll;
	int    changes=0;

	/* radius of concentric regions */
	int    nradius = (int)( log( 1.0 / ((double)w) ) / log(alpha) ); 
	int   *radius = calloc(nradius,sizeof(int));
	double radius_tmp = w;
	radius[0] = (int)w;
	for (ii = 1; ii < nradius; ii++) {
		radius_tmp = alpha*radius_tmp;
		radius[ii] = (int)(radius_tmp);
	}

	list candidate_list = NewList();

	/* go through indices in msk1 */
	idxelement i2rnd;
	int samples_per_radius = LIST_MAX/nradius + 1;
	for (ii = 0; ii < lidx1; ii++) {

		idxelement *i1 = idx1 + ii;

		/* coordinates of current destination pixel */
		list current_list = &nnf[ i1->pos ];

		for (ll = 0; ll < current_list->quant; ll++)
		{
			idxelement i2;
			i2.x = i1->x + current_list->l[ll].ox;
			i2.y = i1->y + current_list->l[ll].oy;

			/* for each radius */
			candidate_list->quant 		= 0;
			candidate_list->bestScore  = 0;
			candidate_list->worstScore = 0;
			for (kk = 0; kk < nradius ; kk++) 
			for (jj = 0; jj < samples_per_radius; jj++) {

				i2rnd = sample_inside_circle (i2,radius[kk],msk2);
				if (i2rnd.pos < 0) continue;

				listElement tmp_el;
				tmp_el.ox = i2rnd.x - i1->x;
				tmp_el.oy = i2rnd.y - i1->y;
				tmp_el.score = PATCH_DISTANCE(i1->pos, msk2->val[i2rnd.pos], globaldata);
				AddToList(candidate_list,tmp_el);

			}

			/*SortList(candidate_list);*/
			changes += MergeListsIntoA(current_list,candidate_list);
		}
	}

	free(radius);
	DelList(candidate_list);
	return changes;

}


/* Main function ------------------------------------------------------------------------ */

/* ONLY COMPLETE PATCHES */
extern void NNfield_basic(int *mask1, int ncol1, int nrow1,
                          int *mask2, int ncol2, int nrow2,
                          int  maxit, list nnf, void *input_data,
                          DOUBLE (*patch_dist)(int, int , void *))
{
	/* Attention if using from matlab, ncol = nrow and viceversa */

	int i,j,it;

	double alpha = 0.5;
	int w = max(ncol2, nrow2);

	/* global variables */
	globaldata = input_data;
	globaldistfunc= patch_dist;

	/* random seed */
	srand(time(NULL));

#ifdef DEBUG_STATS
	{
		srand(0);
		printf("\nWARNING NNF NOT RANDOM \n");
		FAILED_RANDOM_SEARCH=0;
		SUCCESFUL_RANDOM_SEARCH=0;
	}
#endif

	/* copy the images to a framed image */
	iImage msk1 = copy_data_to_iimage_with_frame ( ncol1, nrow1, 1, 1, 1, OUT_OF_BOUND, mask1);
	iImage msk2 = copy_data_to_iimage_with_frame ( ncol2, nrow2, 1, 1, 1, OUT_OF_BOUND, mask2);
	/* NOTE: we add only a one pixel frame. This is only used in forward and backward prop.
	 * to avoid verifying if the neighbors are in the image. */

	/* fill the index arrays */ 
	int lidx1=0; 
	int lidx2=0;
	idxelement * idx1 = compute_idxelement( msk1 , &lidx1 );
	idxelement * idx2 = compute_idxelement( msk2 , &lidx2 );
	/* NOTE: idx1.pos is an index of a position in the image without the reference frame. Where as 
	 *       idx1.x and idx1.y are coordinates of positions with the reference frame.  */

#ifdef DEBUG_STATS
	clock_t start_time = clock();
	FILE *nnfstats    = fopen("stats.nnf","w");
	FILE *nnfchanges  = fopen("changes.nnf","w");
	{
		/* file headers */
		fprintf(nnfstats,"%% Patch-Match with lists.\n"); 
		fprintf(nnfstats,"%% 	list length = %d \n",LIST_MAX);
		fprintf(nnfstats,"%% 	iterations  = %d \n%%\n",maxit);
	}
#endif

	/* retrieve the input field and update scores */
	int valid_nnf = update_nnf_score (idx1, lidx1, idx2, lidx2, msk2, nnf);

	/* if that failed, initialize randomly */
	if (!valid_nnf)
	{
		if (VERBOSE) printf("found invalid nnf ~ init!\n");
		random_initialize_nnf (idx1, lidx1, idx2, lidx2, msk2, nnf);
	}
//	else
//		printf("valid nnf found - no random init\n");

#ifdef DEBUG_STATS
	if (!valid_nnf) {
		/* log time taken by random init and global nnf energy */
		fprintf(nnfstats  ,"%f ",((double)clock() - start_time)/ CLOCKS_PER_SEC);
		for (j = 0; j < LIST_MAX; j++) {
			double score = 0;
			for (i = 0; i < lidx1; i++) score += nnf[i].l[j].score;
			fprintf(nnfstats  ,"%f "   ,score);
		}
		fprintf(nnfstats,"\n");
	}
#endif

	/* main loop ------------------------------------------------------------------- */
	int changes;
	int phase = (maxit < 0);
	maxit = phase ? -maxit : maxit;
	for(it = 0; it < maxit; it++) {

		changes = 0;

		if ( (it + phase) % 2 ) {
			changes = forward_pass (msk1, msk2, idx1, lidx1, idx2, lidx2, alpha, w, nnf);
			if (VERBOSE) printf("% 3d : forward\n", it);
		}
		else {
			changes = backward_pass(msk1, msk2, idx1, lidx1, idx2, lidx2, alpha, w, nnf);
			if (VERBOSE) printf("% 3d : backward\n", it);
		}

#ifdef DEBUG_STATS
		{
			fprintf(nnfstats  ,"%f ",((double)clock() - start_time)/ CLOCKS_PER_SEC);
			for (j = 0; j < LIST_MAX; j++) 
			{
				double score = 0;
				for (i = 0; i < lidx1; i++) score += nnf[i].l[j].score;
				fprintf(nnfstats  ,"%f "   ,score);
			}
			fprintf(nnfstats,"\n");
			fprintf(nnfchanges,"%d\n",changes);
		}
#endif

	}

#ifdef DEBUG_STATS
	{
		fclose(nnfstats);
		fclose(nnfchanges);
		double TOTAL = (double)(FAILED_RANDOM_SEARCH + SUCCESFUL_RANDOM_SEARCH);
		printf("       -> rnd search trials: failed: %g%% succesful: %g%% \n",
				(double)FAILED_RANDOM_SEARCH   /TOTAL*100,
				(double)SUCCESFUL_RANDOM_SEARCH/TOTAL*100);
		}
#endif

	free(idx1);
	free(idx2);
	del_iimage(msk1);
	del_iimage(msk2);

}


/* ONLY COMPLETE PATCHES */
extern void NNfield_decoupled(int *mask1, int ncol1, int nrow1,
                              int *mask2, int ncol2, int nrow2,
                              int  maxit, list nnfinout, void *input_data,
                              DOUBLE (*patch_dist)(int, int , void *))
{
	/* Attention if using from matlab, ncol = nrow and viceversa */

	int i,j,it;

	double alpha = 0.5;
	int w = max(ncol2, nrow2);

	/* global variables */
	globaldata = input_data;
	globaldistfunc= patch_dist;


#ifndef DEBUG_STATS
	srand(time(NULL));
#else
	srand(0);
	printf("\nWARNING NNF NOT RANDOM \n");
	FAILED_RANDOM_SEARCH=0;
	SUCCESFUL_RANDOM_SEARCH=0;
#endif


	/* copy the images to a framed image */
	iImage msk1 = copy_data_to_iimage_with_frame ( ncol1, nrow1, 1, 1, 1, OUT_OF_BOUND, mask1);
	iImage msk2 = copy_data_to_iimage_with_frame ( ncol2, nrow2, 1, 1, 1, OUT_OF_BOUND, mask2);
	/* NOTE: we add only a one pixel frame. This is only used in forward and backward prop.
	 * to avoid verifying if the neighbors are in the image. */

	/* fill the index arrays */ 
	int lidx1=0; 
	int lidx2=0;
	idxelement * idx1 = compute_idxelement( msk1 , &lidx1 );
	idxelement * idx2 = compute_idxelement( msk2 , &lidx2 );
	/* NOTE: idx1.pos is an index of a position in the image without the reference frame. Where as 
	 *       idx1.x and idx1.y are coordinates of positions with the reference frame.  */

	/* retrieve the input field */
	list nnf = nnfinout;
	list nnfPREV = (list) calloc(lidx1,sizeof(struct list_s));
	for (i = 0; i < lidx1; i++) 
	{
		nnfPREV[i].l = (listElement*)calloc(LIST_MAX,sizeof(listElement));
		CopyListToA(&nnfPREV[i],&nnf[i]);
	}

	int valid_nnfPREV = update_nnf_score (idx1, lidx1, idx2, lidx2, msk2, nnfPREV);


#ifdef DEBUG_STATS /* -------------------------------------------- */
	clock_t start_time = clock();
	FILE *nnfstats    = fopen("stats.nnf","w");
	FILE *nnfchanges  = fopen("changes.nnf","w");

	fprintf(nnfstats,"%% Patch-Match with lists.\n"); 
	fprintf(nnfstats,"%% 		list length = %d \n",LIST_MAX);
	fprintf(nnfstats,"%% 		iterations  = %d \n%%\n",maxit);
#endif             /* -------------------------------------------- */


	/* initialize the field with random offsets and compute the scores */
	random_initialize_nnf (idx1, lidx1, idx2, lidx2, msk2, nnf);

	/*Merge the random field with the initialization */
	if (valid_nnfPREV) 
		for (i=0;i<lidx1;i++) MergeListsIntoA(&nnf[i],&nnfPREV[i]); 

	for (i = 0; i< lidx1; i++) free(nnfPREV[i].l);
	free(nnfPREV);

#ifdef DEBUG_STATS
		fprintf(nnfstats  ,"%f ",((double)clock() - start_time)/ CLOCKS_PER_SEC);

		for (j = 0; j < LIST_MAX; j++) 
		{
			double score = 0;
			for (i = 0; i < lidx1; i++) score += nnf[i].l[j].score;
			fprintf(nnfstats  ,"%f "   ,score);
		}
		fprintf(nnfstats,"\n");
#endif

	/* main loop ------------------------------------------------------------------- */
	int changes_pro,changes_sch;
	for(it = 0 ;it < maxit; it++)  
	{
		changes_pro = 0;

		if ( it % 2 ) changes_pro = forward_propagation  (msk1, msk2, idx1, lidx1, idx2, lidx2, nnf);
		else          changes_pro = backward_propagation (msk1, msk2, idx1, lidx1, idx2, lidx2, nnf);

#ifdef DEBUG_STATS
		fprintf(nnfstats  ,"%f ",((double)clock() - start_time)/ CLOCKS_PER_SEC);

		for (j = 0; j < LIST_MAX; j++) 
		{
			double score = 0;
			for (i = 0; i < lidx1; i++) score += nnf[i].l[j].score;
			fprintf(nnfstats  ,"%f "   ,score);
		}
		fprintf(nnfstats,"\n");
#endif

		changes_sch = random_search_all (idx1,lidx1, idx2,lidx2, msk2, alpha, w, nnf);


#ifdef DEBUG_STATS
		fprintf(nnfchanges,"%d %d\n",changes_pro,changes_sch);
		fprintf(nnfstats  ,"%f ",((double)clock() - start_time)/ CLOCKS_PER_SEC);
		for (j = 0; j < nnf[0].quant; j++) 
		{
			double score = 0;
			for (i = 0; i < lidx1; i++) score += nnf[i].l[j].score;
			fprintf(nnfstats  ,"%f "   ,score);
		}
		fprintf(nnfstats,"\n");
#endif

	}


#ifdef DEBUG_STATS
	fclose(nnfstats);
	fclose(nnfchanges);
	double TOTAL = (double)(FAILED_RANDOM_SEARCH + SUCCESFUL_RANDOM_SEARCH);
	printf("       -> rnd search trials: failed: %g%% succesful: %g%% \n", (double)FAILED_RANDOM_SEARCH   /TOTAL*100, 
			                                                                 (double)SUCCESFUL_RANDOM_SEARCH/TOTAL*100);
#endif

	free(idx1);
	free(idx2);
	del_iimage(msk1);
	del_iimage(msk2);

}


/* ONLY COMPLETE PATCHES */
extern void NNfield_exhaustive(int *mask1, int ncol1, int nrow1, 
                               int *mask2, int ncol2, int nrow2,
                               int  maxit, list nnfinout, void *input_data,
                               DOUBLE (*patch_dist)(int, int , void *))
{
	int i,j;

	/* global variables */
	globaldata = input_data;
	globaldistfunc= patch_dist;

	/* copy the images to a framed image */
	iImage msk1 = copy_data_to_iimage_with_frame ( ncol1, nrow1, 1, 1 , 1, OUT_OF_BOUND, mask1);
	iImage msk2 = copy_data_to_iimage_with_frame ( ncol2, nrow2, 1, 1 , 1, OUT_OF_BOUND, mask2);

	list nnf = nnfinout;

	/* fill the index arrays */ 
	int lidx1=0;
	int lidx2=0;
	idxelement * idx1 = compute_idxelement( msk1 , &lidx1 );
	idxelement * idx2 = compute_idxelement( msk2 , &lidx2 );
	/* NOTE: idx1.pos is an index of a position in the image without the reference frame. Where as 
	 *       idx1.x and idx1.y are coordinates of positions with the reference frame.  */

	random_initialize_nnf (idx1, lidx1, idx2, lidx2, msk2, nnf);

	listElement tmp_el;
//	#pragma omp parallel for
	for (i = 0; i < lidx1; i++)
	{
		idxelement i1 = idx1[i];
		if (i % 1000 == 0)
			printf("    [%04d,%04d]\n",i,lidx1);

		/* search for minimum */
		for (j = 0; j < lidx2; j++)
		{
			idxelement i2 = idx2[j];

			tmp_el.ox = i2.x - i1.x;
			tmp_el.oy = i2.y - i1.y;
			tmp_el.score = PATCH_DISTANCE(i1.pos,i2.pos,globaldata);
			AddToList(&nnf[i1.pos],tmp_el);

		}
	}

#ifdef DEBUG_STATS /* -------------------------------------------- */
	clock_t start_time = clock();
	FILE *nnfstats = fopen("stats.exhnnf","w");

	fprintf(nnfstats,"%% Exhaustive NNF with lists.\n"); 
	fprintf(nnfstats,"%% 		list length = %d \n",LIST_MAX);

	for (j = 0; j < nnf[0].quant; j++) 
	{
		double score = 0;
		for (i = 0; i < lidx1; i++) score += nnf[i].l[j].score;
		fprintf(nnfstats,"%f ",score);
	}
	fprintf(nnfstats,"\n");

	fclose(nnfstats);
#endif             /* -------------------------------------------- */

	free(idx1);
	free(idx2);
	del_iimage(msk1);
	del_iimage(msk2);
}
