typedef double DOUBLE;   /* scores, weights and stuff */

#include <time.h>
#include <limits.h>
#include <math.h>
#include "list.h" 

#define OUT_OF_BOUND -1

#define PATCH_DISTANCE (*globaldistfunc)


/*#define DEBUG_STATS*/


/* GLOBAL VARIABLES */
void * globaldata;
DOUBLE (*globaldistfunc)(int, int, void *);

extern list NNfield_create_nnflists ( int size );
extern void NNfield_destroy_nnflists ( list nnf, int size );

extern void NNfield_basic(int *mask1, int ncol1, int nrow1,
                          int *mask2, int ncol2, int nrow2,
                          int  maxit, list nnfinout, void *input_data,
                          DOUBLE (*patch_dist)(int, int , void *));

extern void NNfield_exhaustive(int *mask1, int ncol1, int nrow1,
                               int *mask2, int ncol2, int nrow2,
                               int  maxit, list nnfinout, void *input_data,
                               DOUBLE (*patch_dist)(int, int , void *));
