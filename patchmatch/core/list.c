#include <stdio.h>
#include <stdlib.h>
#include "list.h"

/* THIS IS A GLOBAL VARIABLE IT ESTABLISHES THE MAXIMUM SIZE OF THE LISTS */
int list_max_global = 10;
#define LIST_MAX list_max_global

#define SWAPel(a,b) {listElement swap=a;a=b;b=swap;}


/* create a single list */
list NewList ( void ) { 
	list a = calloc( 1 , sizeof(struct list_s));
	a->l = (listElement*) calloc( LIST_MAX , sizeof(listElement));
	a->quant=0;
	a->bestScore=0;
	a->worstScore=0;

	return a;
}

/* free list's memory */
void DelList ( list a) {
	free (a->l);
	free (a);
}

/* creates a new list copied from s */
list CopyList(list s) {
   int i;

   list d = NewList();

   for (i = 0; i < s->quant; i++) { d->l[i] = s->l[i]; }
	/* TODO try this: memcpy ( d->l, a->l, a->quant*sizeof(labelListElement) ); */

   d->bestScore  = s->bestScore;
   d->worstScore = s->worstScore;
   d->quant = s->quant;

   return d;
}

/* copies contents of b into a */
void CopyListToA(list a,list b) {
   int i;

   for (i = 0; i < b->quant; i++) { a->l[i] = b->l[i]; }
	/* TODO memcpy ( a->l, c->l, mx_length*sizeof(labelListElement) ); */

   a->bestScore  = b->bestScore;
   a->worstScore = b->worstScore;
   a->quant = b->quant;
}

/* comparison function for qsort ~ ascending order */
int cmp_listElement(const void *a, const void *b) {
   const listElement *da = ((const listElement*) a);
   const listElement *db = ((const listElement*) b);
   return (da->score > db->score) - (da->score < db->score);
}

/* sorts the list a without repetitions */
void SortList(list a) {
	/* trivial case */
	if(a->quant == 0) 
		return;

	if(a->quant == 1) 
	{
		a->bestScore  = a->l[0].score;
		a->worstScore = a->l[0].score;
		return;
	}

	/* allocate a new listElement e */
	listElement *e = (listElement*) calloc( LIST_MAX, sizeof(listElement));

	/* quicksort a->l */
	qsort(a->l,a->quant,sizeof(listElement),*cmp_listElement);
	
	/* scan a->l for repetitions */
	double prev = -10000000000.;
	int ie=0, i=0;
	while (i < a->quant) {
		if((prev == a->l[i].score) && (e[ie-1].ox == a->l[i].ox) && (e[ie-1].oy == a->l[i].oy) )  {
			i++;
		} else {
			e[ie] = a->l[i];
			prev = e[ie].score;
			ie++;
			i++;
		}
	}

	/* dealloc a->l */
	free(a->l);	
	/* replace a->l = e */
	a->l = e;
	a->quant = ie;
	a->bestScore = e[0].score;
	a->worstScore = e[ie-1].score;
}

/* adds an element to the list with bubble sort - w/out repetitions */
int AddToList(list a, listElement e) {						
	/*
	 * It checks if e is in the list. If it is, it does not add it. 
	 * Returns 1 if e has been added to the list.
	 *
	 */

	/* bubble sort the element e into a->l */
	if (e.score < a->worstScore || a->quant < LIST_MAX) { 

		/* determine if e is not already in a->l */ /* TODO: DYADIC SEARCH IS FASTER */
		int i=a->quant-1;
		while (i>=0) {
			if( ( a->l[i].score == e.score ) && ( a->l[i].ox == e.ox ) && ( a->l[i].oy == e.oy ) ) 
				return 0;
			i--;
		}

		/* insert e into a->l */
		if(a->quant == LIST_MAX) {
			a->l[a->quant-1] = e;
		} else {
			a->l[a->quant] = e;
			a->quant++;
		}

		/* REORDER THE LIST */
    	i=a->quant-1;
    	while (i>=1 && a->l[i].score < a->l[i-1].score) {
      	SWAPel(a->l[i], a->l[i-1]);
      	i--;
    	}
		a->bestScore  = a->l[    0     ].score;
		a->worstScore = a->l[a->quant-1].score;
		return 1;
	}
	return 0;
}

/* merges list a and b into a ~ w/out repetitions */
int MergeListsIntoA(list a, list b) {
	/* allocate a new list e */
	listElement *e  = (listElement*) calloc( LIST_MAX, sizeof(listElement));

	/* merge all the elements a->l, b->l into e*/
	int ia=0,ib=0,ie=0,changes=0; 								                         /* indices a, b and e*/
	if(a->quant == b->quant == LIST_MAX) {
		while(ie < LIST_MAX) {
			if(a->l[ia].score < b->l[ib].score) {                                    /* a less than b */
				e[ie++] = a->l[ia++];
			} else { if(a->l[ia].score > b->l[ib].score) {                           /* b less than a */
				e[ie++] = b->l[ib++];
				changes++;
			} else {                                                                 /* same score    */
				if( (a->l[ia].ox == b->l[ib].ox) && (a->l[ia].oy == b->l[ib].oy) ) {  /* same element  */ 
					e[ie++] = a->l[ia++];
					ib++;
				} else {                                                              /* different element */
					e[ie++] = a->l[ia++];
					if(ie < LIST_MAX) {
					   e[ie++] = b->l[ib++];
					   changes++;
					}
					
				}
			}
			}
		}
	} else {                                                                         /* at least one of the lists is not full */
		while(ie < LIST_MAX) {
			if(ia == a->quant) {                                                       /* done with list a */
				if(ib == b->quant){break;}                                              /* done also with list b: finished! */ 
				e[ie++] = b->l[ib++];
				changes++;
			} else if(ib == b->quant) {                                                /* done with list b but not with list a */
				e[ie++] = a->l[ia++];
			} else {                                                                   /* not done with a nor with b */		
				if(a->l[ia].score < b->l[ib].score) {                                   /* a less than b */
					e[ie++] = a->l[ia++];
				} else { if(a->l[ia].score > b->l[ib].score) {	                        /* b less than a */
					e[ie++] = b->l[ib++];
					changes++;
				} else {																				      /* same score */
					if( (a->l[ia].ox == b->l[ib].ox) && (a->l[ia].oy == b->l[ib].oy) ) { /* same element */ 
						e[ie++] = a->l[ia++];
						ib++;
					} else {										                                 /* different element */
						e[ie++] = a->l[ia++];
						if(ie < LIST_MAX) {
						   e[ie++] = b->l[ib++];
						   changes++;
						}
					}
				}
				}
			}
		}
	}
	/* dealloc a->l */
	free(a->l);	
	/* replace a->l = e */
	a->l = e;
	a->quant = ie;
	a->bestScore = e[0].score;
	a->worstScore = e[ie-1].score;
	return changes;
}






