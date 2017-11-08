#ifndef LIST_ORDERED_BY_SCORE_H
#define LIST_ORDERED_BY_SCORE_H


/* ACCES TO THE GLOBAL VARIABLE INSIDE THE OBJECT FILE */
/* IT INDICATES THE MAXIMUM SIZE OF THE IMAGE */
extern int list_max_global;
#define LIST_MAX list_max_global

/* offset element */
typedef struct {
   int ox;              /* offset in the x direction */
   int oy;
   int pos;					/* this field is NOT USED */
   double score;			/* score (i.e. patch distance) associated to offset */
} listElement;

/* offset list, ordered by offset score */
typedef struct list_s {
   listElement* l;      /* list of offsets   */
   double bestScore;		/* l[0]->score       */
   double worstScore;	/* l[quant-1]->score */
   int quant;				/* quant number elements in the list */
} *list;


/* create a single list */
list NewList ( void );

/* free list's memory */
void DelList ( list a);

/* creates a new list copied from s */
list CopyList(list s);

/* copies contents of b into a */
void CopyListToA(list a,list b);

/* comparison function for qsort ~ ascending order by score */
int cmp_listElement(const void *a, const void *b);

/* sorts list w/out repetitions */
void SortList(list a) ;

/* adds an element to the list with bubble sort - w/out repetitions */
int AddToList(list a, listElement e);

/* merges list a and b into a ~ w/out repetitions */
int MergeListsIntoA(list a, list b);


/* funciones de debuggeo en matlab */
#ifdef MATLAB_MEX_FILE

/* |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| */
/* funcion de debuggeo: permite analizar una variable en matlab *
 *                      escupe la lista en una estructura       */
int mexPutnBrake_list(
   list   pdata,          /* pointer to the data                */
   char   *name)          /* variable name                      */
{
   /* create array of offsets */
   const char *listElement_fields[] = {"ox","oy","score","pos"};
   mxArray *offset_list = mxCreateStructMatrix(pdata->quant,1,4,listElement_fields);

   /* fill offset_list */
   int i;
   for (i = 0; i < pdata->quant; i++) 
   {
       mxSetField(offset_list,i,"ox"   ,mxCreateDoubleScalar((double)(pdata->l[i].ox   )));
       mxSetField(offset_list,i,"oy"   ,mxCreateDoubleScalar((double)(pdata->l[i].oy   )));
       mxSetField(offset_list,i,"score",mxCreateDoubleScalar((double)(pdata->l[i].score)));
       mxSetField(offset_list,i,"pos"  ,mxCreateDoubleScalar((double)(pdata->l[i].pos  ))); 
   }

   /* create list struct */
   const char *list_fields[] = {"quant","bestScore","worstScore","l"};
   mxArray *ma_list = mxCreateStructMatrix(1,1,4,list_fields);

   mxSetField(ma_list,0,"l"         ,offset_list);
   mxSetField(ma_list,0,"quant"     ,mxCreateDoubleScalar((double)(pdata->quant     )));
   mxSetField(ma_list,0,"bestScore" ,mxCreateDoubleScalar((double)(pdata->bestScore )));
   mxSetField(ma_list,0,"worstScore",mxCreateDoubleScalar((double)(pdata->worstScore))); 

   /* put it in matlab workspace */
   int status = mexPutVariable("caller",name,ma_list);

   /* go to matlab prompt */
   mexPrintf("mexPutnBrake -> %s \n",name);
   mexEvalString("keyboard");

   /* remove variable from workspace */
   char clear_str[100] = "clear ";
   mexEvalString(strcat(clear_str,name));

   /* destroy allocated mxArray structs */
   /* POR ALGUN MOTIVO ESTO PETA!!!
    *
    * mxDestroyArray(offset_list);
   mxDestroyArray(ma_list);*/

   return status;
}

/* |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| */
/* funcion de debuggeo: permite analizar una variable en matlab. para un array de listas de 
 *                      offsets, asociados a posiciones en una imagen, genera dos imagenes 
 *                      con los offsets obtenidos del indice off_idx (argumento de entrada) 
 *                      de cada una de las listas                                         */
int mexPutnBrake_nnf(
   list   pdata,          /* pointer to the data                */
   int    *pos,           /* positions on image                 */
   int    lpos,
   int    off_idx,        /* offset list index                  */ 
   int    *sz,            /* image size                         */
   char   *name)          /* variable name                      */
{
   /* create output mxArrays */
   /*int dims[3] = {sz[0],sz[1],2};
   mxArray *mxnnf = mxCreateDoubleArray(3,dims,mxDOUBLE_CLASS,mxREAL);*/
   mxArray *mxnnfx = mxCreateDoubleMatrix(sz[0],sz[1],mxREAL);
   mxArray *mxnnfy = mxCreateDoubleMatrix(sz[0],sz[1],mxREAL);
   mxArray *mxnnfs = mxCreateDoubleMatrix(sz[0],sz[1],mxREAL);
  /* mxArray *mxnnfx = mxCreateDoubleMatrix(3,dims,mxDOUBLE_CLASS,mxREAL);
   mxArray *mxnnfy = mxCreateDoubleMatrix(3,dims,mxDOUBLE_CLASS,mxREAL);*/

   /* fill output mxArrays with data */
   int i;
   double *p_mxnnfx = mxGetPr(mxnnfx);
   double *p_mxnnfy = mxGetPr(mxnnfy);
   double *p_mxnnfs = mxGetPr(mxnnfs);

   /*for (i = 0; i < lpos; i++) p_mxnnf[        pos[i]] = pdata[i].l[off_idx].ox;
   for (i = 0; i < lpos; i++) p_mxnnf[imdim + pos[i]] = pdata[i].l[off_idx].oy;*/
   for (i = 0; i < lpos; i++) p_mxnnfx[pos[i]] = pdata[i].l[off_idx].ox;
   for (i = 0; i < lpos; i++) p_mxnnfy[pos[i]] = pdata[i].l[off_idx].oy;
   for (i = 0; i < lpos; i++) p_mxnnfs[pos[i]] = pdata[i].l[off_idx].score;
      
   /* put it in matlab workspace */
   /*int status = mexPutVariable("caller",name,mxnnf);*/
   char namex[100]; strcpy(namex,name); strcat(namex,"x");
   char namey[100]; strcpy(namey,name); strcat(namey,"y");
   char names[100]; strcpy(names,name); strcat(names,"s");
   int status  = mexPutVariable("caller",namex,mxnnfx);
       status *= mexPutVariable("caller",namey,mxnnfy);
       status *= mexPutVariable("caller",names,mxnnfs);

   /* go to matlab prompt */
   mexPrintf("mexPutnBrake_nnf -> %s :: list idx = %i\n",name,off_idx);
   mexEvalString("keyboard");

   /* remove variable from workspace */
   char clearx_str[100] = "clear ";
   strcat(clearx_str,namex);
   mexPrintf("mexPutnBrake -> %s \n",clearx_str);
   mexEvalString(clearx_str);

   char cleary_str[100] = "clear ";
   strcat(cleary_str,namey);
   mexPrintf("mexPutnBrake -> %s \n",cleary_str);
   mexEvalString(strcat(cleary_str,name));

   char clears_str[100] = "clear ";
   strcat(clears_str,names);
   mexPrintf("mexPutnBrake -> %s \n",clears_str);
   mexEvalString(strcat(clears_str,name));

   /* destroy allocated mxArray structs */
   /*mxDestroyArray(mxnnf);*/
   /*mxDestroyArray(mxnnfx);
   mxDestroyArray(mxnnfy);*/

   return status;
}

#endif
#endif
