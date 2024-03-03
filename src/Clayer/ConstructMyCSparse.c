# include <cs.h>
# include "clayer.h"

MyCSparse *ConstructMyCSparse(int nrow, int ncol, int nval, int row [], int col [], 
    double val [])
{
    /* Define sparse matrix */
    MyCSparse *A = calloc(1, sizeof(MyCSparse));

    /* Allocate */
    A->row = (int*) calloc(nval, sizeof(int));
    A->col = (int*) calloc(nval, sizeof(int));
    A->val = (double*) calloc(nval, sizeof(double));

    /* Add data */
    A->row = row; 
    A->col = col; 
    A->val = val;
    A->ncol = ncol;
    A->nval = nval;
    A->nrow = nrow;

    return A;

}