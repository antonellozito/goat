# include <cs.h>
# include "clayer.h"

cs *ConvertMyCSparseToCS(MyCSparse A)
{
    /* Define variables*/
    cs *B;
    int i, j; 

    /* Allocate */
    B = cs_spalloc(A.nrow, A.ncol, A.nval, A.nval, A.nval);

    /* Add entries */
    for (i = 0; i < A.nval; i++)
    {
        cs_entry(B, A.row[i], A.col[i], A.val[i]);
    };

    return B; 
}