# include <cs.h> 
# include "clayer.h"
# include <stdlib.h> 
# include <string.h>
MyCSparse SpMM(MyCSparse *A, MyCSparse *B)
{
    /*
    Sparse matrix-matrix multiplication wrapper for interfacing to 
    fortran. The structure 'MyCSparse' should be mimicked in fortran 
    for easier manipulation. Under the hood, this routine simply calls 
    the sparse matrix-matrix multiply routine 'cs_multiply' from the 
    SuiteSparse library. Here, the matrices are in CSC format, so the
    resulting matrix must be converted back into triplet format to 
    comply with the format of GOAT. This may lead to some overhead and
    may be revisited in the future. 
    */

    
    /* Local variables */
    cs *csA, *csB, *csC, *csA0, *csB0;
    int i, p1, p2, p;
    MyCSparse C;

    /* Convert */
    csA0 = ConvertMyCSparseToCS(A);
    csB0 = ConvertMyCSparseToCS(B);

    /* Compress to column storage*/
    csA = cs_compress(csA0);
    csB = cs_compress(csB0);

    /* Compute */
    csC = cs_multiply(csA, csB);

    /* Reconvert */
    C.nrow = csC->m; 
    C.ncol = csC->n; 
    
    C.nval = csC->p[csC->n];
    C.val = (double *) calloc(C.nval, sizeof(double));
    C.row = (int *) calloc(C.nval, sizeof(int));
    memcpy(C.val, csC->x, sizeof(double)*C.nval);
    memcpy(C.row, csC->i, sizeof(int)*C.nval);
    C.col = (int *) calloc(C.nval, sizeof(int));
    
    for (i = 0; i < C.ncol; i++)
    {
        p1 = csC->p[i];
        p2 = csC->p[i+1];

        for(p = p1; p < p2; p++)
        {
            C.col[p] = i;
        }

    }
    cs_di_spfree(csA);
    cs_di_spfree(csB);
    cs_di_spfree(csA0);
    cs_di_spfree(csB0);
    cs_di_spfree(csC);

    return (C);


    
    
}

