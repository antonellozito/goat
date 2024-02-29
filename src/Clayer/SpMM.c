# include <cs.h> 
# include "clayer.h"
MyCSparse SpMM(const MyCSparse A, const MyCSparse B)
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
    cs *csA, *csB, *csC;
    int i, p1, p2, p;
    MyCSparse C; 

    /* Convert */
    csA = ConvertMyCSparseToCS(A);
    csB = ConvertMyCSparseToCS(B);

    /* Compress to column storage*/
    csA = cs_compress(csA);
    csB = cs_compress(csB);

    /* Compute */
    csC = cs_multiply(csA, csB);

    /* Reconvert */
    C.nrow = csC->m; 
    C.ncol = csC->n; 
    
    C.nval = csC->p[csC->n];
    C.val = csC->x;
    C.row = csC->i; 
    C.col = malloc(C.nval);
    
    for (i = 0; i < C.ncol; i++)
    {
        p1 = csC->p[i];
        p2 = csC->p[i+1];

        for(p = p1; p < p2; p++)
        {
            C.col[p] = i;
        }

    }

    return C;


    
    
}

