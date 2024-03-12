# include <cs.h> 
# include "clayer.h"
MyCSparse SpMMBasic(const int nr1, const int nc1, const int nv1, 
const int nr2, const int nc2, const int nv2, int *r1, int *c1, double *v1, int *r2, int *c2, double *v2)
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
    cs *csA, *csB, *csC, *csA2, *csB2;
    MyCSparse *A, *B;
    int i, p1, p2, p;
    MyCSparse C; 

    /* Construct */
    A = ConstructMyCSparse(nr1, nc1, nv1, r1, c1, v1);
    B = ConstructMyCSparse(nr2, nc2, nv2, r2, c2, v2);
    csA = ConvertMyCSparseToCS(A);
    csB = ConvertMyCSparseToCS(B);

    /* Compress to column storage*/
    csA2 = cs_compress(csA);
    csB2 = cs_compress(csB);

    /* Compute */
    csC = cs_multiply(csA2, csB2);

    /* Reconvert */
    C.nrow = csC->m; 
    C.ncol = csC->n;     
    C.nval = csC->p[csC->n];

    C.val = (double* ) calloc(C.nval, sizeof(double));
    C.row = (int*) calloc(C.nval, sizeof(int));
    C.col = (int*) calloc(C.nval, sizeof(int));
    C.val = csC->x;
    C.row = csC->i; 
    C.col = (int*) calloc(C.nval, sizeof(int));
    
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

