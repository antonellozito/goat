# include <stdio.h>
# include <cs.h>
# include "clayer.h"

/*
Just a small test routine for compilation etc 
*/
main() 
{
    /* Construct a matrix */
    cs *T, *A, *B, *C, *thisCS; 
    MyCSparse this, thismm;
    int rows[]    = {1, 2, 3, 4, 5};
    int cols[]    = {1, 2, 3, 4, 5};
    double vals[] = {0.22, 1.23, 4.56, 4, 5.0}; 
    int i;
    int ncol;
    int nz; 
    int p1, p2, length, p; 

    /* Allocate */
    this = ConstructMyCSparse(10, 10, 5, rows, cols, vals);
    //this.nrow = 10;
    //this.ncol = 10;
    //this.nval = 5;
    //this.row = rows;
    //this.col = cols;
    //this.val = vals;

    T = cs_spalloc(10, 10, 1, 1, 1);
    thismm = SpMM(this, this);



    /* Add entries */
    for (i = 0; i < 5; i++)
    {
        cs_entry(T, rows[i], cols[i], vals[i]);
    };

    /* Print */
    cs_print(T, i);

    /* Try some multiplications*/
    A = cs_compress(T); // convert to csc format for multiplication
    B = cs_multiply(A, A);

    /* Convert to triplet format again - fully based on 
    umfpack_col_to_triplet
     */
    C = B;
    ncol = B->n;
    C->nz = B->p[ncol];
    for (i = 0; i < ncol; i++)
    {
        p1 = B->p[i];
        p2 = B->p[i+1];
        length = p2-p1;

        for(p = p1; p < p2; p++)
        {
            C->p[p] = i;
        }

    }
    cs_print(C, i);


};



