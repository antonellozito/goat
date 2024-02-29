# include <cs.h>
# include "clayer.h"

MyCSparse ConstructMyCSparse(int nrow, int ncol, int nval, int row [], int col [], 
    double val [])
{
    /* Define sparse matrix */
    MyCSparse A;

    /* Add data */
    A.row = row; 
    A.col = col; 
    A.val = val;
    A.ncol = ncol;
    A.nval = nval;
    A.nrow = nrow;

    return A;

}