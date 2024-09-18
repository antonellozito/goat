# include <cs.h>
/*
Type definitions
*/

/* Own C sparse type for interfacing with Fortran*/
typedef struct MyCSparse
{
    int *row;
    int *col;
    double *val; 
    int nrow;
    int ncol;
    int nval;
} MyCSparse;
cs *ConvertMyCSparseToCS(MyCSparse *);
MyCSparse ParseTest(MyCSparse);
MyCSparse *ParseTestP2(MyCSparse *);
MyCSparse SpMM(MyCSparse *, MyCSparse *);
MyCSparse *ConstructMyCSparse(int, int, int, int [], int [], 
    double []);
MyCSparse SpMMV(MyCSparse, MyCSparse);
