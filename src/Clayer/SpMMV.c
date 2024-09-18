# include <cs.h> 
# include "clayer.h"

MyCSparse SpMMV(MyCSparse A, MyCSparse B)
{
    MyCSparse C;
    MyCSparse *D, *E, *F;
    D=&A;
    E=&B;

    F = SpMM(D, E);
    C = *F;
    free(F->row);
    free(F->col);
    free(F->val);
    free(F);
    return C;

}
