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

    return C;

}