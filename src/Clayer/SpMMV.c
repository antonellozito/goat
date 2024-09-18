# include <cs.h> 
# include "clayer.h"

MyCSparse SpMMV(MyCSparse A, MyCSparse B)
{
    MyCSparse C;
    MyCSparse *D, *E;
    D=&A;
    E=&B;

    C = SpMM(D, E);

    return C;

}
