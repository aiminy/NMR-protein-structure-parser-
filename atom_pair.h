#ifndef ATOM_PAIR_H
#define ATOM_PAIR_H

#include "linked_list.h"


class Atom_Pair
{
 public:
    Atom_Pair();
    void SetFirstAtom(Atom T);
    void SetSecondAtom(Atom N);
    Atom GetFirstAtom();
    Atom GetSecondAtom();
    double Distance_Between_Two_Atom(Atom T, Atom N);

    
 private:

   Atom First;
   Atom Second;
   double d;
   
};
#endif
